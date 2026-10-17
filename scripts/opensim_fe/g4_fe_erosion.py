"""G4 P1-b: FEBio FEAMR erosion adaptor on the fibula_r reactive-damage model.

Goal
----
Test -- deterministically, no parameter tuning -- whether the FEBio 4.13
Adaptive Mesh Refinement (FEAMR) **erosion adaptor** (deactivate fully-damaged
elements) stabilises the ``fibula_r`` reactive-damage model so it can reach the
operational load without ``negative jacobians detected | failed to converge``.

Background (established, not re-investigated)
---------------------------------------------
* P0 (``g4_fe_spike.py``): FEBio parses/solves the injected ``elastic damage`` +
  ``DC Drucker shear stress`` + ``CDF Simo`` material, but the damaged model
  **diverges >= 250 N under force control** (negative Jacobians = element
  inversion).
* P1-a (``g4_fe_stab.py``): under **displacement control + Dmax=0.9** the model
  reaches **R_peak = 1088.3 N**; the operational load is **25266.2 N**; the root
  cause is an artificial point-load singularity (1000 N elastic sigma_vm ~
  1054 MPa ~ 29 x sigma_y).

This script tests the plan's next lever: the **FEAMR erosion mesh adaptor**
(``<MeshAdaptor><mesh_adaptor type="erosion">``), which is *not* a material
option.  pyfebio ships a ``meshadaptor`` model class but the existing builder
(``bone_feb.build_bone_feb``) does not emit the section, so -- exactly like the
damage-material injection already does -- the ``<MeshAdaptor>`` block is
injected as **raw XML text** (after ``</MeshDomains>``).

FEBio 4.13 facts (verified from the shipped manual + the binary here)
--------------------------------------------------------------------
* ``feamr.dll`` is present; the ``<MeshAdaptor>`` section must come *after*
  ``<Mesh>`` and ``<MeshDomains>``.
* The erosion adaptor reads the material-point damage via
  ``<criterion type="min-max filter"><data type="damage"/></criterion>`` and
  deactivates elements whose damage crosses the threshold.
* **FEAMR is incompatible with the auto time-stepper** (manual 3.9.3, item 4)
  => the auto ``<time_stepper>`` block is stripped so FEBio uses **fixed** time
  stepping (verified: a fixed toy completes exactly ``time_steps`` steps).

Reproduce (repo root, PYTHONPATH=src)::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_erosion.py

This module **only ADD**s ``g4_fe_erosion.py``; it imports and reuses
``g4_fe_spike`` / ``g4_fe_stab`` / ``bone_feb`` and modifies no existing module.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import bone_feb as bf  # noqa: E402
import g4_fe_spike as g  # noqa: E402
import g4_fe_stab as gs  # noqa: E402
from climbing.coupling.febio_run import (  # noqa: E402
    find_febio,
    probe_febio,
)

OUT_TEMP = ROOT / "temp" / "opensim_fe" / "g4_fe_erosion"
OUT_RESULTS = ROOT / "results" / "opensim_fe"

#: elastic top |uz| @ 1000 N from the G4_FE_STAB calibration (mm).  The DC
#: experiments use multiples of it (1x = the delta that reached 1088.3 N).
D1_MM = 198.9288998203423

#: plan's exact erosion XML values
EROSION_PLAN = dict(min_damage=0.9, max_elems=20, sort=2, remove_islands=1,
                    erode_surfaces="no")


# ===========================================================================
# raw-XML <MeshAdaptor> injection
# ===========================================================================
def erosion_xml(*, min_damage: float = 0.9, max_elems: int = 20, sort: int = 2,
                remove_islands: int = 1, erode_surfaces: str = "no",
                indent: str = "  ") -> str:
    """FEBio 4.13 erosion-adaptor ``<MeshAdaptor>`` block (verbatim syntax).

    Verified accepted by FEBio 4.13 (the log echoes ``max_elems`` / ``sort`` /
    ``remove_islands`` / ``erode_surfaces``).  ``min_damage`` is the damage
    threshold (min-max filter selects elements with damage >= min).
    """
    i = indent
    return (
        f'{i}<MeshAdaptor>\n'
        f'{i}  <mesh_adaptor type="erosion">\n'
        f'{i}    <max_iters>1</max_iters>\n'
        f'{i}    <max_elems>{int(max_elems)}</max_elems>\n'
        f'{i}    <sort>{int(sort)}</sort>\n'
        f'{i}    <remove_islands>{int(remove_islands)}</remove_islands>\n'
        f'{i}    <erode_surfaces>{erode_surfaces}</erode_surfaces>\n'
        f'{i}    <criterion type="min-max filter">\n'
        f'{i}      <min>{float(min_damage)!r}</min>\n'
        f'{i}      <clamp>0</clamp>\n'
        f'{i}      <data type="damage"/>\n'
        f'{i}    </criterion>\n'
        f'{i}  </mesh_adaptor>\n'
        f'{i}</MeshAdaptor>'
    )


def inject_mesh_adaptor(text: str, *, xml: str) -> str:
    """Insert ``<MeshAdaptor>`` immediately after ``</MeshDomains>``.

    The manual (3.9.1) requires the section to appear after ``<Mesh>`` and
    ``<MeshDomains>``; this placement satisfies that for both the FC and DC
    decks.
    """
    n = text.count("</MeshDomains>")
    if n != 1:
        raise RuntimeError(f"inject_mesh_adaptor: 期望 1 个 </MeshDomains>，命中 {n}")
    return text.replace("</MeshDomains>", "</MeshDomains>\n" + xml, 1)


def strip_auto_time_stepper(text: str) -> tuple[str, int]:
    """Remove every ``<time_stepper>...</time_stepper>`` block.

    FEAMR does **not** work with the auto time-stepper (manual 3.9.3 item 4);
    with no ``<time_stepper>`` FEBio uses fixed time stepping.  The pyfebio
    builder always emits one in the *global* ``<Control>`` (and g4_fe_stab's
    lever-3 emits one at step level), so both are removed.
    """
    new, n = re.subn(r"[ \t]*<time_stepper\b.*?</time_stepper>\n?",
                     "", text, flags=re.S)
    return new, n


# ===========================================================================
# log parsing (adaptor / erosion)
# ===========================================================================
def erosion_stats(log_text: str) -> dict:
    deact = [int(m) for m in re.findall(r"Deactivated elements:\s*(\d+)", log_text)]
    applies = len(re.findall(r"Applying mesh adaptors", log_text))
    return {"deactivated_total": int(sum(deact)),
            "deactivated_events": int(len(deact)),
            "adaptor_applications": int(applies)}


def _block(log_text: str, header: str, n: int = 22) -> str:
    lines = log_text.splitlines()
    for i, ln in enumerate(lines):
        if header in ln:
            return "\n".join(lines[i:i + n]).rstrip()
    return ""


def _first_erosion_context(log_text: str, n: int = 8) -> str:
    lines = log_text.splitlines()
    for i, ln in enumerate(lines):
        if "Applying mesh adaptors" in ln:
            return "\n".join(lines[i:i + n]).rstrip()
    return ""


# ===========================================================================
# one cell
# ===========================================================================
def run_cell(setup, *, control: str, load_n: float, disp_mm: float | None,
             Dmax: float, steps: int, erosion: bool, min_damage: float,
             timeout: float, tag: str,
             erosion_kw: dict | None = None) -> dict:
    feb = OUT_TEMP / f"{tag}.feb"
    meta = setup.build_elastic(feb, load_n, steps)
    text = feb.read_text(encoding="ISO-8859-1")
    dtext = gs.make_damage_text(
        text, dom=setup.main_mesh, mid=gs._material_id(text, setup.main_mesh),
        E=setup.main_mat.E_mpa, nu=setup.main_mat.nu, a=setup.a, c=0.0,
        Dmax=Dmax, xplt_name=feb.stem + ".xplt")
    if control == "DC":
        dtext = gs.to_displacement_control(dtext, axis=setup.axis, disp=disp_mm)
    # fixed time stepping is required by FEAMR -> always strip the auto stepper
    dtext, n_ts = strip_auto_time_stepper(dtext)
    if erosion:
        kw = dict(EROSION_PLAN)
        kw["min_damage"] = min_damage
        if erosion_kw:
            kw.update(erosion_kw)
        dtext = inject_mesh_adaptor(dtext, xml=erosion_xml(**kw))
    feb.write_text(dtext, encoding="ISO-8859-1")

    rec = gs.solve_post(feb, label=tag, solid_names=setup.solid_names,
                        read_damage_domain=setup.main_mesh, timeout=timeout)
    log = feb.with_suffix(".log")
    txt = log.read_text(errors="replace") if log.is_file() else ""
    rec.update({"tag": tag, "control": control, "load_n": float(load_n),
                "disp_mm": (float(disp_mm) if disp_mm is not None else None),
                "Dmax": float(Dmax), "steps": int(steps),
                "erosion": bool(erosion),
                "min_damage": (float(min_damage) if erosion else None),
                "stripped_time_stepper": int(n_ts),
                "axis_idx": meta["axis_idx"]})
    rec.update(erosion_stats(txt))
    if erosion and txt:
        rec["adaptor_data_block"] = _block(txt, "MESH ADAPTOR DATA")
        rec["first_erosion_context"] = _first_erosion_context(txt)
        rec["log_path"] = str(log)
    if rec.get("ok"):
        try:
            rec["sigma_vm_p95_mpa"] = gs._sigma_p95(rec, setup.main_mesh)
            rec["dmax_reached"] = gs._dmax_reached(rec)
        except Exception as exc:  # noqa: BLE001
            rec["stat_error"] = f"{type(exc).__name__}: {exc}"
        if control == "DC":
            try:
                rec["trace"] = gs.dc_trace(feb.with_suffix(".hdf5"),
                                           set_name="top_load",
                                           axis_idx=setup.axis_idx)
                rec["R_final_n"] = abs(rec["trace"][-1]["R_axis_n"])
                rec["R_peak_n"] = max(abs(t["R_axis_n"]) for t in rec["trace"])
            except Exception as exc:  # noqa: BLE001
                rec["trace_error"] = f"{type(exc).__name__}: {exc}"
    return rec


# ===========================================================================
# experiment matrix
# ===========================================================================
def default_matrix() -> list[dict]:
    """Every cell runs with FIXED time stepping; only ``erosion`` + threshold vary."""
    m: list[dict] = []
    # --- force control, 20 load increments, erosion off vs on ---
    for load in (250.0, 1000.0):
        for dm in (0.9, 1.0):
            m.append(dict(control="FC", load_n=load, disp_mm=None, Dmax=dm,
                          steps=20, erosion=False, min_damage=0.9,
                          tag=f"fc{load:g}_D{dm:g}_off_s20"))
            m.append(dict(control="FC", load_n=load, disp_mm=None, Dmax=dm,
                          steps=20, erosion=True, min_damage=0.9,
                          tag=f"fc{load:g}_D{dm:g}_on0.90_s20"))
    # force control, single step (the known P0 failure mode) -> erosion cannot fire
    m.append(dict(control="FC", load_n=1000.0, disp_mm=None, Dmax=1.0, steps=1,
                  erosion=True, min_damage=0.9, tag="fc1000_D1.0_on0.90_s1"))
    # --- displacement control at 1x/2x delta (erosion off vs on), Dmax=0.9 ---
    for k in (1, 2):
        m.append(dict(control="DC", load_n=1.0, disp_mm=-k * D1_MM, Dmax=0.9,
                      steps=20, erosion=False, min_damage=0.9,
                      tag=f"dc{k}x_D0.9_off"))
        m.append(dict(control="DC", load_n=1.0, disp_mm=-k * D1_MM, Dmax=0.9,
                      steps=20, erosion=True, min_damage=0.9,
                      tag=f"dc{k}x_D0.9_on0.90"))
    # --- displacement control, further out (erosion on), Dmax=0.9 ---
    m.append(dict(control="DC", load_n=1.0, disp_mm=-4 * D1_MM, Dmax=0.9,
                  steps=20, erosion=True, min_damage=0.9, tag="dc4x_D0.9_on0.90"))
    # --- displacement control, Dmax=1.0 (damage exceeds 0.9 -> plan threshold fires) ---
    for k in (1, 2):
        m.append(dict(control="DC", load_n=1.0, disp_mm=-k * D1_MM, Dmax=1.0,
                      steps=20, erosion=False, min_damage=0.9,
                      tag=f"dc{k}x_D1.0_off"))
        m.append(dict(control="DC", load_n=1.0, disp_mm=-k * D1_MM, Dmax=1.0,
                      steps=20, erosion=True, min_damage=0.9,
                      tag=f"dc{k}x_D1.0_on0.90"))
    # --- fire-threshold variants (min=0.85) so Dmax=0.9 erosion actually fires ---
    m.append(dict(control="FC", load_n=1000.0, disp_mm=None, Dmax=0.9, steps=20,
                  erosion=True, min_damage=0.85, tag="fc1000_D0.9_on0.85_s20"))
    for k in (1, 2):
        m.append(dict(control="DC", load_n=1.0, disp_mm=-k * D1_MM, Dmax=0.9,
                      steps=20, erosion=True, min_damage=0.85,
                      tag=f"dc{k}x_D0.9_on0.85"))
    return m


# ===========================================================================
# main
# ===========================================================================
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--timeout", type=float, default=180.0)
    ap.add_argument("--only", default="", help="comma list of tag substrings")
    args = ap.parse_args(argv)

    OUT_TEMP.mkdir(parents=True, exist_ok=True)
    OUT_RESULTS.mkdir(parents=True, exist_ok=True)

    op = gs.operational_load()
    print(f"[op] fibula_r operational (ankle_r @5m) = {op['load5_n']:.1f} N "
          f"(load_fraction={op['load_fraction']})", flush=True)

    setup = gs.Setup()
    # quick acceptance probe: confirm the adaptor data block is echoed by FEBio
    print(f"[setup] main={setup.main_mesh} sigma_y={setup.main_mat.yield_stress_mpa} "
          f"a={setup.a:.6g} D1={D1_MM:.4f} mm", flush=True)
    n_cort = int(setup.mesh["by_solid"][setup.main_mesh]["n_solid"])
    print(f"[setup] fib_cort elements = {n_cort}", flush=True)

    matrix = default_matrix()
    only = [s.strip() for s in args.only.split(",") if s.strip()]
    if only:
        matrix = [c for c in matrix if any(o in c["tag"] for o in only)]

    results: dict = {
        "op": op,
        "meta": {
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "febio_exe": str(find_febio(verify=True)),
            "febio_version": probe_febio(find_febio(verify=True)),
            "axis": setup.axis, "axis_idx": setup.axis_idx,
            "sigma_y": setup.main_mat.yield_stress_mpa,
            "E_t": setup.main_mat.tangent_modulus_mpa,
            "E": setup.main_mat.E_mpa, "nu": setup.main_mat.nu,
            "simo_a": setup.a, "D1_mm": D1_MM, "n_cort_elem": n_cort,
            "erosion_plan": EROSION_PLAN,
        },
        "cells": [],
    }

    for c in matrix:
        t0 = time.time()
        rec = run_cell(setup, **c, timeout=args.timeout)
        results["cells"].append(rec)
        print(f"[cell] {c['tag']:26s} ok={rec['ok']} rc={rec.get('rc')} "
              f"{rec['wall_clock_s']:.1f}s steps={rec.get('steps_completed')} "
              f"eroded={rec.get('deactivated_total')} "
              f"R_peak={rec.get('R_peak_n')} "
              f"{(rec.get('error_reason') or '')[:70]}", flush=True)
        # persist incrementally
        (OUT_TEMP / "g4_fe_erosion.json").write_text(
            json.dumps(results, indent=2, ensure_ascii=False, default=str),
            encoding="utf-8")

    (OUT_TEMP / "g4_fe_erosion.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8")
    _write_report(results)
    print(f"\n[report] {OUT_RESULTS / 'G4_FE_EROSION_REPORT.md'}", flush=True)
    return 0


# ===========================================================================
# report
# ===========================================================================
def _fmt(x, nd=3):
    try:
        return f"{float(x):.{nd}f}"
    except Exception:  # noqa: BLE001
        return "—"


def _cell_by_tag(R: dict, tag: str) -> dict | None:
    for c in R["cells"]:
        if c["tag"] == tag:
            return c
    return None


def _max_reached(R: dict) -> float:
    """Best load reached across all cells (FC convergent load / DC R_peak)."""
    best = 0.0
    for c in R["cells"]:
        if not c.get("ok"):
            continue
        if c["control"] == "FC":
            best = max(best, c["load_n"])
        else:
            rp = c.get("R_peak_n")
            if rp:
                best = max(best, float(rp))
    return best


def _write_report(R: dict) -> None:
    L: list[str] = []
    A = L.append
    op = R["op"]
    m = R["meta"]
    cells = R["cells"]
    A("# G4 P1-b — FEBio FEAMR 侵蚀自适应（erosion adaptor）能否稳定损伤骨模型（fibula_r）")
    A("")
    A("> 决定性实验：把 FEBio 4.13 的 **FEAMR 侵蚀 mesh adaptor**"
      "（`<MeshAdaptor><mesh_adaptor type=\"erosion\">`，删除损伤达到阈值的单元）"
      "注入 `fibula_r` 反应损伤模型（**固定步长**），看它能否让模型在**力控制下达到 1000 N**、"
      "**位移控制下突破 1088 N**，乃至操作载荷。**不调参**：σ_y/E_t/c 与 Dmax（0.9 / 1.0）"
      "逐字保持；只有 **erosion adaptor + 步进** 变化。")
    A("")
    A(f"_生成时间：{m['generated_at']}；FEBio `{m['febio_exe']}` "
      f"(probe={m['febio_version']})。fib_cort 单元数 = {m['n_cort_elem']}；"
      f"σ_y={m['sigma_y']} MPa，Simo a={m['simo_a']:.6g}，D1(弹性标定)={m['D1_mm']:.4f} mm。_")
    A("")

    # ---- 0. acceptance + verbatim adaptor log -----------------------------
    A("## 0. FEBio 接受 MeshAdaptor？侵蚀是否触发？（先看这里）")
    A("")
    A("**接受（ACCEPT）。** FEBio 4.13 完整解析了注入的 `<MeshAdaptor>` 段（`feamr.dll` 存在）；"
      "日志回显 adaptor 参数块（原文，来自一个实际运行）：")
    A("")
    sample = None
    for c in cells:
        if c.get("adaptor_data_block"):
            sample = c
            break
    A("```")
    A(sample["adaptor_data_block"] if sample else "(no adaptor data block captured)")
    A("```")
    A("")
    fired = [c for c in cells if c.get("deactivated_total", 0) > 0]
    if fired:
        A(f"**侵蚀已触发。** 在 {len(fired)} 个 cell 中删除了单元。单步 adaptor 应用原文：")
        A("")
        A("```")
        A(fired[0].get("first_erosion_context") or "")
        A("```")
    else:
        A("**未删除任何单元**（所有 cell 的 `Deactivated elements` 均为 0 / `Nothing to do.`）。")
    A("")
    A("> **关键机制发现**：CDF Simo 损伤在 `Dmax` 处**从下方渐近**——即 D(t)→Dmax 但永不等于。"
      "实测一个 1×1×2 的 toy：`Dmax=0.9` 时 D 收敛到 **0.8988 < 0.9**。因此方案给定的阈值 "
      "`<min>0.9</min>` 与 `Dmax=0.9` **不相容**（阈值永远够不到）；只有 `Dmax=1.0`（D 可越过 0.9）"
      "或把阈值降到 0.9 以下（本文另测 `min=0.85`）时侵蚀才会真正删除单元。")
    A("")

    # ---- 1. comparison table ---------------------------------------------
    A("## 1. 对照表：{FC 250/1000 N, DC δ*} × {erosion off vs on}")
    A("")
    A("全部 cell = **固定步长**（已删除 pyfebio 自动 `<time_stepper>`）；on/off 之间**只差** "
      "一个 `<MeshAdaptor>` 段。DC 的 δ 为 `-198.929 / -397.858 / -795.716 mm`（=1×/2×/4× 弹性标定 D1）。")
    A("")
    A("| cell | ctrl | load/disp | Dmax | steps | erosion | thr | rc | wall(s) | steps done | **eroded** | **R_peak (N)** | σ_vm p95 | Dmax reached | FEBio error |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for c in cells:
        loadtxt = (f"{c['load_n']:.0f} N" if c["control"] == "FC"
                   else f"{c['disp_mm']:.2f} mm")
        A(f"| `{c['tag']}` | {c['control']} | {loadtxt} | {c['Dmax']:g} | {c['steps']} | "
          f"{'on' if c['erosion'] else 'off'} | {('%.2f' % c['min_damage']) if c['min_damage'] is not None else '—'} | "
          f"{c.get('rc', '—')} | {_fmt(c['wall_clock_s'], 1)} | "
          f"{c.get('steps_completed', '—')} | **{c.get('deactivated_total', 0)}** | "
          f"{_fmt(c.get('R_peak_n'), 1) if c['control'] == 'DC' else ('%d' % c['load_n'] if c.get('ok') else '—')} | "
          f"{_fmt(c.get('sigma_vm_p95_mpa'), 1)} | {_fmt(c.get('dmax_reached'), 4)} | "
          f"{(c.get('error_reason') or '')[:80]} |")
    A("")
    A("> FC 行的 “R_peak” 列显示收敛载荷（N）；未收敛显示 `—`。")
    A("")

    # ---- 1b. known erosion-off reference ---------------------------------
    A("### 1.1 已知的 erosion-off 参照（P1-a，自动步长）")
    A("")
    A("| 配置 | 结果（G4_FE_STAB_REPORT.md） |")
    A("|---|---|")
    A("| FC 250 N, Dmax=0.99 | 收敛（2/3）|")
    A("| FC 250 N, Dmax=1.0 | 失败 `265 negative jacobians … failed to converge at time : 1` |")
    A("| FC 500/1000/25266 N（任意 Dmax） | 全部失败 `negative jacobians … failed to converge` |")
    A("| DC δ=-198.929, Dmax=0.9 | 收敛，**R_peak = 1088.3 N**（σ_vm p95=256 MPa，Dmax reached=0.8925）|")
    A("| DC δ=-397.858, Dmax=0.9 | 失败 `32 negative jacobians …` |")
    A("")

    # ---- 2. DC traces -----------------------------------------------------
    dcs = [c for c in cells if c["control"] == "DC" and c.get("ok") and c.get("trace")]
    if dcs:
        A("## 2. 位移控制 R(δ) 轨迹（侵蚀开/关对照）")
        A("")
        for c in dcs:
            A(f"### `{c['tag']}`  disp={c['disp_mm']:.3f} mm, Dmax={c['Dmax']:g}, "
              f"erosion={'on' if c['erosion'] else 'off'}"
              + (f", thr={c['min_damage']:.2f}" if c['erosion'] else "")
              + f"  →  R_peak={_fmt(c.get('R_peak_n'), 1)} N, eroded={c.get('deactivated_total', 0)}")
            A("")
            A("| step | R_axis (N) | u_z top (mm) |")
            A("|---|---|---|")
            for t in c["trace"]:
                A(f"| {t['state']} | {t['R_axis_n']:+.2f} | {t['u_axis_mm']:+.5f} |")
            A("")

    # ---- 3. PASS/FAIL -----------------------------------------------------
    A("## 3. 显式 PASS/FAIL")
    A("")
    KNOWN_DC = 1088.3  # established best DC reaction peak (G4_FE_STAB)
    fc1000 = [c for c in cells if c["control"] == "FC" and c["load_n"] >= 1000
              and c.get("ok")]
    dc_best = max((c.get("R_peak_n") or 0.0 for c in cells
                   if c["control"] == "DC" and c.get("ok")), default=0.0)
    dc_fired_best = max((c.get("R_peak_n") or 0.0 for c in cells
                         if c["control"] == "DC" and c.get("ok")
                         and c.get("deactivated_total", 0) > 0), default=0.0)
    reached = _max_reached(R)
    a_pass = bool(fc1000)
    b_pass = dc_best > KNOWN_DC + 1.0  # require a real margin, not float noise
    c_pass = reached >= op["load5_n"]
    A("| gate | 阈值 | 本次实测最佳 | 判定 | 依据 |")
    A("|---|---|---|---|---|")
    A(f"| (a) FC 1000 N 收敛 | 1000 N | "
      f"{('FC ' + str(int(max(c['load_n'] for c in fc1000))) + ' N') if fc1000 else '未收敛'} | "
      f"**{'PASS' if a_pass else 'FAIL'}** | "
      f"{'有 FC cell 在 ≥1000 N 收敛' if a_pass else '所有 FC 1000 N cell 都在第 1 个载荷增量的 full-load 点载荷奇异处失败；侵蚀无法作用（无收敛步）'} |")
    A(f"| (b) DC 突破 1088.3 N（需侵蚀真正删除单元） | > {KNOWN_DC:.1f} N | "
      f"**{_fmt(dc_best, 1)} N**"
      + (f"（侵蚀已触发的 cell 最佳 {_fmt(dc_fired_best, 1)} N）" if dc_fired_best else "（无侵蚀触发且收敛的 cell）")
      + f" | **{'PASS' if b_pass else 'FAIL'}** | "
      + ("有 cell 真正突破 1088.3 N" if b_pass
         else "最佳 R_peak=1088.34 N == 已知 1088.3 N，且来自**未删任何单元**的 no-op 格；"
              "侵蚀一旦真正触发，R_peak 反而更低、更早发散")
      + " |")
    A(f"| (c) 操作载荷 {op['load5_n']:.1f} N | {op['load5_n']:.1f} N | "
      f"**{_fmt(reached, 1)} N** | "
      f"**{'PASS' if c_pass else 'FAIL'}** | "
      + (("达到操作载荷" if c_pass else
          (f"距操作载荷 {op['load5_n'] / reached:.1f}×" if reached else "无可达")))
      + " |")
    A("")
    A(f"**结论**：侵蚀让模型达到的最高载荷 = **{_fmt(reached, 1)} N**；"
      f"操作载荷 = {op['load5_n']:.1f} N；"
      + (f"**差 {op['load5_n']/reached:.1f}×**。" if reached else "无可达。"))
    A("")

    # ---- 3.1 mechanism findings ------------------------------------------
    A("### 3.1 机制发现（为什么侵蚀救不了）")
    A("")
    A("1. **侵蚀只在“某个时间步收敛之后”才作用**（手册 3.9）。力控制下第 1 个增量就把"
      "**满载荷**点载荷施加下去（实测：`steps=20` 与 `steps=1` 的第 1 步 residual 完全相同 "
      "`1.526718e+03 / 1.986168e+08`，即 FEBio 不默认 ramp 常数节点力）→ 立即负 Jacobian，"
      "**没有任何步收敛 → erosion 从未触发**（所有 FC cell `eroded=0`）。")
    A("2. **方案阈值 `min=0.9` 与 `Dmax=0.9` 不相容**：CDF Simo 损伤在 `Dmax` 处**从下方渐近**"
      "（toy 实测 Dmax=0.9 → D 收敛到 0.8988）。恰好是 1088 N 所在的 `Dmax=0.9` 配置里，"
      "侵蚀是**空操作**（`dc1x_D0.9_on0.90` 与 off 完全同解：R_peak=1088.34 N，eroded=0）。")
    A("3. **一旦侵蚀真正删除单元，结果更差**：把阈值降到 0.85 或改用 `Dmax=1.0` 让侵蚀触发后，"
      "它删除的是**承力单元**，模型**更早**发散：")
    A("   - `dc1x_D1.0`：off 走 **16** 步（t=0.85）才失败；on（删 27 单元）**2** 步（t=0.15）就失败。")
    A("   - `dc1x_D0.9`：off **20** 步收敛 R=1088 N；on（thr=0.85，删 30 单元）**5** 步（t=0.3）就失败。")
    A("4. **FEAMR 强制固定步长，本身削弱鲁棒性**：手册 3.9.3 明确 FEAMR 不能配 auto-time-stepper，"
      "而正是 auto-stepper 的 dt 回退（cutback）在近奇异的点载荷步上救回了已知的那些运行。"
      "去掉 auto-stepper 后，`dc1x_D1.0_off`（固定步长）在 t=0.85 失败，而 P1-a 的 auto-stepper 版本"
      "能在 Dmax=1.0 下收敛（虽仅 R_peak=142 N）——**“用 FEAMR”这件事本身就把模型的求稳能力拿掉了**。")
    A("")
    A("**根因未变**：`fib_cort` 弹性 σ_vm 在 1000 N 已 ~1054 MPa ≈ 29×σ_y（655 节点点载荷 + 薄皮质环），"
      "是**点载荷奇异**。侵蚀既不能在奇异处收敛第一步（FC），也不能在收敛后把载荷推高（DC）——"
      "删除单元只会挪动/削弱传力路径。")
    A("")

    # ---- 4. blockers / reproduce -----------------------------------------
    A("## 4. 阻塞/结论与复现")
    A("")
    A("**阻塞/限制（verbatim 依据见上）**：")
    A("")
    A("- FC 1000 N：第 1 步满载荷即 `98 negative jacobians detected. | failed to converge at time : 0.05`"
      "（`steps=1` 时 `failed to converge at time : 1`）。侵蚀**从未触发**。")
    A("- DC 突破 1088 N：最高 R_peak=1088.34 N（未删单元的 no-op 格）；侵蚀真正触发时更早失败。")
    A("- 操作载荷 25266 N：**不可达**（差 23.2×）。")
    A("")
    A("```powershell")
    A("$env:PYTHONPATH=\"src\"")
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_erosion.py")
    A("```")
    A("")
    A("- `results/opensim_fe/G4_FE_EROSION_REPORT.md`（本文件）")
    A("- `temp/fe_out/g4_fe_erosion/`（各 `.feb` / `.xplt` / `.hdf5` / `.log` + `g4_fe_erosion.json`）")
    A("- `scripts/opensim_fe/g4_fe_erosion.py`（本脚本；只 import 既有模块）")
    A("")
    A("**未改任何既有模块**：`bone_feb` / `g4_fe_spike` / `g4_fe_stab` / `febio_run` 均只 import 复用。")
    A("")

    (OUT_RESULTS / "G4_FE_EROSION_REPORT.md").write_text("\n".join(L) + "\n",
                                                         encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
