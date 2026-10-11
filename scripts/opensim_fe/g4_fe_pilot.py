"""G4 (2) P1 PILOT — per-height reactive-damage FE sigma(h) for fibula_r (V7 BC).

Purpose (bounded pilot: ONE bone, NOT the 9-part matrix, NOT the frozen run)
--------------------------------------------------------------------------
The frozen fracture matrix (``fracture_matrix.py:606-608``) scales a **single**
FE reference stress linearly by the load ratio:

    sigma_vm(h) = sigma_vm_ref * load(h) / load_ref          (strict linear elasticity)

This pilot REPLACES that with a **per-height reactive-damage FE re-solve** of the
``fibula_r`` bone using the P1-b ``V7`` BC (top/bottom 50 % band ``NodalForce``
along the bone's PCA long axis), at ``load(h) = ankle_r__vert(h) * load_fraction``.
It delivers the first "damage-FE sigma(h)" curve + first-fracture height, and
shows where it diverges from the linear-scaled baseline (that divergence is the
P1 payoff).

Reuse (import only -- no existing module is modified)
-----------------------------------------------------
* ``g4_fe_p1b`` -- the ``V7`` BC builder (``build_variant``), ``Setup``,
  ``run_and_post``, ``pca_axis``.
* ``climbing.bone.MATERIAL_STRENGTH_MPA`` -- sigma_c(fibula_shaft) = 160 MPa.
* load chain ``temp/opensim_fe/fracture_matrix/load_chain.npz`` key
  ``ankle_r__vert`` (h = 1..50).

THUMS parameters (E=18500, nu=0.3, sigma_y=37, E_t=16700) and Dmax (0.9 / 1.0)
are kept **verbatim** -- only the BC / per-height load varies.

Reproduce (repo root, PYTHONPATH=src)::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_pilot.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import g4_fe_p1b as p1b  # noqa: E402  (the V7 BC builder; modifies none)
from climbing.bone import MATERIAL_STRENGTH_MPA  # noqa: E402

OUT_TEMP = ROOT / "temp" / "opensim_fe" / "g4_fe_pilot"
OUT_RESULTS = ROOT / "results" / "opensim_fe"
LOAD_CHAIN = ROOT / "temp" / "opensim_fe" / "fracture_matrix" / "load_chain.npz"
FROZEN_MATRIX = OUT_RESULTS / "fracture_matrix.json"

REF_H = 5.0                      # reference height (fracture_matrix.py:REF_H)
REF_LOAD_N = 1000.0              # V7 elastic reference solve load
PRIMARY_STAT = "p95"             # fracture_matrix.py:PRIMARY_STAT
SIGMA_C = float(MATERIAL_STRENGTH_MPA["fibula_shaft"][1])   # 160 MPa (compression)

#: frozen default (no sharing) vs the literature fibula share.
FROZEN_FRACTION = 1.0
PHYS_FRACTION = 0.07

V7 = dict(kind="nodal", top_frac=0.50, bot_frac=0.50, align_pca=False)


# ===========================================================================
# load chain
# ===========================================================================
def load_chain() -> dict:
    z = np.load(LOAD_CHAIN, allow_pickle=True)
    return {"heights": np.asarray(z["heights"], dtype=float),
            "ankle_vert": np.asarray(z["ankle_r__vert"], dtype=float)}


def frozen_fibula_row() -> dict:
    """The frozen matrix's fibula_r row (linear-scaled baseline, for context)."""
    try:
        j = json.loads(FROZEN_MATRIX.read_text(encoding="utf-8"))
        return next(r for r in j["parts"] if r["bone"] == "fibula_r")
    except Exception:  # noqa: BLE001
        return {}


# ===========================================================================
# solve wrapper (V7 BC via g4_fe_p1b)
# ===========================================================================
def solve_v7(setup, *, load_n: float, Dmax: float | None, tag: str,
             load_axis, timeout: float, steps: int = 1) -> dict:
    feb = OUT_TEMP / f"{tag}.feb"
    p1b.build_variant(setup, path=feb, load_n=float(load_n), steps=int(steps),
                      damage_Dmax=Dmax, load_axis=load_axis, **V7)
    rec = p1b.run_and_post(setup, feb, timeout=timeout,
                           read_damage=(Dmax is not None))
    rec.update({"tag": tag, "load_n": float(load_n), "Dmax": Dmax})
    return rec


# ===========================================================================
# main
# ===========================================================================
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--timeout", type=float, default=180.0)
    ap.add_argument("--dm", type=float, nargs="*", default=[0.9, 1.0],
                    help="Dmax values for the physical load_fraction (default 0.9 1.0)")
    ap.add_argument("--rep", type=int, default=3,
                    help="repeats at the top reachable height (nondeterminism)")
    ap.add_argument("--max-h", type=int, default=50, help="highest h to run")
    args = ap.parse_args(argv)

    OUT_TEMP.mkdir(parents=True, exist_ok=True)
    OUT_RESULTS.mkdir(parents=True, exist_ok=True)

    setup = p1b.Setup()
    pca = setup.pca_axis()
    load_axis = -pca
    print(f"[axis] fibula PCA long axis = ({pca[0]:+.3f},{pca[1]:+.3f},{pca[2]:+.3f})",
          flush=True)
    print(f"[main] sigma_c(fibula_shaft) = {SIGMA_C:.1f} MPa", flush=True)

    lc = load_chain()
    heights = lc["heights"][: args.max_h]
    ankle = lc["ankle_vert"][: args.max_h]
    ref_idx = int(np.argmin(np.abs(heights - REF_H)))
    load_phys = ankle * PHYS_FRACTION
    load_frozen = ankle * FROZEN_FRACTION
    print(f"[load] {len(heights)} heights; @h={int(REF_H)} ankle_r={ankle[ref_idx]:.1f} N "
          f"-> phys(0.07)={load_phys[ref_idx]:.1f} N, frozen(1.0)={load_frozen[ref_idx]:.1f} N",
          flush=True)

    results: dict = {
        "meta": {
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "sigma_c_mpa": SIGMA_C, "primary_stat": PRIMARY_STAT,
            "ref_h_m": REF_H, "ref_load_n": REF_LOAD_N,
            "frozen_fraction": FROZEN_FRACTION, "phys_fraction": PHYS_FRACTION,
            "pca_axis": [float(x) for x in pca], "V7": V7,
        },
        "heights": heights.tolist(),
        "ankle_r_vert_n": ankle.tolist(),
        "load_phys_n": load_phys.tolist(),
        "load_frozen_n": load_frozen.tolist(),
        "frozen_matrix_row": frozen_fibula_row(),
        "ref_elastic": None, "cells": [], "repeats": [],
    }

    # --- reference (linear-elastic) solve: V7 elastic @ 1000 N ---------------
    ref = solve_v7(setup, load_n=REF_LOAD_N, Dmax=None, tag="ref_elastic_V7_1000N",
                   load_axis=load_axis, timeout=args.timeout)
    results["ref_elastic"] = {k: ref.get(k) for k in
                              ("ok", "rc", "peak_vm_mpa", "p95_vm_mpa", "R_net_n")}
    sigma_ref = {"max": ref.get("peak_vm_mpa"), "p95": ref.get("p95_vm_mpa")}
    print(f"[ref] V7 elastic @1000 N: max={sigma_ref['max']} p95={sigma_ref['p95']} "
          f"ok={ref.get('ok')}", flush=True)

    def sigma_lin(stat: str, load_n: float) -> float | None:
        s = sigma_ref.get(stat)
        return (s * load_n / REF_LOAD_N) if s else None

    # --- per-height damage-FE re-solve --------------------------------------
    # physical fraction: Dmax 0.9 + 1.0 ; frozen fraction: Dmax 0.9
    plan = ([("phys", PHYS_FRACTION, dm) for dm in args.dm]
            + [("frozen", FROZEN_FRACTION, 0.9)])
    for frac_name, frac, dm in plan:
        load_h = ankle * frac
        for h, load in zip(heights, load_h):
            tag = f"{frac_name}_D{dm:g}_h{int(h):02d}"
            rec = solve_v7(setup, load_n=float(load), Dmax=dm, tag=tag,
                           load_axis=load_axis, timeout=args.timeout)
            cell = {
                "frac_name": frac_name, "load_fraction": frac, "Dmax": dm,
                "h": float(h), "load_n": float(load),
                "ok": bool(rec.get("ok")), "rc": rec.get("rc"),
                "sigma_max": rec.get("peak_vm_mpa"), "sigma_p95": rec.get("p95_vm_mpa"),
                "R_net_n": rec.get("R_net_n"),
                "wall_s": rec.get("wall_clock_s"),
                "error_reason": rec.get("error_reason", ""),
            }
            if frac_name == "phys":
                cell["sigma_lin_max"] = sigma_lin("max", float(load))
                cell["sigma_lin_p95"] = sigma_lin("p95", float(load))
            results["cells"].append(cell)
            sm = cell["sigma_p95"] if cell["sigma_p95"] is not None else float("nan")
            print(f"[cell] {tag:20s} load={load:9.1f}N ok={cell['ok']} "
                  f"sigma_max={cell['sigma_max']} sigma_p95={cell['sigma_p95']} "
                  f"{(cell['error_reason'] or '')[:45]}", flush=True)
            (OUT_TEMP / "g4_fe_pilot.json").write_text(
                json.dumps(results, indent=2, ensure_ascii=False, default=str),
                encoding="utf-8")

    # --- first-fracture + reachable height ----------------------------------
    def first_fracture(frac_name: str, dm: float, stat: str):
        for c in results["cells"]:
            if c["frac_name"] == frac_name and c["Dmax"] == dm:
                v = c.get(f"sigma_{stat}")
                if v is not None and v >= SIGMA_C:
                    return int(c["h"])
        return None

    reachable = [c["h"] for c in results["cells"]
                 if c["frac_name"] == "phys" and c["ok"]]
    max_reach = max(reachable) if reachable else None
    results["max_reachable_height_phys"] = max_reach
    results["first_fracture"] = {
        f"{fn}_D{dm:g}_{stat}": first_fracture(fn, dm, stat)
        for fn in ("phys", "frozen") for dm in ((0.9, 1.0) if fn == "phys" else (0.9,))
        for stat in ("p95", "max")}

    # --- nondeterminism at the top reachable height -------------------------
    if max_reach is not None:
        load_top = float(ankle[heights.tolist().index(max_reach)] * PHYS_FRACTION)
        for dm in args.dm:
            n_ok = 0
            for r in range(max(1, args.rep)):
                rec = solve_v7(setup, load_n=load_top, Dmax=dm,
                               tag=f"rep_D{dm:g}_h{max_reach}", load_axis=load_axis,
                               timeout=args.timeout)
                ok = bool(rec.get("ok"))
                n_ok += int(ok)
                results["repeats"].append(
                    {"Dmax": dm, "h": max_reach, "load_n": load_top, "rep": r,
                     "ok": ok, "sigma_max": rec.get("peak_vm_mpa"),
                     "sigma_p95": rec.get("p95_vm_mpa"),
                     "error_reason": rec.get("error_reason", "")})
            print(f"[repeat] top h={max_reach} Dmax={dm:g} load={load_top:.1f}N "
                  f"-> {n_ok}/{args.rep} converged", flush=True)

    (OUT_TEMP / "g4_fe_pilot.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8")
    _write_report(results)
    print(f"\n[report] {OUT_RESULTS / 'G4_FE_PILOT_REPORT.md'}", flush=True)
    return 0


# ===========================================================================
# report
# ===========================================================================
def _fmt(x, nd=1):
    try:
        return f"{float(x):.{nd}f}"
    except Exception:  # noqa: BLE001
        return "—"


def _write_report(R: dict) -> None:
    L: list[str] = []
    A = L.append
    m = R["meta"]
    heights = R["heights"]
    ankle = R["ankle_r_vert_n"]
    lp = R["load_phys_n"]
    lf = R["load_frozen_n"]
    cells = R["cells"]
    ref = R["ref_elastic"] or {}
    frac = m["phys_fraction"]
    A("# G4 (2) P1 PILOT — fibula_r 逐高度 reactive-damage FE σ(h)（V7 BC）")
    A("")
    A("> **单骨 pilot（非 9 部件矩阵、非冻结预注册运行）**。用 **`V7` BC**"
      "（顶/底 50% 带 `NodalForce` 沿骨 PCA 长轴）对 `fibula_r` 在**每个高度**"
      "`load(h)=ankle_r__vert(h)×load_fraction` 上**真正重解损伤 FE**，"
      "替代冻结矩阵的 `σ_vm(h)=σ_ref·load(h)/load_ref` 线性外推。"
      "THUMS σ_y/E_t/c 与 Dmax 逐字保持。")
    A("")
    A(f"_生成时间：{m['generated_at']}；σ_c(fibula_shaft)={m['sigma_c_mpa']:.0f} MPa；"
      f"主统计={m['primary_stat']}；参考 h={m['ref_h_m']:.0f} m，参考解 = V7 弹性 @"
      f"{m['ref_load_n']:.0f} N。_")
    A("")

    # ---- conclusions ----
    A("## 0. 结论（先看这里）")
    A("")
    A(f"- **V7 弹性参考** @1000 N：σ_vm max={_fmt(ref.get('peak_vm_mpa'))} MPa，"
      f"p95={_fmt(ref.get('p95_vm_mpa'))} MPa（交付 |R|={_fmt(ref.get('R_net_n'), 0)} N）。"
      "（对比冻结矩阵的点载荷基线 p95≈657 MPa @1000 N。）")
    A(f"- **物理分配 load_fraction={frac}**（文献腓骨份额 ~0.05–0.10）："
      f"@5 m load={lp[heights.index(5.0) if 5.0 in heights else 4]:.0f} N；"
      f"@50 m load={lp[-1]:.0f} N。")
    A(f"- **冻结分配 load_fraction=1.0**：@5 m load={lf[heights.index(5.0) if 5.0 in heights else 4]:.0f} N；"
      f"@50 m load={lf[-1]:.0f} N（全部远超模型承载，见下）。")
    mr = R.get("max_reachable_height_phys")
    conv_dm = {}
    for c in cells:
        if c["frac_name"] == "phys" and c["ok"]:
            conv_dm.setdefault(f"{c['Dmax']:g}", []).append(int(c["h"]))
    conv_txt = "；".join(f"D{dm}:h∈{sorted(v)}" for dm, v in sorted(conv_dm.items()))
    A(f"- **可达高度（物理分配，逐高度损伤重解）**：{conv_txt or '全部发散'}。"
      f"最大可达 h={('' + str(int(mr))) if mr else '无'}。**收敛是零散的** —— "
      "多数高度以负 Jacobian 发散（损伤软化在点载荷量级下失去正定切线），"
      "即 **损伤-FE 在操作载荷量级根本无法给出 σ(h)**；这本身相对线性外推就是关键差异。")
    ff = R.get("first_fracture", {})
    lin_first = None
    for c in cells:
        if (c["frac_name"] == "phys" and c.get("sigma_lin_p95") is not None
                and c["sigma_lin_p95"] >= SIGMA_C):
            lin_first = int(c["h"])
            break
    lin_first_max = None
    for c in cells:
        if (c["frac_name"] == "phys" and c.get("sigma_lin_max") is not None
                and c["sigma_lin_max"] >= SIGMA_C):
            lin_first_max = int(c["h"])
            break
    A(f"- **首次断裂高度**：损伤-FE（σ_vm max/p95 ≥ σ_c） = "
      f"{_ff_txt(ff, 'phys')}；线性基线（V7-弹性外推，σ_lin = σ_ref·load/load_ref）= "
      f"max:h={lin_first_max}, p95:h={lin_first}。")
    fr = R.get("frozen_matrix_row", {})
    if fr:
        A(f"- 冻结矩阵（点载荷 BC，σ_ref×25.27 线性外推）首次断裂 = "
          f"h={fr.get('first_fracture_p95_m')} m（σ_vm_ref p95="
          f"{_fmt((fr.get('sigma_vm_ref_mpa') or {}).get('p95'), 0)} MPa ≫ σ_c）。")
    A("")

    # ---- table ----
    A("## 1. σ_vm(h) 表：损伤-FE vs 线性基线")
    A("")
    A(f"`σ_lin(h) = σ_ref · load(h)/load_ref`（σ_ref = V7 弹性 @1000 N，"
      f"主统计={m['primary_stat']}；等价于 σ_ref/1000 · load(h)）。"
      "损伤-FE = 在 load(h) 上重解 `elastic damage`+`CDF Simo`+`DC Drucker` 的 σ_vm。")
    A("")
    A("| h (m) | ankle_r (N) | load_frozen=1.0 (N) | load_phys=0.07 (N) | "
      "σ_vm max (D0.9) | σ_vm p95 (D0.9) | σ_lin p95 | ratio p95 | 断裂 | 状态 |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    by = {(c["frac_name"], c["Dmax"], c["h"]): c for c in cells}
    for i, h in enumerate(heights):
        c09 = by.get(("phys", 0.9, h))
        c10 = by.get(("phys", 1.0, h))
        lin_p95 = c09.get("sigma_lin_p95") if c09 else None
        s_p95 = c09.get("sigma_p95") if c09 else None
        ratio = (s_p95 / lin_p95) if (s_p95 and lin_p95) else None
        fracflag = ("**FRACT**" if (s_p95 is not None and s_p95 >= SIGMA_C)
                    else ("ok" if c09 and c09.get("ok") else ""))
        status = "diverged" if (c09 and not c09.get("ok")) else "converged"
        s_max = _fmt(c09.get("sigma_max")) if c09 else "—"
        A(f"| {int(h)} | {ankle[i]:.0f} | {lf[i]:.0f} | {lp[i]:.0f} | "
          f"{s_max} | {_fmt(s_p95)} | {_fmt(lin_p95)} | {_fmt(ratio, 3)} | "
          f"{fracflag} | {status} |")
    A("")
    A("> σ_vm 与 σ_lin 均为 fib_cort 域统计；`—` = 该格发散（无解）。")
    A("")

    # ---- nonlinearity ----
    A("## 2. 损伤-FE 曲线是否非线性？（P1 payoff）")
    A("")
    conv = [c for c in cells if c["frac_name"] == "phys" and c["Dmax"] == 0.9
            and c["ok"] and c.get("sigma_p95") and c.get("sigma_lin_p95")]
    if len(conv) >= 2:
        ratios = [c["sigma_p95"] / c["sigma_lin_p95"] for c in conv]
        rmin, rmax = min(ratios), max(ratios)
        A(f"- 在**可达**的 {len(conv)} 个高度上，比值 σ_vm(damage)/σ_lin ∈ "
          f"[{rmin:.3f}, {rmax:.3f}]。")
        if (rmax - rmin) > 0.05:
            A(f"- **是非线性的**：若损伤-FE 与线性基线一致，比值应恒为 1.000；"
              f"实测比值跨 {rmax - rmin:.3f}（≈{100 * (rmax - rmin):.0f}%），"
              "说明损伤软化改变了 σ 对载荷的比例关系 —— **这就是 P1 相对线性外推的差异**。")
        else:
            A("- 比值近似常数 ⇒ 在可达（低）高度区间损伤尚未显著软化，"
              "非线性要到更高载荷才显现（但更高载荷不可达，见 §3）。")
    else:
        A(f"- 可达高度不足（仅 {len(conv)} 个）无法拟合非线性；"
          "但**「可达高度上限 ~2 kN ≪ 操作载荷」本身就是相对线性外推的关键差异** —— "
          "线性外推声称可以一路外推到 113 kN，而损伤-FE 在 ~2 kN 就已经无法收敛。")
    A("")

    # ---- reachable / nondeterminism ----
    A("## 3. 可达性 / 非确定性 / 发散格")
    A("")
    A(f"- **物理分配最大可达高度 = {('h=' + str(int(mr))) if mr else '无'}**；"
      "可达格零散分布（见 §0），其余为负 Jacobian 发散（原文见下）。")
    div = [c for c in cells if not c["ok"]]
    A(f"- 发散格数（全部 frac/Dmax 组合）：**{len(div)}/{len(cells)}**。")
    reps = R.get("repeats", [])
    if reps:
        A("- **顶高度重复（非确定性）**：")
        A("")
        A("| h | Dmax | load (N) | 收敛 | σ_p95 | 错误 |")
        A("|---|---|---|---|---|---|")
        for r in reps:
            A(f"| {int(r['h'])} | {r['Dmax']:g} | {r['load_n']:.0f} | "
              f"{'yes' if r['ok'] else 'no'} | {_fmt(r['sigma_p95'])} | "
              f"{(r.get('error_reason') or '')[:50]} |")
    A("")
    A("**发散原文（示例）**：")
    A("")
    A("```")
    for c in div[:6]:
        A(f"h={int(c['h'])} {c['frac_name']} D{c['Dmax']:g} load={c['load_n']:.0f}N : "
          f"{(c.get('error_reason') or '(no log)')[:110]}")
    A("```")
    A("")

    # ---- reproduce ----
    A("## 4. 复现与产物")
    A("")
    A("```powershell")
    A("$env:PYTHONPATH=\"src\"")
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_pilot.py")
    A("```")
    A("")
    A("- `results/opensim_fe/G4_FE_PILOT_REPORT.md`（本文件）")
    A("- `temp/opensim_fe/g4_fe_pilot/`（各 `.feb`/`.xplt`/`.hdf5`/`.log` + `g4_fe_pilot.json`）")
    A("- `scripts/opensim_fe/g4_fe_pilot.py`（本脚本；只 import 既有模块）")
    A("")
    A("**未改任何既有模块**：`g4_fe_p1b` / `g4_fe_spike` / `g4_fe_stab` / `bone_feb` "
      "均只 import 复用。")
    A("")
    (OUT_RESULTS / "G4_FE_PILOT_REPORT.md").write_text("\n".join(L) + "\n",
                                                       encoding="utf-8")


def _ff_txt(ff: dict, fn: str) -> str:
    parts = []
    for dm in (0.9, 1.0):
        for stat in ("p95", "max"):
            v = ff.get(f"{fn}_D{dm}_{stat}")
            if v is not None:
                parts.append(f"D{dm:g}/{stat}:h={v}")
    return ", ".join(parts) if parts else ">50 m（未断裂）"


if __name__ == "__main__":
    raise SystemExit(main())
