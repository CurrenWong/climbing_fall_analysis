"""G4 — FAITHFUL first-failure + damage-onset exports (viz2).

The P1-a/P1-b ``deformed_first_failure.vtk`` came from the run with the fixed
preset ``disp@-397.858`` whose adaptive stepper jumped across the inversion
onset -> a post-buckling fold (max|u|=217.8 mm).  This script produces a
**minimal-overshoot** first-failure state by (a) using a **fine displacement
ramp** (N small increments) and (b) **bisecting the preset delta** to the
smallest delta* that still diverges (between the fully-converging elastic-1000 N
delta = 198.929 mm and the known-failing 397.858 mm).  The last converged state
of that delta* run is the closest exportable state to element inversion.

It also exports the **damage-onset** state: the first increment at which any
``fib_cort`` element has D>0.

Reuses ``g4_fe_viz`` (geometry/VTK/plot/Jacobian helpers) and ``g4_fe_stab``
(``Setup`` + ``cell_disp``); modifies no existing module (only an opt-in
``suptitle`` kwarg was added to two ``g4_fe_viz`` plot helpers).

Reproduce (repo root, PYTHONPATH=src)::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_viz2.py
    # --steps N  (fine-ramp increments, default 100)   --iters K (bisection, default 6)
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import g4_fe_viz as gv  # noqa: E402
import g4_fe_stab as stab  # noqa: E402
from climbing.coupling.fe_post import build_hdf5  # noqa: E402

OUT = ROOT / "results" / "opensim_fe" / "g4_fe_viz"
VIZ_TMP = ROOT / "temp" / "opensim_fe" / "g4_fe_viz"

#: elastic top |uz| at 1000 N (P1-a calibration)
ELASTIC_DELTA = 198.9288998203423
#: the known-failing preset delta (P1-a `disp@-397.858`)
KNOWN_FAIL = 397.8577996406846


# ===========================================================================
# state readers
# ===========================================================================
def _read_state_at(h5p: Path, sid: int) -> dict:
    import h5py
    with h5py.File(h5p, "r") as h5:
        st = h5[f"states/{sid}"]
        try:
            t = float(np.asarray(st.attrs.get("time", np.nan)).reshape(-1)[0])
        except Exception:  # noqa: BLE001
            t = float("nan")
        nd = st["node_data"]
        u = None
        for nm in nd:
            if "displacement" in nm.lower():
                it = nd[nm]
                u = np.asarray(it if isinstance(it, h5py.Dataset) else it[list(it)[0]],
                               dtype=np.float64)
                break
        damage, vm = {}, {}
        ed = st.get("element_data", {})
        if "damage" in ed:
            for dom in ed["damage"]:
                damage[dom] = np.asarray(ed["damage"][dom], dtype=np.float64).reshape(-1)
        if "stress" in ed:
            for dom in ed["stress"]:
                vm[dom] = gv.von_mises_from_voigt(np.asarray(ed["stress"][dom]))
        sids = sorted(h5["states"].keys(), key=int)
    return {"sid": int(sid), "time": t, "u": u, "damage": damage, "vm": vm,
            "n_states": len(sids), "last_sid": int(sids[-1])}


def find_damage_onset(h5p: Path, domain: str = "fib_cort") -> dict | None:
    """First state (smallest displacement) where any `domain` element has D>0."""
    import h5py
    with h5py.File(h5p, "r") as h5:
        sids = sorted(h5["states"].keys(), key=int)
        for sid in sids:
            ed = h5[f"states/{sid}"].get("element_data", {})
            if "damage" not in ed or domain not in ed["damage"]:
                continue
            d = np.asarray(ed["damage"][domain], dtype=np.float64).reshape(-1)
            if d.max() > 0.0:
                nd = h5[f"states/{sid}/node_data"]
                u = None
                for nm in nd:
                    if "displacement" in nm.lower():
                        it = nd[nm]
                        u = np.asarray(it if isinstance(it, h5py.Dataset)
                                       else it[list(it)[0]], dtype=np.float64)
                        break
                return {"sid": int(sid), "d_max": float(d.max()),
                        "n_damaged": int((d > 1e-12).sum()),
                        "argmax_elem": int(np.argmax(d)),
                        "max_abs_u_mm": float(np.linalg.norm(u, axis=1).max()),
                        "u": u}
    return None


# ===========================================================================
# run / export
# ===========================================================================
def run_dc(setup, delta: float, steps: int, timeout: float):
    tag = f"faith{abs(delta):.3f}".replace(".", "p")
    rec = stab.cell_disp(setup, disp=-abs(delta), Dmax=1.0, steps=steps,
                         timeout=timeout, tag=tag)
    feb = stab.OUT_TEMP / f"disp_{tag}.feb"
    return feb, rec


def _hdf5_for(feb: Path) -> Path:
    h5p = VIZ_TMP / f"{feb.stem}.hdf5"
    if not h5p.is_file():
        build_hdf5(feb.with_suffix(".xplt"), hdf5_path=h5p)
    return h5p


def _jacobian_stats(deformed, blocks, pts):
    jd = gv.element_jacobians(deformed, blocks)
    ji = gv._initial_jacobians(pts, blocks)
    neg_total, min_ratio, by_dom = 0, np.inf, {}
    for b in blocks:
        denom = np.where(np.abs(ji[b["name"]]) > 1e-30, np.abs(ji[b["name"]]), np.nan)
        r = jd[b["name"]] / denom
        n = int((jd[b["name"]] < 0).sum())
        neg_total += n
        by_dom[b["name"]] = n
        if len(r):
            min_ratio = min(min_ratio, float(np.nanmin(r)))
    return {"n_neg_jacobian": neg_total, "min_j_ratio": float(min_ratio),
            "neg_by_domain": by_dom}


def export(feb: Path, blocks, pts, sets, sid: int, stem_out: str, *,
           png_kind: str, suptitle: str) -> dict:
    h5p = _hdf5_for(feb)
    st = _read_state_at(h5p, sid)
    u = st["u"]
    deformed = pts + u
    damage = {b["name"]: st["damage"].get(b["name"], np.zeros(len(b["conn"])))
              for b in blocks}
    vm = {b["name"]: st["vm"].get(b["name"], np.zeros(len(b["conn"])))
          for b in blocks}
    m = gv._meshio_mesh(deformed, blocks, damage, vm, u=u, sets=sets)
    gv.write_vtk(m, OUT / f"{stem_out}.vtk")
    jac = _jacobian_stats(deformed, blocks, pts)
    if png_kind == "damage":
        gv.plot_damage(pts, blocks, u, damage, sets, OUT / f"{stem_out}.png",
                       suptitle=suptitle)
    else:
        gv.plot_overlay(pts, blocks, u, sets, OUT / f"{stem_out}.png",
                        suptitle=suptitle)
    return {"feb": str(feb), "state_index": sid, "n_states": st["n_states"],
            "time": st["time"], "max_abs_u_mm": float(np.linalg.norm(u, axis=1).max()),
            "max_damage": float(max((a.max() if a.size else 0.0)
                                    for a in damage.values())),
            **jac}


# ===========================================================================
# main
# ===========================================================================
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--steps", type=int, default=100, help="fine-ramp increments")
    ap.add_argument("--iters", type=int, default=6, help="bisection iterations")
    ap.add_argument("--timeout", type=float, default=300.0)
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)

    print("[viz2] setup ...", flush=True)
    setup = stab.Setup()
    setup.fill_materials()
    probe = VIZ_TMP / "probe_geom.feb"
    setup.build_elastic(probe, 1.0, 1)          # identical mesh, tiny load
    pts, nids, blocks, sets = gv.parse_feb_geometry(probe)
    print(f"[viz2] mesh {pts.shape} "
          f"blocks={[(b['name'], b['type'], len(b['conn'])) for b in blocks]}",
          flush=True)

    # ---- bisection: smallest preset delta that still diverges --------------
    print(f"\n[viz2] bisect preset delta in ({ELASTIC_DELTA:.3f}, {KNOWN_FAIL:.3f}), "
          f"N={args.steps} ...", flush=True)
    lo, hi = ELASTIC_DELTA, KNOWN_FAIL
    best_fail: tuple[Path, float] | None = None
    bisect_log = []
    for it in range(args.iters):
        mid = 0.5 * (lo + hi)
        feb, rec = run_dc(setup, mid, args.steps, args.timeout)
        ok = bool(rec.get("ok"))
        last_t = np.nan
        if ok:
            try:
                h5p = _hdf5_for(feb)
                last_t = _read_state_at(h5p, _read_state_at(h5p, 0)["last_sid"])["time"]
            except Exception:  # noqa: BLE001
                last_t = np.nan
        diverged = (not ok) or (np.isfinite(last_t) and last_t < 0.999)
        bisect_log.append({"delta": float(mid), "rc": rec.get("rc"),
                           "ok": ok, "last_time": float(last_t) if np.isfinite(last_t) else None,
                           "diverged": bool(diverged),
                           "wall_clock_s": rec.get("wall_clock_s"),
                           "error": (rec.get("error_reason") or "")[:120]})
        print(f"[viz2] it{it}: delta={mid:8.3f} rc={rec.get('rc')} "
              f"last_t={last_t if np.isfinite(last_t) else float('nan'):.4f} "
              f"-> {'DIVERGED' if diverged else 'converged'} "
              f"({rec.get('wall_clock_s', 0):.1f}s)", flush=True)
        if diverged:
            hi = mid
            if best_fail is None or mid < best_fail[1]:
                best_fail = (feb, mid)
        else:
            lo = mid
    delta_star = hi
    print(f"[viz2] delta* = {delta_star:.4f} mm "
          f"(overshoot vs elastic delta {ELASTIC_DELTA:.3f} mm = "
          f"{delta_star / ELASTIC_DELTA:.3f}x)", flush=True)

    results: dict = {"generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                     "elastic_delta_mm": ELASTIC_DELTA,
                     "known_fail_delta_mm": KNOWN_FAIL,
                     "delta_star_mm": float(delta_star),
                     "overshoot_factor": float(delta_star / ELASTIC_DELTA),
                     "steps": args.steps, "bisect_log": bisect_log}

    # ---- faithful first failure (last converged state of the delta* run) ----
    if best_fail is None:
        print("[viz2] NO diverging run found in the bracket -> cannot export "
              "a faithful first-failure state.", flush=True)
        results["faithful"] = {"error": "no diverging run in bracket"}
    else:
        feb_f, d_f = best_fail
        h5p = _hdf5_for(feb_f)
        last_sid = _read_state_at(h5p, 0)["last_sid"]
        last_t = _read_state_at(h5p, last_sid)["time"]
        results["faithful"] = export(
            feb_f, blocks, pts, sets, last_sid, "deformed_first_failure_faithful",
            png_kind="overlay",
            suptitle=(f"fibula_r 忠实首失败：preset δ*={d_f:.1f} mm，末收敛 state "
                      f"{last_sid}（t={last_t:.4f}，min detJ/|detJ0|→0）"))
        results["faithful"]["delta_star_mm"] = float(d_f)
        results["faithful"]["overshoot_factor"] = float(d_f / ELASTIC_DELTA)
        print(f"[viz2] faithful: run δ*={d_f:.2f}, state {last_sid}, "
              f"max|u|={results['faithful']['max_abs_u_mm']:.2f} mm, "
              f"negJ={results['faithful']['n_neg_jacobian']}, "
              f"minJ={results['faithful']['min_j_ratio']:.5f}", flush=True)

        # ---- damage onset (scan all states of the same run) ----------------
        onset = find_damage_onset(h5p)
        if onset:
            # need the onset state's full fields for the VTK
            results["damage_onset"] = export(
                feb_f, blocks, pts, sets, onset["sid"], "damage_onset",
                png_kind="damage",
                suptitle=(f"fibula_r 损伤起始：state {onset['sid']} "
                          f"(Dmax={onset['d_max']:.3g}, {onset['n_damaged']} elems, "
                          f"max|u|={onset['max_abs_u_mm']:.2f} mm)"))
            results["damage_onset"].update({
                "d_max": onset["d_max"], "n_damaged": onset["n_damaged"],
                "argmax_elem": onset["argmax_elem"], "max_abs_u_mm": onset["max_abs_u_mm"]})
            print(f"[viz2] onset: state {onset['sid']}, max|u|={onset['max_abs_u_mm']:.2f} mm, "
                  f"Dmax={onset['d_max']:.4g}, n_damaged={onset['n_damaged']}", flush=True)
        else:
            print("[viz2] no damage-onset state found in this run.", flush=True)
            results["damage_onset"] = {"error": "no D>0 state found"}

    (VIZ_TMP / "g4_fe_viz2.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    _update_readme(results)
    print(f"\n[viz2] done -> {OUT}", flush=True)
    return 0


def _update_readme(R: dict) -> None:
    p = OUT / "README.md"
    base = p.read_text(encoding="utf-8") if p.is_file() else "# g4_fe_viz\n"
    base = re.split(r"\n## 忠实首失败", base)[0].rstrip() + "\n"
    A = base.split("\n")
    A.append("")
    A.append("## 忠实首失败（viz2：最小过冲）+ 损伤起始")
    A.append("")
    A.append(f"_生成时间：{R['generated_at']}。脚本：`scripts/opensim_fe/g4_fe_viz2.py`。_")
    A.append("")
    A.append(f"- **δ\\*** = **{R['delta_star_mm']:.3f} mm**"
             f"（在弹性 δ={R['elastic_delta_mm']:.3f} mm 与已知发散 δ="
             f"{R['known_fail_delta_mm']:.1f} mm 之间二分得到的最小发散 preset δ；"
             f"**过冲 = {R['overshoot_factor']:.3f}×**，即仅比弹性-1000 N 位移大 "
             f"{(R['overshoot_factor']-1)*100:.1f}%）。")
    fa = R.get("faithful") or {}
    if fa.get("state_index") is not None:
        A.append(f"- **忠实首失败态**：preset δ*={fa.get('delta_star_mm'):.2f} mm 运行的"
                 f"**末收敛 state {fa['state_index']}**（共 {fa['n_states']}；"
                 f"t={fa['time']:.4f}）→ `deformed_first_failure_faithful.vtk`。")
        A.append(f"  - `max|u| = {fa['max_abs_u_mm']:.2f} mm`；"
                 f"**负 Jacobian 单元数 = {fa['n_neg_jacobian']}**"
                 f"（min detJ/|detJ0| = {fa['min_j_ratio']:.5f}，即压到 "
                 f"{100*fa['min_j_ratio']:.2f}% 体积、即将反转）。")
    do = R.get("damage_onset") or {}
    if do.get("state_index") is not None:
        A.append(f"- **损伤起始态**：state {do['state_index']}（`fib_cort` 首个 D>0；"
                 f"`max|u| = {do['max_abs_u_mm']:.2f} mm`，`Dmax = {do['d_max']:.3g}`，"
                 f"{do['n_damaged']} 个单元）→ `damage_onset.vtk`。")
    A.append("")
    A.append("### 新增图（各 1–2 行说明）")
    A.append("")
    A.append("- `overlay_first_failure_faithful.png`：δ* 运行**末收敛态**——初始灰线框 vs "
             "变形网格（按 |u| 着色）；折叠幅度最小、且含最接近反转的压溃单元（min detJ→0）。")
    A.append("- `damage_onset.png`：**损伤起始**态（`fib_cort` 首次 D>0）按单元 D 填充着色；"
             "此时损伤仅出现在点载荷热点，网格几乎未折叠。")
    A.append("")
    A.append("### δ* 二分记录")
    A.append("")
    A.append("| preset δ (mm) | rc | last t | 结果 | wall (s) |")
    A.append("|---|---|---|---|---|")
    for b in R.get("bisect_log", []):
        lt = "—" if b.get("last_time") is None else f"{b['last_time']:.4f}"
        A.append(f"| {b['delta']:.3f} | {b.get('rc')} | {lt} | "
                 f"{'发散' if b['diverged'] else '收敛'} | "
                 f"{(b.get('wall_clock_s') or 0):.1f} |")
    A.append("")
    p.write_text("\n".join(A) + "\n", encoding="utf-8")
    print(f"[readme] updated {p}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
