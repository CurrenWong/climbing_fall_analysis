# -*- coding: utf-8 -*-
"""T6 · 读全 50 高度的壳域（判据域 shell_t6_cort = h5 的 __unnamed_2）+ 判 0/1 形态。"""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import h5py  # noqa: E402

BASE = ROOT / "temp" / "opensim_fe" / "fe_bc_stab" / "T6"
RES = ROOT / "results" / "opensim_fe"
SIG_C = 150.0

fj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
heights = fj["heights_m"]
loads = fj["loads_n"]["lumbar"]

rows, missing = [], []
for i in range(1, 51):
    h5p = BASE / f"A2_axial_lat_90_h{i:02d}" / "T6.hdf5"
    if not h5p.is_file():
        missing.append(i)
        continue
    with h5py.File(h5p, "r") as h:
        doms = list(h["meshes/0/domains"].keys())
        tgt = next((d for d in doms if "t6_cort" in d), None) or \
            next((d for d in doms if d.startswith("__unnamed")), None)
        sid = sorted(h["states"].keys(), key=int)[-1]
        s = np.asarray(h[f"states/{sid}/element_data/stress/{tgt}"][:], dtype=float)
    vm = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
                        + (s[:, 2] - s[:, 0]) ** 2)
                 + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))
    med = float(np.median(vm))
    rows.append({"h_no": i, "h_m": float(heights[i - 1]),
                 "load_n": float(loads[i - 1]),
                 "shell_median_mpa": med, "ratio": med / SIG_C,
                 "frac": bool(med >= SIG_C)})

print(f"{'h':>3} {'h_m':>5} {'load_N':>9} {'shell_med':>10} {'r/150':>7} frac")
for r in rows:
    print(f"{r['h_no']:>3} {r['h_m']:>5.0f} {r['load_n']:>9.0f} "
          f"{r['shell_median_mpa']:>10.3f} {r['ratio']:>7.3f} "
          f"{'1' if r['frac'] else '0'}")

m = np.array([r["frac"] for r in rows], dtype=int)
n = m.sum()
print(f"\n已解 {len(rows)}/50，缺 {missing}")
print(f"骨折高度数 = {n}/{len(rows)}")
first = next((r["h_m"] for r in rows if r["frac"]), None)
print(f"首次骨折高度 = {first} m")
if n == 0:
    behav = "all_zero"
elif n == len(rows) and rows[0]["frac"]:
    behav = "all_one"
elif n > 0 and not rows[0]["frac"] and all(rows[k]["frac"] for k in range(
        next(k for k, r in enumerate(rows) if r["frac"]), len(rows))):
    behav = "complete_separation"
else:
    behav = "ok(非单调/混合)"
print(f"形态 = {behav}   （常值 = all_zero/all_one；阶跃 = complete_separation）")
print(f"论文靶（胸椎 ∈ PAPER_STEP）= 阶跃")

(RES / "T6_PERHEIGHT.json").write_text(json.dumps({
    "bone": "T6", "cn": "胸椎", "domain": "shell_t6_cort",
    "bc": {"variant": "A2_axial_lat_90", "bottom_fix_axial_only": True,
           "lateral_support": [0.05, 0.95],
           "scaling": "per-height re-solve (no fallback, no extrapolation)"},
    "sigma_c_mpa": SIG_C, "rows": rows, "missing_h": missing,
    "n_solved": len(rows), "n_fracture": n, "first_fracture_h_m": first,
    "behav": behav, "paper_target": "阶跃",
}, ensure_ascii=False, indent=2,
        default=lambda o: o.item() if hasattr(o, "item") else str(o)), encoding="utf-8")
print(f"\n写出 {RES / 'T6_PERHEIGHT.json'}")
