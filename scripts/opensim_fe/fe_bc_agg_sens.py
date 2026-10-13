# -*- coding: utf-8 -*-
"""FE_BC · 任务①②：口径敏感性 + 修好的「侧向约束」诊断（**零新增解算**）。

读 `temp/opensim_fe/fe_bc_perheight/femur_r/h*/femur_r.hdf5`（G9 口径侧向 `[(0.31,0.69)]`
的逐高度真解），对每个高度算**多种聚合口径**，看 A2-股骨的判定对口径有多敏感。

① 口径敏感性（候选 1：`fem_cort` shaft-median 对**弯曲主导**部位是否系统性偏低）
     - shaft_median（现行）｜shaft_p95｜shaft_max
     - whole_median｜whole_p95｜whole_max（整个 `fem_cort`）
     - 阈值 220 MPa 不变 ⇒ 看**哪个口径在哪个高度跨越阈值** ⇒ 判定是否翻成「阶跃」

② 修好的诊断 `vm_over_abs_szz`：**限定在 mask 内**（此前全域算得 1.97，越界不可用）
     - 对照：≈1.0 自由（单轴应力）／≈0.571 全约束（单轴应变 ν=0.3）
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

import paper_alignment as pa  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
FE_DIR = ROOT / "temp" / "opensim_fe" / "fe_bc_perheight" / "femur_r"
OUT = RES / "FE_BC_AGG_SENS.json"
DOM = "fem_cort"
SIG_S = 220.0
SIG_E = 80.0


def load(h5: Path):
    import h5py

    with h5py.File(h5, "r") as h:
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
        conn = h["meshes/0/domains"][DOM][:]
        sid = sorted(h["states"].keys(), key=int)[-1]
        stress = h[f"states/{sid}/element_data/stress/{DOM}"][:].astype(float)
    return coords, conn, stress


def main() -> int:
    fmj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    heights = np.asarray(fmj["heights_m"], dtype=float)
    rows = []
    for d in sorted(FE_DIR.glob("h*")):
        h5 = d / "femur_r.hdf5"
        if not h5.is_file():
            continue
        i = int(d.name[1:]) - 1
        coords, conn, stress = load(h5)
        idx = conn[:, 1:] - 1
        cent = coords[idx].mean(axis=1)
        p = pa.partition_ends_shaft(coords, cent)
        em, sm = p["ends_mask"], p["shaft_mask"]
        vm = pa.vm_voigt(stress)
        az = np.abs(stress[:, 2])
        with np.errstate(invalid="ignore", divide="ignore"):
            r_shaft = (np.median(vm[sm]) / np.median(az[sm])
                       if np.median(az[sm]) > 0 else float("nan"))
            r_ends = (np.median(vm[em]) / np.median(az[em])
                      if np.median(az[em]) > 0 else float("nan"))
            r_whole = (np.median(vm) / np.median(az)
                       if np.median(az) > 0 else float("nan"))
        rows.append({
            "i": i, "h_m": float(heights[i]),
            "load_n": float(fmj["loads_n"]["hip_r"][i]),
            "shaft_median": float(np.median(vm[sm])),
            "shaft_p95": float(np.percentile(vm[sm], 95)),
            "shaft_max": float(vm[sm].max()),
            "whole_median": float(np.median(vm)),
            "whole_p95": float(np.percentile(vm, 95)),
            "whole_max": float(vm.max()),
            "ends_median": float(np.median(vm[em])),
            "r_vm_szz_shaft": float(r_shaft),
            "r_vm_szz_ends": float(r_ends),
            "r_vm_szz_whole": float(r_whole),
        })

    print(f"{'h':>5} {'load_N':>9} | {'shaft_med':>9} {'shaft_p95':>9} {'shaft_max':>9} "
          f"| {'whole_p95':>9} {'whole_max':>9} | {'vm/szz':>7} {'vm/szz_s':>8}")
    for r in rows:
        print(f"{r['h_m']:>5.0f} {r['load_n']:>9.0f} | {r['shaft_median']:>9.2f} "
              f"{r['shaft_p95']:>9.2f} {r['shaft_max']:>9.2f} | "
              f"{r['whole_p95']:>9.2f} {r['whole_max']:>9.2f} | "
              f"{r['r_vm_szz_whole']:>7.3f} {r['r_vm_szz_shaft']:>8.3f}")

    print("\n=== ① 口径敏感性：哪个口径在哪个高度跨 220 MPa ===")
    verdicts = {}
    for key, thr in (("shaft_median", SIG_S), ("shaft_p95", SIG_S),
                     ("shaft_max", SIG_S), ("whole_p95", SIG_S),
                     ("whole_max", SIG_S)):
        cross = [r["h_m"] for r in rows if r[key] >= thr]
        verdicts[key] = {"threshold": thr, "n_cross": len(cross),
                         "first_cross_h_m": cross[0] if cross else None,
                         "max_value": max(r[key] for r in rows)}
        print(f"  {key:14s} (阈{thr:g}) 跨越 {len(cross):2d}/{len(rows)} 个高度"
              f"  首次={cross[0] if cross else None}  最大={verdicts[key]['max_value']:.1f}")

    print("\n=== ② 修好的诊断（限定 mask）：≈1.0 自由 / ≈0.571 全约束 ===")
    for key in ("r_vm_szz_shaft", "r_vm_szz_ends", "r_vm_szz_whole"):
        v = [r[key] for r in rows if np.isfinite(r[key])]
        print(f"  {key:18s} min={min(v):.4f}  median={np.median(v):.4f}  "
              f"max={max(v):.4f}")

    OUT.write_text(json.dumps({
        "bone": "femur_r", "domain": DOM, "bc": {
            "lateral_support": [0.31, 0.69], "note": "G9-faithful free-segment band",
            "scaling": "per-height re-solve"},
        "thresholds": {"shaft": SIG_S, "ends": SIG_E},
        "aggregation_verdicts": verdicts,
        "rows": rows, "reads_frozen_only": False,
        "note_2": "r_vm_szz_* = median(VM)/median(|szz|) inside each mask",
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
