# -*- coding: utf-8 -*-
"""T6 注册修正后：**指定高度**的 FE 实测（把首折从两点外推变实测）。

用途：A2 的 5/5 目前建立在「V1 两点（h=10/30）+ 严格线弹性缩放」之上。
本脚本补 h=5 / h=6 / h=17 三档实测：
  h=5  = 官方矩阵的参考高度（`fracture_matrix.py` 用 h=5 单点解 + 缩放）
  h=6  = 修正后的首折点（相变窗）
  h=17 = A1 的目标窗中心（用于确认「靠载荷降/升能否动」以外的实测锚）
变量与 V1 逐字一致：只注入压板 z 注册修正，其余不动。
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import g4_fe_plate as plate            # noqa: E402
import p2_head_contact as p2           # noqa: E402
import p5_head_contact as p5           # noqa: E402

_ORIG = p5.plate_geom_from_mesh


def regfix(nodes_rot, geom):
    pg = _ORIG(nodes_rot, geom)
    z_top = float(max(nodes_rot[np.asarray(f["poly"], int) - 1, 2].max()
                      for f in geom["cap"]))
    zw = z_top + p5.GAP_MM
    pg["z_wall"] = zw
    pg["z_bot"] = zw - plate.PLATE_THICK_MM
    pg["z_mid"] = zw - 0.5 * plate.PLATE_THICK_MM
    pg["regfix"] = {"a_max": float(geom["a_max"]), "z_cap_top": z_top,
                    "shift_mm": float(_ORIG(nodes_rot, geom)["z_wall"] - zw)}
    return pg


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--heights", nargs="*", type=int, default=[5, 6, 17])
    ap.add_argument("--workdir", default="p2_t6_heights")
    args = ap.parse_args()

    p5.plate_geom_from_mesh = regfix
    p2.BASE = ROOT / "temp" / "opensim_fe" / args.workdir
    SIGMA_C = 150.0

    print("=" * 92)
    print("T6 注册修正 · 指定高度实测   heights=%s" % args.heights)
    print("=" * 92, flush=True)
    t0 = time.time()
    res = p2.run_bone("T6", args.heights, force_loads=None, with_shell=True)
    print("  (run_bone 耗时 %.1fs)" % (time.time() - t0), flush=True)

    out = {"probe": "T6 regfix 指定高度实测", "heights": args.heights,
           "sigma_c_mpa": SIGMA_C, "rows": []}
    for r in res["rows"]:
        d = p2.BASE / "T6" / ("h%02d" % r["h_no"])
        h5 = d / "T6.hdf5"
        rec = {"h_no": r["h_no"], "load_n": r["load_n"],
               "penalty": r["penalty"], "converged": bool(r.get("converged")),
               "error": r.get("error")}
        if h5.is_file():
            with h5py.File(h5, "r") as h:
                doms = {k: int(h[f"meshes/0/domains/{k}"].shape[0])
                        for k in h["meshes/0/domains"].keys()}
            key = next((k for k, n in doms.items() if n == 712), None)
            if key:
                st = p2.domain_stat(h5, key)
                rec.update({"domain_h5": key,
                            "sigma_median": st.get("median_mpa"),
                            "sigma_p95": st.get("p95_mpa"),
                            "sigma_max": st.get("max_mpa"),
                            "szz_min": st.get("szz_min_mpa"),
                            "szz_max": st.get("szz_max_mpa")})
                rec["fracture"] = bool((st.get("median_mpa") or 0) > SIGMA_C)
            try:
                bs = p2.load_bear_stats(h5, r["load_n"], "T6")
                rec.update({"bear_frac": bs["bear_frac"],
                            "R_axis_check": bs["R_axis_check"]})
            except Exception as e:                    # noqa: BLE001
                rec["bear_error"] = str(e)[:100]
        out["rows"].append(rec)
        print("  h=%2d F=%9.1fN conv=%-5s σ_med=%-9s (>σ_c=%s)  bear=%s"
              % (rec["h_no"], rec["load_n"], rec["converged"],
                 ("%.3f" % rec["sigma_median"]) if rec.get("sigma_median") else "—",
                 rec.get("fracture"), rec.get("bear_frac")), flush=True)

    dst = ROOT / "results" / "opensim_fe" / "P2_T6_HEIGHTS_MEASURED.json"
    dst.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=str),
                   encoding="utf-8")
    print("\n产物:", dst.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
