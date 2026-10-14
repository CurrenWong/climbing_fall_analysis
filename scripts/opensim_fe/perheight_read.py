# -*- coding: utf-8 -*-
"""E′ · 通用逐高度读数：按各骨**自己的 main_domain**（含壳域 `__unnamed_*` 兜底）读应力并判形态。

用法::

    python scripts/opensim_fe/perheight_read.py --bone fibula_r \
        --variant A2_axial_lat_90

口径（与 `s4_alignment_bc2.py:159-192` 同源，仅把"解算载荷用哪个域"换成 main_domain）：
* **长骨**：`ends`/`shaft` mask 上的 **median of von Mises**，OR 判据（`s4` 逐字）
* **非长骨**：main_domain 全域 **median of von Mises** ≥ σ_c
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import h5py  # noqa: E402
import paper_alignment as pa  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
BASE = ROOT / "temp" / "opensim_fe" / "fe_bc_stab"


def vm_voigt(s):
    return np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
                          + (s[:, 2] - s[:, 0]) ** 2)
                   + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bone", required=True)
    ap.add_argument("--variant", default="A2_axial_lat_90")
    args = ap.parse_args()
    bone = args.bone

    fj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    heights = fj["heights_m"]
    pr = next(p for p in fj["parts"] if p["bone"] == bone)
    joint = pr["joint"]
    loads = fj["loads_n"][joint]
    sig_c = float(pr["sigma_c_mpa"])
    main_dom = pr["main_domain"]
    reg = next(r for r in pa.REGIONS if r[1] == bone)
    cn, subs = reg[0], reg[3]
    is_long = bone in pa.LONG_BONES
    sig_e = pa._sc(subs[0][1]) if is_long else None
    sig_s = pa._sc(subs[1][1]) if is_long else None

    rows, missing = [], []
    for i in range(1, 51):
        h5p = BASE / bone / f"{args.variant}_h{i:02d}" / f"{bone}.hdf5"
        if not h5p.is_file():
            missing.append(i)
            continue
        with h5py.File(h5p, "r") as h:
            doms = list(h["meshes/0/domains"].keys())
            tgt = (main_dom if main_dom in doms else None) or \
                next((d for d in doms if main_dom and main_dom.split("_")[-1] in d),
                     None) or \
                next((d for d in doms if d.startswith("__unnamed")), None) or doms[0]
            raw = h["meshes/0/nodes"][:]
            coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
            conn = h["meshes/0/domains"][tgt][:]
            sid = sorted(h["states"].keys(), key=int)[-1]
            s = np.asarray(h[f"states/{sid}/element_data/stress/{tgt}"][:], dtype=float)
        vm = vm_voigt(s)
        r = {"h_no": i, "h_m": float(heights[i - 1]), "load_n": float(loads[i - 1]),
             "domain_used": tgt}
        if is_long:
            idx = conn[:, 1:] - 1
            p = pa.partition_ends_shaft(coords, coords[idx].mean(axis=1))
            e, sh = float(np.median(vm[p["ends_mask"]])), float(np.median(vm[p["shaft_mask"]]))
            r.update(ends_stat_mpa=e, shaft_stat_mpa=sh,
                     frac=bool(e >= sig_e or sh >= sig_s),
                     ratio_e=e / sig_e, ratio_s=sh / sig_s)
        else:
            med = float(np.median(vm))
            r.update(whole_median_mpa=med, frac=bool(med >= sig_c),
                     ratio=med / sig_c)
        rows.append(r)

    print(f"=== {bone} ({cn}) · variant={args.variant} · main_domain={main_dom} "
          f"· σ_c={sig_c} ===")
    print(f"{'h':>3} {'load_N':>9} " + ("{:>10} {:>10} {:>7}".format(
        "ends", "shaft", "r_s") if is_long else "{:>10} {:>7}".format("median", "ratio"))
        + "  frac")
    for r in rows:
        if is_long:
            print(f"{r['h_no']:>3} {r['load_n']:>9.0f} {r['ends_stat_mpa']:>10.2f} "
                  f"{r['shaft_stat_mpa']:>10.2f} {r['ratio_s']:>7.3f}  "
                  f"{'1' if r['frac'] else '0'}")
        else:
            print(f"{r['h_no']:>3} {r['load_n']:>9.0f} {r['whole_median_mpa']:>10.3f} "
                  f"{r['ratio']:>7.3f}  {'1' if r['frac'] else '0'}")
    print(f"\n已解 {len(rows)}/50，缺 {missing}")
    if rows:
        m = [r["frac"] for r in rows]
        n = sum(m)
        first = next((r["h_m"] for r in rows if r["frac"]), None)
        if n == 0:
            behav = "all_zero"
        elif n == len(rows) and m[0]:
            behav = "all_one"
        elif n and not m[0]:
            k = next(k for k, x in enumerate(m) if x)
            behav = "complete_separation" if all(m[k:]) else "mixed(非单调)"
        else:
            behav = "mixed(非单调)"
        bmap = {"all_zero": "常值", "all_one": "常值",
                "complete_separation": "阶跃", "mixed(非单调)": "渐变(ok?)"}
        tgt = "阶跃" if cn in pa.PAPER_STEP else "渐进"
        print(f"骨折数={n}/{len(rows)}  首次={first} m  形态={behav}"
              f"（{bmap[behav]}）  论文靶={tgt}"
              f"  => {'✓ 命中' if bmap[behav] == tgt else '✗ 未命中'}")
        out = RES / f"PERHEIGHT_{bone}.json"
        out.write_text(json.dumps({
            "bone": bone, "cn": cn, "variant": args.variant,
            "main_domain": main_dom, "domain_used": rows[0]["domain_used"],
            "sigma_c_mpa": sig_c, "is_long": is_long,
            "rows": rows, "missing_h": missing, "n_solved": len(rows),
            "n_fracture": n, "first_fracture_h_m": first, "behav": behav,
            "target": tgt, "hit": bool(bmap[behav] == tgt),
        }, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"写出 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
