# -*- coding: utf-8 -*-
"""FE_BC · 任务A 前置取证：**换聚合口径会怎样影响全部 9 部位**（零新增解算）。

红线（`results/opensim_fe/FE_BC_STAB.md` §4-ter.3）：
  选哪个口径必须由独立论证决定，**且须同时说明其他 8 部位会怎样**。
本脚本把 A2 的逐部位判定在 {max, p95, mean, median} 四种口径下**并列算出**。

口径来源
--------
* 非长骨：`fracture_matrix_bc2.json["sigma_vm_mpa"][stat][bone]`（50 长度，已冻结四种）
* 长骨  ：现 FE 只用 `parts[*]["fe_solve_load_n"]` + **冻结 hdf5** 重算 ends/shaft
          （`s4_alignment_bc2.py:159-192` 的同一路径；`sel` 按口径取 median / p95 / max / mean）
判据/强度：`pa.REGIONS` + `pa._sc` + `pa.LONG_BONES`（**与 S4 同源，不改**）
行为分类：`pa._behav(pa.logistic_fit(...))`（**参考实现原封**）
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
BC2_DIR = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2"
OUT = RES / "FE_BC_AGG_SENS_9BONES.json"

STATS = ("median", "mean", "p95", "max")
SEL = {"median": np.median, "mean": np.mean,
       "p95": lambda a: np.percentile(a, 95), "max": np.max}


def _behav(rec) -> str:
    """`paper_alignment.py:465-473` 的嵌套 `_behav` 逐字复刻。"""
    st = rec["status"]
    if st in ("all_zero", "all_one"):
        return "常值"
    if st == "complete_separation":
        return "阶跃"
    if st == "ok":
        return "渐进"
    return st


def long_bone_stats(stat: str) -> dict:
    """按口径重算三个长骨的 ends/shaft（读冻结 hdf5，纯只读）。"""
    import h5py

    fmj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    r1 = json.loads((RES / "bc_robust_metric_route1.json").read_text(encoding="utf-8"))
    dom = {r["bone"]: r for r in r1["rows"]}
    heights = np.asarray(fmj["heights_m"], float)
    loads = fmj["loads_n"]
    out = {}
    for _cn, bone, joint, subs in pa.REGIONS:
        if bone not in pa.LONG_BONES:
            continue
        h5p = BC2_DIR / bone / f"{bone}.hdf5"
        if not h5p.is_file():
            out[bone] = None
            continue
        d = dom[bone]["solid_dom"]
        with h5py.File(h5p, "r") as h:
            raw = h["meshes/0/nodes"][:]
            coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
            conn = h["meshes/0/domains"][d][:]
            sid = sorted(h["states"].keys(), key=int)[-1]
            stress = h[f"states/{sid}/element_data/stress/{d}"][:].astype(float)
        idx = conn[:, 1:] - 1
        cent = coords[idx].mean(axis=1)
        p = pa.partition_ends_shaft(coords, cent)
        vm = pa.vm_voigt(stress)
        e = float(SEL[stat](vm[p["ends_mask"]]))
        s = float(SEL[stat](vm[p["shaft_mask"]]))
        fe = float(next(x for x in fmj["parts"] if x["bone"] == bone)["fe_solve_load_n"])
        scale = np.asarray(loads[joint], float) / fe
        sig_e, sig_s = pa._sc(subs[0][1]), pa._sc(subs[1][1])
        m = (e * scale >= sig_e) | (s * scale >= sig_s)
        out[bone] = {"ends_stat": e, "shaft_stat": s, "fe_solve_n": fe,
                     "n_ends": int(p["n_ends"]), "mask": m}
    return out, heights


def main() -> int:
    fmj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    heights = np.asarray(fmj["heights_m"], float)

    result = {}
    for stat in STATS:
        lb, _ = long_bone_stats(stat)
        rows = {}
        for cn, bone, joint, subs in pa.REGIONS:
            if bone in pa.LONG_BONES:
                if lb.get(bone) is None:
                    continue
                m = lb[bone]["mask"]
            else:
                vm_h = np.asarray(fmj["sigma_vm_mpa"][stat][bone], float)
                m = vm_h >= pa._sc(subs[0][1])
            fit = pa.logistic_fit(heights, m.astype(int))
            b = _behav(fit)
            n = int(np.asarray(m).sum())
            target = "阶跃" if cn in pa.PAPER_STEP else "渐进"
            ok = (b == target)
            rows[cn] = {"bone": bone, "behav": b, "status": fit["status"], "n_fracture": n,
                        "target": target, "hit": bool(ok)}
        result[stat] = rows

    print("口径敏感性 · A2 逐部位（`PAPER_STEP` = 6 个体部位；strict/逐部位口径）")
    hdr = f"{'部位':<6}{'骨':<14}" + "".join(f"{s:>10}" for s in STATS)
    print(hdr)
    for cn, bone, _j, _s in pa.REGIONS:
        line = f"{cn:<6}{bone:<14}"
        for s in STATS:
            r = result[s].get(cn)
            line += f"{(r['behav'] if r else '-'):>10}"
        print(line)

    print("\n逐口径：阶跃 6 靶命中数（部位口径）")
    summary = {}
    for s in STATS:
        rows = result[s]
        n_step = sum(1 for cn in pa.PAPER_STEP if rows.get(cn, {}).get("hit"))
        n_grad = sum(1 for cn in pa.PAPER_GRADUAL if rows.get(cn, {}).get("hit"))
        summary[s] = {"step_hits": n_step, "step_total": len(pa.PAPER_STEP),
                      "gradual_hits": n_grad, "gradual_total": len(pa.PAPER_GRADUAL),
                      "note": "6 靶逐部位"}
        print(f"  {s:>7}: 阶跃 {n_step}/{len(pa.PAPER_STEP)}   "
              f"渐进 {n_grad}/{len(pa.PAPER_GRADUAL)}")
        miss = [cn for cn in pa.PAPER_STEP if not rows.get(cn, {}).get("hit")]
        print(f"           阶跃漏判: {miss}")

    OUT.write_text(json.dumps({
        "purpose": "A2 口径敏感性（零新增解算；读冻结 json + 冻结 hdf5）",
        "stats": list(STATS), "per_stat": result, "summary": summary,
        "target": {"PAPER_STEP": sorted(pa.PAPER_STEP),
                   "PAPER_GRADUAL": sorted(pa.PAPER_GRADUAL)},
        "warning": "本件只测敏感性，不裁定口径；裁定须独立论证",
    }, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(f"\n写出 {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
