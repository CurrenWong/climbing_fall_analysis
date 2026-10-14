# -*- coding: utf-8 -*-
"""E · 口径一致重算 A2 的**范围取证**（零新增解算）。

问题：frozen BC2 是「无侧向 + 线性外推（含降载回退）」。
新 BC 族是「侧向支撑 + 逐高度真解」。要「口径一致」，**所有部位必须在同一族下**。

判据：某部位**必须重算** ⟺ 下列任一成立
  (a) **降载回退**（`load_reduction_factor != 1`）—— 应力含"解在低载、线性放回"的误差
  (b) **非线性显著**（冻结记录自带线性自检 `ratio_median > 2.5`，期望 2.0）
  两项都不成立 ⇒ 冻结数据可沿用（线性 + 无回退）。
"""
from __future__ import annotations

import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import paper_alignment as pa  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
BC2 = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2"

fj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
STAT = "median"
SIG = {p["bone"]: p["sigma_c_mpa"] for p in fj["parts"]}

print(f"{'部位':<6}{'骨':<12}{'lrf':>8}{'ratio_med':>11}{'需重算?':>10}  原因")
need = []
for cn, bone, joint, subs in pa.REGIONS:
    rec = json.loads((BC2 / bone / f"{bone}_fe.json").read_text(encoding="utf-8"))
    lrf = float(rec.get("load_reduction_factor") or 1.0)
    lin = rec.get("linearity") or {}
    rm = lin.get("ratio_median")
    a = abs(lrf - 1.0) > 1e-6
    b = (rm is not None and rm > 2.5)
    why = []
    if a:
        why.append(f"降载回退 lrf={lrf:.3f}")
    if b:
        why.append(f"非线性 ratio_med={rm:.2f}")
    if rm is None and not a:
        why.append("自检缺失")
    print(f"{cn:<6}{bone:<12}{lrf:>8.3f}"
          f"{(f'{rm:.2f}' if rm is not None else 'n/a'):>11}"
          f"{('是' if (a or b) else '否'):>10}  {'; '.join(why) or '线性且无回退 ⇒ 可沿用冻结'}")
    if a or b:
        need.append(bone)

print(f"\n必须重算：{len(need)} 骨 -> {need}")
skip = [p["bone"] for p in fj["parts"] if p["bone"] not in need]
print(f"可沿用冻结：{len(skip)} 骨 -> {skip}")

# ---- 现有新数据能给什么计数 -------------------------------------------------
print("\n=== 现有新数据（本会话已跑）===")
FE_DIR = ROOT / "temp" / "opensim_fe" / "fe_bc_perheight" / "femur_r"
fem = sorted(FE_DIR.glob("h*"))
print(f"  股骨 femur_r : {len(fem)}/50 高度（侧向 [(0.31,0.69)] + 逐高度）")
t6 = json.loads((RES / "T6_PERHEIGHT.json").read_text(encoding="utf-8"))
print(f"  T6           : {t6['n_solved']}/50 高度（侧向 [(0.05,0.95)] + 逐高度）"
      f"  behav={t6['behav']}")

print("\n=== A2 逐部位：现状（冻结） vs 已有新数据 ===")
rows = []
for cn, bone, joint, subs in pa.REGIONS:
    if bone in pa.LONG_BONES:
        h5 = BC2 / bone / f"{bone}.hdf5"
        d = {r["bone"]: r for r in json.loads(
            (RES / "bc_robust_metric_route1.json").read_text(encoding="utf-8"))["rows"]}
        with __import__("h5py").File(h5, "r") as h:
            raw = h["meshes/0/nodes"][:]
            coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
            conn = h["meshes/0/domains"][d[bone]["solid_dom"]][:]
            sid = sorted(h["states"].keys(), key=int)[-1]
            st = h[f"states/{sid}/element_data/stress/{d[bone]['solid_dom']}"][:].astype(float)
        idx = conn[:, 1:] - 1
        p = pa.partition_ends_shaft(coords, coords[idx].mean(axis=1))
        vm = pa.vm_voigt(st)
        e, s_ = float(np.median(vm[p["ends_mask"]])), float(np.median(vm[p["shaft_mask"]]))
        fe = float(next(x for x in fj["parts"] if x["bone"] == bone)["fe_solve_load_n"])
        sc = np.asarray(fj["loads_n"][joint], float) / fe
        se, ss = pa._sc(subs[0][1]), pa._sc(subs[1][1])
        m = (e * sc >= se) | (s_ * sc >= ss)
    else:
        m = np.asarray(fj["sigma_vm_mpa"][STAT][bone], float) >= pa._sc(subs[0][1])
    b = {"all_zero": "常值", "all_one": "常值",
         "complete_separation": "阶跃", "ok": "渐进"}[\
        pa.logistic_fit(np.asarray(fj["heights_m"], float), m.astype(int))["status"]]
    tgt = "阶跃" if cn in pa.PAPER_STEP else "渐进"
    rows.append((cn, bone, b, tgt, b == tgt))
    print(f"  {cn:<6}{bone:<12} 现状={b:<4} 靶={tgt:<4} {'✓' if b==tgt else '✗'}"
          f"   {'【已有新数据】' if bone in ('femur_r','T6') else ''}")

n_hit = sum(1 for r in rows if r[4])
print(f"\n现状逐部位总命中 = {n_hit}/9  "
      f"（阶跃 {sum(1 for r in rows if r[3]=='阶跃' and r[4])}/{len(pa.PAPER_STEP)}"
      f"，渐进 {sum(1 for r in rows if r[3]=='渐进' and r[4])}/{len(pa.PAPER_GRADUAL)}）")
print("⚠️ 本件**不**给「新官方计数」—— 口径未一致（只有 2 个部位有新数据）。")

(RES / "A2_RECOMPUTE_SCOPE.json").write_text(json.dumps({
    "must_rerun": need, "can_keep_frozen": skip,
    "new_data_available": {"femur_r": "50/50 stabilized per-height",
                           "T6": f"{t6['n_solved']}/50 stabilized per-height"},
    "current_per_part": [{"cn": a, "bone": b, "behav": c, "target": d, "hit": e}
                         for a, b, c, d, e in rows],
    "note": "口径一致重算需先把 must_rerun 的 5 骨在新 BC 族下跑完",
}, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"\n写出 {RES / 'A2_RECOMPUTE_SCOPE.json'}")
