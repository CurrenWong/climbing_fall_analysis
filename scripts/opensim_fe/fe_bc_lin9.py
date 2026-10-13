# -*- coding: utf-8 -*-
"""W3-b 第一步：读冻结 9 骨 FE 记录里**自带**的 linearity 自检（零解算成本）。

`solve_bone_fe` 在 ref_load 与 ref_load/2 各解一次并记 ratio_*（期望 2.0 = 线性弹性）。
> 2.0 越多 ⇒ 该骨的「线性外推」假设越不可靠。

⚠️ 注意：冻结记录**不可复现**（FE_BC_PROBE.md §3，差 6–22%）⇒ 本表只能作**定性指示**，
   最终须用**当前代码重解**复核（W3-b 第二步）。
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BC2 = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2"
RES = ROOT / "results" / "opensim_fe"

frozen = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
r1 = json.loads((RES / "bc_robust_metric_route1.json").read_text(encoding="utf-8"))
dom = {r["bone"]: r for r in r1["rows"]}

print(f"{'bone':16s} {'joint':10s} {'ref_load_N':>11s} {'status':6s} "
      f"{'ratio_max':>9s} {'ratio_p95':>9s} {'ratio_med':>9s} {'ratio_mean':>10s} "
      f"{'main_dom':14s} {'vm_med@ref':>10s}")
rows = []
for p in frozen["parts"]:
    b = p["bone"]
    f = BC2 / b / f"{b}_fe.json"
    if not f.is_file():
        print(f"{b:16s} MISSING {f}")
        continue
    d = json.loads(f.read_text(encoding="utf-8"))
    lin = d.get("linearity") or {}
    post = d.get("post", {}).get("domains", {})
    dname = dom.get(b, {}).get("solid_dom")
    vm_med = None
    for k, v in post.items():
        if k.endswith(dname or "\x00"):
            vm_med = v.get("vm_median_mpa")
    ratio = {k: lin.get(k) for k in
             ("ratio_max", "ratio_p95", "ratio_median", "ratio_mean")}
    print(f"{b:16s} {p['joint']:10s} {float(p['fe_solve_load_n']):11.1f} "
          f"{str(d.get('status')):6s} "
          + " ".join(f"{(ratio[k] if ratio[k] is not None else float('nan')):9.2f}"
                     for k in ("ratio_max", "ratio_p95", "ratio_median", "ratio_mean"))
          + f" {(dname or '-'):14s} "
          + f"{(vm_med if vm_med is not None else float('nan')):10.3f}")
    rows.append({"bone": b, "joint": p["joint"],
                 "ref_load_n": float(p["fe_solve_load_n"]),
                 "status": d.get("status"), "domain": dname,
                 "vm_median_at_ref_mpa": vm_med, **ratio})

print(f"\n共 {len(rows)} 骨；期望 ratio = 2.0（线性弹性）")
bad = [r for r in rows if r.get("ratio_median") and r["ratio_median"] > 2.5]
print(f"ratio_median > 2.5 的骨（线性外推明显失真）：{len(bad)} 个")
for r in sorted(bad, key=lambda x: -x["ratio_median"]):
    print(f"  {r['bone']:16s} ratio_med={r['ratio_median']:6.2f} "
          f"ratio_max={r['ratio_max']:6.2f}")

OUT = RES / "FE_BC_LINEARITY_9BONES.json"
OUT.write_text(json.dumps(
    {"source": "frozen BC2 *_fe.json linearity (NOT reproducible, qualitative only)",
     "expected_ratio": 2.0, "rows": rows}, ensure_ascii=False, indent=2),
    encoding="utf-8")
print(f"\n写出 {OUT}")
