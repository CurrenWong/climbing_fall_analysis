# -*- coding: utf-8 -*-
"""A1 边界锚定诊断（零 FE）—— #1(5.0) / #2(12.0) 由哪些部位的骨折行决定？

做法：逐部位把该行替换为「常值0（从不折）」或「常值1（总是折）」，
看 A1 的两个边界怎么动 ⇒ 识别锚定部位。
再对锚定部位扫「阶跃@h*」，看 #1/#2 能否被合法地推到目标 (7, 16)。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import paper_alignment as pa            # noqa: E402
import s4_alignment_bc2 as s4           # noqa: E402

yi = s4.route_yi("median")
names = pa.REGION_NAMES
base_rows = {cn: np.asarray(yi["matrix_by_part"][cn], dtype=int) for cn in names}

print("=" * 92)
print("① 现役相段与各段骨折部位组成")
print("=" * 92)
print("  边界 = %s   Δ = %s   -> %s"
      % (s4.a1_verdict(yi["phases"])["ours_boundaries_m"],
         s4.a1_verdict(yi["phases"])["abs_delta"],
         s4.a1_verdict(yi["phases"])["verdict"]))
for p in yi["phases"]:
    print("   相 %s  骨折部位=%s" % (p["range_m"], p.get("fractured_regions")))
print("\n  逐部位首折：")
for cn in names:
    r = base_rows[cn]
    fh = next((i + 1 for i, v in enumerate(r) if v), None)
    print("    %-8s first=%-6s n_frac=%d" % (cn, fh if fh else ">50", int(r.sum())))


def a1_of(over):
    r2 = {cn: over.get(cn, base_rows[cn]) for cn in names}
    M = np.array([[int(r2[cn][i]) for cn in names] for i in range(50)], dtype=int)
    clus = pa.jaccard_clusters(M, 5)
    ph = s4.three_phases_seq(clus["labels"], M)
    v = s4.a1_verdict(ph)
    return v, ph


print()
print("=" * 92)
print("② 锚定诊断：逐部位置「常值0（从不折）」/「常值1（总是折）」")
print("=" * 92)
z = np.zeros(50, dtype=int)
o = np.ones(50, dtype=int)
print("  %-10s %-24s %-24s" % ("部位", "置 常值0 → 边界/Δ", "置 常值1 → 边界/Δ"))
for cn in names:
    v0, _ = a1_of({cn: z})
    v1, _ = a1_of({cn: o})
    m0 = "" if v0["ours_boundaries_m"] == [5.0, 12.0] else "  ←★"
    m1 = "" if v1["ours_boundaries_m"] == [5.0, 12.0] else "  ←★"
    print("  %-10s %-24s %-24s" % (cn,
          "%s %s%s" % (v0["ours_boundaries_m"], v0["abs_delta"], m0),
          "%s %s%s" % (v1["ours_boundaries_m"], v1["abs_delta"], m1)))


def step_at(h):
    return np.array([1 if (i + 1) >= h else 0 for i in range(50)], dtype=int)


print()
print("=" * 92)
print("③ 对「敏感部位」扫阶跃 h*（找能否把 #1→7 / #2→16）")
print("=" * 92)
cands = []
for cn in names:
    v0, _ = a1_of({cn: z})
    v1, _ = a1_of({cn: o})
    if v0["ours_boundaries_m"] != [5.0, 12.0] or v1["ours_boundaries_m"] != [5.0, 12.0]:
        cands.append(cn)
print("  敏感部位（置常值0/1 会改变边界）：%s" % cands)

for cn in cands:
    print("\n  【%s】h* → 边界/Δ" % cn)
    for h in (4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 25, 33):
        v, ph = a1_of({cn: step_at(h)})
        ok = ""
        if v["abs_delta"][0] < 2.0 and v["abs_delta"][1] < 2.0:
            ok = "  ← ✅ 双边界过"
        elif v["abs_delta"][1] < 2.0:
            ok = "  ← #2 过"
        elif v["abs_delta"][0] < 2.0:
            ok = "  ← #1 过"
        print("     h*=%2d  边界=%-14s Δ=%-14s%s"
              % (h, v["ours_boundaries_m"], v["abs_delta"], ok))
