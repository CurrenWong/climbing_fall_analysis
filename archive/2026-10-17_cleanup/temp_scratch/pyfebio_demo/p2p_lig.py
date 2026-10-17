"""P2' 第一步：韧带本构参数 + 与解剖对照 + 逻辑自查。"""
import io
import os
import re

import numpy as np

R = r"D:\Project\climbing_fall_analysis"

print("=" * 92)
print("① lig_read_WK.npz 里的韧带参数（模型自己算的）")
print("=" * 92)
p = os.path.join(R, "temp", "thums", "febio_lig", "lig_read_WK.npz")
d = np.load(p, allow_pickle=True)
names = d["lig_names"]
ea = d["lig_ea"]
npair = d["lig_npair"]
print(f"{'韧带':<36}{'EA (N)':>12}{'配对数':>8}")
for n, e, k in zip(names, ea, npair):
    print(f"{str(n):<36}{float(e):>12.2f}{int(k):>8}")
print(f"\n  EA 合计 = {float(np.sum(ea)):.1f} N")

print("\n" + "=" * 92)
print("② 模型里 shell 韧带的材料/截面参数")
print("=" * 92)
for fn in ("ankle_lig_feb.py", "ligaments.py", "build_feb.py"):
    t = io.open(os.path.join(R, "scripts", "ankle_fe", fn),
                encoding="utf-8").read()
    for i, l in enumerate(t.splitlines()):
        if re.search(r"LIG_E|LIG_MODE|shell|厚度|thickness|E\s*=\s*260|260|"
                     r"a_cross|EA|prescrib|truss", l, re.I):
            s = l.strip()
            if len(s) > 4 and not s.startswith("#!"):
                print(f"  [{fn}:{i+1}] {s[:128]}")

print("\n" + "=" * 92)
print("③ 解剖对照（文献量级）")
print("=" * 92)
LIT = {
    "ATFL": dict(A=(5, 12), L=(15, 20), E=(100, 300), note="最弱、最先断"),
    "CFL": dict(A=(10, 30), L=(20, 30), E=(100, 300), note="内翻时主要抵抗"),
    "PTFL": dict(A=(10, 25), L=(15, 25), E=(100, 300), note="最强"),
    "MD":  dict(A=(10, 40), L=(20, 40), E=(100, 300), note="三角韧带"),
}
for k, v in LIT.items():
    lo, hi = v["A"]
    elo, ehi = v["E"]
    print(f"  {k:<6} A {lo}-{hi} mm²  L {v['L'][0]}-{v['L'][1]} mm  "
          f"E {elo}-{ehi} MPa ⇒ EA ≈ {lo*elo}-{hi*ehi} N   ({v['note']})")

print("\n" + "=" * 92)
print("④ 逻辑自查：两个测量互不相容")
print("=" * 92)
k1 = 592.602 / 0.0142
print(f"  (a) 由「关节只闭 0.0142 mm」反推  k_joint ≈ {k1:,.0f} N/mm")
print(f"  (b) 由软骨材料反推 E·A/t = 12 MPa × 300 mm² / 2.6 mm "
      f"= {12*300/2.6:,.0f} N/mm")
print(f"  (c) 由「总刚度差 2.2×」反推：k_lig + k_joint(tet4) = "
      f"{259.3/0.1472:,.0f}，k_lig + k_joint(penta6) = {116.9/0.1472:,.0f}")
print(f"      ⇒ Δk_joint = {259.3/0.1472 - 116.9/0.1472:,.0f} N/mm")
print()
print("  ⇒ (a) 与 (c) 不可能同时成立：(a) 说关节极刚，则改单元类型不该有 2.2× 差异。")
print("  ⇒ 必须用**无符号的间隙法**重测关节闭合（不靠轴向投影，避免符号约定风险）。")
