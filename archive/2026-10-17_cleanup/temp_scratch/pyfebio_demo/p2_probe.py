"""P2 第一步：定位退化软骨单元。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"

print("=" * 86)
print("① temp/thums/cart_m1 目录")
print("=" * 86)
d = os.path.join(R, "temp", "thums", "cart_m1")
if os.path.isdir(d):
    for f in sorted(os.listdir(d)):
        print(f"  {f:<44} {os.path.getsize(os.path.join(d,f)):>10} B")
else:
    print("  ✗ 不存在")

print("\n" + "=" * 86)
print("② build_feb.py 的 DET_MIN / det(J) 逻辑")
print("=" * 86)
b = io.open(os.path.join(R, "scripts", "ankle_fe", "build_feb.py"),
            encoding="utf-8").read()
for i, l in enumerate(b.splitlines()):
    if re.search(r"DET_MIN|det.?\(J\)|jacob|退化|degener", l, re.I):
        print(f"  {i+1:>5}: {l.strip()[:140]}")

print("\n" + "=" * 86)
print("③ cart_patch.py 的结构（怎么建软骨）")
print("=" * 86)
c = io.open(os.path.join(R, "scripts", "ankle_fe", "cart_patch.py"),
            encoding="utf-8").read()
print(f"  总 {len(c.encode('utf-8'))} B / {len(c.splitlines())} 行")
for i, l in enumerate(c.splitlines()):
    if re.match(r"^\s*(def |class )", l) or re.search(
            r"det|jacob|min_thick|厚度|退化|_lift|抬回|重心", l, re.I):
        print(f"  {i+1:>5}: {l.strip()[:136]}")

print("\n" + "=" * 86)
print("④ cartilage.py 里的 det / 质量检查")
print("=" * 86)
g = io.open(os.path.join(R, "scripts", "ankle_fe", "cartilage.py"),
            encoding="utf-8").read()
for i, l in enumerate(g.splitlines()):
    if re.search(r"det|jacob|min_thick|厚度|退化|quality|质量", l, re.I):
        print(f"  {i+1:>5}: {l.strip()[:136]}")
