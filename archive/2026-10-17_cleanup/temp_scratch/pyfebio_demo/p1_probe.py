"""P1 现状排查：骨应力读数链路。
① 踝模型写出的输出里有没有应力 ② 宏观线怎么读应力 ③ 骨强度判据在哪
"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
A = os.path.join(R, "scripts", "ankle_fe", "ankle_lig_feb.py")
BF = os.path.join(R, "scripts", "opensim_fe", "bone_feb.py")
BONE = os.path.join(R, "src", "climbing", "bone.py")

t = io.open(A, encoding="utf-8").read()
print("=" * 88)
print("① 踝模型的 <Output> 段（写什么）")
print("=" * 88)
m = re.search(r"<Output>.*?</Output>", t, re.S)
print(m.group(0)[:900] if m else "✗ 没有 <Output> 段")
print("\n--- 提到 logfile/xplt/plotfile/stress 的行 ---")
for i, l in enumerate(t.splitlines()):
    if re.search(r"logfile|xplt|plotfile|LOG_DATA|stress|element_data", l, re.I):
        print(f"  {i+1:>5}: {l.strip()[:140]}")

print("\n" + "=" * 88)
print("② 宏观线 bond_feb.py 怎么读应力")
print("=" * 88)
if os.path.exists(BF):
    b = io.open(BF, encoding="utf-8").read()
    for i, l in enumerate(b.splitlines()):
        if re.search(r"von|mises|stress|h5|xplt|logfile|npz", l, re.I):
            print(f"  {i+1:>5}: {l.strip()[:140]}")
else:
    print("✗ 不存在")

print("\n" + "=" * 88)
print("③ 骨强度判据 src/climbing/bone.py")
print("=" * 88)
if os.path.exists(BONE):
    bo = io.open(BONE, encoding="utf-8").read()
    for i, l in enumerate(bo.splitlines()):
        if re.search(r"^(def |class )|strength|强度|fail|骨折|MPa|yield|ultimate",
                     l, re.I):
            print(f"  {i+1:>5}: {l.strip()[:140]}")
else:
    print("✗ 不存在")

print("\n" + "=" * 88)
print("④ 踝模型刚跑完的 npz 里有哪些键")
print("=" * 88)
import numpy as np
for p in (os.path.join(R, "temp", "thums", "febio_lig", "lig_read_WK.npz"),):
    if os.path.exists(p):
        d = np.load(p, allow_pickle=True)
        print(f"  {os.path.basename(p)}: {len(d.files)} 键 / "
              f"{os.path.getsize(p)} B")
        for k in sorted(d.files):
            v = d[k]
            print(f"    {k:<24} {getattr(v, 'shape', '-')} {getattr(v, 'dtype', '')}")
    else:
        print(f"  ✗ 缺 {p}")
print("\n--- temp/thums/febio_lig 目录下的产物 ---")
d = os.path.join(R, "temp", "thums", "febio_lig")
if os.path.isdir(d):
    for f in sorted(os.listdir(d)):
        print(f"  {f:<40} {os.path.getsize(os.path.join(d, f)):>10} B")
