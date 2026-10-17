"""P2b：厚度分布 —— 验证"超薄单元"新假设。"""
import re
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"
txt = FEB.read_text(encoding="latin-1", errors="replace")

XYZ = {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,'
                     r'\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])

DOMS = {}
for m in re.finditer(r'<Elements\s+type="(hex8|penta6)"[^>]*?name="([^"]*)"'
                     r'[^>]*>(.*?)</Elements>', txt, re.S):
    kind, nm, body = m.group(1), m.group(2), m.group(3)
    conn = [[int(x) for x in c.split(",")] for c in
            re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*</elem>', body)]
    DOMS[nm] = (kind, [c for c in conn if len(c) >= 4])

print("=" * 100)
print("软骨单元厚度（底面三角形质心 → 顶面三角形质心；wedge 用 0-1-2/3-4-5，hex 用 0-3-4-7/1-2-6-5）")
print("=" * 100)
allt = {}
for nm, (kind, conn) in DOMS.items():
    ts = []
    for c in conn:
        try:
            P = np.array([XYZ[n] for n in c], float)
        except KeyError:
            continue
        if kind == "penta6":
            bot = P[[0, 1, 2]].mean(0)
            top = P[[3, 4, 5]].mean(0)
        else:
            bot = P[[0, 3, 4, 7]].mean(0)
            top = P[[1, 2, 6, 5]].mean(0)
        ts.append(np.linalg.norm(top - bot))
    ts = np.array(ts)
    allt[nm] = ts
    if not len(ts):
        continue
    print(f"{nm:<34}{kind:<8}{len(ts):>5}  "
          f"min {ts.min():8.4f}  p1 {np.percentile(ts,1):8.4f}  "
          f"p5 {np.percentile(ts,5):8.4f}  中位 {np.median(ts):8.4f}  "
          f"max {ts.max():8.4f}")
    for th in (0.01, 0.05, 0.1, 0.3, 0.5, 1.0):
        n = (ts < th).sum()
        if n:
            print(f"      < {th:<5} mm : {n:>4} 个 ({n/len(ts)*100:5.1f}%)")

print("\n" + "=" * 100)
print("软骨（两域合计）厚度分位")
print("=" * 100)
cart = np.concatenate([allt[n] for n in allt if "tibiotalar" in n])
cart.sort()
for q in (0, 0.5, 1, 2, 5, 10, 25, 50, 75, 90, 100):
    print(f"  p{q:<4} = {np.percentile(cart, q):8.4f} mm")
print(f"\n  < 0.10 mm 的单元：{(cart<0.10).sum():>4} / {len(cart)} "
      f"({(cart<0.10).sum()/len(cart)*100:.1f}%)")
print(f"  < 0.05 mm 的单元：{(cart<0.05).sum():>4} / {len(cart)} "
      f"({(cart<0.05).sum()/len(cart)*100:.1f}%)")
print(f"  < 0.02 mm 的单元：{(cart<0.02).sum():>4} / {len(cart)}")

print("\n" + "=" * 100)
print("👉 判定")
print("=" * 100)
print(f"  最薄单元 {cart[0]:.5f} mm ⇒ 局部压缩超过 {cart[0]:.5f} mm 即翻转")
print(f"  当前压缩上限 +1.6~1.7 mm ⇒ 若最薄单元只有 {cart[0]:.4f} mm ⇒ 必然翻")
print("  ⇒ 假设成立与否：见上面'<0.1mm'的占比。")
