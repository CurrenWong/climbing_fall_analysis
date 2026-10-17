"""P2'a 复核：用 Gauss 点（而非中心）算 penta6 的 det(J)，确认近奇异单元。"""
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
for m in re.finditer(r'<Elements\s+type="(penta6|hex8)"[^>]*?name="([^"]*)"'
                     r'[^>]*>(.*?)</Elements>', txt, re.S):
    nm, kind, body = m.group(2), m.group(1), m.group(3)
    conn = [[int(x) for x in c.split(",")] for c in
            re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*</elem>', body)]
    DOMS[nm] = (kind, [c for c in conn if len(c) >= 4])


def dN_wedge(r, s, t):
    return np.array([
        [-(1 - t), -(1 - t), -(1 - r - s)],
        [(1 - t), 0, -r],
        [0, (1 - t), -s],
        [-t, -t, (1 - r - s)],
        [t, 0, r],
        [0, t, s],
    ])


def dN_hex(xi, et, ze):
    sg = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                   [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    dN = np.zeros((8, 3))
    for i in range(8):
        a, b, c = sg[i]
        dN[i] = [a * (1 + b * et) * (1 + c * ze) / 8,
                 b * (1 + a * xi) * (1 + c * ze) / 8,
                 c * (1 + a * xi) * (1 + b * et) / 8]
    return dN


# wedge 的 3 点 Gauss（FEBio 常用），另加角点
WEDGE_PTS = [(1/6, 1/6, 1/2), (2/3, 1/6, 1/2), (1/6, 2/3, 1/2),
             (1/3, 1/3, 0.0), (1/3, 1/3, 1.0),
             (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
             (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0)]
HEX_PTS = [tuple(g) for g in
           np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                     [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]) *
           (1 / np.sqrt(3))]

print("=" * 100)
print("penta6/hex8 det(J)：中心 vs Gauss/角点最小值")
print("=" * 100)
print(f"{'域':<30}{'单元':>6}{'中心 min':>13}{'Gauss min':>13}"
      f"{'<1e-2':>7}{'<0.1':>7}{'<0.5':>7}")
for nm, (kind, conn) in sorted(DOMS.items()):
    ctr, gss = [], []
    for c in conn:
        try:
            P = np.array([XYZ[n] for n in c], float)
        except KeyError:
            continue
        if kind == "penta6":
            ctr.append(np.linalg.det(dN_wedge(1/3, 1/3, 1/2).T @ P))
            gss.append(min(np.linalg.det(dN_wedge(*p).T @ P)
                           for p in WEDGE_PTS))
        else:
            ctr.append(min(np.linalg.det(dN_hex(*p).T @ P) for p in HEX_PTS))
            gss.append(ctr[-1])
    ctr, gss = np.array(ctr), np.array(gss)
    if not len(ctr):
        continue
    print(f"{nm:<30}{len(ctr):>6}{ctr.min():>13.4e}{gss.min():>13.4e}"
          f"{(gss<1e-2).sum():>7}{(gss<0.1).sum():>7}{(gss<0.5).sum():>7}")

print("\n" + "=" * 100)
print("builder 自报（同一次运行）")
print("=" * 100)
print("  tibiotalar_tibia  penta6 保留 145 剔除  9  det(J) min=1.9085e-01")
print("  tibiotalar_talus  penta6 保留 182 剔除 22  det(J) min=2.3914e-03")
print("  ⇒ 若本表 Gauss min 与之同量级 ⇒ 近奇异单元确实存在于建好的网格里")
