"""P2：逐单元算 det(J)，定位退化软骨单元（按域分组 + 找最差）。"""
import re
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"
txt = FEB.read_text(encoding="latin-1", errors="replace")

# 节点
XYZ = {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,'
                     r'\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
print(f"节点 {len(XYZ):,}")

# 域：<Elements type="penta6" ... name="XXX"> <elem id="1">a,b,c,d,e,f</elem> ...
DOMS = {}
for m in re.finditer(r'<Elements\s+type="(hex8|penta6|tet4|tri3|quad4)"'
                     r'[^>]*?name="([^"]*)"[^>]*>(.*?)</Elements>', txt, re.S):
    kind, name, body = m.group(1), m.group(2), m.group(3)
    conn = [[int(x) for x in c.split(",")] for c in
            re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*</elem>', body)]
    conn = [c for c in conn if len(c) >= 4]
    DOMS.setdefault(kind, []).append((name, conn))
for k, v in DOMS.items():
    print(f"  {k:<8} {len(v)} 组，共 {sum(len(c) for _, c in v):,} 单元")
    for nm, c in v[:30]:
        print(f"      {nm:<42} {len(c)}")


def detJ_wedge(P):
    """P:(6,3) 一点/多点 det(J)（penta6，取中心 + 2 点 Gauss）。"""
    r, s, t = 1/3, 1/3, 1/2
    dN = np.array([
        [-(1-t), -(1-t), -(1-r-s)],
        [(1-t), 0, -r],
        [0, (1-t), -s],
        [-t, -t, (1-r-s)],
        [t, 0, r],
        [0, t, s],
    ])
    J = dN.T @ P
    return np.linalg.det(J)


def detJ_hex(P):
    g = 1/np.sqrt(3)
    ds = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                   [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    out = []
    for xi, et, ze in ds * g:
        dN = np.zeros((8, 3))
        sg = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
        for i in range(8):
            a, b, c = sg[i]
            dN[i] = [a*(1+b*et)*(1+c*ze)/8, b*(1+a*xi)*(1+c*ze)/8,
                     c*(1+a*xi)*(1+b*et)/8]
        out.append(np.linalg.det(dN.T @ P))
    return min(out)


print("\n" + "=" * 96)
print("逐域 det(J) 统计（penta6 取中心；hex8 取 8 Gauss 最小）")
print("=" * 96)
print(f"{'域(id)':<12}{'单元':>6}{'min':>13}{'p1':>13}{'中位':>13}"
      f"{'<1e-2':>7}{'<5e-2':>7}{'<0.1':>7}")
rows = []
for kind, groups in DOMS.items():
    if kind not in ("penta6", "hex8"):
        continue
    f = detJ_wedge if kind == "penta6" else detJ_hex
    for did, conn in groups:
        if len(conn) < 3:
            continue
        ds = []
        for c in conn:
            try:
                P = np.array([XYZ[n] for n in c], float)
            except KeyError:
                continue
            ds.append(f(P))
        ds = np.array(ds)
        if not len(ds):
            continue
        rows.append((kind, did, len(ds), ds))
        print(f"{kind+':'+str(did):<12}{len(ds):>6}{ds.min():>13.4e}"
              f"{np.percentile(ds,1):>13.4e}{np.median(ds):>13.4e}"
              f"{(ds<1e-2).sum():>7}{(ds<5e-2).sum():>7}{(ds<0.1).sum():>7}")

print("\n" + "=" * 96)
print("最差的 15 个单元（全局）")
print("=" * 96)
allrows = []
for kind, did, n, ds in rows:
    for i, d in enumerate(ds):
        allrows.append((d, kind, did, i))
allrows.sort()
for d, kind, did, i in allrows[:15]:
    print(f"  {kind}:{did:<6} 单元#{i:<5} det(J)={d:>12.5e}")

print("\n" + "=" * 96)
print("👉 结论要点")
print("=" * 96)
worst = allrows[0]
print(f"  全局最小 det(J) = {worst[0]:.5e}（{worst[1]}:{worst[2]}）")
neg = [r for r in allrows if r[0] <= 0]
print(f"  det(J) <= 0 的单元数 = {len(neg)}")
small = [r for r in allrows if 0 < r[0] < 1e-2]
print(f"  0 < det(J) < 1e-2 的单元数 = {len(small)}")
