"""TIE 诊断 ①：从生成的 .feb 取地面真值 —— tie contact 块 + 两个 surface 的实际节点。"""
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
FEB = W / "ankle_tibiotalar.feb"
txt = FEB.read_text(encoding="latin-1")

print("=" * 72)
print(f"① tie contact 块（原样）  {FEB.name} {FEB.stat().st_size/1024:.0f} KB")
print("=" * 72)
for m in re.finditer(
        r'<contact\b[^>]*type="tied-elastic".*?</contact>', txt, re.S):
    print("  " + re.sub(r"\n\s*", "\n  ", m.group(0).strip()))

print("\n" + "=" * 72)
print("② 所有 SurfacePair / Surface 的节点数")
print("=" * 72)
pair2surf = {}
for m in re.finditer(r'<SurfacePair\b[^>]*name="([^"]+)"[^>]*>(.*?)</SurfacePair>', txt, re.S):
    nm, body = m.group(1), m.group(2)
    a = re.search(r'<primary>([^<]+)</primary>', body)
    b = re.search(r'<secondary>([^<]+)</secondary>', body)
    pair2surf[nm] = (a.group(1) if a else "?", b.group(1) if b else "?")

surf_nodes, surf_faces = {}, {}
node_ids = np.array([int(x) for x in re.findall(r'<node id="(\d+)"', txt)])
for m in re.finditer(r'<Surface\b[^>]*name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    nm, body = m.group(1), m.group(2)
    surf_nodes[nm] = np.array([int(x) for x in re.findall(r"<(?:tri3|quad4|quad8)\b[^>]*>\s*([\d,]+)\s*</", body)])
    surf_faces[nm] = [tuple(int(y) for y in f.split(","))
                      for f in re.findall(r"<(?:tri3|quad4|quad8)\b[^>]*>\s*([\d,]+)\s*</", body)]

for nm in sorted(surf_nodes):
    print(f"  {nm:28s} 节点 {len(surf_nodes[nm]):>6d} 面 {len(surf_faces[nm]):>6d}")
print("\n  SurfacePair:")
for nm, (a, b) in sorted(pair2surf.items()):
    na = len(surf_nodes.get(a, []))
    nb = len(surf_nodes.get(b, []))
    print(f"  {nm:28s} primary={a}({na})  secondary={b}({nb})")

# ── ③ 节点坐标 ──
XYZ = np.zeros((int(node_ids.max()) + 1, 3))
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)',
                     txt):
    i = int(m.group(1))
    XYZ[i] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
print(f"\n  节点总数 {len(node_ids)}  坐标解析 {int((np.abs(XYZ).sum(1) > 0).sum())}")


def pt2tri(P, A, B, C):
    """点到三角形的最短距离（向量化）。"""
    AB, AC, AP = B - A, C - A, P - A
    d1, d2 = AB @ AP.T, AC @ AP.T
    BP = P - B
    d3, d4 = AB @ BP.T, AC @ BP.T
    CP = P - C
    d5, d6 = AB @ CP.T, AC @ CP.T
    VC = np.cross(AB, AC)
    l2 = VC @ VC.T
    W = A.T + ((d1[:, None] * AC - d2[:, None] * AB) / l2).T
    U = B.T + ((d3[:, None] * AC - d4[:, None] * AB) / l2).T
    V = C.T + ((d5[:, None] * AC - d6[:, None] * AB) / l2).T
    Pn = np.repeat(P, A.shape[0], 0)
    cand = np.stack([P[:, None, :] - W.T, P[:, None, :] - U.T, P[:, None, :] - V.T], 0)
    dd = np.linalg.norm(cand, axis=3)
    dd[0][(d1 < 0) | (d1 > 1)] = np.inf
    dd[1][(d3 > 0) | (d4 > 1)] = np.inf
    dd[2][(d5 < 0) | (d6 > 1)] = np.inf
    n = np.cross(VC, AP)
    d = (n @ VC.T) / (l2 ** 1.5)
    ok = (d1 >= 0) & (d2 >= 0) & (d1 + d2 <= 1)
    dd.append(np.where(ok, np.abs(d), np.inf))
    return dd.min(0)


print("\n" + "=" * 72)
print("③ primary 节点 → secondary 真实面 的最短距离（核心判据）")
print("=" * 72)
for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    pr, se = f"tiesurf_{j}", f"bonesurf_{j}"
    if pr not in surf_nodes or se not in surf_nodes:
        print(f"  {j}: 缺 surface，跳过")
        continue
    P = XYZ[surf_nodes[pr]]
    Fc = np.array([[XYZ[k] for k in f] for f in surf_faces[se]])
    A, B, C = Fc[:, 0], Fc[:, 1], Fc[:, 2]
    # 分块避免内存爆
    best = np.full(len(P), np.inf)
    for s in range(0, len(P), 256):
        e = min(s + 256, len(P))
        best[s:e] = pt2tri(P[s:e], A, B, C).min(1)
    print(f"\n  【{j}】primary 节点 {len(P)}  secondary 面 {len(A)}")
    for t in (1e-9, 1e-6, 1e-4, 0.01, 0.05, 0.1, 0.2, 0.5):
        n = int((best > t).sum())
        print(f"    距 secondary 面 > {t:<8g}: {n:>5d} 个 ({n/len(P)*100:5.1f}%)")
    fin = best[np.isfinite(best)]
    if len(fin):
        print(f"    有限距离: min {fin.min():.2e}  中位 {np.median(fin):.2e}  "
              f"p95 {np.percentile(fin,95):.2e}  max {fin.max():.2e}")

print("\n" + "=" * 72)
print("④ 软骨底面节点 ↔ 骨面节点的实际间隙（为什么 TIE_SEC_TOL=0.05 会漏）")
print("=" * 72)
from scipy.spatial import cKDTree
for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    pr, se = f"tiesurf_{j}", f"bonesurf_{j}"
    if pr not in surf_nodes or se not in surf_nodes:
        continue
    P = XYZ[surf_nodes[pr]]
    Q = XYZ[surf_nodes[se]]
    d = cKDTree(Q).query(P)[0]
    print(f"  【{j}】primary↔secondary 节点距: min {d.min():.2e} 中位 {np.median(d):.2e} "
          f"p95 {np.percentile(d,95):.2e} max {d.max():.2e}  >0.05 的 {int((d>0.05).sum())}/{len(P)}")

# ── ⑤ 整块骨面（未筛）有多少面 / 命中多少 ──
print("\n" + "=" * 72)
print("⑤ 关节区 bbox")
print("=" * 72)
for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    pr = f"tiesurf_{j}"
    if pr in surf_nodes:
        P = XYZ[surf_nodes[pr]]
        print(f"  {j:20s} x[{P[:,0].min():7.2f},{P[:,0].max():7.2f}] "
              f"y[{P[:,1].min():7.2f},{P[:,1].max():7.2f}] "
              f"z[{P[:,2].min():7.2f},{P[:,2].max():7.2f}]")
