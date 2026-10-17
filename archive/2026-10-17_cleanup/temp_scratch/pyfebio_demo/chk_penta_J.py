"""多点采样的 penta6 雅可比检查：找出 FEBio 认定为负、但质心采样查不出的单元。

penta6 自然坐标 (r,s,t): r,s>=0, r+s<=1, t∈[-1,1]
  节点序 0..2 = 底面三角形(r,s) 顶点, 3..5 = 顶面
  N_i = (1-t)/2 * L_i(r,s)  (i=0,1,2)   N_{i+3} = (1+t)/2 * L_i(r,s)
  L_0=1-r-s, L_1=r, L_2=s
质心 (1/3,1/3,0) 处若单元"扭转"，正负贡献会抵消 ⇒ 必须多点采样。
"""
from pathlib import Path

import numpy as np

P = Path(r"D:\Project\climbing_fall_analysis\temp\thums\cart_new\cartilage_mesh.npz")
d = np.load(P)

# 采样点：6 个顶点 + 8 个内部点
PTS = [(1.0, 0.0, -1.0), (0.0, 1.0, -1.0), (0.0, 0.0, -1.0),
       (1.0, 0.0, 1.0), (0.0, 1.0, 1.0), (0.0, 0.0, 1.0),
       (1 / 3, 1 / 3, 0.0), (0.25, 0.25, -0.5), (0.25, 0.25, 0.5),
       (0.5, 0.25, -0.5), (0.25, 0.5, -0.5), (0.5, 0.25, 0.5), (0.25, 0.5, 0.5),
       (0.4, 0.4, 0.0)]


def dN(r, s, t):
    """∂N/∂(r,s,t)，形状 (6,3)。"""
    L0, L1, L2 = 1 - r - s, r, s
    dL = [(-1.0, -1.0), (1.0, 0.0), (0.0, 1.0)]
    L = [L0, L1, L2]
    G = np.zeros((6, 3))
    for i in range(3):
        G[i] = np.array([(1 - t) / 2 * dL[i][0], (1 - t) / 2 * dL[i][1], -L[i] / 2])
        G[i + 3] = np.array([(1 + t) / 2 * dL[i][0], (1 + t) / 2 * dL[i][1], L[i] / 2])
    return G


def bad_prisms(N, E, swap=None):
    out = []
    for k, row in enumerate(E):
        c = [int(x) for x in row if int(x) != 0]
        if len(c) != 6:
            continue
        idx = [(-x - 1) if x < 0 else (x - 1) for x in c]
        if swap:
            idx = list(idx)
            idx[4], idx[5] = idx[5], idx[4]
        X = N[idx]
        worst = np.inf
        for (r, s, t) in PTS:
            J = X.T @ dN(r, s, t)
            worst = min(worst, float(np.linalg.det(J)))
        if worst <= 0:
            out.append((k, worst))
    return out


for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    N = np.asarray(d[f"{j}_nodes"], float)
    E = np.asarray(d[f"{j}_elems"])
    b0 = bad_prisms(N, E)
    b1 = bad_prisms(N, E, swap=True)
    print(f"\n===== {j} ({len(E)} prism) =====")
    print(f"  原样:     负雅可比 {len(b0):>4} 个  {[x[0] for x in b0][:15]}")
    print(f"  对调 4/5: 负雅可比 {len(b1):>4} 个  {[x[0] for x in b1][:15]}")
    if b0:
        print("  原样最差 5 个 (prism#, min det(J)):")
        for k, w in sorted(b0, key=lambda z: z[1])[:5]:
            print(f"    #{k:<5} min detJ = {w:.4e}")
