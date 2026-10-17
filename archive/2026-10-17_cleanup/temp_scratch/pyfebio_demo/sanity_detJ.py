"""★ 前置体检：my detJ_hex8 对吗？基线拆分后的**每个 tet** 符号如何？"""
import numpy as np
from pathlib import Path

G = 1 / np.sqrt(3)
GP = [(-G, -G, -G), (G, -G, -G), (G, G, -G), (-G, G, -G),
      (-G, -G, G), (G, -G, G), (G, G, G), (-G, G, G)]
REF = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
       (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]


def detJ_min(P):
    w = 1e30
    for s in GP:
        dN = np.zeros((8, 3))
        for i in range(8):
            a, b, c = REF[i]
            dN[i, 0] = a * (1 + b * s[1]) * (1 + c * s[2]) / 8
            dN[i, 1] = b * (1 + a * s[0]) * (1 + c * s[2]) / 8
            dN[i, 2] = c * (1 + a * s[0]) * (1 + b * s[1]) / 8
        w = min(w, float(np.linalg.det(P.T @ dN)))
    return w


print("=" * 74)
print("A. detJ 函数体检（已知答案）")
print("=" * 74)
cube = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                 [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
print(f"  单位立方体        detJ = {detJ_min(cube):+.6f}   期望 +0.125000  "
      f"{'✅' if abs(detJ_min(cube)-0.125) < 1e-9 else '❌'}")
slab = cube.copy()
slab[4:, 2] = 0.7
print(f"  薄板 1×1×0.7      detJ = {detJ_min(slab):+.6f}   期望 +0.087500  "
      f"{'✅' if abs(detJ_min(slab)-0.0875) < 1e-9 else '❌'}")
inv = cube[[0, 1, 2, 3, 4, 5, 6, 7]][[3, 2, 1, 0, 7, 6, 5, 4]]
print(f"  底面倒序(应反号)  detJ = {detJ_min(inv):+.6f}   期望 −0.125000  "
      f"{'✅' if detJ_min(inv) < 0 else '❌'}")
# 故意扭曲：把顶面平移使它折叠
tw = cube.copy()
tw[[4, 5], 0] += 1.2
print(f"  顶面平移致折叠    detJ = {detJ_min(tw):+.6f}   "
      f"{'(应为负)' if detJ_min(tw) < 0 else '(仍正)'}")

HEX_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
            (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
print("\n  同一立方体的 6-tet 分解，每个 tet 的 det：")
for t in HEX_TETS:
    P = cube[list(t)]
    print(f"    {t}  det={np.linalg.det(np.stack([P[1]-P[0],P[2]-P[0],P[3]-P[0]])):+.6f}")

print()
print("=" * 74)
print("B. 基线软骨：**拆分后用到的每个 tet** 的符号（这才是 build_feb 看到的）")
print("=" * 74)
THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}
WEDGE_TETS = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]

print(f"{'关节':22s} {'格':>5s} {'拆出tet':>8s} {'负tet':>7s} {'零tet':>6s} "
      f"{'6tet和<0':>9s} {'Gauss<0':>8s}")
for j in ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]:
    E = np.asarray(src[f"{j}_elems"], np.int64)
    off = src[f"{j}_offsets"]

    def P_of(r):
        return off[-r - 1] if r < 0 else bpos[r]

    nt = nneg = nzero = nsum = ngauss = 0
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        P = np.array([P_of(r) for r in v])
        tabs = HEX_TETS if m == 4 else WEDGE_TETS
        s = 0.0
        for t in tabs:
            nt += 1
            d = float(np.linalg.det(np.stack([P[t[1]] - P[t[0]],
                                              P[t[2]] - P[t[0]],
                                              P[t[3]] - P[t[0]]])))
            s += d
            if d < -1e-12:
                nneg += 1
            elif abs(d) <= 1e-12:
                nzero += 1
        if s <= 0:
            nsum += 1
        if m == 4 and detJ_min(P) <= 0:
            ngauss += 1
    print(f"{j:22s} {len(E):5d} {nt:8d} {nneg:7d} {nzero:6d} {nsum:9d} {ngauss:8d}")
