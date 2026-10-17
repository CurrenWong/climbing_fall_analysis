"""对照实验：同一「负体积」判据套到**现有基线软骨**上。
若基线也大面积报负 ⇒ 判据本身错（6-tet 分解对扁 hex 不可靠），不是网格错。
"""
import numpy as np
from pathlib import Path

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}

HEX_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
            (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
WEDGE_TETS = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]


def tvol(P):
    return float(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0])) / 6.0


def cvol(P, m):
    return sum(tvol(P[list(t)]) for t in (HEX_TETS if m == 4 else WEDGE_TETS))


# ── hex8 在 8 个 Gauss 点的 det(J)（FEBio 的真正判据）──
G = 1 / np.sqrt(3)
GP = [(-G, -G, -G), (G, -G, -G), (G, G, -G), (-G, G, -G),
      (-G, -G, G), (G, -G, G), (G, G, G), (-G, G, G)]
HN = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
      (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]


def detJ_hex8(P, xi, eta, ze):
    J = np.zeros((3, 3))
    for k in range(8):
        sx, sy, sz = HN[k]
        dx = sx * (1 + sy * eta) * (1 + sz * ze) / 8
        dy = sy * (1 + sx * xi) * (1 + sz * ze) / 8
        dz = sz * (1 + sx * xi) * (1 + sy * eta) / 8
        J[0] += dx * P[k]
        J[1] += dy * P[k]
        J[2] += dz * P[k]
    return float(np.linalg.det(J))


print(f"{'关节':22s} {'格数':>5s} {'tri%':>6s} | "
      f"{'6tet负':>8s} | {'Gauss detJ<0 的格数':>20s}")
print("-" * 78)
TOT = {}
for j in ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]:
    E = np.asarray(src[f"{j}_elems"], np.int64)
    off = src[f"{j}_offsets"]

    def P_of(r):
        return off[-r - 1] if r < 0 else bpos[r]

    nneg6, ngauss_bad, ntri, n = 0, 0, 0, 0
    minj = 1e9
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        if not v:
            continue
        m = len(v) // 2
        n += 1
        if m == 3:
            ntri += 1
        P = np.array([P_of(r) for r in v])
        if cvol(P, m) <= 0:
            nneg6 += 1
        if m == 4:
            ds = [detJ_hex8(P, *g) for g in GP]
            if min(ds) <= 0:
                ngauss_bad += 1
            minj = min(minj, min(ds))
    TOT[j] = (n, nneg6, ngauss_bad)
    print(f"{j:22s} {n:5d} {100*ntri/n:5.1f}% | {nneg6:8d} | {ngauss_bad:20d}"
          f"   (min detJ={minj:+.4f})")

print("-" * 78)
s6 = sum(v[1] for v in TOT.values())
sg = sum(v[2] for v in TOT.values())
sn = sum(v[0] for v in TOT.values())
print(f"基线合计 {sn} 格：6-tet 判据报负 {s6} ({100*s6/sn:.1f}%)；"
      f"Gauss detJ<0 只有 {sg} ({100*sg/sn:.1f}%)")
print()
print("★ 若基线也大面积报负 ⇒ **6-tet 分解判据对扁 hex 不可靠**，")
print("  应以 Gauss 点 det(J) 为准（FEBio 实际用的就是这个）。")
