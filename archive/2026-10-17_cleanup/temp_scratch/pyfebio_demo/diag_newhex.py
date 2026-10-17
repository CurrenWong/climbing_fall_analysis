"""我的新 hex 网格在「真实判据」下表现如何？+ 底面 quads 本身是不是瘦长条？"""
import numpy as np
from pathlib import Path

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
new = np.load(THUMS / "cartilage_mesh_hex.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}

G = 1 / np.sqrt(3)
GP = [(-G, -G, -G), (G, -G, -G), (G, G, -G), (-G, G, -G),
      (-G, -G, G), (G, -G, G), (G, G, G), (-G, G, G)]
REF = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
       (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
HEX_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
            (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]


def detJ_hex8_min(P):
    worst = 1e30
    for s in GP:
        dN = np.zeros((8, 3))
        for i in range(8):
            a, b, c = REF[i]
            dN[i, 0] = a * (1 + b * s[1]) * (1 + c * s[2]) / 8
            dN[i, 1] = b * (1 + a * s[0]) * (1 + c * s[2]) / 8
            dN[i, 2] = c * (1 + a * s[0]) * (1 + b * s[1]) / 8
        worst = min(worst, float(np.linalg.det(P.T @ dN)))
    return worst


def tdet(P):
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0], P[3] - P[0]])))


print("=" * 80)
print("A. 我的新网格：Gauss det(J) 判据 + 6-tet 内部各 tet 符号")
print("=" * 80)
for j in ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]:
    nd = new[f"{j}_nodes"]
    E = new[f"{j}_elems"]

    def P_of(r):
        return nd[-r - 1] if r < 0 else bpos[r]

    nhex, nwedge = 0, 0
    ngauss_bad, nminj = 0, 1e9
    ntneg = 0
    ntet_in_neg_hex = 0
    for row in E:
        blk = [int(x) for x in row if int(x) != 0]
        m = len(blk) // 2
        P = np.array([P_of(r) for r in blk])
        if m == 4:
            nhex += 1
            dj = detJ_hex8_min(P)
            nminj = min(nminj, dj)
            if dj <= 0:
                ngauss_bad += 1
            for t in HEX_TETS:
                if tdet(P[list(t)]) <= 0:
                    ntneg += 1
                    ntet_in_neg_hex += 1
        else:
            nwedge += 1
    print(f"{j:22s} hex {nhex:5d} wedge {nwedge:5d} | "
          f"Gauss detJ<0: {ngauss_bad:5d} ({100*ngauss_bad/max(nhex,1):5.1f}% of hex)"
          f" min={nminj:+.4e} | 6-tet 中负 tet: {ntneg}")

print()
print("=" * 80)
print("B. 底面 quads 本身的质量（这是根因候选）")
print("=" * 80)
print(f"{'关节':22s} {'底面quad':>8s} {'面积中位':>9s} {'最长棱中位':>10s} "
      f"{'长宽比中位':>10s} {'最小内角':>9s}")
for j in ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]:
    E = np.asarray(src[f"{j}_elems"], np.int64)
    ar, lr, amin = [], [], []
    nq = 0
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        if m != 4:
            continue
        nq += 1
        P = np.array([bpos[r] for r in v[:4]])
        e = [np.linalg.norm(P[(k + 1) % 4] - P[k]) for k in range(4)]
        a = 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0])) + \
            0.5 * np.linalg.norm(np.cross(P[2] - P[0], P[3] - P[0]))
        ar.append(a)
        lr.append(max(e) / max(min(e), 1e-9))
        # 最小内角
        angs = []
        for k in range(4):
            u = P[(k - 1) % 4] - P[k]
            w = P[(k + 1) % 4] - P[k]
            c = np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w) + 1e-15)
            angs.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
        amin.append(min(angs))
    if nq:
        print(f"{j:22s} {nq:8d} {np.median(ar):9.3f} {np.median(lr):10.2f} "
              f"{np.median(amin):9.1f}°   [最小内角 p5={np.percentile(amin,5):.1f}°, "
              f"<20° 的占比 {100*np.mean(np.array(amin)<20):.0f}%]")
