"""按 FEBio 的判据检查：hex8 在 8 个高斯点的 detJ、tet4 体积、
以及各部件「唯一节点数」直方图（找被静默丢弃的单元）。"""
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict

OUT = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
db = np.load(OUT / "thums_lowerlimb.npz", allow_pickle=True)
nid_all, xyz_all = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid_all)}

PARTS = {
    "tibia_cort": "81000700_solid", "tibia_spon": "81000600_solid",
    "tibia_marrow": "81000601_solid", "fibula_cort": "81000900_solid",
    "fibula_spon": "81000800_solid", "talus_cort": "81001100_solid",
    "talus_spon": "81001000_solid", "calcaneus_cort": "81001300_solid",
    "calcaneus_spon": "81001200_solid", "navicular_spon": "81001400_solid",
}

G = [-0.577350269189626, 0.577350269189626]


def detJ_hex8(P):
    """返回 8 个高斯点的 detJ 最小值。"""
    H = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                  [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    worst = 1e30
    for g in G:
        for h in G:
            for r in G:
                dN = np.zeros((3, 8))
                for i in range(8):
                    xi, et, ze = H[i]
                    dN[0, i] = 0.125 * xi * (1 + et * h) * (1 + ze * r)
                    dN[1, i] = 0.125 * et * (1 + xi * g) * (1 + ze * r)
                    dN[2, i] = 0.125 * ze * (1 + xi * g) * (1 + et * h)
                J = dN @ P
                worst = min(worst, float(np.linalg.det(J)))
    return worst


print(f"{'部件':18s} {'唯一节点数直方图':40s} {'hex8':>6} {'minDetJ<0':>10}")
print("-" * 96)
bad_elems = defaultdict(list)
for dom, key in PARTS.items():
    hist = Counter()
    rows = []
    for row in db[key]:
        conn = []
        for x in row:
            xi = int(x)
            if xi not in (0, -1) and xi not in conn:
                conn.append(xi)
        hist[len(conn)] += 1
        rows.append(conn)
    nhex = nbad = 0
    for ci, conn in enumerate(rows):
        if len(conn) == 8:
            nhex += 1
            P = xyz_all[[idx[c] for c in conn]]
            dj = detJ_hex8(P)
            if dj <= 0:
                nbad += 1
                if len(bad_elems[dom]) < 6:
                    bad_elems[dom].append((ci, round(dj, 6), conn))
    flag = "  <<<" if nbad else ""
    print(f"{dom:18s} {str(dict(hist)):40s} {nhex:>6} {nbad:>10}{flag}")

print("\n=== 负 detJ 的 hex8 明细 ===")
for dom, lst in bad_elems.items():
    print(f"\n{dom}:")
    for ci, dj, conn in lst:
        P = xyz_all[[idx[c] for c in conn]]
        print(f"  行{ci} detJ={dj} 节点={conn}")
        print(f"        bbox 对角 {np.linalg.norm(P.max(0)-P.min(0)):.3f} mm")
