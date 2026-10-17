"""定位负 Jacobian：逐单元算体积，按域汇总。"""
import numpy as np
from pathlib import Path
from collections import defaultdict

OUT = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
db = np.load(OUT / "thums_lowerlimb.npz", allow_pickle=True)
cg = np.load(OUT / "cartilage_mesh.npz")

BONE_PARTS = {
    "tibia_cort": "81000700_solid", "tibia_spon": "81000600_solid",
    "tibia_marrow": "81000601_solid", "fibula_cort": "81000900_solid",
    "fibula_spon": "81000800_solid", "talus_cort": "81001100_solid",
    "talus_spon": "81001000_solid", "calcaneus_cort": "81001300_solid",
    "calcaneus_spon": "81001200_solid", "navicular_spon": "81001400_solid",
}
WEDGE_T = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
HEX_T = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
         (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]

XYZ = None
stat = defaultdict(lambda: [0, 0, 0])   # dom -> [n, neg, zero]
rows = []
for dom, key in BONE_PARTS.items():
    for row in db[key]:
        conn = []
        for x in row:
            xi = int(x)
            if xi not in (0, -1) and xi not in conn:
                conn.append(xi)
        u = len(conn)
        if u == 8:
            for t in HEX_T:
                rows.append((dom, [conn[z] for z in t]))
        elif u == 4:
            rows.append((dom, conn))
        elif u == 6:
            for t in WEDGE_T:
                rows.append((dom, [conn[z] for z in t]))
        else:
            stat[dom][2] += 1
print(f"分解后 tet4 总数 {len(rows)}")

# 坐标
nid_all, xyz_all = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid_all)}


def vol(n):
    P = xyz_all[[idx[c] for c in n]]
    return float(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0])) / 6.0


for dom, n in rows:
    v = vol(n)
    st = stat[dom]
    st[0] += 1
    if v < 0:
        st[1] += 1
    if abs(v) < 1e-12:
        st[2] += 1

print(f"\n{'域':18s} {'tet数':>7} {'负体积':>7} {'零体积':>7}")
for dom in BONE_PARTS:
    st = stat[dom]
    flag = "  <<<" if st[1] or st[2] else ""
    print(f"{dom:18s} {st[0]:>7} {st[1]:>7} {st[2]:>7}{flag}")

# 韧带壳的绕向一致性（相邻面共享边应反向遍历）
print("\n=== 韧带壳绕向一致性 ===")
LIGP = [81100701, 81100801, 81100901, 81101001, 81101101,
        81101201, 81101202, 81101203, 81101204, 81101409]
for pid in LIGP:
    e = db[f"{pid:08d}_shell"]
    # 每条有向边出现次数：一致的绕向 -> 内部边正反各一次
    dc = defaultdict(int)
    for q in e:
        n = [int(v) for v in q]
        for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
            dc[(n[a], n[b])] += 1
    same = sum(1 for (a, b), c in dc.items() if c > 1 and dc.get((b, a), 0) == 0)
    opp = sum(1 for (a, b), c in dc.items() if dc.get((b, a), 0) > 0)
    # 面积法向 vs 质心外指（判断整体是否朝外）
    P = xyz_all[[idx[int(v)] for v in e.ravel()]]
    ctr = P.mean(0)
    nrm = np.zeros(3)
    for q in e:
        p = xyz_all[[idx[int(v)] for v in q]]
        nrm += np.cross(p[1] - p[0], p[2] - p[0]) / 2
    print(f"  {pid}: 同向重复边 {same}, 反向配对边 {opp}, "
          f"净法向模 {np.linalg.norm(nrm):8.1f}")
