"""量化「薄片尾巴」：我的新网格 vs 基线，拆出 tet 的归一化体积分布。
归一化体积 nv = V / lmax^3；正四面体 ≈ 0.1179。nv < 1e-3 即极端薄片。
"""
import numpy as np
from pathlib import Path

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
new = np.load(THUMS / "cartilage_mesh_hex.npz", allow_pickle=True)

HEX_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
            (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
WEDGE_TETS = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}


def tets_of(j, NPZ, bone_ref_positive):
    E = np.asarray(NPZ[f"{j}_elems"], np.int64)
    if f"{j}_nodes" in NPZ.files:
        nd = NPZ[f"{j}_nodes"]

        def P_of(r):
            return nd[-r - 1] if r < 0 else bpos[r]
    else:
        off = NPZ[f"{j}_offsets"]

        def P_of(r):
            return off[-r - 1] if r < 0 else bpos[r]

    out = []
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        P = np.array([P_of(r) for r in v])
        tabs = HEX_TETS if m == 4 else WEDGE_TETS
        for t in tabs:
            out.append(P[list(t)])
    return out


print("=" * 88)
print("拆出 tet 的质量分布   nv = V/lmax³（正四面体 0.1179；<1e-3 = 极端薄片）")
print("=" * 88)
print(f"{'来源':22s} {'tet数':>7s} {'nv中位':>9s} {'nv p1':>10s} {'nv p0.1':>10s} "
      f"{'nv<1e-3':>9s} {'nv<1e-2':>9s} {'V<1e-4':>8s}")
for tag, NPZ in (("基线 cartilage_mesh", src), ("新 hex(细)", new)):
    for j in ["tibiotalar_tibia", "subtalar_calcaneus"]:
        ts = tets_of(j, NPZ, True)
        nv, V = [], []
        for P in ts:
            vol = abs(float(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0])) / 6.0)
            lm = max(np.linalg.norm(P[(i + 1) % 4] - P[i]) for i in range(4))
            l2 = max(np.linalg.norm(P[j2] - P[i]) for i in range(4) for j2 in range(i + 1, 4))
            nv.append(vol / (l2 ** 3 + 1e-30))
            V.append(vol)
        nv = np.array(nv)
        V = np.array(V)
        print(f"{tag[:11]+' '+j[:10]:22s} {len(ts):7d} {np.median(nv):9.2e} "
              f"{np.percentile(nv,1):10.2e} {np.percentile(nv,0.1):10.2e} "
              f"{100*np.mean(nv<1e-3):8.2f}% {100*np.mean(nv<1e-2):8.2f}% "
              f"{np.sum(V<1e-4):8d}")
print()
print("★ 若新网格的 nv<1e-3 尾巴显著大于基线 ⇒ 细分把「母单元的翘曲」放大成了薄片")
print("★ 对策候选：① 按 nv 剔除薄片(会留洞) ② 细分后局部优化节点 ③ 重新分网底面")
