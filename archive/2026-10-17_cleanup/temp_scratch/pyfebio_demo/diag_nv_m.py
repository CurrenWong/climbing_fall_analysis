"""m=1/2/3 的软骨 tet 归一化体积分布 —— 判断「薄片是少数还是普遍」。"""
import numpy as np
from pathlib import Path

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}
WEDGE_TETS = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]

print(f"{'来源':16s} {'关节':20s} {'tet':>7s} {'nv中位':>9s} {'nv p1':>9s} "
      f"{'nv<1e-4':>8s} {'nv<1e-3':>8s} {'nv<1e-2':>8s}")
for tag, fn in (("基线(t0d)", "cartilage_mesh_t0d.npz"),
                ("m=1 ✔", "_pd02.npz"),
                ("m=2 ✘", "cartilage_mesh_patch_m2.npz"),
                ("m=3 ✘", "cartilage_mesh_patch_m3.npz")):
    p = THUMS / fn
    if not p.exists():
        continue
    npz = np.load(p, allow_pickle=True)
    for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
        if f"{j}_elems" not in npz.files:
            continue
        nd = npz[f"{j}_nodes"]
        E = npz[f"{j}_elems"]

        def P_of(r):
            return nd[-r - 1] if r < 0 else bpos[r]

        nv = []
        for row in E:
            v = [int(x) for x in row if int(x) != 0]
            m = len(v) // 2
            P = np.array([P_of(r) for r in v])
            for t in WEDGE_TETS:
                q = P[list(t)]
                vol = abs(float(np.dot(np.cross(q[1] - q[0], q[2] - q[0]),
                                       q[3] - q[0])) / 6.0)
                lm = max(np.linalg.norm(q[i] - q[k])
                         for i in range(4) for k in range(i + 1, 4))
                nv.append(vol / (lm ** 3 + 1e-30))
        nv = np.array(nv)
        print(f"{tag:16s} {j:20s} {len(nv):7d} {np.median(nv):9.2e} "
              f"{np.percentile(nv,1):9.2e} {100*np.mean(nv<1e-4):7.2f}% "
              f"{100*np.mean(nv<1e-3):7.2f}% {100*np.mean(nv<1e-2):7.2f}%")
print()
print("★ 正四面体 nv=0.1179。若 nv<1e-4 只有极少数 ⇒ 局部修即可；")
print("  若普遍很低 ⇒ 母三角本身瘦长，必须换三角化。")
