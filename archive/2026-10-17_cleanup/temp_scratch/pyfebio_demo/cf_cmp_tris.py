"""对比两个 cart_patch 输出目录的「底面三角集合」是否一致（判断基线用了哪个 CP_DIAG）。"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
S0 = np.load(ROOT / "temp/thums/cartilage_mesh.npz", allow_pickle=True)
bids, bxyz = S0["bone_node_ids"], S0["bone_xyz"]
key = {tuple(np.round(bxyz[i], 9)): int(v) for i, v in enumerate(bids)}


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


def tri_sets(p, j):
    d = np.load(p, allow_pickle=True)
    N = np.asarray(d[f"{j}_nodes"], float)
    E = np.asarray(d[f"{j}_elems"], np.int64)
    C = np.array([[ix(x) for x in row] for row in E])
    out = set()
    for e in range(len(E)):
        t = []
        for a in range(3):
            pt = N[int(C[e, a])]
            k = key.get(tuple(np.round(pt, 9)))
            if k is None:
                dd = np.linalg.norm(bxyz - pt, axis=1)
                k = int(bids[int(np.argmin(dd))])
            t.append(k)
        out.add(frozenset(t))
    return out


A = ROOT / "temp/thums/cart_m1/cartilage_mesh.npz"
B = ROOT / "temp/thums/cart_scratch_d02/cartilage_mesh.npz"
C = ROOT / "temp/thums/cart_m1_fix/cartilage_mesh.npz"
for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    a = tri_sets(A, j)
    b = tri_sets(B, j)
    c = tri_sets(C, j)
    print(f"{j}: cart_m1 {len(a)}  d02 {len(b)}  fix(默认DIAG=best) {len(c)}")
    print(f"   cart_m1 vs d02      : 相同 {len(a & b)} / 仅m1 {len(a - b)} / 仅d02 {len(b - a)}")
    print(f"   cart_m1 vs fix(best): 相同 {len(a & c)} / 仅m1 {len(a - c)} / 仅fix {len(c - a)}")
