"""CARTFIX：把 cart_patch.py 的解缠函数单独抽出来扫参（不写文件）。

做法：读 cart_patch.py 源码，只 exec 到 `src = np.load(` 之前的纯函数部分。
"""
import io
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
SRC = io.open(ROOT / "scripts/ankle_fe/cart_patch.py", encoding="utf-8").read()
head = SRC.split("\nsrc = np.load(")[0]
ns = {"__name__": "cp_head", "__file__": str(ROOT / "scripts/ankle_fe/cart_patch.py")}
exec(compile(head, "cart_patch_head", "exec"), ns)
untangle_disp = ns["untangle_disp"]
wedge_dets = ns["wedge_dets"]

CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
S0 = np.load(ROOT / "temp/thums/cartilage_mesh.npz", allow_pickle=True)
bids, bxyz = S0["bone_node_ids"], S0["bone_xyz"]
bpos = {int(v): bxyz[i] for i, v in enumerate(bids)}


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


def prep(j):
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    # 底面节点位置 -> 骨 id
    key = {}
    for i, v in enumerate(bids):
        key[tuple(np.round(bxyz[i], 9))] = int(v)
    tris, disp = [], {}
    for e in range(len(E)):
        t = []
        for a in range(3):
            p = N[int(B[e, a])]
            k = tuple(np.round(p, 9))
            bid = key.get(k)
            if bid is None:
                d = np.linalg.norm(bxyz - p, axis=1)
                bid = int(bids[int(np.argmin(d))])
            t.append(bid)
            disp[bid] = N[int(T[e, a])] - N[int(B[e, a])]
        tris.append(tuple(t))
    return tris, disp, N, E, C


for j in ["tibiotalar_tibia", "tibiotalar_talus", "subtalar_talus", "subtalar_calcaneus"]:
    tris, disp, N, E, C = prep(j)
    t0 = time.time()
    dnew, st = untangle_disp(dict(disp), tris, bpos)
    dt = st["after"]
    print(f"[{j:22s}] {len(tris):5d} 三角  worst {st['before']:+9.4f} → {st['after']:+9.4f}  "
          f"≤0 {st['nneg_before']:3d}→{st['nneg_after']:3d}  ≤1e-2 {st['nbad_after']:3d}  "
          f"Δd {st['dmove']:.3f}mm  {time.time()-t0:.1f}s", flush=True)
