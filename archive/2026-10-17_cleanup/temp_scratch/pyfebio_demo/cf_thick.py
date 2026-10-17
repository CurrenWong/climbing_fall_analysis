"""CARTFIX 诊断 6：cart_m1 的「厚度」到底怎么量？为什么逐顶点对差是 [0.107, 2.33]？"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    off = np.asarray(CM1[f"{j}_is_offset"], bool)
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    print("=" * 76)
    print(f"{j}: 底 {int((~off).sum())} 顶 {int(off.sum())} elems {len(E)}")
    print(f"  C[:,3:] 中不是 offset 的个数: {int((~off[T]).sum())}")
    print(f"  C[:,:3] 中是 offset 的个数  : {int(off[B].sum())}")
    # 每列两两配对是否唯一（同一底节点是否只有一个顶）
    pairmap = {}
    for r in range(len(E)):
        for a in range(3):
            pairmap.setdefault(int(B[r, a]), set()).add(int(T[r, a]))
    bad = {k: v for k, v in pairmap.items() if len(v) > 1}
    print(f"  底->顶 多对一违反: {len(bad)}  (样例 {list(bad.items())[:3]})")
    vp = np.linalg.norm(N[T] - N[B], axis=2)          # (n,3)
    print(f"  逐顶点对 |top-bottom|: min {vp.min():.4f} 中位 {np.median(vp):.4f} max {vp.max():.4f}")
    # 逐底节点的 disp（唯一化）
    pm = {}
    for r in range(len(E)):
        for a in range(3):
            pm[int(B[r, a])] = int(T[r, a])
    dm = {b: N[t] - N[b] for b, t in pm.items()}
    mags = np.array([np.linalg.norm(v) for v in dm.values()])
    print(f"  底节点 |disp|: min {mags.min():.4f} 中位 {np.median(mags):.4f} max {mags.max():.4f}")
    print(f"  ★ 注意：单元内 3 个顶点对的 |t-b| 与「节点 disp」不同 ⇒ 说明 T/B 列配对错位")
    # 检查配对：单元内 B[a] 的 disp（按 pm） vs 实际 T[r,a] 的位置差
    mis = 0
    for r in range(len(E)):
        for a in range(3):
            b_, t_ = int(B[r, a]), int(T[r, a])
            if not np.allclose(N[t_] - N[b_], dm[b_], atol=1e-9):
                mis += 1
    print(f"  配对不一致的 (elem,col) 数: {mis}")
    # 单元「真厚度」：底面质心到顶面质心的距离，及顶面质心到底面的投影
    cb, ct = N[B].mean(1), N[T].mean(1)
    d = ct - cb
    nb = np.cross(N[B[:, 1]] - N[B[:, 0]], N[B[:, 2]] - N[B[:, 0]])
    nb = nb / (np.linalg.norm(nb, axis=1, keepdims=True) + 1e-15)
    perp = np.abs(np.einsum("ij,ij->i", d, nb))
    print(f"  质心距离: min {np.linalg.norm(d,axis=1).min():.4f} max {np.linalg.norm(d,axis=1).max():.4f}")
    print(f"  质心距离·底面法向(垂直厚度): min {perp.min():.4f} 中位 {np.median(perp):.4f} max {perp.max():.4f}")
    i = int(np.argmax(vp) // 3)
    print(f"  最大逐顶点差的单元 #{i}: T={T[i]} B={B[i]}")
    for a in range(3):
        print(f"     col{a}: B={B[i,a]} T={T[i,a]}  disp={np.round(N[T[i,a]]-N[B[i,a]],5).tolist()}  "
              f"|d|={np.linalg.norm(N[T[i,a]]-N[B[i,a]]):.4f}  offB={off[B[i,a]]} offT={off[T[i,a]]}")
    # 是否顶/底都被标 True 等
    print(f"  节点表前 6 / 后 6 is_offset: {off[:6].tolist()} ... {off[-6:].tolist()}")
