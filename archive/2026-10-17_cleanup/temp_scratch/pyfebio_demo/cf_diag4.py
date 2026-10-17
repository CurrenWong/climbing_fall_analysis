"""CARTFIX 诊断 4：是「骨面本身折了」还是「offset 场噪声」？

A. 底面上每个节点：相邻底面三角法向最大夹角（曲面折角），对比 disp 方向扩散角
B. 坏单元所在区域的底面法向：是否 up-facing / side-facing 混合（补丁边缘翻卷）
C. 坏单元是否靠近补丁边界（边界节点的 star 不闭合）
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
SRC = np.load(ROOT / "temp/thums/cartilage_mesh.npz", allow_pickle=True)
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
bids, bxyz = SRC["bone_node_ids"], SRC["bone_xyz"]
bpos = {int(v): bxyz[i] for i, v in enumerate(bids)}


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    off = np.asarray(CM1[f"{j}_is_offset"], bool)
    pair = {}
    for row in E:
        c = [ix(x) for x in row]
        for b, t in zip(c[:3], c[3:]):
            pair[b] = t
    disp = {b: N[pair[b]] - N[b] for b in pair}
    # 底面三角（每单元的底 3 点）法向
    star = {}
    for row in E:
        c = [ix(x) for x in row][:3]
        P = N[c]
        nv = np.cross(P[1] - P[0], P[2] - P[0])
        nv = nv / (np.linalg.norm(nv) + 1e-15)
        for q in c:
            star.setdefault(q, []).append(nv)
    # 边计数（判断是否补丁边界）
    ecnt = {}
    for row in E:
        c = [ix(x) for x in row][:3]
        for a in range(3):
            e = tuple(sorted((c[a], c[(a + 1) % 3])))
            ecnt[e] = ecnt.get(e, 0) + 1

    print("=" * 84)
    print(f"{j}: 底节点 {len(star)}")
    rows = []
    for b, ns in star.items():
        ns = np.array(ns)
        sp_face = max((np.degrees(np.arccos(np.clip(ns[i] @ ns[k], -1, 1)))
                       for i in range(len(ns)) for k in range(i + 1, len(ns))), default=0.0)
        ds = list(disp[b])
        dsp = np.degrees(np.arccos(np.clip(np.dot(ds / np.linalg.norm(ds), ns[0]), -1, 1)))
        rows.append((b, len(ns), sp_face, dsp, np.linalg.norm(disp[b]),
                     bpos.get(None) if False else None))
    # disp 方向在相邻节点之间的最大夹角（这才是 fold 驱动量）
    adj = {}
    for row in E:
        c = [ix(x) for x in row][:3]
        for a in range(3):
            for bb in range(3):
                if a != bb:
                    adj.setdefault(c[a], set()).add(c[bb])
    out = []
    for b, nface, sp, dsp, mag, _ in rows:
        nb = list(adj.get(b, []))
        dd = 0.0
        u = disp[b] / np.linalg.norm(disp[b])
        for x in nb:
            v = disp[x] / np.linalg.norm(disp[x])
            dd = max(dd, np.degrees(np.arccos(np.clip(np.dot(u, v), -1, 1))))
        # 该节点周围最短边长
        mind = min((np.linalg.norm(N[b] - N[x]) for x in nb), default=np.inf)
        # 是否边界节点（某条边只出现 1 次）
        bnd = any(ecnt[tuple(sorted((b, x)))] == 1 for x in nb)
        out.append((b, nface, sp, dd, mind, bnd, mag))
    out.sort(key=lambda z: -z[3])
    print("  ── 相邻节点 disp 方向夹角最大的 12 个 ──")
    print(f"     {'node':>5} {'面数':>4} {'底面折角°':>9} {'disp方向最大夹角°':>16} {'最短边':>7} {'边界?':>5} {'|disp|':>7}")
    for b, nf, sp, dd, mind, bnd, mag in out[:12]:
        print(f"     {b:>5} {nf:>4} {sp:9.2f} {dd:16.2f} {mind:7.3f} {str(bnd):>5} {mag:7.4f}")
    print("  ── 中位数对比 ──")
    print(f"     底面折角  中位={np.median([z[2] for z in out]):.2f}°   "
          f"disp方向夹角 中位={np.median([z[3] for z in out]):.2f}°  "
          f"  90分位={np.percentile([z[3] for z in out], 90):.2f}°")
    print(f"     边界节点数 {sum(1 for z in out if z[5])} / {len(out)}")

    # 坏单元区域：底面法向的"混合程度"（up vs side）
    print("  ── 坏单元附近底面法向分布 ──")
    ecnt2 = {}
    for row in E:
        c = [ix(x) for x in row][:3]
        for a in range(3):
            e = tuple(sorted((c[a], c[(a + 1) % 3])))
            ecnt2[e] = ecnt2.get(e, 0) + 1
    ctr = N[[ix(x) for x in E[0]][:3]].mean(0)
    zs = []
    for row in E:
        c = [ix(x) for x in row]
        P = N[c[:3]]
        nv = np.cross(P[1] - P[0], P[2] - P[0])
        nv = nv / (np.linalg.norm(nv) + 1e-15)
        zs.append(nv)
    zs = np.array(zs)
    print(f"     底面法向 z 分量: 中位={np.median(zs[:,2]):+.3f}  "
          f"z<0 的三角 {int((zs[:,2]<0).sum())}/{len(zs)}  "
          f"|z|<0.5（近竖直）的 {int((np.abs(zs[:,2])<0.5).sum())}")
