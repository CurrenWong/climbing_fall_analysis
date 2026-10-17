"""CARTFIX 诊断 2：扭转的真实来源。

思路：底面节点位置 == 原骨节点位置（精确抬回）⇒ 可反查骨节点 id。
  disp(node) = 顶面位置 − 底面位置。
若 disp 场在某个骨节点上"跳变"（spike），则该点周围所有楔形都被剪切/旋转
⇒ 小底三角 + 大 |∇disp| = 翻转。

打印：
  A. 坏单元 vs 邻居：底/顶三角法向、逐顶点 disp、扭角、min-detJ 的 (r,s,t) 位置
  B. disp 场的局部光滑性：per-node disp 与其邻居中位数的偏差 / 该点最短边长
  C. 顶面相对底面的极分解（旋转角）
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
SRC = np.load(ROOT / "temp/thums/cartilage_mesh.npz", allow_pickle=True)
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
bids, bxyz = SRC["bone_node_ids"], SRC["bone_xyz"]
bpos = {int(v): bxyz[i] for i, v in enumerate(bids)}
# 位置 -> 骨节点 id（精确）
key2id = {}
for i, v in enumerate(bids):
    key2id[tuple(np.round(bxyz[i], 9))] = int(v)


def bone_id(p, tol=1e-9):
    k = tuple(np.round(p, 9))
    if k in key2id:
        return key2id[k]
    # 容差匹配
    d = np.linalg.norm(bxyz - p, axis=1)
    j = int(np.argmin(d))
    return int(bids[j]) if d[j] < 1e-6 else None


def det_pt(P, r, s, t):
    G = np.array([
        [-(1 - t), -(1 - t), -(1 - r - s)],
        [1 - t, 0.0, -r],
        [0.0, 1 - t, -s],
        [-t, -t, (1 - r - s)],
        [t, 0.0, r],
        [0.0, t, s],
    ])
    return float(np.linalg.det(P.T @ G))


RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)


def worst(P):
    w, arg = None, None
    for r in RS:
        for s in RS:
            if r + s > 1 + 1e-9:
                continue
            for t in TS:
                v = det_pt(P, r, s, t)
                if w is None or v < w:
                    w, arg = v, (r, s, t)
    return w, arg


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    off = np.asarray(CM1[f"{j}_is_offset"], bool)
    bt = [i for i in range(len(N)) if not off[i]]
    tp = [i for i in range(len(N)) if off[i]]
    print("=" * 80)
    print(f"{j}: nodes {len(N)} (底 {len(bt)} / 顶 {len(tp)})  elems {len(E)}")

    # ── 反查：底节点 -> 骨 id；顶节点 -> 其"配对"底节点（同 disp）──
    B2BONE, disp = {}, {}
    for i in bt:
        bid = bone_id(N[i])
        B2BONE[i] = bid
    # 逐单元配对（按 tie_faces / elems）
    pair = {}
    for row in E:
        c = [ix(x) for x in row]
        for b, t in zip(c[:3], c[3:]):
            pair.setdefault(b, set()).add(t)
    multi = {b: len(s) for b, s in pair.items() if len(s) > 1}
    print(f"  底节点有多个顶配对: {len(multi)} 个（应为 0，否则共享不一致）")
    for b, ts in pair.items():
        for t in ts:
            disp[b] = N[t] - N[b]

    # ── B. disp 场光滑性 ──
    adj = {b: set() for b in bt}
    for row in E:
        c = [ix(x) for x in row]
        for a in range(3):
            for bb in range(3):
                if a != bb:
                    adj[c[a]].add(c[bb])
    devs, elen = [], []
    for b in bt:
        nb = [x for x in adj[b]]
        if not nb:
            continue
        med = np.median(np.array([disp[x] for x in nb]), axis=0)
        dev = np.linalg.norm(disp[b] - med)
        mind = min(np.linalg.norm(N[b] - N[x]) for x in nb)
        devs.append((b, dev, mind, dev / max(mind, 1e-9)))
    devs.sort(key=lambda z: -z[3])
    print("  ── disp 场局部跳变最猛的 10 个底节点（dev/最短边）──")
    print(f"     {'node':>5} {'骨id':>11} {'|disp|':>8} {'dev':>8} {'最短边':>8} {'dev/边':>8}")
    for b, dev, mind, ratio in devs[:10]:
        print(f"     {b:>5} {str(B2BONE[b]):>11} {np.linalg.norm(disp[b]):8.4f} "
              f"{dev:8.4f} {mind:8.4f} {ratio:8.3f}")

    # ── 坏单元 vs 邻居 ──
    print("  ── 坏单元逐条 ──")
    rows = []
    for k, row in enumerate(E):
        c = [ix(x) for x in row]
        B = c[:3]
        P = N[c]
        w, arg = worst(P)
        nb = np.cross(N[B[1]] - N[B[0]], N[B[2]] - N[B[0]])
        rows.append((k, w, arg, B, nb))
    bad = [r for r in rows if r[1] <= 1e-2]
    print(f"     坏单元 {len(bad)} 个；最差 4 个详情：")
    for k, w, arg, B, nb in sorted(bad, key=lambda z: z[1])[:4]:
        row = E[k]
        c = [ix(x) for x in row]
        T = c[3:]
        nt = np.cross(N[T[1]] - N[T[0]], N[T[2]] - N[T[0]])
        ab, at = np.linalg.norm(nb), np.linalg.norm(nt)
        cb, ct = N[B].mean(0), N[T].mean(0)
        ax = ct - cb
        print(f"\n     #{k}  worst detJ={w:+.4f} @ (r,s,t)={arg}  底面法向·轴={np.dot(nb,ax):+.3f} "
              f"顶面法向·轴={np.dot(nt,ax):+.3f}")
        print(f"        底面面积={ab/2:.3f}mm² 顶面面积={at/2:.3f}mm² "
              f"扭角={np.degrees(np.arccos(np.clip(np.dot(nb,nt)/(ab*at+1e-15),-1,1))):.2f}°")
        print(f"        顶点  骨id          |disp|     disp 向量")
        for q in range(3):
            b_ = B[q]
            print(f"          {b_:>3}  {str(B2BONE[b_]):>10}  {np.linalg.norm(disp[b_]):7.4f}   "
                  f"({disp[b_][0]:+.4f},{disp[b_][1]:+.4f},{disp[b_][2]:+.4f})")
        # 邻居
        nbr = []
        for m, row2 in enumerate(E):
            if m == k:
                continue
            c2 = [ix(x) for x in row2]
            if len(set(c2[:3]) & set(B)) >= 2:
                w2, _ = worst(N[c2])
                nbr.append((m, w2, len(set(c2[:3]) & set(B))))
        nbr.sort(key=lambda z: -z[2])
        print(f"        共享边的邻居(前 4): " + ", ".join(
            f"#{m}(共{s}点,detJ={w2:+.3f})" for m, w2, s in nbr[:4]))
        for m, w2, s in nbr[:3]:
            c2 = [ix(x) for x in E[m]]
            B2 = c2[:3]
            print(f"          #{m} 顶点骨id={[B2BONE[x] for x in B2]}  "
                  f"|disp|={[round(float(np.linalg.norm(disp[x])),4) for x in B2]}")

    # ── C. 极分解：顶面相对底面的面内旋转 ──
    print("  ── 顶面/底面 面内极分解（旋转角）──")
    rot_bad, rot_good = [], []
    for k, row in enumerate(E):
        c = [ix(x) for x in row]
        B, T = np.array(c[:3]), np.array(c[3:])
        P = N[c]
        nb = np.cross(N[B[1]] - N[B[0]], N[B[2]] - N[B[0]])
        nt = np.cross(N[T[1]] - N[T[0]], N[T[2]] - N[T[0]])
        if np.linalg.norm(nb) < 1e-12 or np.linalg.norm(nt) < 1e-12:
            continue
        e3 = nb / np.linalg.norm(nb)
        e1 = (N[B[1]] - N[B[0]])
        e1 = e1 / np.linalg.norm(e1)
        e2 = np.cross(e3, e1)
        # 底面 2D
        qb = np.array([[np.dot(N[B[m]] - N[B[0]], e1), np.dot(N[B[m]] - N[B[0]], e2)] for m in range(3)])
        qt = np.array([[np.dot(N[T[m]] - N[T[0]], e1), np.dot(N[T[m]] - N[T[0]], e2)] for m in range(3)])
        M = np.linalg.lstsq(qb[1:], qt[1:], rcond=None)[0].T
        u, sv, vt = np.linalg.svd(M)
        R = u @ vt
        ang = np.degrees(np.arctan2(R[1, 0], R[0, 0]))
        w, _ = worst(P)
        (rot_bad if w <= 1e-2 else rot_good).append((k, ang, sv[0] / max(sv[1], 1e-12), w))
    for nm, lst in (("坏", rot_bad), ("好", rot_good)):
        a = np.array([abs(z[1]) for z in lst])
        s = np.array([z[2] for z in lst])
        print(f"     {nm}: n={len(lst)}  |面内旋转| 中位={np.median(a):.2f}° max={a.max():.2f}°   "
              f"奇异值比 中位={np.median(s):.2f} max={s.max():.2f}")
