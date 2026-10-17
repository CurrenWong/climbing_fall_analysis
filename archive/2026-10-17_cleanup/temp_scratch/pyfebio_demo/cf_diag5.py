"""CARTFIX 诊断 5：底面朝向一致性 + 真实折角 + 补丁形状。

A. 存储顺序给出的底面法向 nb（有向）在单元间的分布：相对"参考法向"翻转的比例
B. 无向折角：相邻底面 facet 法向取 |dot| 后的夹角（真实曲率）
C. 补丁形状：连通分量 / 边界边数 / bbox / PCA 展布
D. disp 方向夹角（无向）与坏单元位置的关系
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)


def worst(P):
    w = None
    for r in RS:
        for s in RS:
            if r + s > 1 + 1e-9:
                continue
            for t in TS:
                G = np.array([
                    [-(1 - t), -(1 - t), -(1 - r - s)],
                    [1 - t, 0.0, -r],
                    [0.0, 1 - t, -s],
                    [-t, -t, (1 - r - s)],
                    [t, 0.0, r],
                    [0.0, t, s],
                ])
                v = float(np.linalg.det(P.T @ G))
                if w is None or v < w:
                    w = v
    return w


for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    off = np.asarray(CM1[f"{j}_is_offset"], bool)
    tris = [[ix(x) for x in row][:3] for row in E]
    tops = [[ix(x) for x in row][3:] for row in E]
    print("=" * 84)
    print(f"{j}")
    C = np.array([N[t].mean(0) for t in tris])
    c0 = C.mean(0)
    Z = N[[q for t in tris for q in t]].mean(0)
    print(f"  底面节点 bbox  {N[[q for t in tris for q in t]].min(0)} .. "
          f"{N[[q for t in tris for q in t]].max(0)}")
    print(f"  底面三角质心 bbox {C.min(0)} .. {C.max(0)}")

    # ── 有向法向 vs 参考（用质心相对模型中心的方向做参考）──
    nvs = []
    for t in tris:
        P = N[t]
        nv = np.cross(P[1] - P[0], P[2] - P[0])
        nvs.append(nv / (np.linalg.norm(nv) + 1e-15))
    nvs = np.array(nvs)
    ref = Z - c0
    ref = ref / np.linalg.norm(ref)
    flip = (nvs @ ref) < 0
    print(f"  A. 有向法向 相对(补丁质心→全体质心) 为负的单元: {int(flip.sum())}/{len(nvs)}")
    print(f"     底面法向 z: 中位 {np.median(nvs[:,2]):+.3f}  |z|<0.5 的 {int((np.abs(nvs[:,2])<0.5).sum())}")
    # 邻接一致性：共享边的两三角，有向法向是否同号
    e2t = {}
    for i, t in enumerate(tris):
        for a in range(3):
            e = tuple(sorted((t[a], t[(a + 1) % 3])))
            e2t.setdefault(e, []).append(i)
    ncons, ntwin = 0, 0
    for e, ts in e2t.items():
        if len(ts) == 2:
            ntwin += 1
            if (nvs[ts[0]] @ nvs[ts[1]]) < 0:
                ncons += 1
            elif np.degrees(np.arccos(np.clip(nvs[ts[0]] @ nvs[ts[1]], -1, 1))) > 100:
                ncons += 1
    print(f"  A2. 共享边的三角对 {ntwin}，其中有向法向夹角 >100° 的 {ncons} "
          f"（= 翻转/极锐折）")

    # ── B. 无向折角（|dot|）──
    angs = []
    for e, ts in e2t.items():
        if len(ts) == 2:
            a = np.degrees(np.arccos(np.clip(abs(nvs[ts[0]] @ nvs[ts[1]]), -1, 1)))
            angs.append(a)
    angs = np.array(angs)
    print(f"  B. 无向折角(|dot|，90°=正交): 中位 {np.median(angs):.1f}°  "
          f"p90 {np.percentile(angs,90):.1f}°  max {angs.max():.1f}°  "
          f">60° 的 {int((angs>60).sum())}/{len(angs)}")

    # ── C. 连通分量 / 边界边 ──
    parent = list(range(len(tris)))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for e, ts in e2t.items():
        if len(ts) == 2:
            ra, rb = find(ts[0]), find(ts[1])
            if ra != rb:
                parent[ra] = rb
    comps = {}
    for i in range(len(tris)):
        comps.setdefault(find(i), []).append(i)
    nb_edge = sum(1 for e, ts in e2t.items() if len(ts) == 1)
    print(f"  C. 连通分量 {len(comps)} 个（大小 {sorted((len(v) for v in comps.values()), reverse=True)[:10]}）"
          f"  边界边 {nb_edge}  非流形边 {sum(1 for e,ts in e2t.items() if len(ts)>2)}")
    V = len(set(q for t in tris for q in t))
    Ee = len(e2t)
    print(f"     V={V} E={Ee} F={len(tris)}  Euler V-E+F={V-Ee+len(tris)} (单连通盘=1)")

    # ── D. disp 方向（无向）夹角 vs 坏单元 ──
    pair = {}
    for row in E:
        c = [ix(x) for x in row]
        for b, t in zip(c[:3], c[3:]):
            pair[b] = t
    disp = {b: N[pair[b]] - N[b] for b in pair}
    udir = {b: disp[b] / np.linalg.norm(disp[b]) for b in disp}
    adj = {}
    for t in tris:
        for a in range(3):
            for bb in range(3):
                if a != bb:
                    adj.setdefault(t[a], set()).add(t[bb])
    W = np.array([worst(N[[ix(x) for x in row]]) for row in E])
    print(f"  D. 坏单元 {int((W<=1e-2).sum())} / {len(W)}")
    print(f"     {'elem':>5} {'detJ':>9} {'disp方向在3顶点间的最大夹角(无向)°':>34} {'最短边':>7} {'底面积':>8}")
    for k in np.argsort(W)[:8]:
        t = tris[k]
        ds = [udir[q] for q in t]
        ang = max(np.degrees(np.arccos(np.clip(abs(np.dot(ds[a], ds[b])), -1, 1)))
                  for a in range(3) for b in range(a + 1, 3))
        ed = min(np.linalg.norm(N[t[a]] - N[t[b]]) for a in range(3) for b in range(a + 1, 3))
        ar = 0.5 * np.linalg.norm(np.cross(N[t[1]] - N[t[0]], N[t[2]] - N[t[0]]))
        print(f"     {k:>5} {W[k]:9.3f} {ang:34.2f} {ed:7.3f} {ar:8.3f}")
    # 全局：单元内 disp 方向最大夹角 vs detJ 的单调性
    allang = []
    for k, t in enumerate(tris):
        ds = [udir[q] for q in t]
        ang = max(np.degrees(np.arccos(np.clip(abs(np.dot(ds[a], ds[b])), -1, 1)))
                  for a in range(3) for b in range(a + 1, 3))
        allang.append(ang)
    allang = np.array(allang)
    gd = W > 1e-2
    print(f"     单元内 disp 方向最大夹角（无向）: 坏 {np.median(allang[~gd]):.2f}° vs 好 {np.median(allang[gd]):.2f}°")
    for thr in (20, 30, 40, 60):
        m = allang > thr
        if m.sum():
            print(f"      阈值 {thr}°: 命中 {int(m.sum())} 个，其中真坏 {int((W[m]<=1e-2).sum())}"
                  f"（漏 {int(((W<=1e-2) & ~m).sum())}）")
