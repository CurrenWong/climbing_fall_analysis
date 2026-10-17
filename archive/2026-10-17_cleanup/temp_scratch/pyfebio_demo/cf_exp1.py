"""CARTFIX 实验 1：偏置向量场的「方向平滑」能否消掉翻转？

底面保持完全不动（⇒ tie/接触面几何不变），只改顶面节点的位置：
  d(r) = t(r) * u(r)   ->   d'(r) = t(r) * u_smooth(r)
扫描平滑次数 K，报告：bad 数 / 厚度范围 / 顶面最大位移。
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


def worst_all(P6):
    """批量：P6 (n,6,3) -> (n,) worst detJ"""
    out = np.full(len(P6), np.inf)
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
                d = np.einsum("nia,ic->nac", P6, G)
                out = np.minimum(out, np.linalg.det(d))
    return out


THRESH = float(__import__("os").environ.get("CF_THRESH", "1e-2"))
for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    off = np.asarray(CM1[f"{j}_is_offset"], bool)
    bt = np.array([i for i in range(len(N)) if not off[i]])
    tp = np.array([i for i in range(len(N)) if off[i]])
    idx_of = {int(v): i for i, v in enumerate(bt)}
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    # 逐底节点 disp（取该底的顶配对）
    nb_of_top = {}
    for row in C:
        for b, t in zip(row[:3], row[3:]):
            nb_of_top.setdefault(int(b), []).append(int(t))
    print("=" * 84)
    print(f"{j}: 底 {len(bt)} 顶 {len(tp)} elems {len(E)}   "
          f"底的顶配对数唯一? {all(len(v)==1 for v in nb_of_top.values())}")
    disp = np.zeros((len(bt), 3))
    multi = {b: ts for b, ts in nb_of_top.items() if len(ts) > 1}
    if multi:
        print(f"  ⚠ 底节点有多个顶配对的: {len(multi)} 个，样例 "
              f"{list(multi.items())[:3]}")
        for b, ts in list(multi.items())[:2]:
            print(f"     底节点 {b} -> 顶 {ts}；各自 disp = "
                  f"{[np.round(N[t]-N[b], 4).tolist() for t in ts]}")
    for b, ts in nb_of_top.items():
        disp[idx_of[b]] = N[ts[0]] - N[b]
    D = {i: disp[idx_of[i]] for i in range(len(N)) if i in idx_of}
    # 邻接（底节点）
    adj = {i: set() for i in bt}
    for row in C:
        for a in range(3):
            for b2 in range(3):
                if a != b2:
                    adj[int(row[a])].add(int(row[b2]))

    def build(dnew):
        """用新 disp 重建顶面位置 + 返回 worst detJ 与厚度。"""
        N2 = N.copy()
        for i in bt:
            for t in nb_of_top[i]:
                N2[t] = N[i] + dnew[i]
        P6 = N2[C]
        w = worst_all(P6)
        th = np.linalg.norm(N2[T] - N2[B], axis=1)
        return w, th, N2

    disp0 = {i: D[i] for i in D}
    w0, th0, _ = build(disp0)
    print(f"  K=0(原): bad(≤{THRESH}) {int((w0<=THRESH).sum()):>3}  min detJ {w0.min():+9.4f}  "
          f"厚度 [{th0.min():.4f}, {th0.max():.4f}] 中位 {np.median(th0):.4f}")
    for K in (1, 2, 4, 8, 16, 32, 64):
        cur = {i: D[i].copy() for i in D}
        for _ in range(K):
            nxt = {}
            for i in cur:
                u = cur[i] / (np.linalg.norm(cur[i]) + 1e-15)
                nbu = [cur[x] / (np.linalg.norm(cur[x]) + 1e-15) for x in adj[i]]
                if nbu:
                    v = u + np.mean(nbu, axis=0)
                else:
                    v = u
                nxt[i] = v / (np.linalg.norm(v) + 1e-15)
            cur = {i: nxt[i] * np.linalg.norm(D[i]) for i in D}   # 模长保持
        w, th, N2 = build(cur)
        moved = max(np.linalg.norm(N2[t] - N[t]) for t in tp)
        print(f"  K={K:<3}: bad(≤{THRESH}) {int((w<=THRESH).sum()):>3}  min detJ {w.min():+9.4f}  "
              f"厚度 [{th.min():.4f}, {th.max():.4f}] 中位 {np.median(th):.4f}  "
              f"顶面最大位移 {moved:.4f}mm")
    # 只平滑方向、不平滑模长做对照：把模长也一起平滑
    for K in (16, 32):
        cur = {i: D[i].copy() for i in D}
        for _ in range(K):
            nxt = {}
            for i in cur:
                nb = [cur[x] for x in adj[i]]
                v = cur[i] + (np.mean(nb, axis=0) if nb else 0)
                nxt[i] = v
            cur = nxt
        w, th, _ = build(cur)
        print(f"  K={K}(向量全平滑,含模长): bad {int((w<=THRESH).sum()):>3}  min detJ {w.min():+9.4f}  "
              f"厚度 [{th.min():.4f}, {th.max():.4f}]")
