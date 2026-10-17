"""A1 失败分析（只测不改）：subtalar 软骨接入后 FEBio 负 Jacobian 的几何/拓扑嫌疑点。

1) 孤立节点：既不在存活单元里、又不在 tie 面里 ⇒ **完全自由** ⇒ 位移跑飞
2) 存活单元连通分量（岛）：是否有岛既不 tie、又只有接触（法向约束）⇒ 漂浮
3) 存活单元的完整采样 det(J)：≤0 的先翻（这就是 FEBio 报的负 Jacobian 候选）
4) tie 面覆盖：tie primary 里有多少面引用了「已被闸门剔除」的单元
用法: python temp/pyfebio_demo/a1_cart_bisect.py
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
cg = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
DET_MIN = 0.0
J = ["tibiotalar_tibia", "tibiotalar_talus", "subtalar_talus",
     "subtalar_calcaneus"]

PENTA_FACES = [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]


def cnidx(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


def penta6_dN(r, s, t):
    L = (1 - r - s, r, s)
    dL = ((-1.0, -1.0), (1.0, 0.0), (0.0, 1.0))
    G = np.zeros((6, 3))
    for i in range(3):
        G[i] = ((1 - t) / 2 * dL[i][0], (1 - t) / 2 * dL[i][1], -L[i] / 2)
        G[i + 3] = ((1 + t) / 2 * dL[i][0], (1 + t) / 2 * dL[i][1], L[i] / 2)
    return G


def detJ(X, idx, full=True):
    P = X[list(idx)]
    RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3) + ((1.0,) if full else ())
    TS = (0.0, 0.25, 0.5, 0.75, 1.0)
    worst = None
    for r in RS:
        for s in RS:
            if r + s > 1.0 + 1e-9:
                continue
            for t in TS:
                d = float(np.linalg.det(P.T @ penta6_dN(r, s, t)))
                worst = d if worst is None else min(worst, d)
    return worst


def sv_penta6(P):
    dt = (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3.0
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0], dt])))


def gate(X, idx):
    """返回 (存活?, detmin原序, 用的排列)。与 ankle_lig_feb 一致（DET_MIN=0.0）。"""
    d = detJ(X, idx, full=False)
    if sv_penta6(X[idx]) > 0 and d > DET_MIN:
        return True, d, idx
    alt = idx[3:] + idx[:3]
    d2 = detJ(X, alt, full=False)
    if sv_penta6(X[alt]) > 0 and d2 > DET_MIN:
        return True, d2, alt
    return False, min(d, d2), idx


for j in J:
    X = np.asarray(cg[f"{j}_nodes"], float)
    isoff = np.asarray(cg[f"{j}_is_offset"], bool)
    allnodes = set(range(len(X)))
    alive, dead = [], []
    for row in cg[f"{j}_elems"]:
        raw = [int(x) for x in row if int(x) != 0]
        if len(raw) != 6:
            continue
        idx = [cnidx(x) for x in raw]
        ok, d, used = gate(X, idx)
        (alive if ok else dead).append((tuple(used), d))
    alive_ids = [a for a, _ in alive]
    tie_faces = []
    for row in cg[f"{j}_tie_faces"]:
        fs = [cnidx(x) for x in row if int(x) != 0]
        if len(fs) == 3:
            tie_faces.append(tuple(fs))
    tie_nodes = {x for f in tie_faces for x in f}
    elem_nodes = {x for a in alive_ids for x in a}

    free = allnodes - elem_nodes - tie_nodes
    dead_nodes = allnodes - elem_nodes          # 不在存活单元里的节点
    # 孤立节点里「既不在存活单元也不在 tie」= 真·自由
    # 连通分量（按共享节点）
    parent = {n: n for n in elem_nodes}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def uni(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a in alive_ids:
        for n in a[1:]:
            uni(a[0], n)
    comps = {}
    for n in elem_nodes:
        comps.setdefault(find(n), set()).add(n)
    # 每个连通分量：大小、是否含 tie 节点、是否含 offset（接触）节点
    isl = []
    for root, mem in comps.items():
        isl.append(dict(n=len(mem), tie=len(mem & tie_nodes),
                        off=len(mem & set(np.where(isoff)[0]))))
    isl.sort(key=lambda d: d["n"])
    # 死单元的 tie 面（FEBio 报 invalid 的来源）
    dead_set = set()
    for a, _ in dead:
        dead_set.add(a)
    tie_from_dead = sum(1 for f in tie_faces
                        if tuple(sorted(f)) not in {tuple(sorted(a[:3])) for a in alive_ids}
                        and all(x in dead_nodes or True for x in f))
    # 完整采样：存活单元里 ≤0 的（FEBio 负 Jacobian 的直接候选）
    fullneg = [(a, d) for a, d in alive if detJ(X, a, full=True) <= 0]
    print(f"\n=== {j} ===")
    print(f"  单元: 存活 {len(alive)} / 剔除 {len(dead)}   节点总数 {len(allnodes)}")
    print(f"  节点: 在存活单元 {len(elem_nodes)} | 在 tie 面 {len(tie_nodes)} | "
          f"**两者都不在(真自由) {len(free)}** {sorted(free)[:12]}")
    print(f"  连通分量 {len(isl)} 个；最小的 5 个: {isl[:5]}")
    print(f"  存活但**完整采样 det(J)<=0** 的单元: {len(fullneg)} 个 "
          f"(min={min((d for _, d in fullneg), default=float('nan')):.3e})")
    print(f"  tie 面总数 {len(tie_faces)}（FEBio 报 invalid 的应来自被剔除单元的底面）")
