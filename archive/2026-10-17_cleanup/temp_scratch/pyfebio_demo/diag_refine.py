"""单单元级诊断：红加密到底哪里错了。

对每个粗单元做一级红加密，比较：
  * fan 体积（6 四面体绕 0-6 对角线）
  * 精确体积（把 6 个面三角化后用散度定理）
并找出第一个出错的单元，打印其底面/外面的坐标。
"""
import numpy as np
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\temp\thums\cartilage_mesh.npz")
d = np.load(P, allow_pickle=True)
bpos = {int(v): d["bone_xyz"][i] for i, v in enumerate(d["bone_node_ids"])}

FAN8 = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
        (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
FAN6 = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]


def fan_vol(P8, n):
    t = 0.0
    for q in (FAN8 if n == 4 else FAN6):
        a, b, c, e = (P8[i] for i in q)
        t += np.dot(np.cross(b - a, c - a), e - a) / 6.0
    return t


def exact_vol(P8, n):
    """散度定理：把多面体表面三角化。base 在 z=0 附近、outer 在另一侧。
    对 hex8：6 个四边形面各拆 2 个三角形；对 wedge6：2 个三角面 + 3 个四边形。
    """
    if n == 4:
        quads = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
                 (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
        tris = []
        for q in quads:
            tris += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
    else:
        tris = [(0, 2, 1), (3, 4, 5)]
        for q in [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]:
            tris += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
    v = 0.0
    for t in tris:
        a, b, c = (P8[i] for i in t)
        v += np.dot(a, np.cross(b, c)) / 6.0
    return v


def refine(P8, n):
    b = [np.asarray(P8[i], float) for i in range(n)]
    o = [np.asarray(P8[n + i], float) for i in range(n)]
    m = [(b[i] + b[(i + 1) % n]) / 2 for i in range(n)]
    mo = [(o[i] + o[(i + 1) % n]) / 2 for i in range(n)]
    if n == 4:
        cb, co = np.mean(b, 0), np.mean(o, 0)
        subs = [[b[0], m[0], cb, m[3], o[0], mo[0], co, mo[3]],
                [m[0], b[1], m[1], cb, mo[0], o[1], mo[1], co],
                [cb, m[1], b[2], m[2], co, mo[1], o[2], mo[2]],
                [m[3], cb, m[2], b[3], mo[3], co, mo[2], o[3]]]
    else:
        subs = [[b[0], m[0], m[2], o[0], mo[0], mo[2]],
                [m[0], b[1], m[1], mo[0], o[1], mo[1]],
                [m[2], m[1], b[2], mo[2], mo[1], o[2]]]
    return subs


for j in ("tibiotalar_tibia",):
    E = np.asarray(d[f"{j}_elems"], dtype=np.int64)
    off = np.asarray(d[f"{j}_offsets"], float)
    print(f"=== {j} ===")
    n_hex = n_wedge = 0
    bad = []
    for ri, row in enumerate(E):
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        P = np.array([bpos[x] if x > 0 else off[-x - 1] for x in v])
        vf, ve = fan_vol(P, m), exact_vol(P, m)
        subs = refine(P, m)
        sf = sum(fan_vol(np.array(s), m) for s in subs)
        se = sum(exact_vol(np.array(s), m) for s in subs)
        nsub = len(subs)
        if m == 4:
            n_hex += 1
        else:
            n_wedge += 1
        if min(fan_vol(np.array(s), m) for s in subs) <= 0 or abs(se - ve) > 1e-6:
            bad.append((ri, m, vf, ve, sf, se, nsub))
    print(f"hex8 {n_hex} / wedge6 {n_wedge}；可疑 {len(bad)}")
    for ri, m, vf, ve, sf, se, ns in bad[:5]:
        print(f"  #{ri} n={m} 子单元 {ns}  fan {vf:+.4f}->{sf:+.4f}  "
              f"exact {ve:+.4f}->{se:+.4f}")
    if bad:
        ri = bad[0][0]
        v = [int(x) for x in E[ri] if int(x) != 0]
        m = len(v) // 2
        P = np.array([bpos[x] if x > 0 else off[-x - 1] for x in v])
        print(f"\n  首个可疑单元 #{ri} (n={m}) 坐标:")
        for i, p in enumerate(P):
            tag = "base" if i < m else "out"
            print(f"    {i} {tag}: {np.round(p, 3)}")
        print(f"  底面法向(叉积) = {np.round(np.cross(P[1]-P[0], P[2]-P[0]), 4)}")
        print(f"  外面向量      = {np.round(P[m]-P[0], 4)}")
        print(f"  dot(叉积, 外向) = {np.dot(np.cross(P[1]-P[0], P[2]-P[0]), P[m]-P[0]):+.4f}")
