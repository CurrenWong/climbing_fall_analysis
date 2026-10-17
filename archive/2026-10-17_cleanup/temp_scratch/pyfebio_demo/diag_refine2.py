"""逐类型统计：父单元 vs 子单元的 tet 体积与负体积分布。"""
import numpy as np
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\temp\thums\cartilage_mesh.npz")
d = np.load(P, allow_pickle=True)
bpos = {int(v): d["bone_xyz"][i] for i, v in enumerate(d["bone_node_ids"])}
FAN8 = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
        (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
FAN6 = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]


def tets(P8, n):
    return np.array([
        float(np.dot(np.cross(P8[q[1]] - P8[q[0]], P8[q[2]] - P8[q[0]]),
                     P8[q[3]] - P8[q[0]])) / 6.0
        for q in (FAN8 if n == 4 else FAN6)])


def refine(P8, n):
    b = [np.asarray(P8[i], float) for i in range(n)]
    o = [np.asarray(P8[n + i], float) for i in range(n)]
    m = [(b[i] + b[(i + 1) % n]) / 2 for i in range(n)]
    mo = [(o[i] + o[(i + 1) % n]) / 2 for i in range(n)]
    if n == 4:
        cb, co = np.mean(b, 0), np.mean(o, 0)
        return [[b[0], m[0], cb, m[3], o[0], mo[0], co, mo[3]],
                [m[0], b[1], m[1], cb, mo[0], o[1], mo[1], co],
                [cb, m[1], b[2], m[2], co, mo[1], o[2], mo[2]],
                [m[3], cb, m[2], b[3], mo[3], co, mo[2], o[3]]]
    return [[b[0], m[0], m[2], o[0], mo[0], mo[2]],
            [m[0], b[1], m[1], mo[0], o[1], mo[1]],
            [m[2], m[1], b[2], mo[2], mo[1], o[2]]]


for j in ("tibiotalar_tibia", "tibiotalar_talus", "subtalar_talus"):
    E = np.asarray(d[f"{j}_elems"], dtype=np.int64)
    off = np.asarray(d[f"{j}_offsets"], float)
    for n in (4, 3):
        rows = [[int(x) for x in r if int(x) != 0] for r in E]
        rows = [r for r in rows if len(r) // 2 == n]
        if not rows:
            continue
        pt_v, ch_v = [], []
        neg_c = 0
        nsub = 0
        for v in rows:
            P = np.array([bpos[x] if x > 0 else off[-x - 1] for x in v])
            pt_v.append(tets(P, n))
            for s in refine(P, n):
                t = tets(np.array(s), n)
                ch_v.append(t)
                nsub += 1
                if t.sum() <= 0:
                    neg_c += 1
        pv = np.abs(np.concatenate(pt_v))
        cv = np.abs(np.concatenate(ch_v))
        print(f"[{j}] {'hex8' if n == 4 else 'wedge6'}  父 {len(rows)} -> 子 {nsub}")
        print(f"   父 |V| min {pv.min():.5f}  中位 {np.median(pv):.5f}")
        print(f"   子 |V| min {cv.min():.5f}  中位 {np.median(cv):.5f}   "
              f"改善 {(pv.min()/max(cv.min(),1e-12)):.1f}x")
        print(f"   子单元总 tet 体积 <=0 的：{neg_c}/{nsub}")
        print()
