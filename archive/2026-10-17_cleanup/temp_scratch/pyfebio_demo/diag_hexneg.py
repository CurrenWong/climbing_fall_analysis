"""诊断：负体积单元长什么样？法向场有什么问题？"""
import numpy as np
from pathlib import Path
from collections import Counter

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
hexf = np.load(THUMS / "cartilage_mesh_hex.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}

WEDGE_TETS = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
HEX_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
            (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]


def tvol(P):
    return float(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0])) / 6.0


def cvol(P, m):
    return sum(tvol(P[list(t)]) for t in (HEX_TETS if m == 4 else WEDGE_TETS))


for j in ["tibiotalar_tibia", "subtalar_calcaneus"]:
    nd = hexf[f"{j}_nodes"]
    tg = hexf[f"{j}_is_offset"]
    E = hexf[f"{j}_elems"]

    def P_of(r):
        return nd[-r - 1] if r < 0 else bpos[r]

    # 重建底面面法向（按旧 offset 定向）
    old_off = src[f"{j}_offsets"]
    srcE = np.asarray(src[f"{j}_elems"], np.int64)
    fdir = {}          # 骨节点对 -> 旧偏移方向（用于定向）
    for row in srcE:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        for r, o in zip(v[:m], v[m:]):
            d = old_off[-o - 1] - bpos[r]
            if r not in fdir:
                fdir[r] = d

    neg, pos = [], []
    for ci, row in enumerate(E):
        blk = [int(x) for x in row[:8] if int(x) != 0]
        m = len(blk) // 2
        P = np.array([P_of(r) for r in blk])
        V = cvol(P, m)
        bot = newell = None
        b4 = P[:m]
        t4 = P[m:]
        nb = np.cross(b4[1] - b4[0], b4[2] - b4[0])
        nt = np.cross(t4[1] - t4[0], t4[2] - t4[0])
        nb /= (np.linalg.norm(nb) + 1e-15)
        nt /= (np.linalg.norm(nt) + 1e-15)
        rec = dict(ci=ci, V=V, dot=float(np.dot(nb, nt)),
                   area=float(np.linalg.norm(np.cross(b4[1] - b4[0], b4[2] - b4[0])) * .5),
                   h=float(np.mean([np.linalg.norm(t4[k] - b4[k]) for k in range(m)])),
                   edge=float(max(np.linalg.norm(b4[(k + 1) % m] - b4[k]) for k in range(m))))
        (neg if V <= 0 else pos).append(rec)

    print(f"\n{'='*78}\n[{j}]  负 {len(neg)} / 总 {len(E)}"
          f"  ({100*len(neg)/len(E):.1f}%)\n{'='*78}")
    for tag, arr in (("负", neg), ("正", pos)):
        if not arr:
            continue
        d = np.array([r["dot"] for r in arr])
        a = np.array([r["area"] for r in arr])
        h = np.array([r["h"] for r in arr])
        e = np.array([r["edge"] for r in arr])
        print(f"  {tag}: dot(底法向,顶法向) 中位={np.median(d):+.3f} "
              f"分位5%={np.percentile(d,5):+.3f} 95%={np.percentile(d,95):+.3f}")
        print(f"      底面面积 中位={np.median(a):9.4f} mm² "
              f"({np.sqrt(np.median(a)*2):.2f} mm 等效边长)  "
              f"高度 中位={np.median(h):.3f} mm  最长棱 中位={np.median(e):.3f}")
        print(f"      dot<0 的比例 = {100*np.mean(d<0):.1f}%")
    # γ 的真实分布
    g = [abs(r["V"]) / max(r["edge"], 1e-9) ** 3 for r in neg + pos]
    print(f"  γ(按最长底棱归一) 中位={np.median(g):.3f}  → 说明单元有多扁")
    print(f"  长宽比 等效边长/高度 中位 = "
          f"{np.median([np.sqrt(r['area']*2)/max(r['h'],1e-9) for r in neg+pos]):.2f}")
