"""软骨重做 v2：**保留原 `toward` 位移场，插值「位移」而非绝对坐标**。

── 为什么是 v2（v1 的假设被自己的实验证伪）────────────────────────────
v1 假设「28% 翻负 = 偏移方向用 toward 导致顶面扭转」，于是把方向改成骨面法向。
**实测证伪**：
  * 基线（toward，不细分）        6-tet 判据负体积 **0%**
  * v1（法向，K=1 不细分）        负体积 **38%**，FEBio 直接
    `Negative jacobian detected during mesh initialization`（0.2 s 就挂）
  * v1（法向，K=2 细分）          负体积 20%，rc=1
⇒ **法向偏移比 toward 差得多，v1 方向错了。**

── v2 的真机制 ─────────────────────────────────────────────────────────
`cart_refine.py` 的做法是：底面与偏移面**各自独立**做双线性细分。
问题在于偏移面是一个**扭曲的四边形**（4 个角被 toward 不规则地推开），
对扭曲四边形做双线性插值 ⇒ 内部点飞离真实位置 ⇒ 子单元扭转 ⇒ det<0。

**正解：插值「位移场」而不是插值绝对坐标**
    外层面(格点) = 底面(格点) + 插值后的位移
其中位移 d_i = old_offset_i − base_i （**原 toward 位移，逐角保留**）。
这样外层面的到达位置由「底面 + 光滑位移场」决定，位移场的变化幅度有限
⇒ 单元不可能比母单元扭得更多。而且沿共享棱两侧都是「同一棱的线性插值」
⇒ **天然共形**。

── 输出编码（与 `cart_refine.py` 兼容）────────────────────────────────
``{j}_nodes`` (n,3) / ``{j}_is_offset`` (n,) / ``{j}_elems`` (M,2m)
``{j}_tie_faces`` (F,m) 底面 tie primary / ``{j}_bone_faces`` (F,4)

── 开关 ────────────────────────────────────────────────────────────────
CH_K 面内倍数(quad->k×k, tri->1->4)   CH_N 层数   CH_OUT 输出名
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
THUMS = ROOT / "temp/thums"

K = int(os.environ.get("CH_K", "2") or 2)
N = int(os.environ.get("CH_N", "2") or 2)
OUTNPZ = os.environ.get("CH_OUT", "cartilage_mesh_hex.npz")

JOINTS = ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]

HEX_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
            (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
WEDGE_TETS = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]


def tet_vol(P):
    return float(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0])) / 6.0


def cell_vol(P, m):
    return sum(tet_vol(P[list(t)]) for t in (HEX_TETS if m == 4 else WEDGE_TETS))


src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bone_ids, bone_xyz = src["bone_node_ids"], src["bone_xyz"]
bpos = {int(v): bone_xyz[i] for i, v in enumerate(bone_ids)}

out: dict[str, np.ndarray] = {"bone_node_ids": bone_ids, "bone_xyz": bone_xyz}
rep: dict[str, dict] = {}

for j in JOINTS:
    elems = np.asarray(src[f"{j}_elems"], dtype=np.int64)
    old_off = np.asarray(src[f"{j}_offsets"], float)

    # ── 1. 解码：(base_refs, outer_refs) ──
    cells = []
    for row in elems:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        assert m in (3, 4), f"{j}: 只处理 tri3/quad4，实际 {len(v)} 节点"
        cells.append((v[:m], v[m:]))

    # ── 2. 节点表 ──
    nodes: list[np.ndarray] = []
    tag: list[bool] = []

    def add(p, is_off):
        nodes.append(np.asarray(p, float))
        tag.append(bool(is_off))
        return -(len(nodes))

    # 底面：复制成软骨专属节点（与骨解耦；共形交给 tied-elastic）
    bnodes = sorted({r for b, _ in cells for r in b})
    base_ref = {r: add(bpos[r], False) for r in bnodes}
    # 逐节点位移（同一骨节点被多面共用 → 取中位数，稳健）
    dacc: dict[int, list[np.ndarray]] = {r: [] for r in bnodes}
    tacc: dict[int, list[float]] = {r: [] for r in bnodes}
    for base, outer in cells:
        for r, o in zip(base, outer):
            d = old_off[-o - 1] - bpos[r]
            dacc[r].append(d)
            tacc[r].append(float(np.linalg.norm(d)))
    disp = {r: np.median(np.array(v), axis=0) for r, v in dacc.items() if v}
    for r in bnodes:
        disp.setdefault(r, np.zeros(3))

    def Pt0(r):
        """层 0（底面）位置。"""
        return nodes[-base_ref[r] - 1]

    # ── 3. 网格点：位置 = 底面的双线性插值 + (l/N)·位移场插值 ──
    emap: dict[tuple, int] = {}
    fmap: dict[tuple, int] = {}

    def edge_pt(r_from, r_to, l, i):
        """棱上第 i 个内分点（i=0..K），参数从 r_from 出发。"""
        if i == 0:
            return base_ref[r_from] if l == 0 else _lay(r_from, l)
        if i == K:
            return base_ref[r_to] if l == 0 else _lay(r_to, l)
        # ★ 规范化方向，保证相邻单元拿到同一个点
        if r_from < r_to:
            a, b, ii = r_from, r_to, i
        else:
            a, b, ii = r_to, r_from, K - i
        key = (a, b, l, ii)
        if key not in emap:
            t = ii / K
            p = (1 - t) * Pt0(a) + t * Pt0(b) + (l / N) * ((1 - t) * disp[a] + t * disp[b])
            emap[key] = add(p, l == N)
        return emap[key]

    def interior_pt(base, l, i, kk):
        key = (tuple(base), l, i, kk)
        if key not in fmap:
            u, v = i / K, kk / K
            w = np.array([(1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v])
            Pc = np.array([Pt0(base[0]), Pt0(base[1]), Pt0(base[2]), Pt0(base[3])])
            Dc = np.array([disp[base[0]], disp[base[1]], disp[base[2]], disp[base[3]]])
            p = w @ Pc + (l / N) * (w @ Dc)
            fmap[key] = add(p, l == N)
        return fmap[key]

    def qgrid(base, l, i, kk):
        if (i, kk) == (0, 0):
            return base_ref[base[0]] if l == 0 else _lay(base[0], l)
        if (i, kk) == (K, 0):
            return base_ref[base[1]] if l == 0 else _lay(base[1], l)
        if (i, kk) == (K, K):
            return base_ref[base[2]] if l == 0 else _lay(base[2], l)
        if (i, kk) == (0, K):
            return base_ref[base[3]] if l == 0 else _lay(base[3], l)
        if kk == 0:
            return edge_pt(base[0], base[1], l, i)
        if i == K:
            return edge_pt(base[1], base[2], l, kk)
        if kk == K:
            return edge_pt(base[3], base[2], l, K - i)
        if i == 0:
            return edge_pt(base[0], base[3], l, kk)
        return interior_pt(base, l, i, kk)

    # 「角点在各层的副本」：角点必须全局唯一（相邻单元共用）⇒ 按 (r,l) 缓存
    _layc: dict[tuple[int, int], int] = {}

    def _lay(r, l):
        key = (r, l)
        if key not in _layc:
            _layc[key] = add(Pt0(r) + (l / N) * disp[r], l == N)
        return _layc[key]

    # ── 4. 生成单元 ──
    new_cells: list[tuple[list[int], list[int]]] = []
    tie_faces: list[list[int]] = []
    bone_faces: list[list[int]] = []

    for base, outer in cells:
        if len(base) == 4:
            for a in range(K):
                for b in range(K):
                    cr = [(a, b), (a + 1, b), (a + 1, b + 1), (a, b + 1)]
                    col = [[qgrid(base, l, i, kk) for (i, kk) in cr]
                           for l in range(N + 1)]
                    tie_faces.append(list(col[0]))
                    for l in range(N):
                        new_cells.append((col[l], col[l + 1]))
            bone_faces.append(list(base))
        else:
            c0 = [_lay(base[0], l) if l else base_ref[base[0]] for l in range(N + 1)]
            c1 = [_lay(base[1], l) if l else base_ref[base[1]] for l in range(N + 1)]
            c2 = [_lay(base[2], l) if l else base_ref[base[2]] for l in range(N + 1)]
            if K == 1:
                # ★ K=1 表示「完全不细分」——三角也保持 1 块，才能与基线逐格可比
                subs = [(c0, c1, c2)]
            else:
                m01 = [edge_pt(base[0], base[1], l, 1) for l in range(N + 1)]
                m12 = [edge_pt(base[1], base[2], l, 1) for l in range(N + 1)]
                m20 = [edge_pt(base[2], base[0], l, 1) for l in range(N + 1)]
                subs = [(c0, m01, m20), (m01, c1, m12), (m20, m12, c2),
                        (m01, m12, m20)]       # ★ 末项 = 中央（LRN-051）
            for A, B, C in subs:
                tie_faces.append([A[0], B[0], C[0]])
                for l in range(N):
                    new_cells.append(([A[l], B[l], C[l]],
                                      [A[l + 1], B[l + 1], C[l + 1]]))
            bone_faces.append(list(base))

    # ── 5. 自检 ──
    def P_of(r):
        return nodes[-r - 1]

    nneg, vols = 0, []
    for blk, top in new_cells:
        P = np.array([P_of(r) for r in blk + top])
        v = cell_vol(P, len(blk))
        vols.append(v)
        if v <= 0:
            nneg += 1

    nq = sum(1 for b, _ in cells if len(b) == 4)
    nt = sum(1 for b, _ in cells if len(b) == 3)
    exp = (nq * K * K + nt * (4 if K > 1 else 1)) * N
    n_face = nq * K * K + nt * (4 if K > 1 else 1)
    tk = np.array([np.linalg.norm(Pt0(r) + disp[r] - Pt0(r)) for r in bnodes])

    rep[j] = dict(cells=len(new_cells), exp=exp, nneg=nneg, nface=n_face,
                  tmin=float(tk.min()), tmed=float(np.median(tk)),
                  tmax=float(tk.max()), nnode=len(nodes),
                  vol=float(sum(vols)))

    out[f"{j}_nodes"] = np.array(nodes, float)
    out[f"{j}_is_offset"] = np.array(tag, bool)
    E = np.zeros((len(new_cells), 8), np.int64)
    for i, (b_, o_) in enumerate(new_cells):
        E[i, :len(b_)] = b_
        E[i, len(b_):len(b_) + len(o_)] = o_
    out[f"{j}_elems"] = E
    Ft = max((len(f) for f in tie_faces), default=3)
    T = np.zeros((len(tie_faces), Ft), np.int64)
    for i, f in enumerate(tie_faces):
        T[i, :len(f)] = f
    out[f"{j}_tie_faces"] = T
    BF = np.zeros((len(bone_faces), 4), np.int64)
    for i, f in enumerate(bone_faces):
        BF[i, :len(f)] = f
    out[f"{j}_bone_faces"] = BF

    print(f"[{j:20s}] 单元 {len(new_cells):6d}(期望{exp:6d}) 负体积 {nneg:5d} "
          f" 厚度 {rep[j]['tmin']:.3f}/{rep[j]['tmed']:.3f}/{rep[j]['tmax']:.3f}"
          f"  节点 {len(nodes):6d}  底面格 {n_face:5d}")

np.savez_compressed(THUMS / OUTNPZ, **out)
print(f"\n写出 {OUTNPZ}  K={K} N={N}  "
      f"({(THUMS / OUTNPZ).stat().st_size / 1e6:.2f} MB)")
tn = sum(r["nneg"] for r in rep.values())
tc = sum(r["cells"] for r in rep.values())
tf = sum(r["nface"] for r in rep.values())
print(f"① 6-tet 判据负体积 : {tn:6d} / {tc} ({100*tn/tc:.1f}%)   "
      f"（基线 0%；v1 法向 K=1 时 38%）")
print(f"② 关节面底面格数   : {tf:6d}   （基线 1119）")
