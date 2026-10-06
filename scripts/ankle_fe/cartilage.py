"""B3-B5: 从 THUMS 骨面生成关节软骨层。

方法（贴合吴恺"在各骨面加厚生成软骨"）：
  1. 由 CORT ∪ SPON 的单元集提取**带朝外定向的**外表面（面出现 1 次 = 外表面；
     朝向由「面法向 vs 面心−单元心」判定，不依赖模板绕向）
  2. 面积加权求每个表面节点的**外法向**
  3. 对每个节点测「到对侧骨表面的最近距离」= 局部关节间隙 g
  4. 软骨厚度 t = min(t_target, fill_frac × g)   ← 用间隙做上限，**保证两侧不穿透**
  5. 沿法向 offset 生成外表面；原面 + offset 面 之间建 hex(四边形面) / wedge(三角面)
  6. 校验：体积、厚度分布、**与对侧骨的最小距离（必须 ≥0）**

★ 为什么必须测间隙而不是直接取 2 mm：
  实测踝关节骨间最小距离只有 1.406 mm。若两侧各取 2 mm，软骨会互相穿透 2.6 mm，
  FE 里表现为初始过盈 → 虚假的高接触应力。所以厚度必须由几何约束。
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import os as _env
NPZ = Path(_env.environ.get("THUMS_NPZ", "").strip()
           or r"D:\Project\climbing_fall_analysis\temp\thums\thums_lowerlimb.npz")
OUT = Path(_env.environ.get("CART_OUT", "").strip()
           or r"D:\Project\climbing_fall_analysis\temp\thums")
OUT.mkdir(parents=True, exist_ok=True)

TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
WEDGE_FACES = [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]

GROUPS = {
    "tibia":     ["81000700_solid", "81000600_solid", "81000601_solid"],
    "fibula":    ["81000900_solid", "81000800_solid", "81000801_solid"],
    "talus":     ["81001100_solid", "81001000_solid"],
    "calcaneus": ["81001300_solid", "81001200_solid"],
}
JOINTS = [("tibia", "talus", "tibiotalar_tibia"),
          ("talus", "tibia", "tibiotalar_talus"),
          ("talus", "calcaneus", "subtalar_talus"),
          ("calcaneus", "talus", "subtalar_calcaneus")]
# 双向都要生成：真实关节是「骨-软骨-软骨-骨」四层。
# 只做单侧会让关节一侧留空，接触应力分布不对称。

# 参数（可调；每个都写清来源）
# ★ T_TARGET 已 env 门控（2026-10-05）。动机：thickness = min(T_TARGET, 0.5×gap)
#   在 gap > 2×T_TARGET 处**封顶生效** ⇒ 两层合计 2×T_TARGET < gap ⇒ **留下间隙**。
#   实测中位间隙 0.9732mm（对应局部 gap≈3.9mm），而该间隙是「600 N 到不了」的直接原因
#   ——模型在 ~0.62mm 位移 / ~340 N 处撞天花板（DET_MIN=0.95 实测），
#   而间隙消失后位移天花板应随之消失。默认 1.5 = 逐字保持既有行为。
T_TARGET = float(_env.environ.get("T_TARGET", "").strip() or 1.5)
GAP_MAX = 6.0       # 关节面判定：节点距对侧骨 < 此值才算关节面
# ★ FILL_FRAC = 0.5: the two layers exactly fill the gap and are tangent,
#   NOT interpenetrating.
#   (0.6 / 0.52 were tried to make contact active at t=0; initial
#    interpenetration fires a contact force at t=0 and the first Newton step
#    crushes 4-5 elements. The real fix for rigid-body modes is clamping the
#    loaded cross-section -- see the disp driver in build_feb.py.)
FILL_FRAC = 0.5
# ---- mesh quality gates (switchable; default OFF) ----
# EXPERIMENT RESULT (2026-10-04, do not enable casually): these gates DO improve
#   tet quality a lot (|V| min 0.0013 -> 0.053 mm^3, 42x; penetration still 0),
#   but FEBio convergence gets WORSE: same 0.30 mm / 60-step case, progress
#   drops from 0.965 to 0.642; the known working point 0.1472 mm goes 17 s -> 60 s.
#   => the convergence bottleneck is the CONTACT, not sliver tets. Filtering faces
#     also punches holes / fragments into the contact surface, which makes the
#   penalty contact even harder to converge. Kept for the record; enable with
#   CART_QUALITY=1
import os as _os
_Q = _os.environ.get("CART_QUALITY", "0") == "1"
SMOOTH_ITERS = 4 if _Q else 0   # thickness-field smoothing passes (down-only)
T_MIN = 0.40 if _Q else 0.0     # mm: reject faces whose min thickness is below
T_RATIO = 2.5 if _Q else 1e9    # reject faces whose max/min thickness exceeds this
MIN_TET_V = 0.05 if _Q else 1e-9  # mm^3: drop whole cells with a sliver tet

# ---- cap each layer by the ALREADY-BUILT partner offset surface ----------
# Each side used t = 0.5 * (its own local gap) and both offset toward the
# *nearest* point of the opposite bone. Where the two sides are not
# geometrically symmetric, one layer pokes into the other at t=0; penalty
# contact then fires a huge force and crushes elements in step 1
# (reproduced: the 2-joint model dies with '2 negative jacobians' at ~zero
# load). Capping by the partner's offset surface removes it by construction.
# Order matters: the partner must be generated first (joint list is A,B,A,B).
PARTNER = {"tibiotalar_tibia": "tibiotalar_talus",
           "tibiotalar_talus": "tibiotalar_tibia",
           "subtalar_talus": "subtalar_calcaneus",
           "subtalar_calcaneus": "subtalar_talus"}
CAP_SAFETY = float(_os.environ.get("CAP_SAFETY", "0.97") or 0.97)
_OFFSET_PTS: dict[str, np.ndarray] = {}

d = np.load(NPZ)
nid = d["node_ids"]
xyz = d["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}


def oriented_boundary(keys: list[str]):
    """返回 {face_key: (verts tuple 朝外, normal, area, count)}，只保留外表面。"""
    acc: dict[tuple, dict] = {}
    for k in keys:
        arr = d[k]
        for row in arr:
            nds = [int(x) for x in row]
            uniq = sorted(set(nds))
            if len(uniq) == 4:
                tmpl = TET_FACES
            elif len(uniq) == 8:
                tmpl = HEX_FACES
            elif len(uniq) == 6:
                tmpl = WEDGE_FACES
            else:
                continue
            P = np.array([[xyz[idx[v]] for v in uniq]], float)[0]
            cen = P.mean(0)
            for f in tmpl:
                if max(f) >= len(uniq):
                    continue
                verts = [uniq[z] for z in f]
                pts = np.array([xyz[idx[v]] for v in verts], float)
                # Newell 法向
                n = np.zeros(3)
                m = len(pts)
                for i in range(m):
                    a, b = pts[i], pts[(i + 1) % m]
                    n += np.cross(a, b)
                fc = pts.mean(0)
                if np.dot(n, fc - cen) < 0:          # 朝内则翻转
                    n = -n
                    verts = verts[::-1]
                    # ★ 必须同步重排 pts！否则 f["pts"] 与 f["verts"] 顺序不一致，
                    #   下游 node_normals 用 pts 算出的法向会**反号** ——
                    #   这正是"软骨单元 Jacobian 全负"的根因。
                    pts = pts[::-1]
                key = tuple(sorted(verts))
                rec = acc.get(key)
                if rec is None:
                    acc[key] = {"verts": tuple(verts), "n": n, "cnt": 1,
                                "area": 0.0, "pts": pts}
                else:
                    rec["cnt"] += 1
    return [v for v in acc.values() if v["cnt"] == 1]


print("=== 提取带定向的外表面 ===")
surf = {}
for g, keys in GROUPS.items():
    faces = oriented_boundary(keys)
    surf[g] = faces
    nt = sum(1 for f in faces if len(f["verts"]) == 3)
    nq = sum(1 for f in faces if len(f["verts"]) == 4)
    print(f"  {g:10s} 外表面 {len(faces):>6,} 面 (tri {nt:>6,}/quad {nq:>6,})")


def node_normals(faces, nodes: set[int]) -> dict[int, np.ndarray]:
    """面积加权节点法向（Newell 法向模长 = 2×面积）。"""
    acc = {v: np.zeros(3) for v in nodes}
    for f in faces:
        pts = f["pts"]
        if len(pts) >= 3:
            a = pts[1] - pts[0]
            b = pts[2] - pts[0]
            fn = np.cross(a, b) * 0.5
            nrm = np.linalg.norm(fn)
            if nrm > 1e-12:
                unit = fn / nrm
            else:
                unit = f["n"] / max(np.linalg.norm(f["n"]), 1e-12)
        for v in f["verts"]:
            if v in acc:
                acc[v] = acc[v] + unit * max(nrm, 1e-12)
    out = {}
    for v, s in acc.items():
        ln = np.linalg.norm(s)
        out[v] = s / ln if ln > 1e-12 else np.array([0.0, 0.0, 1.0])
    return out


print("\n=== 逐关节生成软骨 ===")
cart = {}          # joint -> {"nodes": [...], "faces": [...], "elems": [...]}
report = {}
for A, B, jname in JOINTS:
    fA, fB = surf[A], surf[B]
    nA = set().union(*[set(f["verts"]) for f in fA])
    nB = set().union(*[set(f["verts"]) for f in fB])
    # ★ 必须由 sorted 后的 id 数组构造坐标，否则 PA 顺序与 idsA 不一致，
    #   gap 会配到错误的节点上（曾因此出现"关节面节点距对侧 146 mm"的荒谬值）
    idsA = np.array(sorted(nA), dtype=np.int64)
    idsB = np.array(sorted(nB), dtype=np.int64)
    PA = xyz[np.fromiter((idx[int(v)] for v in idsA), np.int64, len(idsA))]
    PB = xyz[np.fromiter((idx[int(v)] for v in idsB), np.int64, len(idsB))]
    treeB = cKDTree(PB)
    # ★ 只用「最近节点」估方向太粗糙（粗网格下最近节点常偏在侧面），
    #   会把关节面误判掉。改用「节点 + 面心」加密点云，方向估计稳得多。
    extra = np.array([f["pts"].mean(0) for f in fB])
    QB = np.vstack([PB, extra]) if len(extra) else PB
    treeQ = cKDTree(QB)
    gapA, nq = treeQ.query(PA, k=1)             # 到对侧**表面**的距离

    normals = node_normals(fA, nA)
    Narr = np.array([normals[int(v)] for v in idsA])
    dirA = QB[nq] - PA
    dl = np.linalg.norm(dirA, axis=1, keepdims=True)
    dirA = dirA / np.maximum(dl, 1e-12)
    align = np.einsum("ij,ij->i", Narr, dirA)   # 法向 · 指向对侧骨的单位向量

    ART_COS = np.cos(np.radians(60.0))
    # ---- 偏移方向 / 关节面筛选（2026-10-04 新增开关；默认保持旧行为）----
    # 探针实测（temp/pyfebio_demo/norm_probe.py，tibia 侧 73 个关节面节点）：
    #   · 合成模长 ≈ Σ面面积  ⇒ **顶点法向对齐良好，没有"相加抵消"**
    #     （旧注释说"严重不一致、几乎抵消"经实测为**误判**）
    #   · 但 `gap<6mm` 是**纯距离**判据 ⇒ 混入 36% 反向(p10 cos=-0.91)、
    #     55% 夹角>60° 的侧壁/边缘节点 ⇒ 沿法向挤会把软骨甩向侧面
    # ⇒ 正确做法（对齐 Anderson 2006）：**距离 + 法向对齐双重判据**，
    #   再沿**顶点法向**挤出 ⇒ 得到不歪斜的 hex8 层。
    CART_OFFSET = (_os.environ.get("CART_OFFSET", "toward") or "toward").strip()
    ART_ALIGN = float(_os.environ.get("ART_ALIGN", "0") or 0)
    art = {int(v): float(g) for v, g in zip(idsA, gapA) if g < GAP_MAX}
    if ART_ALIGN > 0:
        _al = {int(v): float(c) for v, c in zip(idsA, align) if c > ART_ALIGN}
        _n0 = len(art)
        art = {v: g for v, g in art.items() if v in _al}
        print(f"    法向对齐过滤 (cos>{ART_ALIGN:+.2f}): {_n0} -> {len(art)} 个关节面节点")
    print(f"    偏移方向 = {'顶点法向' if CART_OFFSET == 'normal' else '指向对侧'}"
          f"（CART_OFFSET={CART_OFFSET}）")
    aligned = int((align[gapA < GAP_MAX] > ART_COS).sum())
    print(f"\n  【{jname}】{A} -> {B}")
    print(f"    {A} 面节点 {len(nA):,}；关节面节点（距 {B} < {GAP_MAX} mm）= "
          f"{len(art):,}  其中法向对齐(>cos60°)仅 {aligned:,} 个"
          f"  <- 边缘脊线法向不良定义，故不过滤")
    if not art:
        print("    ⚠️ 无关节面节点，跳过")
        continue

    # ---- ① 厚度场平滑：只降不升，保证零穿透 ----
    #   为什么要平滑：t = FILL_FRAC×gap，而 gap 在关节面上从 0.16 mm 到 6 mm 剧烈变化，
    #   于是同一个面片 4 个节点的 t 能差 7 倍 → 六面体拆出的四面体里出**薄片单元**
    #   （实测 det(J) 最小 0.0048，体积仅 0.0008 mm³，比中位数小 3 个数量级），
    #   载荷下第一批翻负，是收敛包络上不去的主因之一。
    #   平滑后取 min(原值, 邻居均值) ⇒ 峰值被削平、谷值保持（谷处本来就薄），
    #   且**每个节点只会变薄** ⇒ 绝不产生穿透。
    adj: dict[int, set[int]] = {}
    for f in fA:
        vs = [v for v in f["verts"] if v in art]
        for a in vs:
            adj.setdefault(a, set()).update(b for b in vs if b != a)
    g = dict(art)
    for _ in range(SMOOTH_ITERS):
        nxt = {}
        for v, gv in g.items():
            nb = adj.get(v)
            nxt[v] = min(gv, float(np.mean([g[n] for n in nb]))) if nb else gv
        g = nxt

    # 厚度：min(t_target, FILL_FRAC × gap)
    thick = {v: min(T_TARGET, FILL_FRAC * gv) for v, gv in g.items()}
    # ★ 追加「对侧 offset 面」上限，避免初始过盈（见上方 PARTNER 注释）
    _pj = PARTNER.get(jname)
    if CAP_SAFETY > 0 and _pj in _OFFSET_PTS:
        from scipy.spatial import cKDTree as _KD
        _pts = np.array([xyz[idx[v]] for v in art])
        _d, _ = _KD(_OFFSET_PTS[_pj]).query(_pts)
        _ncap = 0
        for _i, v in enumerate(art):
            _lim = CAP_SAFETY * float(_d[_i])
            if _lim < thick[v]:
                thick[v] = _lim
                _ncap += 1
        print(f"    对侧({_pj})面限制: {_ncap}/{len(art)} 节点被压薄 "
              f"(到对侧 offset 面最近 {float(_d.min()):.3f} mm)")
    # ★ 偏移方向用「指向对侧骨表面」的单位向量（dirA），**不用面积加权节点法向**。
    #   实测：粗网格上相邻节点的面积加权法向严重不一致（相加几乎抵消），
    #   会让 offset 面几乎与 base 面重合（25 mm² 的面只移动 0.375 mm）→ 单元退化、
    #   Jacobian 全负。指向对侧骨的方向既稳健，也符合"软骨朝关节对侧生长"的物理。
    # ★ 偏移方向：默认「指向对侧骨表面」（旧行为）；CART_OFFSET=normal 时用
    #   **面积加权顶点法向** —— 实测该法向对齐良好（mag ≈ Σ面积），
    #   配合 ART_ALIGN 过滤掉侧壁/反向节点后，可得到不歪斜的 hex8 层。
    if CART_OFFSET == "normal":
        dir_by_id = {int(v): normals[int(v)] for v in art}
    else:
        dir_by_id = {int(v): dirA[i] for i, v in enumerate(idsA)}
    # 反向兜底：若某个节点的法向与指向对侧方向相反，翻正它（避免软骨长进骨头里）
    if CART_OFFSET == "normal":
        _rev = 0
        for i, v in enumerate(idsA):
            if int(v) in dir_by_id and float(np.dot(dir_by_id[int(v)], dirA[i])) < 0:
                dir_by_id[int(v)] = -dir_by_id[int(v)]
                _rev += 1
        if _rev:
            print(f"    法向反向翻转 {_rev} 个节点（原法向朝骨内）")
    off = {v: xyz[idx[v]] + dir_by_id[v] * thick[v] for v in art}
    tv = np.array(list(thick.values()))
    print(f"    厚度 t: min {tv.min():.3f} / 中位 {np.median(tv):.3f} / "
          f"max {tv.max():.3f} mm（上限 {T_TARGET} mm）")

    # 只用「整面都在关节面内」的表面面片建单元（保证贴合）
    # ★ 再加两道**逐面质量闸门**：太薄 / 面内厚度落差太大的面片直接不要。
    #   关节在某些位置本来就几乎闭合（gap 最小 0.16 mm），那里物理上放不下软骨层，
    #   硬生成只会得到必然翻负的薄片单元 —— 与其留个隐患，不如这段不建软骨
    #   （软骨是独立域，缺一小块不影响骨网格，也不影响接触面定义）。
    new_nodes: dict[int, int] = {}          # 原节点 -> 新 offset 节点索引（本地）
    local_pts: list[np.ndarray] = []
    base_faces: list[list[int]] = []
    n_layer = 0
    drop_thin = drop_ratio = 0
    for f in fA:
        vs = list(f["verts"])
        if not all(v in art for v in vs):
            continue
        tf = [thick[v] for v in vs]
        if min(tf) < T_MIN:                 # 层太薄 -> 必出薄片
            drop_thin += 1
            continue
        if max(tf) / min(tf) > T_RATIO:      # 面内落差太大 -> 棱台畸形
            drop_ratio += 1
            continue
        n_layer += 1
        for v in vs:
            if v not in new_nodes:
                new_nodes[v] = len(local_pts)
                local_pts.append(off[v])
        base_faces.append(vs)
    if drop_thin or drop_ratio:
        print(f"    质量闸门：剔除过薄面 {drop_thin}、厚度落差过大面 {drop_ratio} "
              f"（T_MIN={T_MIN} mm, T_RATIO={T_RATIO}）")

    if not base_faces:
        print("    ⚠️ 没有完整落在关节面内的面片，跳过")
        continue
    _OFFSET_PTS[jname] = np.asarray(local_pts, float)

    # ---- 组装：骨面节点沿用全局索引；offset 节点追加 ----
    n_existing = len(xyz)
    all_pts = np.asarray(list(xyz) + local_pts, float)
    off_index = {v: n_existing + new_nodes[v] for v in new_nodes}

    cells = []
    for vs in base_faces:
        # ★ 必须把 node ID 映射成数组下标（idx），否则 all_pts 索引越界
        cells.append([idx[v] for v in vs] + [off_index[v] for v in vs])

    # 用已验证的分解求体积，并保证为正（朝向不对就翻转）
    HEX_T = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
             (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
    WEDGE_T = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]

    def _tet6(a, b, c, e):
        return float(np.dot(np.cross(b - a, c - a), e - a)) / 6.0

    def vol_of(cell):
        pts = all_pts[cell]
        m = len(cell) // 2
        ts = HEX_T if m == 4 else WEDGE_T
        tot = 0.0
        for t in ts:
            tot += _tet6(pts[t[0]], pts[t[1]], pts[t[2]], pts[t[3]])
        return tot

    def tet_vols(cell):
        pts = all_pts[cell]
        m = len(cell) // 2
        ts = HEX_T if m == 4 else WEDGE_T
        return np.array([_tet6(pts[t[0]], pts[t[1]], pts[t[2]], pts[t[3]]) for t in ts])

    good = []
    skipped = bad_dir = sliver = 0
    alltv: list[np.ndarray] = []
    for cell in cells:
        # ★ 不要再"体积为负就翻转"——那样只反转 base 半边会破坏 base↔offset
        #   的一一对应，反而制造出负 Jacobian。
        #   base 面已是朝外 CCW、offset 沿外法向，构造出的单元天然正定向；
        #   体积非正的只可能是退化单元，直接剔除即可。
        tvol = tet_vols(cell)
        # 朝向交给 build_feb 修（它现在用的是**正确的奇置换**），
        # 这里只用整格**总体积**判退化；朝向反了不是网格缺陷。
        if tvol.sum() <= 1e-9:
            bad_dir += 1
            continue
        # ★ 逐四面体查薄片：整格体积正常 ≠ 每个四面体都正常。
        #   六面体拆 6 个四面体时，厚度落差大的面会产出极扁的四面体
        #   （实测 |V| 小到 0.0008 mm³，中位数的 1/2500），载荷下第一批翻负。
        #   一格里有一个薄片就**整格丢弃**（留半格会出内部自由面，更糟）。
        #   用绝对值判：朝向反的交给 build_feb 用**奇置换**修，不是网格缺陷。
        if np.abs(tvol).min() < MIN_TET_V:
            sliver += 1
            continue
        good.append(cell)
        alltv.append(tvol)

    flat = np.concatenate(alltv) if alltv else np.array([0.0])
    favol = np.abs(flat)
    print(f"    软骨单元 {len(good):,}（退化丢弃 {bad_dir}，含薄片整格丢弃 {sliver}）"
          f" | 体积 {flat.sum():,.1f} mm³ = {flat.sum()/1000:.3f} cm³")
    print(f"    四面体 |V| 分布: min {favol.min():.5f} / p1 {np.percentile(favol,1):.4f} "
          f"/ 中位 {np.median(favol):.4f} mm³   "
          f"(MIN_TET_V={MIN_TET_V}, 需 build_feb 翻正的 {int((flat<0).sum()):,})")

    # ★ 编码成「可直接与全局合并」的形式：
    #   骨节点 -> 正的 node ID；offset 节点 -> -(k+1)，k 为该关节 offsets 的下标
    #   （各关节的 offset 基址必须分开，否则全局编号会撞车）
    conv = []
    for cell in good:
        m = len(cell) // 2
        row = []
        for e in cell:
            if e >= n_existing:
                row.append(-(e - n_existing + 1))
            else:
                row.append(int(nid[e]))
        conv.append(row)
    cart[jname] = {"elems": conv,
                   "offsets": np.asarray(local_pts, float),
                   "n_new": len(local_pts)}

    # ---- 校验：offset 后是否穿透对侧骨 ----
    newP = all_pts[n_existing:]
    dd, _ = treeB.query(newP, k=1) if len(newP) else (np.array([]), None)
    report[jname] = {
        "A": A, "B": B,
        "n_surface_nodes_A": len(nA),
        "n_articular_nodes": len(art),
        "n_aligned_nodes": aligned,
        "normal_alignment_note": "边缘脊线节点法向不良定义；Newell 与三点法可能符号相反",
        "n_cartilage_elements": len(good),
        "n_new_nodes": len(local_pts),
        "volume_mm3": float(flat.sum()),
        "thickness_mm": {"min": float(tv.min()), "median": float(np.median(tv)),
                         "max": float(tv.max()), "target": T_TARGET},
        "penetration_check": {
            "min_dist_offset_to_B_mm": float(dd.min()) if len(dd) else None,
            "n_penetrating": int((dd < 0).sum()) if len(dd) else None,
            "n_below_0p05mm": int((dd < 0.05).sum()) if len(dd) else None,
        },
        "params": {"T_TARGET": T_TARGET, "GAP_MAX": GAP_MAX,
                   "FILL_FRAC": FILL_FRAC},
    }

print("\n=== 校验汇总 ===")
print(f"{'joint':12s} {'elem':>7s} {'vol(cm³)':>10s} {'t_min':>7s} {'t_med':>7s} "
      f"{'t_max':>7s} {'offset↔B min':>13s} {'穿透过':>7s}")
for j, r in report.items():
    pc = r["penetration_check"]
    print(f"{j:12s} {r['n_cartilage_elements']:>7,} "
          f"{r['volume_mm3']/1000:>10.3f} {r['thickness_mm']['min']:>7.3f} "
          f"{r['thickness_mm']['median']:>7.3f} {r['thickness_mm']['max']:>7.3f} "
          f"{(pc['min_dist_offset_to_B_mm'] if pc['min_dist_offset_to_B_mm'] is not None else float('nan')):>13.3f} "
          f"{pc['n_penetrating']:>7}")

(OUT / "cartilage.json").write_text(
    json.dumps({"params": {"T_TARGET": T_TARGET, "GAP_MAX": GAP_MAX,
                           "FILL_FRAC": FILL_FRAC},
                "joints": report}, ensure_ascii=False, indent=2), encoding="utf-8")


def pad_elems(elems):
    """六面体(8)与楔(6)混存 -> 补齐到同宽。
    ★ 填充值用 0，不能用 -1 —— 因为 offset 节点编码就是 -(k+1)，
      第 0 个 offset 恰好也是 -1，会与填充混淆。节点 id 恒为正，故 0 无歧义。
    """
    if not elems:
        return np.zeros((0, 8), np.int64)
    w = max(len(e) for e in elems)
    a = np.zeros((len(elems), w), dtype=np.int64)
    for i, e in enumerate(elems):
        a[i, : len(e)] = e
    return a


# 落盘网格供 A 用
if cart:
    arrays = {"bone_node_ids": nid, "bone_xyz": xyz}
    for j, v in cart.items():
        arrays[f"{j}_offsets"] = v["offsets"]
        arrays[f"{j}_elems"] = pad_elems(v["elems"])
    np.savez_compressed(OUT / "cartilage_mesh.npz", **arrays)
    print(f"\n写出 {OUT/'cartilage.json'} 与 {OUT/'cartilage_mesh.npz'}")
else:
    print("\n⚠️ 未生成任何软骨")

# ── GAP_MAX 敏感性扫描：关节面大小会不会被阈值截断？ ──
print("\n=== GAP_MAX 敏感性（关节面规模 vs 阈值）===")
print(f"{'joint':12s} {'bone':10s} " +
      " ".join(f"{f'<{g}mm':>10s}" for g in (3, 4, 6, 8, 10, 15, 20)))
for A, B, jname in JOINTS:
    nA = set().union(*[set(f["verts"]) for f in surf[A]])
    idsA = np.array(sorted(nA), dtype=np.int64)
    PA = xyz[np.fromiter((idx[int(v)] for v in idsA), np.int64, len(idsA))]
    nB = set().union(*[set(f["verts"]) for f in surf[B]])
    idsB = np.array(sorted(nB), dtype=np.int64)
    PB = xyz[np.fromiter((idx[int(v)] for v in idsB), np.int64, len(idsB))]
    g, _ = cKDTree(PB).query(PA, k=1)
    counts = [int((g < t).sum()) for t in (3, 4, 6, 8, 10, 15, 20)]
    print(f"{jname:12s} {A:10s} " + " ".join(f"{c:>10,}" for c in counts))

