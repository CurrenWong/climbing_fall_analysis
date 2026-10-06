"""A1-A4: 把 THUMS 骨 + 生成的软骨组装成 FEBio 模型并求解（冒烟）。

单位：mm - N - MPa - tonne（与 THUMS 一致；E 已是 MPa，零换算）
材料：皮质/松质/骨髓取自 THUMS 材料卡；**软骨 THUMS 没有** -> 用吴恺 2012 的 E=260, ν=0.4

阶段化：先跑「胫骨 + 距骨 + 双侧软骨」的 tibiotalar 单关节（接触对最少，最易收敛），
       跑通再上 fibula / calcaneus。
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# 二分定位用：BONES_ONLY=1 时只跑骨头（无软骨、无接触、用指定位移加载）
BONES_ONLY = os.environ.get("BONES_ONLY") == "1"
# 驱动方式：force=600N 面力（默认）；disp=顶面全约束+指定位移
# ★ 为什么需要 disp：胫骨只靠接触与距骨相连，接触未咬合时胫骨存在**刚体模态**
#   → 刚度矩阵奇异 → 位移爆掉 → 单元大批翻转（实测 2200 个）。
#   disp 模式把胫骨顶面三个自由度全部约束住（=夹持端），刚体模态彻底消失，
#   不再依赖接触咬合来稳定，是标准的"压缩试验"配置。
DRIVER = os.environ.get("DRIVER", "force").strip() or "force"
DISP_MM = float(os.environ.get("DISP_MM", "0.05") or 0.05)
# 载荷步数：大变形/接触问题收敛困难时，细化载荷步通常最有效
TIME_STEPS = int(os.environ.get("TIME_STEPS", "20") or 20)
# contact knobs: PENALTY 收敛差, AUGLAG 是标准解法; node_reloc 把接触节点投回从面
# （软骨被压溃时很关键）。实测收敛瓶颈在接触，不在软骨网格质量。
LAUGON = (os.environ.get("LAUGON", "PENALTY") or "PENALTY").upper()
NODE_RELOC = os.environ.get("NODE_RELOC", "0") == "1"
TWO_PASS = os.environ.get("TWO_PASS", "0") == "1"
PENALTY = float(os.environ.get("PENALTY", "1") or 1)

OUT = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
WORK = OUT / "febio"
WORK.mkdir(parents=True, exist_ok=True)

# ── 参与的骨部件（pid_solid 在 npz 里的键）与材料 ──
BONE_PARTS = {
    "tibia_cort": ("81000700_solid", "mat_cort_tibia"),
    "tibia_end_spon": ("81000600_solid", "mat_spon_long"),
    "tibia_marrow": ("81000601_solid", "mat_marrow"),
    "talus_cort": ("81001100_solid", "mat_cort_tarsal"),
    "talus_spon": ("81001000_solid", "mat_spon_tarsal"),
}
# FULL_ANKLE=1 -> add the subtalar joint: calcaneus bone + both cartilage layers
# + a second contact pair; the fixed face moves from the talus bottom to the
# calcaneus bottom (anatomically the most distal fix of the ankle complex).
FULL_ANKLE = os.environ.get("FULL_ANKLE", "0") == "1" and not BONES_ONLY
if FULL_ANKLE:
    BONE_PARTS["calcaneus_cort"] = ("81001300_solid", "mat_cort_tarsal")
    BONE_PARTS["calcaneus_spon"] = ("81001200_solid", "mat_spon_tarsal")
if BONES_ONLY:
    PAIRS = []
else:
    PAIRS = [("tibiotalar_pair", "tibiotalar_tibia", "tibiotalar_talus")]
    if FULL_ANKLE:
        PAIRS.append(("subtalar_pair", "subtalar_talus", "subtalar_calcaneus"))
CARTS = [j for _, a, b in PAIRS for j in (a, b)]
FIX_BONE = "calcaneus_cort" if FULL_ANKLE else "talus_cort"
_mode = "BONES_ONLY" if BONES_ONLY else (
    "full-ankle(2 joints)" if FULL_ANKLE else "tibiotalar-only")
print(f"[mode] {_mode}; pairs={[p[0] for p in PAIRS]}; fixed on {FIX_BONE}")
if BONES_ONLY:
    W = None  # 标记：无软骨
print(f"[mode] {'BONES_ONLY' if BONES_ONLY else 'full(含软骨+接触)'}")

MATERIALS = {
    # 名称: (type, 参数)  —— 来源见 C 阶段；软骨来自吴恺 2012 表 1
    "mat_cort_tibia": ("isotropic elastic",
                       {"density": 2.0e-9, "E": 18000.0, "v": 0.30}),
    "mat_cort_tarsal": ("isotropic elastic",
                        {"density": 2.0e-9, "E": 15000.0, "v": 0.30}),
    "mat_spon_long": ("isotropic elastic",
                      {"density": 8.615e-10, "E": 160.0, "v": 0.45}),
    "mat_spon_tarsal": ("isotropic elastic",
                        {"density": 1.0e-9, "E": 73.4, "v": 0.45}),
    "mat_marrow": ("isotropic elastic",
                   {"density": 1.0e-9, "E": 12.0, "v": 0.499}),
    "mat_cartilage": ("isotropic elastic",
                      {"density": 1.1e-9, "E": 260.0, "v": 0.40}),
}

FACES = {4: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
         8: [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
         6: [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]}
ELTYPE = {4: "tet4", 8: "hex8", 6: "penta6"}
NFACE = {4: "tri3", 8: "quad4", 6: "quad4"}
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
MIN_FACET_AREA = 1e-8   # mm^2: facets smaller than this are degenerate


def facet_ok(f):
    """Facet must have 3 distinct nodes and non-zero area.
    FEBio warns 'The surface X has N invalid facets' for these, and a penalty
    contact built on them blows up at t=0."""
    if len(set(f)) != len(f):
        return False
    p = XYZ[[i - 1 for i in f]]
    return float(np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))) * 0.5 \
        > MIN_FACET_AREA

# 软骨：按「半边节点数」映射 —— 3 节点底 -> 楔(penta6)，4 节点底 -> 六面体(hex8)
CART_ELTYPE = {3: "penta6", 4: "hex8"}

# ── 1. 读取并建立紧凑节点编号 ──
db = np.load(OUT / "thums_lowerlimb.npz")
nid_all, xyz_all = db["node_ids"], db["node_xyz"]
idx_all = {int(v): i for i, v in enumerate(nid_all)}
cg = np.load(OUT / os.environ.get("CART_MESH", "cartilage_mesh.npz"))

mesh_elems: list[tuple[str, str, str, list[int]]] = []   # (domain, etype, mat, conn)
used: set[int] = set()

WEDGE_T = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]   # 楔(6节点) -> 3 四面体
HEX_T = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
         (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]      # 六面体 -> 6 四面体
_CACHE_XYZ = None


def _vol(ids):
    P = XYZ[[i - 1 for i in ids]]
    return float(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0])) / 6.0


for dom, (key, mat) in BONE_PARTS.items():
    arr = db[key]
    for row in arr:
        # ★ 两步都不能少：
        #   ① 去掉 -1 / 0（npz 填充）
        #   ② **按顺序去重** —— 松质/骨髓的四面体是以「8 槽位退化 hex」形式存的
        #      （如 a,b,c,d,d,d,d,d），只过滤填充不去重会把它们整批丢掉
        conn = []
        for x in row:
            xi = int(x)
            if xi not in (0, -1) and xi not in conn:
                conn.append(xi)
        u = len(conn)
        if u == 8:
            mesh_elems.append((dom, "hex8", mat, conn))
            used.update(conn)
        elif u == 4:
            mesh_elems.append((dom, "tet4", mat, conn))
            used.update(conn)
        elif u == 6:
            # ★ 楔不写成 penta6，而是拆成 3 个 tet4。
            #   实测：hex8 混进 penta6 后 FEBio 报 Negative jacobian，
            #   而同样这批 hex8 单独跑 rc=0 —— 嫌疑就是 penta6 的朝向/退化。
            #   tet4 的朝向由体积符号唯一确定，最稳。
            for t in WEDGE_T:
                tet = [conn[z] for z in t]
                mesh_elems.append((dom, "tet4", mat, tet))
                used.update(tet)

for j in CARTS:
    arr = cg[f"{j}_elems"]
    for row in arr:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals:
            continue
        m = len(vals) // 2
        # ★ 负引用一律指向 {j}_nodes（旧文件叫 _offsets，语义相同）。
        #   新格式下底面节点**也在** _nodes 里（tied-elastic 需要底面脱离骨节点），
        #   所以底面**不能**再靠"前 m 个是正数"来判定，正负号本身就是唯一标记。
        #   ★ 也因此不能保留 `vals[0] <= 0: continue` —— 复制后的底面首值必为负，
        #   留着这行会把**整层软骨静默丢掉**（曾经就这样跑出一个"复现基线"的假通过）。
        full = [x if x > 0 else f"off:{j}:{-x - 1}" for x in vals]
        base, offs_ = full[:m], full[m:]
        # ★ 软骨**全部拆成 tet4**：hex8/penta6 在薄层 + 粗网格下极易出现
        #   det(J) ≈ 0.05 的近退化单元，FEBio 会在首个载荷步就报 Negative jacobian。
        #   tet4 的朝向由体积符号唯一确定，最稳。
        for t in (HEX_T if m == 4 else WEDGE_T):
            mesh_elems.append((j, "tet4", "mat_cartilage", [full[z] for z in t]))
        # 只有仍挂在骨上的节点进 used（紧凑编号用）
        used.update(b for b in base if not isinstance(b, str))

print(f"单元总数 {len(mesh_elems):,}；用到的骨节点 {len(used):,}")
from collections import Counter
print("按域统计:", dict(Counter(d for d, _, _, _ in mesh_elems)))
print("按类型统计:", dict(Counter(t for _, t, _, _ in mesh_elems)))

# 紧凑编号：骨节点 id -> 1..N；软骨专属节点紧随其后
# ★ 新格式用 {j}_nodes（旧的 cartilage_mesh.npz 里叫 {j}_offsets，语义相同）；
#   {j}_is_offset 标出哪些落在**外表面**——底面节点现在也在 _nodes 里，
#   若不筛选，cart_outer() 会把底面也当成接触面。
def _cname(j):
    return f"{j}_nodes" if f"{j}_nodes" in cg.files else f"{j}_offsets"


def _cflag(j, k):
    if f"{j}_is_offset" not in cg.files:
        return True                     # 旧文件：全是 offset
    return bool(cg[f"{j}_is_offset"][k])


bone_ids = sorted(used)
node_id_of: dict[int, int] = {v: i + 1 for i, v in enumerate(bone_ids)}
coords: list[tuple[float, float, float]] = [tuple(xyz_all[idx_all[v]]) for v in bone_ids]
for j in CARTS:
    for k, p in enumerate(cg[_cname(j)]):
        node_id_of[f"off:{j}:{k}"] = len(coords) + 1
        coords.append(tuple(p))
print(f"合并后节点 {len(coords):,}")
# offset 节点（软骨外表面节点）的全局编号集合 —— 用于识别软骨外表面面片
OFFSET_IDS = {node_id_of[f"off:{j}:{k}"]
              for j in CARTS for k in range(len(cg[_cname(j)]))
              if _cflag(j, k)}
print(f"其中 offset（软骨外表面）节点 {len(OFFSET_IDS):,}")


def _ref2id(j, x):
    """npz 引用 -> 全局节点号。>0 骨 ID；<0 -> _nodes[-x-1]。"""
    return node_id_of[int(x)] if int(x) > 0 else node_id_of[f"off:{j}:{-int(x) - 1}"]

XYZ = np.asarray(coords, float)


def nid_of(x):
    return node_id_of[x] if not isinstance(x, str) else node_id_of[x]


# ── 2. 归一化连接表 ──
enum = 1
elem_blocks: dict[tuple[str, str], list[tuple[int, list[int]]]] = {}
for dom, etype, mat, conn in mesh_elems:
    ids = [node_id_of[c] for c in conn]
    elem_blocks.setdefault((dom, etype), []).append((enum, ids))
    enum += 1
print(f"写出了 {len(elem_blocks)} 个单元块")

# ── 3. 求主方向（用胫骨→距骨的连线作加载轴），并找加载面 / 固定面 ──
def domain_nodes(dom):
    s = set()
    for (d, _), v in elem_blocks.items():
        if d == dom:
            for _, ids in v:
                s.update(ids)
    return np.array(sorted(s))


tib = XYZ[domain_nodes("tibia_cort") - 1]
tal = XYZ[domain_nodes("talus_cort") - 1]
axis = tib.mean(0) - tal.mean(0)
axis = axis / np.linalg.norm(axis)
print(f"加载轴（距骨→胫骨）= {np.round(axis, 4)}")


def surface_faces(dom, etype):
    """返回该域外表面的 (faces, 朝向法向)"""
    cnt = {}
    info = {}
    for _, ids in elem_blocks[(dom, etype)]:
        P = XYZ[[i - 1 for i in ids]]
        cen = P.mean(0)
        for f in FACES[len(ids)]:
            fv = [ids[z] for z in f]
            pts = XYZ[[i - 1 for i in fv]]
            n = np.zeros(3)
            for i in range(len(pts)):
                n += np.cross(pts[i], pts[(i + 1) % len(pts)])
            if np.dot(n, pts.mean(0) - cen) < 0:
                n = -n
                fv = fv[::-1]
            key = tuple(sorted(fv))
            if not facet_ok(fv):
                continue
            info[key] = (tuple(fv), n)
            cnt[key] = cnt.get(key, 0) + 1
    return [(info[k][0], info[k][1]) for k, c in cnt.items() if c == 1]


tib_faces = []
for et in ("hex8", "tet4", "penta6"):
    if ("tibia_cort", et) in elem_blocks:
        tib_faces += surface_faces("tibia_cort", et)

proj = XYZ[[f[0] - 1 for f, _ in tib_faces]] @ axis
hi = proj.max()
top_sel = [f for f, p in zip([f for f, _ in tib_faces], proj) if p > hi - 3.0]
print(f"胫骨外表面 {len(tib_faces):,} 面；近端加载面 {len(top_sel)} 面 (>95% 高度 -3mm)")

fix_faces = []
for et in ("hex8", "tet4", "penta6"):
    if (FIX_BONE, et) in elem_blocks:
        fix_faces += surface_faces(FIX_BONE, et)
projT = XYZ[[f[0] - 1 for f, _ in fix_faces]] @ axis
lo = projT.min()
bot_sel = [f for f, p in zip([f for f, _ in fix_faces], projT) if p < lo + 3.0]
print(f"{FIX_BONE} 外表面 {len(fix_faces):,} 面；远端固定面 {len(bot_sel)} 面")

# ── 4. 软骨外层接触面（每个软骨单元的最后一半节点即 offset 面） ──
def cart_outer(j):
    """软骨外层接触面 —— ★ 必须从 **tet 单元的真实外表面** 取，不能从 npz 的原始四边形取！

    踩过的坑：软骨单元已拆成 tet4，但接触面若仍用原始 offset 四边形
    （其三角化方式与 tet 分解不一致），面就**与网格实际的面不重合**，
    FEBio 的接触搜索关联不上 ⇒ 接触完全不检测（实测：压缩 0.3 mm 仍零应力，
    两层软骨自由互穿）。正确做法是取「三个节点全为 offset 节点」的 tet 面。
    """
    out, seen = [], set()
    for _, ids in elem_blocks[(j, "tet4")]:
        for f in TET_FACES:
            fv = [ids[z] for z in f]
            if not all(x in OFFSET_IDS for x in fv):
                continue
            key = tuple(sorted(fv))
            if key in seen:
                continue
            if not facet_ok(fv):
                continue
            seen.add(key)
            # 朝外定向：法向应背离该 tet 的对顶点
            opp = ids[[z for z in range(4) if z not in f][0]]
            pts = XYZ[[i - 1 for i in fv]]
            n = np.zeros(3)
            for i in range(3):
                n += np.cross(pts[i], pts[(i + 1) % 3])
            if np.dot(n, pts.mean(0) - XYZ[opp - 1]) < 0:
                fv = fv[::-1]
            out.append(fv)
    return out


contact_surfs = {j: cart_outer(j) for j in CARTS}
for j, v in contact_surfs.items():
    print(f"  {j} 外层接触面 {len(v)} 个")

# ── 4.6 tied-elastic 绑定面（骨 ↔ 软骨底面）────────────────────────────
# 手册 §3.14.5/§3.14.6：tied 用来「连接两个非共形网格」，约束
# primary 的**节点**连到 secondary 的**面**上；tied-elastic 在此基础上
# 强制界面两侧**位移连续**（solid-solid）。这里 primary = 脱开的软骨底面，
# secondary = 粗骨关节面。目的：软骨可独立加密，骨保持已验证的粗网格。
TIE = os.environ.get("TIE", "0") == "1"
TIE_SURF: dict[str, dict[str, list]] = {}
tie_pairs: list[tuple[str, str, str]] = []     # (name, primary, secondary)
if TIE:
    def _faces(key, jj):
        raw = np.asarray(cg[f"{jj}_{key}"], dtype=np.int64)
        return [[int(v) for v in row if int(v) != 0] for row in raw]

    # ★★ 面必须是**单元的真实 facet**，否则 FEBio 报 "N invalid facets" 并静默丢弃
    #   （User Manual §3.6.5 的老坑，cart_outer 踩过一次）。
    #   软骨已全部拆成 tet4 ⇒ facet 只能是**三角形**；npz 里底面却是 quad/tri 混合，
    #   直接照抄会丢掉全部 quad（实测 93 面里 61 个 invalid = 正好是 hex8 的个数）。
    #   所以两边都从**已装配的 elem_blocks** 反查真实 facet，不信任 npz 的四边形。
    FACES_ALL = {"tet4": TET_FACES, "hex8": FACES, "penta6": None}

    def _real_facets(doms):
        """域 -> {sorted(node tuple): 定向面}，直接扫单元。"""
        out = {}
        for (dom, et), v in elem_blocks.items():
            if dom not in doms or et not in ("tet4", "hex8"):
                continue
            for _, ids in v:
                for f in (TET_FACES if et == "tet4" else FACES[len(ids)]):
                    fv = [ids[z] for z in f]
                    out.setdefault(tuple(sorted(fv)), fv)
        return out

    bone_doms = set(BONE_PARTS)
    bone_fac = _real_facets(bone_doms)

    for j in CARTS:
        if f"{j}_tie_faces" not in cg.files or f"{j}_bone_faces" not in cg.files:
            print(f"  [tie] {j}: npz 无 _tie_faces/_bone_faces，跳过")
            continue
        # ---- primary：软骨 tet 的底面（三节点全属于该关节的底面节点）----
        bid = {node_id_of[f"off:{j}:{k}"]
               for k in range(len(cg[_cname(j)])) if not _cflag(j, k)}
        prim, seen = [], set()
        for _, ids in elem_blocks.get((j, "tet4"), []):
            for f in TET_FACES:
                fv = [ids[z] for z in f]
                if not all(x in bid for x in fv):
                    continue
                k = tuple(sorted(fv))
                if k in seen or not facet_ok(fv):
                    continue
                seen.add(k)
                opp = ids[[z for z in range(4) if z not in f][0]]
                pts = XYZ[[i - 1 for i in fv]]
                n = np.zeros(3)
                for i in range(3):
                    n += np.cross(pts[i], pts[(i + 1) % 3])
                if np.dot(n, pts.mean(0) - XYZ[opp - 1]) < 0:
                    fv = fv[::-1]
                prim.append(fv)
        # ---- secondary：该关节「A 骨」外表面中**靠近软骨底面**的局部片 ----
        # ★ 用整块骨外表面太贵（4850 面 -> 395 s）且易误抓远处 facet；
        #   只挑"与底面重合的 facet"又会缺（talus 只匹配 166/200）⇒ 空洞
        #   ⇒ 漏掉的 primary 节点被 search_radius 拉到远处锁死 ⇒ 虚假刚度。
        #   折中：取骨外表面中**质心距任一软骨底面节点 < TOL** 的那些面，
        #   既完整覆盖关节面，又不会引入远处面。
        from scipy.spatial import cKDTree as _KD
        TOL = float(os.environ.get("TIE_SEC_TOL", "0.05"))
        bp_ids = {node_id_of[f"off:{j}:{k}"]
                  for k in range(len(cg[_cname(j)])) if not _cflag(j, k)}
        btree = _KD(XYZ[[i - 1 for i in bp_ids]])
        A_bone = {"tibiotalar_tibia": "tibia_cort",
                  "tibiotalar_talus": "talus_cort",
                  "subtalar_talus": "talus_cort",
                  "subtalar_calcaneus": "calcaneus_cort"}[j]
        sec = []
        for et in ("hex8", "tet4", "penta6"):
            if (A_bone, et) not in elem_blocks:
                continue
            for f, _n in surface_faces(A_bone, et):
                fc = XYZ[[i - 1 for i in f]].mean(0)
                if float(np.linalg.norm(btree.query(fc)[0])) <= 0.0:
                    continue
                d, _ = btree.query(XYZ[[i - 1 for i in f]])
                if float(d.min()) <= TOL:
                    sec.append(f)
        if not prim or not sec:
            print(f"  [tie] {j}: primary {len(prim)} / secondary {len(sec)}，跳过")
            continue
        TIE_SURF[j] = {"primary": prim, "secondary": sec}
        tie_pairs.append((f"tie_{j}", f"tiesurf_{j}", f"bonesurf_{j}"))
        print(f"  [tie] {j}: primary(软骨底面 tri) {len(prim)} 面 / "
              f"secondary(骨关节面, 坐标匹配) {len(sec)} 面")

# ── 5.5 ★ 写入前做 Jacobian 自检（FEBio 会在初始化时拒绝 det(J)<=0 的单元）──
#   不猜 FEBio 的朝向约定，直接按 isoparametric 公式算 det(J) 定位问题域。
GP2 = 1.0 / math.sqrt(3.0)


def detJ_hex8(P):
    """P:(8,3) 在 8 个 Gauss 点算 det(J)，返回最小值。"""
    xi = np.array([[-GP2, -GP2, -GP2], [GP2, -GP2, -GP2], [GP2, GP2, -GP2],
                   [-GP2, GP2, -GP2], [-GP2, -GP2, GP2], [GP2, -GP2, GP2],
                   [GP2, GP2, GP2], [-GP2, GP2, GP2]])
    ref = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    worst = 1e30
    for s in xi:
        dN = np.zeros((8, 3))
        for i in range(8):
            a, b, c = ref[i]
            dN[i, 0] = a * (1 + b * s[1]) * (1 + c * s[2]) / 8
            dN[i, 1] = b * (1 + a * s[0]) * (1 + c * s[2]) / 8
            dN[i, 2] = c * (1 + a * s[0]) * (1 + b * s[1]) / 8
        worst = min(worst, float(np.linalg.det(P.T @ dN)))
    return worst


def detJ_tet4(P):
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0],
                                         P[3] - P[0]])))


def detJ_penta6(P):
    """P:(6,3) 三角面 0-2 与 3-5；返回**采样集上最坏**的 det(J)。

    ⚠️ 2026-10-06 修 `[LRN-20261006-085]`/`[LRN-20261006-086]`：
    原实现只采 **6 个 Gauss 点**（三角 3 点 × z=±1/√3），**不采顶点/棱中点**
    ⇒ 漏掉「顶点附近局部翻转」的楔形单元。实测软骨：
        顶点 worst-det(J) = **−7.89**（9+3 个单元 ≤1e-2，11 个 ≤0）
        而 Gauss 点只报 min = 2.39e-3 ⇒ **静默通过**，坏单元留进模型
        ⇒ 局部应力垃圾（σ_max 尖峰来源之一）。

    判据 = **采样点 + 阈值 + 方向**，三者缺一不可（`[LRN-20261006-085]`）。
    现采：三角 3 Gauss 点 + 3 顶点 + 3 棱中点 + 1 面心，z 取
    ±1/√3（Gauss）/ ±1（上下端面）/ 0（中面）。
    """
    tri = ([(1 / 6, 1 / 6), (2 / 3, 1 / 6), (1 / 6, 2 / 3)]      # Gauss
           + [(1.0, 0.0), (0.0, 1.0), (0.0, 0.0)]                # 顶点
           + [(0.5, 0.5), (0.5, 0.0), (0.0, 0.5)]                # 棱中点
           + [(1 / 3, 1 / 3)])                                    # 面心
    zs = (-GP2, GP2, -1.0, 1.0, 0.0)
    worst = 1e30
    for L1, L2 in tri:
        L3 = 1 - L1 - L2
        for z in zs:
            dN = np.zeros((6, 3))
            # N1..N3 = L_i(1-z)/2 ; N4..N6 = L_i(1+z)/2
            for k, L in enumerate((L1, L2, L3)):
                # d/dL1
                dL1 = 1.0 if k == 0 else 0.0
                dL2 = 1.0 if k == 1 else 0.0
                dN[k, 0] = dL1 * (1 - z) / 2
                dN[k, 1] = dL2 * (1 - z) / 2
                dN[k, 2] = -L / 2
                dN[k + 3, 0] = dL1 * (1 + z) / 2
                dN[k + 3, 1] = dL2 * (1 + z) / 2
                dN[k + 3, 2] = L / 2
            worst = min(worst, float(np.linalg.det(P.T @ dN)))
    return worst


print("\n=== Jacobian 自检 + 自动修正朝向（det(J) 必须 > 0）===")
FN = {"hex8": detJ_hex8, "tet4": detJ_tet4, "penta6": detJ_penta6}
bad_total = 0
fixed_total = 0
small_total = 0
for key in list(elem_blocks):
    dom, et = key
    fn = FN[et]
    v = elem_blocks[key]
    keep, bad_idx = [], []
    dmins = []
    # ★ 近退化「薄片」单元会最先翻转（实测软骨里 det(J) 最小 7.8e-03 ⇒ 体积 0.0013 mm³，
    #   比同类中位数小 2 个数量级）。这类单元在载荷下首批翻负，
    #   所以门槛不能只 >0，要按「同类中位数量级」来设。
    DET_MIN = 0.5
    for e, ids in v:
        d = fn(XYZ[[i - 1 for i in ids]])
        # ★★ 修 LRN-20261004-052 发现的漏洞 ★★
        #   原逻辑：`if d > DET_MIN: keep` 否则**无条件尝试翻转**。
        #   后果：**正但小**的单元（0 < d ≤ DET_MIN）被翻转成负 ⇒ 走进下面的
        #   `if d2 > 0` 判定失败 ⇒ **被当成坏单元剔除** = 网格留洞。
        #   实测：细软骨网格因此被剔除 176 个（基线只剔 9 个）。
        #   正确语义：**正体积就是有效单元**（小 ≠ 无效，只是条件数差）；
        #   只有 d ≤ 0 才需要翻转，而 tet4 的翻转必然变号 ⇒ 不可能两头都负。
        if d > 0:
            keep.append((e, ids))
            dmins.append(d)
            if d <= DET_MIN:
                small_total += 1
            continue
        # ★ 修正朝向必须**按单元类型用正确的置换**：
        #   tet4：交换两个节点（1 次交换 = 奇置换 ⇒ 行列式变号）
        #   hex8/penta6：交换前后两半
        #   ✗ 错误做法：把 4 个节点整体倒序 —— 那是 3 次交换（偶置换），
        #     行列式符号**不变**！结果"翻转后仍为负"→ 被当坏单元剔除，
        #     实测导致 122/462 + 183/618 个软骨单元丢失。
        if et == "tet4":
            alt = [ids[1], ids[0]] + ids[2:]
        else:
            h = len(ids) // 2
            alt = ids[h:] + ids[:h]
        d2 = fn(XYZ[[i - 1 for i in alt]])
        if d2 > 0:
            keep.append((e, alt))
            dmins.append(d2)
            fixed_total += 1
        else:
            bad_idx.append(e)
    elem_blocks[key] = keep
    bad_total += len(bad_idx)
    tag = f"修正 {fixed_total}" if False else ""
    print(f"  {'✅' if not bad_idx else '⚠️'} {dom:22s} {et:7s} 保留 {len(keep):>6,} "
          f"剔除 {len(bad_idx):>4,}  det(J) min={min(dmins) if dmins else float('nan'):.4e}")
print(f"\n  翻转修正 {fixed_total} 个；仍无法修正而剔除 {bad_total} 个；"
      f"保留的「小但正」(≤DET_MIN={DET_MIN}) {small_total} 个")

# ── 6. 写 .feb ──
def area_of(faces):
    a = 0.0
    for f in faces:
        P = XYZ[[i - 1 for i in f]]
        if len(P) == 3:
            a += 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
        else:
            a += 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[3] - P[0])) + \
                 0.5 * np.linalg.norm(np.cross(P[1] - P[2], P[3] - P[2]))
    return a


top_area = area_of(top_sel)
FORCE_N = 600.0
traction = FORCE_N / top_area
print(f"加载面面积 {top_area:.1f} mm² -> 600 N 对应 traction {traction:.4f} MPa")

lines: list[str] = []
w = lines.append
w('<?xml version="1.0" encoding="ISO-8859-1"?>')
w('<febio_spec version="4.0">')
w('  <Module type="solid"/>')
w('  <Globals><Constants><T>0</T><R>0</R><Fc>0</Fc></Constants></Globals>')
w('  <Material>')
mid = {}
for i, (name, (mt, pr)) in enumerate(MATERIALS.items(), 1):
    mid[name] = i
    w(f'    <material id="{i}" name="{name}" type="{mt}">')
    for k, v in pr.items():
        w(f'      <{k}>{v}</{k}>')
    w('    </material>')
w('  </Material>')
w('  <Mesh>')
w('    <Nodes name="all">')
for i, p in enumerate(coords, 1):
    w(f'      <node id="{i}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>')
w('    </Nodes>')
for (dom, et), v in elem_blocks.items():
    # ★ 块名必须唯一：同一域会因混合单元类型(hex8/tet4/penta6)拆成多个 <Elements>，
    #   若都叫 "tibia_cort" FEBio 会报 Duplicate part name。
    w(f'    <Elements type="{et}" name="{dom}__{et}">')
    for e, ids in v:
        w(f'      <elem id="{e}">{",".join(str(x) for x in ids)}</elem>')
    w('    </Elements>')
# ★ NodeSet 的文本**直接就是节点列表**，不能嵌套 <node> 元素（User Manual §3.6.3）
for nm, sel in ((f"{FIX_BONE}_bottom_ns", bot_sel), ("tibia_top_ns", top_sel)):
    w(f'    <NodeSet name="{nm}">'
      f'{",".join(str(x) for x in sorted({x for f in sel for x in f}))}</NodeSet>')
# ★ Surface：子标签**按单元类型命名** <quad4>/<tri3>（User Manual §3.6.5）。
#   既不是 <elem>，也不是 <Elements type="quad4"> —— 用错会报
#   "invalid value for attribute surface"（load 找不到该 surface）
all_surfs = {"tibia_top": top_sel, f"{FIX_BONE}_bottom": bot_sel,
             **{f"surf_{j}": contact_surfs[j] for j in CARTS}}
for _tn, _tp, _ts in tie_pairs:
    all_surfs[_tp] = TIE_SURF[_tn[4:]]["primary"]
    all_surfs[_ts] = TIE_SURF[_tn[4:]]["secondary"]
for nm, sel in all_surfs.items():
    w(f'    <Surface name="{nm}">')
    for et, k in (("quad4", 4), ("tri3", 3)):
        for n, f in enumerate([f for f in sel if len(f) == k], 1):
            w(f'      <{et} id="{n}">{",".join(str(x) for x in f)}</{et}>')
    w('    </Surface>')
# ★ SurfacePair 必须定义在 Surface **之后**（User Manual §3.6.8）
for _tn, _tp, _ts in tie_pairs:
    w(f'    <SurfacePair name="{_tn}">')
    w(f'      <primary>{_tp}</primary>')
    w(f'      <secondary>{_ts}</secondary>')
    w('    </SurfacePair>')
for _pn, _pa, _pb in PAIRS:
    w(f'    <SurfacePair name="{_pn}">')
    w(f'      <primary>surf_{_pa}</primary>')
    w(f'      <secondary>surf_{_pb}</secondary>')
    w('    </SurfacePair>')
w('  </Mesh>')
w('  <MeshDomains>')
_dommat = {**{d: m for d, (k, m) in BONE_PARTS.items()},
           **{j: "mat_cartilage" for j in CARTS}}
for dom, et in elem_blocks:
    # ★ FEBio 4.13 的 SolidDomain.mat 要填**材料名**（不是数字 id）——
    #   写 mat="1" 会报 invalid value for attribute "mat"。
    w(f'    <SolidDomain name="{dom}__{et}" mat="{_dommat[dom]}"/>')
w('  </MeshDomains>')
# ★ 段序必须为 Loads -> Boundary -> Contact（pyfebio 生成的 XML 即此顺序）；
#   顺序错会报 "unrecognized tag" 这种误导性错误。
w('  <LoadData>')
w('    <load_controller id="1" type="loadcurve">')
w('      <interpolate>LINEAR</interpolate><extend>CONSTANT</extend>')
w('      <points><pt>0.0,0.0</pt><pt>1.0,1.0</pt></points>')
w('    </load_controller>')
w('  </LoadData>')
w('  <Loads>')
if DRIVER == "force":
    # ★ 正确写法（实测自 pyfebio 生成的 XML）：
    #   tag 是 <surface_load>，type="traction"，surface 是**属性**，traction 是子元素
    w('    <surface_load type="traction" surface="tibia_top">')
    w(f'      <traction>{-axis[0]*traction:.6f},{-axis[1]*traction:.6f},'
      f'{-axis[2]*traction:.6f}</traction>')
    w('    </surface_load>')
w('  </Loads>')
w('  <Boundary>')
w(f'    <bc name="fix_distal" type="zero displacement" '
  f'node_set="{FIX_BONE}_bottom_ns">')
w('      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>')
w('    </bc>')
if DRIVER == "disp":
    # 顶面三个自由度全部指定位移（=夹持端）⇒ 刚体模态消失，不依赖接触稳定
    kdof = int(np.argmax(np.abs(axis)))
    for di, dn in enumerate("xyz"):
        v = (-DISP_MM if axis[di] > 0 else DISP_MM) if di == kdof else 0.0
        w(f'    <bc name="top_{dn}" type="prescribed displacement" '
          f'node_set="tibia_top_ns">')
        w(f'      <dof>{dn}</dof><value lc="1">{v:.6f}</value><relative>0</relative>')
        w('    </bc>')
w('  </Boundary>')
if PAIRS or tie_pairs:
    w('  <Contact>')
    # ---- 骨 ↔ 软骨底面：tied-elastic（手册 §3.14.6，连接非共形网格）----
    # 参数取手册默认值：penalty=1.0 / auto_penalty=0 / two_pass=0 /
    # tolerance=1.0 / gaptol=0(off) / search_tol=0.01 / search_radius=1.0
    TIE_PEN = float(os.environ.get("TIE_PENALTY", "1.0"))
    TIE_TOL = float(os.environ.get("TIE_TOL", "1.0"))
    TIE_STOL = float(os.environ.get("TIE_SEARCH_TOL", "0.01"))
    TIE_SRAD = float(os.environ.get("TIE_SEARCH_RADIUS", "1.0"))
    TIE_AUTO = os.environ.get("TIE_AUTO_PENALTY", "1") == "1"
    for _tn, _tp, _ts in tie_pairs:
        w(f'    <contact name="{_tn}" surface_pair="{_tn}" '
          'type="tied-elastic">')
        w(f'      <penalty>{TIE_PEN}</penalty>')
        w(f'      <auto_penalty>{1 if TIE_AUTO else 0}</auto_penalty>')
        w('      <two_pass>0</two_pass>')
        w('      <laugon>PENALTY</laugon>')
        w(f'      <tolerance>{TIE_TOL}</tolerance>')
        w('      <symmetric_stiffness>0</symmetric_stiffness>')
        w(f'      <search_tol>{TIE_STOL}</search_tol>')
        w(f'      <search_radius>{TIE_SRAD}</search_radius>')
        w('    </contact>')
    for _pn, _pa, _pb in PAIRS:
        # contact 的 surface_pair 引用的是 **SurfacePair 的名字**（User Manual 3.14）
        w(f'    <contact name="{_pn[:-5]}" surface_pair="{_pn}" '
          'type="sliding-elastic">')
        w('      <laugon>' + LAUGON + '</laugon>')
        w('      <penalty>' + str(PENALTY) + '</penalty>')
        w('      <auto_penalty>' + ('1' if LAUGON == 'PENALTY' else '0') + '</auto_penalty>')
        w('      <two_pass>' + ('1' if TWO_PASS else '0') + '</two_pass>')
        w('      <node_reloc>' + ('1' if NODE_RELOC else '0') + '</node_reloc>')
        w('      <symmetric_stiffness>0</symmetric_stiffness>')
        w('      <tolerance>0.02</tolerance>')
        w('    </contact>')
    w('  </Contact>')
w('  <Output>')
# ★ 不要输出 "contact traction"：它会写出 PLT_FACE_DATA(surface_data)，
#   而 pyfebio.xplt.to_hdf5 的 parse_state 在查 mesh_dict["surfaces"][set_id] 时
#   会 KeyError（它不认我们这些自定义 Surface 集合）→ 整份 .xplt 都读不出来。
w('    <plotfile type="febio"><var type="displacement"/>'
  '<var type="stress"/><var type="reaction forces"/></plotfile>')
w('  </Output>')
# ★ <Control> 是**顶层元素**（单步简写），不包在 <Step> 里 ——
#   见已知可跑通的 temp/pyfebio_demo/my_model.feb
w('  <Control>')
w(f'    <analysis>STATIC</analysis><time_steps>{TIME_STEPS}</time_steps>'
  f'<step_size>{1.0/TIME_STEPS:.6f}</step_size>')
w('    <time_stepper type="default"><max_retries>10</max_retries>'
  '<opt_iter>15</opt_iter><cutback>0.25</cutback></time_stepper>')
w('    <solver type="solid"><linear_solver type="pardiso"/>'
  '<max_refs>25</max_refs><reform_each_time_step>1</reform_each_time_step></solver>')
w('  </Control>')
w('</febio_spec>')

feb = WORK / "ankle_tibiotalar.feb"
feb.write_text("\n".join(lines), encoding="latin-1")
print(f"\n写出 {feb}（{feb.stat().st_size/1024:.0f} KB）")

# ── 6. 求解 ──
if not FEBIO.exists():
    print(f"❌ 找不到 FEBio: {FEBIO}")
    sys.exit(1)
t0 = time.time()
r = subprocess.run([str(FEBIO), "-i", str(feb)], cwd=str(WORK),
                   capture_output=True, text=True, timeout=1800)
print(f"\nFEBio rc={r.returncode}  用时 {time.time()-t0:.1f}s")
out = (r.stdout or "") + (r.stderr or "")
tail = out.strip().split("\n")
print("--- 输出尾部 40 行 ---")
for l in tail[-40:]:
    print("  " + l)
(WORK / "febio.log").write_text(out, encoding="utf-8")
print(f"\n完整日志: {WORK/'febio.log'}")
