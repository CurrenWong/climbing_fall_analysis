"""关节底面补丁重新分网 v2：**2D 参数平面重心网格** + 精确抬回。

── 为什么是 v2 ─────────────────────────────────────────────────────────
v1（提边界环 → 多边形裁剪）**失败**：补丁边界**不是干净闭环**
（tibiotalar_tibia 89 条边界边但只有 69 个节点 ⇒ 有掐点），
多边形追踪碎成 5~23 个环 ⇒ 面积覆盖仅 6~31%。

── v2 的两个关键设计 ───────────────────────────────────────────────────
**① 不做边界裁剪，不做点定位**：直接在**每个原三角内部**打重心网格
   `(i,j,k), i+j+k=m`。新点必然落在原三角内 ⇒ **域 = 原补丁，覆盖率天然 100%**。
**② 用同一组重心坐标精确抬回** ⇒ 新点**精确落在原（分片线性）曲面上**
   ⇒ **底面翘曲恒为 0**（`[LRN-052]` 测出负单元翘曲/对角 0.087~0.123 是根因）。

对比 `cart_hex.py`（v1 细分）失败的原因：它对**扭曲的 quad** 做双线性插值，
中间点飞离两三角形曲面 ⇒ 扭转。这里改成「原 quad 拆 2 三角 + 三角内重心网格 +
精确抬回」，从根上消掉该机制。

**③ 位移场同样用重心坐标插值**（逐角保留原 toward 位移）——
   这是 `[LRN-052]` 里验证过「K=1 逐格复现基线」的正确做法。

── 边的共享 ────────────────────────────────────────────────────────────
三角形 (a,b,c) 的重心网格：
  棱 AB (k=0) 的键 = canonical(a,b,i)    棱 BC (i=0) = (b,c,j)
  棱 CA (j=0) = (c,a,k)                  内部 = (a,b,c,i,j)
同一物理棱两侧算出的**重心坐标相同** ⇒ 位置天然一致；键只用于**共享同一个节点号**。

── v3（2026-10-06）：偏置场解缠，消掉「顶点附近 Jacobian 翻转」───────────
v2 的抬回机制本身没错（底面翘曲恒 0），但**顶面 = 底面 + 逐节点偏置 d 的重心插值**，
当 d 场在单元尺度上变化过大时，楔形被"剪扭"⇒ `det(J)` 在**楔形顶点**上穿零
（完整采样 `detJ_penta6_full`：tibiotalar 两域实测 **42 个 ≤0**）。

**根因（实测，不是"扭角大"或"底三角小"）**：判据
    `min_i ( d_i · n̂_facet ) < 0`  ⟺ 单元翻转
⇒ 某顶点的抬升方向**斜切过它所在底面 facet 的平面**（实测 tibia 12/12 命中、0 误报；
talus 25/31 命中、0 误报；坏单元该夹角中位 **95.9°**、好的 49.5°）。
这正是"楔形三个顶点不都长在底面的同一侧"⇒ 该角 det(J) 翻负。
单指标阈值（扭角 / 面积 / 厚度 / 瘦边）**都会漏**：实测 #30 底面积 17.9mm²（比中位还大）
却 det(J)=−7.89，而扭角 26.7° 的 #32 只有 −0.22。

**做法**：**局部重铺顶面**，不删单元、不动底面（⇒ tie/接触面几何不变）。
把逐节点偏置 `d_i` 重铺成"使所有相邻楔形 det(J) 尽可能正"的场（Adam + 分阶段 τ）。
⚠️ 实测两个重要的**否定结论**：
 ① 严格保模长（只转方向）**做不到**：球面约束 + 线搜索在 talus 上第一步就停滞
    （31/31 坏单元不动）；且 det(J) 的**符号与厚度尺度无关**（对 α 线性）⇒ 缩厚度
    不可能把翻负变正，只会把罚值变小。
 ② 只用"面积加权节点法向"当抬升方向也不行：tibia 仍然 5 个翻负——
    因为单元的**顶点间**方向差才是主因，不是单个顶点的方向选择。

── 开关 ────────────────────────────────────────────────────────────────
CP_M        每三角边的分段数 m（子三角数 = m²），默认 4
CP_OUT      输出文件名，默认 cartilage_mesh.npz
CP_OUTDIR   输出目录（相对 temp/thums），默认 cart_m1_fix
            ⚠️ 硬保护：不允许写到 cart_m1（当前基线），除非 CP_FORCE=1
CP_UNTANGLE 默认 1（=修好的版本）；设 0 关闭解缠，退回 v2 行为
CP_UNT_TAU / CP_UNT_IT / CP_UNT_LR / CP_UNT_SM  解缠超参（一般不用动）
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay

ROOT = Path(r"D:\Project\climbing_fall_analysis")
THUMS = ROOT / "temp/thums"

M = int(os.environ.get("CP_M", "4") or 4)
OUTNPZ = os.environ.get("CP_OUT", "cartilage_mesh.npz")
OUTDIR = THUMS / (os.environ.get("CP_OUTDIR", "cart_m1_fix").strip() or "cart_m1_fix")
CP_FORCE = os.environ.get("CP_FORCE", "").strip() == "1"
if OUTDIR.resolve() == (THUMS / "cart_m1").resolve() and not CP_FORCE:
    raise SystemExit("⛔ 拒写 temp/thums/cart_m1（当前基线）。"
                     "改 CP_OUTDIR，或确实要覆盖时设 CP_FORCE=1。")
# ── 解缠开关（默认开启 = 修好的版本）──
CP_UNT = os.environ.get("CP_UNTANGLE", "1").strip() != "0"
CP_UNT_TAU = float(os.environ.get("CP_UNT_TAU", "0.6") or 0.6)
CP_UNT_IT = int(os.environ.get("CP_UNT_IT", "3000") or 1600)
CP_UNT_LR = float(os.environ.get("CP_UNT_LR", "0.010") or 0.010)
CP_UNT_SM = float(os.environ.get("CP_UNT_SM", "0.001") or 0.001)
CP_UNT_LO = float(os.environ.get("CP_UNT_LO", "0.35") or 0.35)   # 逐节点厚度下界 (mm)
CP_UNT_HI = float(os.environ.get("CP_UNT_HI", "1.58") or 1.58)   # 逐节点厚度上界 (mm)
CP_UNT_TMAX = float(os.environ.get("CP_UNT_TMAX", "300") or 180)  # 单关节秒预算
JOINTS = ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]
WEDGE_TETS = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]

# ── 楔形 det(J) 的完整采样（与 ankle_lig_feb.detJ_penta6_full 同一套）──
_URS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
_UTS = (0.0, 0.25, 0.5, 0.75, 1.0)
_UG = np.array([np.array([
    [-(1 - t), -(1 - t), -(1 - r - s)],
    [1 - t, 0.0, -r],
    [0.0, 1 - t, -s],
    [-t, -t, (1 - r - s)],
    [t, 0.0, r],
    [0.0, t, s],
]) for r in _URS for s in _URS if r + s <= 1.0 + 1e-9 for t in _UTS])   # (nS,6,3)


def wedge_dets(P6):
    """(n,6,3) 楔形角点 -> (n, nS) 逐采样点 det(J)。"""
    return np.linalg.det(np.einsum("nic,mia->nmca", P6, _UG))


def guard_of(jname, bpos, newell):
    """对侧骨的**关节面采样**（= 对侧软骨补丁的底面多边形，顶点+质心）+ 面法向。

    法向定向到"关节腔"一侧（用该面的原始抬升方向判定）⇒ 有符号距 d>0 = 在骨外。
    只覆盖补丁范围（够用：抬升方向只在这附近变化）。
    """
    E = np.asarray(src[f"{jname}_elems"], np.int64)
    OFF = np.asarray(src[f"{jname}_offsets"], float)
    pts, nrm, tri, trin = [], [], [], []
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        base, outer = v[:m], v[m:]
        if len(base) < 3:
            continue
        P = np.array([bpos[r] for r in base])
        D = np.array([OFF[-o - 1] - bpos[r] for r, o in zip(base, outer)])
        nv = newell(P)
        nv = nv / (np.linalg.norm(nv) + 1e-15)
        if float(np.dot(nv, D.mean(axis=0))) < 0:
            nv = -nv                       # 指向关节腔（= 背离本骨）
        for q in P:
            pts.append(q)
            nrm.append(nv)
        pts.append(P.mean(axis=0))
        nrm.append(nv)
        for (i0, i1, i2) in (((0, 1, 2),) if len(base) == 3 else ((0, 1, 2), (0, 2, 3))):
            T = (P[i0], P[i1], P[i2])
            nn = np.cross(T[1] - T[0], T[2] - T[0])
            nn = nn / (np.linalg.norm(nn) + 1e-15)
            if float(np.dot(nn, nv)) < 0:
                nn = -nn                   # 与该面（指向关节腔）同向
            tri.append(T)
            trin.append(nn)
    return np.array(pts), np.array(nrm), np.array(tri), np.array(trin)


def _tri_closest(p, a, b, c):
    """点 p 到三角形 (a,b,c) 的最近点（Ericson）；返回 (q, 是否落在三角内部)。

    ★ 内部性很关键：补丁边界外，最近的三角只能给出"平面延伸"位置，
      那时有符号距离没有意义（会误判成穿透）。
    """
    ab, ac = b - a, c - a
    ap = p - a
    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)
    bp = p - b
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    cp = p - c
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)
    vc = d1 * d4 - d3 * d2
    vb = d5 * d2 - d1 * d6
    va = d3 * d6 - d5 * d4
    den = va + vb + vc
    ok = np.abs(den) > 1e-18
    q = a + ab * np.where(ok, vb / np.where(ok, den, 1.0), 0.0)[:, None] \
        + ac * np.where(ok, vc / np.where(ok, den, 1.0), 0.0)[:, None]
    inner = ok.copy()
    m = (d1 <= 0) & (d2 <= 0)
    q[m] = a[m]
    inner[m] = False
    m = (d3 >= 0) & (d4 <= d3)
    q[m] = b[m]
    inner[m] = False
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    q[m] = a[m] + (d1 / (d1 - d3))[m, None] * ab[m]
    inner[m] = False
    m = (d6 >= 0) & (d5 <= d6)
    q[m] = c[m]
    inner[m] = False
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    q[m] = a[m] + (d2 / (d2 - d6))[m, None] * ac[m]
    inner[m] = False
    m = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    q[m] = b[m] + ((d4 - d3) / ((d4 - d3) + (d5 - d6)))[m, None] * (c[m] - b[m])
    inner[m] = False
    return q, inner


def untangle_disp(disp, tris, pos, guard=None, geps=0.01, grad_rad=2.0):
    """★ 偏置场解缠：重铺逐节点偏置 d_i，使所有楔形 worst-det(J) 远离 0。

    只动 d（顶面），底面 `pos` 完全不动 ⇒ tie/接触面几何不变。
    `tris` = 域三角（底面，顶点=骨节点 id）；每个域三角配一个 M=1 的楔形作为判据
    （单元内偏置场是仿射的 ⇒ M=1 楔形处处 det>0 ⟹ 其任意重心子三角也一样）。

    返回 (新 disp, stats)。
    """
    vlist = sorted({r for t in tris for r in t})
    vid = {r: i for i, r in enumerate(vlist)}
    A = np.array([[vid[r] for r in t] for t in tris], np.int64)      # (nE,3)
    bot = np.array([pos[r] for r in vlist])                         # (nv,3)
    d0 = np.array([disp[r] for r in vlist])                         # (nv,3)
    # 与 build 阶段同一套 newell 定向（先按原始 d 定死，避免目标函数跳变）
    n0 = np.cross(bot[A[:, 1]] - bot[A[:, 0]], bot[A[:, 2]] - bot[A[:, 0]])
    fl = np.einsum("ij,ij->i", n0, d0[A].mean(axis=1)) < 0
    A[fl] = A[fl][:, ::-1]
    B0 = bot[A]                                                     # (nE,3,3)
    edges = np.array(sorted({(min(t[a], t[(a + 1) % 3]), max(t[a], t[(a + 1) % 3]))
                             for t in A for a in range(3)}))
    tt = np.linalg.norm(d0, axis=1)
    lo = min(CP_UNT_LO, float(tt.min()))    # 逐节点厚度只能落在 [下界, 上界] 内
    hi = max(CP_UNT_HI, float(tt.max()))
    # ── 穿透护栏：对侧骨面采样点 + 法向（硬投影，不允许顶节点越过骨面）──
    if guard is not None and len(guard[0]):
        from scipy.spatial import cKDTree
        _Gtri = np.array(guard[2])                            # (m,3,3) 对侧骨面三角
        _Gnrm = np.array(guard[3])                            # (m,3) 指向关节腔
        _gt = cKDTree(_Gtri.mean(axis=1))
        _kq = min(6, len(_Gtri))

    D = d0.copy()

    def _guard_dist(Dv):
        """每个顶节点对「最近的内部三角」的有符号距（近表面才算）。"""
        top = bot + Dv
        _, iq = _gt.query(top, k=_kq)
        R = np.broadcast_to(top[:, None, :], (len(top), iq.shape[1], 3))
        Rt = R.reshape(-1, 3)
        tri = _Gtri[iq.reshape(-1)]
        q, inner = _tri_closest(Rt, tri[:, 0], tri[:, 1], tri[:, 2])
        nr = _Gnrm[iq.reshape(-1)]
        rr = np.linalg.norm(Rt - q, axis=1)
        dd = np.einsum("ij,ij->i", nr, Rt - q)
        dd = np.where(inner & (rr < grad_rad), dd, np.inf).reshape(len(top), -1)
        return dd.min(axis=1)                       # (nv,) 可为 +inf / -inf

    if _gt is not None:
        # ★ 阈值 = min(geps, 原始状态的值) ⇒ 护栏**只会阻止变差，绝不改动原始抬升**
        #   （局部平面在凹面/穹顶旁边会把"贴边"点误判，此规则天然免疫）
        _d0g = _guard_dist(d0)
        _thr = np.minimum(geps, _d0g)

        def _guard_project(Dv):
            dcur = _guard_dist(Dv)
            bad = np.isfinite(dcur) & (dcur < _thr)
            if not bad.any():
                return Dv, 0
            top = bot + Dv
            _, iq = _gt.query(top, k=_kq)
            R = np.broadcast_to(top[:, None, :], (len(top), iq.shape[1], 3))
            Rt = R.reshape(-1, 3)
            tri = _Gtri[iq.reshape(-1)]
            q, inner = _tri_closest(Rt, tri[:, 0], tri[:, 1], tri[:, 2])
            nr = _Gnrm[iq.reshape(-1)]
            rr = np.linalg.norm(Rt - q, axis=1)
            dd = np.einsum("ij,ij->i", nr, Rt - q)
            ddv = np.where((inner & (rr < grad_rad)).reshape(len(top), -1),
                           dd.reshape(len(top), -1), np.inf)
            j = np.argmin(ddv, axis=1)
            idx = np.arange(len(top))
            nj = nr.reshape(len(top), -1, 3)[idx, j]
            push = np.zeros_like(top)
            push[bad] = ((_thr[bad] - dcur[bad])[:, None] * nj[bad])
            return Dv + push, int(bad.sum())
    else:
        def _guard_project(Dv):
            return Dv, 0

    dmin0 = float(wedge_dets(np.concatenate([B0, B0 + D[A]], axis=1)).min())
    nneg0 = int((wedge_dets(np.concatenate([B0, B0 + D[A]], axis=1)).min(axis=1) <= 0).sum())
    best_all = None
    _t_start = time.time()
    for tau, it_n, lr in ((CP_UNT_TAU, CP_UNT_IT, CP_UNT_LR),
                          (0.25, max(CP_UNT_IT // 2, 300), CP_UNT_LR * 0.6),
                          (0.05, max(CP_UNT_IT // 4, 300), CP_UNT_LR * 0.6)):
        m = np.zeros_like(D)
        v = np.zeros_like(D)
        best = None
        if time.time() - _t_start > CP_UNT_TMAX:
            break
        for it in range(it_n):
            if time.time() - _t_start > CP_UNT_TMAX:
                break
            P6 = np.concatenate([B0, B0 + D[A]], axis=1)
            dt = wedge_dets(P6)
            mins = dt.min(axis=1)
            F = float(np.sum(np.maximum(0.0, tau - dt) ** 2))
            if CP_UNT_SM:
                dd = D[edges[:, 0]] - D[edges[:, 1]]
                F += CP_UNT_SM * float(np.sum(dd ** 2))
            cur = (int((mins <= 0).sum()), F)
            if best is None or cur < best[0]:
                best = (cur, D.copy(), float(mins.min()))
            if mins.min() >= tau:
                break
            sel = np.where(mins < tau)[0]
            g = np.zeros_like(D)
            for e in sel:
                Pe = P6[e]
                for mm in np.where(dt[e] < tau)[0]:
                    gm = _UG[mm]
                    J = Pe.T @ gm
                    det = float(np.linalg.det(J))
                    Jinv = np.linalg.inv(J + (1e-9 if det <= 1e-12 else 0.0) * np.eye(3))
                    w = -2.0 * (tau - det) * det
                    for i in range(3, 6):
                        g[A[e, i - 3]] += w * (gm[i] @ Jinv)
            if CP_UNT_SM:
                np.add.at(g, edges[:, 0], 2 * CP_UNT_SM * (D[edges[:, 0]] - D[edges[:, 1]]))
                np.add.at(g, edges[:, 1], 2 * CP_UNT_SM * (D[edges[:, 1]] - D[edges[:, 0]]))
            m = 0.9 * m + 0.1 * g
            v = 0.999 * v + 0.001 * g ** 2
            D = D - lr * (m / (1 - 0.9 ** (it + 1))) / (
                np.sqrt(v / (1 - 0.999 ** (it + 1))) + 1e-12)
            L = np.linalg.norm(D, axis=1)
            D = D * (np.clip(L, lo, hi) / np.maximum(L, 1e-12))[:, None]
            D, _ = _guard_project(D)                       # 穿透护栏（硬投影）
        if best is not None:
            D = best[1]
            if best_all is None or best[0] < best_all[0]:
                best_all = (best[0], best[1].copy(), best[2])
        if best is not None and best[2] >= tau:
            break
    if best_all is not None:
        D = best_all[1]
    dt = wedge_dets(np.concatenate([B0, B0 + D[A]], axis=1))
    mins = dt.min(axis=1)
    stats = dict(before=dmin0, after=float(mins.min()),
                 nneg_before=nneg0,
                 nneg_after=int((mins <= 0).sum()),
                 nbad_after=int((mins <= 1e-2).sum()),
                 dmove=float(np.abs(D - d0).max()))
    return {r: D[vid[r]] for r in vlist}, stats


def tet_vol(P):
    return float(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0])) / 6.0


def min_angle_tri(P):
    a = []
    for k in range(3):
        u = P[(k - 1) % 3] - P[k]
        v = P[(k + 1) % 3] - P[k]
        c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-15)
        a.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return min(a)


def newell(p):
    n = np.zeros(3)
    for i in range(len(p)):
        n += np.cross(p[i], p[(i + 1) % len(p)])
    return n


src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bone_ids, bone_xyz = src["bone_node_ids"], src["bone_xyz"]
bpos = {int(v): bone_xyz[i] for i, v in enumerate(bone_ids)}

# ── ★ 穿透护栏：对侧骨面采样（用于解缠时不让顶节点越过对侧骨）──
OPPOSITE = {"tibiotalar_tibia": "tibiotalar_talus",
            "tibiotalar_talus": "tibiotalar_tibia",
            "subtalar_talus": "subtalar_calcaneus",
            "subtalar_calcaneus": "subtalar_talus"}
GUARD: dict = {}
for _j in JOINTS:
    try:
        GUARD[_j] = guard_of(OPPOSITE[_j], bpos, newell)
        print(f"[guard] {_j:22s} 对侧骨面采样点 {len(GUARD[_j][0])}")
    except Exception as _e:                                    # noqa: BLE001
        GUARD[_j] = None
        print(f"[guard] {_j}: 取不到对侧骨面（{_e}），跳过护栏")

out: dict[str, np.ndarray] = {"bone_node_ids": bone_ids, "bone_xyz": bone_xyz}
rep: dict[str, dict] = {}

for j in JOINTS:
    E = np.asarray(src[f"{j}_elems"], np.int64)
    old_off = np.asarray(src[f"{j}_offsets"], float)
    faces = []
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        faces.append((v[:m], v[m:]))

    nodes = sorted({r for b, _ in faces for r in b})
    # ── 位移场（逐节点中位数，稳健）──
    dacc: dict[int, list[np.ndarray]] = {r: [] for r in nodes}
    for base, outer in faces:
        for r, o in zip(base, outer):
            dacc[r].append(old_off[-o - 1] - bpos[r])
    disp = {r: np.median(np.array(v), axis=0) for r, v in dacc.items() if v}
    for r in nodes:
        disp.setdefault(r, np.zeros(3))

    # ── 域三角 ──
    # CP_DOMAIN=bary（默认）：用**原面的三角化**（quad 选最优对角线）
    # CP_DOMAIN=delaunay  ：用**参数平面 Delaunay**（Delaunay 最大化最小角）
    #   动机：原面三角化里存在瘦长条（实测最小内角最低 8.2°），重心细分会**继承**
    #   父元形状 ⇒ 挤出后出薄片 ⇒ 求解器报 Negative jacobian。
    #   Delaunay 换掉瘦长父元；再按「质心落在原补丁内」过滤掉补丁外的三角。
    #   抬回仍用该域三角的重心坐标 ⇒ **精确落在原曲面**（翘曲 0），无需点定位。
    DOMAIN = os.environ.get("CP_DOMAIN", "bary").strip().lower()
    DIAG = os.environ.get("CP_DIAG", "best").strip().lower()   # best | 02

    def _ma(P):
        a = []
        for k in range(3):
            u = P[(k - 1) % 3] - P[k]
            v = P[(k + 1) % 3] - P[k]
            c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-15)
            a.append(np.arccos(np.clip(c, -1, 1)))
        return min(a)

    # ① 原面三角化（bary 模式用，也是 delaunay 模式的「补丁域」判据）
    tris: list[tuple[int, int, int]] = []
    ndiag_swap = 0
    for base, _ in faces:
        if len(base) == 3:
            tris.append(tuple(base))
            continue
        Q = [np.array(bpos[r]) for r in base]
        m02 = min(_ma(np.array([Q[0], Q[1], Q[2]])), _ma(np.array([Q[0], Q[2], Q[3]])))
        m13 = min(_ma(np.array([Q[1], Q[2], Q[3]])), _ma(np.array([Q[1], Q[3], Q[0]])))
        if DIAG == "02":
            m13 = -1.0        # ★ 强制 0-2 对角（与 build_feb 从 quad 推出的三角一致）
        if m13 > m02:
            tris.append((base[1], base[2], base[3]))
            tris.append((base[1], base[3], base[0]))
            ndiag_swap += 1
        else:
            tris.append((base[0], base[1], base[2]))
            tris.append((base[0], base[2], base[3]))

    if DOMAIN == "delaunay":
        # ② PCA 平面 → 投影（可行性实测：无折叠 / 单连通）
        idx_of = {r: i for i, r in enumerate(nodes)}
        Pp = np.array([bpos[r] for r in nodes])
        cc = Pp.mean(axis=0)
        _, _, Vt = np.linalg.svd(Pp - cc, full_matrices=False)
        nv, e1, e2 = Vt[2], Vt[0], Vt[1]
        if np.dot(np.cross(e1, e2), nv) < 0:
            e2 = -e2
        UVp = np.stack([(Pp - cc) @ e1, (Pp - cc) @ e2], axis=1)
        # 原补丁的 2D 域三角
        DT = UVp[np.array([[idx_of[r] for r in t] for t in tris])]
        # 面积足够大的原三角才作为「域」判据（退化三角会误纳）
        def _ar2(q):
            return 0.5 * abs((q[1, 0] - q[0, 0]) * (q[2, 1] - q[0, 1])
                             - (q[2, 0] - q[0, 0]) * (q[1, 1] - q[0, 1]))
        DT = DT[np.array([_ar2(q) for q in DT]) > 1e-6]

        def in_domain(pts):
            """pts:(n,2) -> 是否落在某个原三角内（含边界）。"""
            hit = np.zeros(len(pts), bool)
            for q in DT:
                if hit.all():
                    break
                d = pts[~hit]
                if not len(d):
                    break
                v0, v1, v2 = q[1] - q[0], q[2] - q[0], d - q[0]
                den = v0[0] * v1[1] - v0[1] * v1[0]
                if abs(den) < 1e-14:
                    continue
                u = (v2[:, 0] * v1[1] - v1[0] * v2[:, 1]) / den
                w = (v0[0] * v2[:, 1] - v2[:, 0] * v0[1]) / den
                ok = (u >= -1e-9) & (w >= -1e-9) & (u + w <= 1 + 1e-9)
                idxs = np.where(~hit)[0][ok]
                hit[idxs] = True
            return hit

        tri2 = Delaunay(UVp, qhull_options="QJ")
        cand = tri2.simplices
        cens = UVp[cand].mean(axis=1)
        keep = in_domain(cens)
        cand = cand[keep]
        tris = [(nodes[c[0]], nodes[c[1]], nodes[c[2]]) for c in cand]
        ndiag_swap = -1   # 标记：用了 delaunay 域

    # ── ★ 偏置场解缠（局部重铺顶面，不动底面）──
    if CP_UNT:
        disp, ust = untangle_disp(disp, tris, bpos, guard=GUARD.get(j))
        print(f"  [{j:18s}] 解缠: worst-det(J) {ust['before']:+.4f} → {ust['after']:+.4f}  "
              f"≤0 单元 {ust['nneg_before']} → {ust['nneg_after']}  "
              f"(≤1e-2: {ust['nbad_after']})  逐节点偏置最大改动 {ust['dmove']:.4f}mm")

    # ── 节点表 ──
    nlist: list[np.ndarray] = []
    tag: list[bool] = []

    def add(p, is_off):
        nlist.append(np.asarray(p, float))
        tag.append(bool(is_off))
        return -(len(nlist))

    _key: dict[tuple, int] = {}

    def node_at(a, b, c, i, kk):
        """三角形 (a,b,c) 的重心网格点（i+j+k=m）——按棱/内部键共享。"""
        jj = M - i - kk
        # 顶点
        if (i, kk) == (M, 0):
            return _vert(a)
        if (i, kk) == (0, 0):
            return _vert(b)
        if (i, kk) == (0, M):
            return _vert(c)
        if kk == 0:                                   # 棱 AB
            x, y = (a, b) if a < b else (b, a)
            ii = i if a < b else M - i
            key = ("E", x, y, ii)
            w = ((ii / M, 1 - ii / M, 0.0) if a < b else ((M - ii) / M, ii / M, 0.0))
        elif i == 0:                                  # 棱 BC
            x, y = (b, c) if b < c else (c, b)
            ii = jj if b < c else M - jj
            key = ("E", x, y, ii)
            w = ((0.0, ii / M, 1 - ii / M) if b < c else (0.0, (M - ii) / M, ii / M))
        elif jj == 0:                                 # 棱 CA
            x, y = (c, a) if c < a else (a, c)
            ii = kk if c < a else M - kk
            key = ("E", x, y, ii)
            w = ((ii / M, 0.0, 1 - ii / M) if c < a else ((M - ii) / M, 0.0, ii / M))
        else:                                         # 内部
            key = ("F", a, b, c, i, kk)
            w = (i / M, jj / M, kk / M)
        if key not in _key:
            P = w[0] * bpos[a] + w[1] * bpos[b] + w[2] * bpos[c]
            D = w[0] * disp[a] + w[1] * disp[b] + w[2] * disp[c]
            _key[key] = (add(P, False), add(P + D, True))
        return _key[key]

    _vertc: dict[int, tuple[int, int]] = {}

    def _vert(r):
        if r not in _vertc:
            _vertc[r] = (add(bpos[r], False), add(bpos[r] + disp[r], True))
        return _vertc[r]

    # ── 生成：每个原三角 → m² 个子三角，每个子三角挤成楔 ──
    new_cells: list[tuple[list[int], list[int]]] = []
    tie_faces: list[list[int]] = []
    for (a, b, c) in tris:
        for i in range(M):
            for k in range(M - i):
                g = [(i, k), (i + 1, k), (i, k + 1)]
                blk = [node_at(a, b, c, u, v) for (u, v) in g]
                bb = [t[0] for t in blk]
                oo = [t[1] for t in blk]
                Pa = np.array([nlist[-r - 1] for r in bb])
                Pb = np.array([nlist[-r - 1] for r in oo])
                if np.dot(newell(Pa), Pb.mean(axis=0) - Pa.mean(axis=0)) < 0:
                    bb, oo = bb[::-1], oo[::-1]
                new_cells.append((bb, oo))
                tie_faces.append(bb)

    # ── 自检 ──
    def P_of(r):
        return nlist[-r - 1]

    nneg, g3, g2, vols = 0, [], [], []
    for blk, top in new_cells:
        Pw = np.array([P_of(r) for r in blk + top])
        v = sum(tet_vol(Pw[list(t)]) for t in WEDGE_TETS)
        vols.append(v)
        if v <= 0:
            nneg += 1
        else:
            g3.append(min_angle_tri(Pw[:3]))
            g2.append(min_angle_tri(Pw[3:]))

    area_new = sum(0.5 * np.linalg.norm(np.cross(
        np.array([bpos[v] for v in t])[1] - np.array([bpos[v] for v in t])[0],
        np.array([bpos[v] for v in t])[2] - np.array([bpos[v] for v in t])[0])) for t in tris)

    rep[j] = dict(ncell=len(new_cells), nnode=len(nlist), nneg=nneg,
                  amin3=min(g3) if g3 else 0,
                  amin3med=float(np.median(g3)) if g3 else 0,
                  amin2med=float(np.median(g2)) if g2 else 0,
                  area=area_new, ntris=len(tris),
                  ntri_ok=sum(1 for t in tris if len(t) == 3),
                  ndiag=ndiag_swap)

    out[f"{j}_nodes"] = np.array(nlist, float)
    out[f"{j}_is_offset"] = np.array(tag, bool)
    Em = np.zeros((len(new_cells), 6), np.int64)
    for ii, (b_, o_) in enumerate(new_cells):
        Em[ii, :3] = b_
        Em[ii, 3:] = o_
    out[f"{j}_elems"] = Em
    Tt = np.zeros((len(tie_faces), 3), np.int64)
    for ii, f in enumerate(tie_faces):
        Tt[ii] = f
    out[f"{j}_tie_faces"] = Tt
    out[f"{j}_bone_faces"] = np.zeros((0, 4), np.int64)

    print(f"[{j:20s}] 面 {len(new_cells):6d}(原 {len(faces):4d}, 三角化 {len(tris):4d}, "
          f"换对角 {ndiag_swap:4d})  负体积 {nneg:5d}  "
          f"底最小内角 {rep[j]['amin3']:5.1f}° / 中位 {rep[j]['amin3med']:5.1f}°  "
          f"节点 {len(nlist):6d}")

OUTDIR.mkdir(parents=True, exist_ok=True)
np.savez_compressed(OUTDIR / OUTNPZ, **out)
print(f"\n写出 {OUTDIR / OUTNPZ}  m={M}  ({(OUTDIR / OUTNPZ).stat().st_size/1e6:.2f} MB)")

# ── 对**成品网格**再做一次完整采样自检（用存下来的连接/顺序 = FEBio 看到的那份）──
print("\n── 成品网格 worst-det(J) 完整采样自检 ──")
_tot_n, _tot_s, _tot = 0, 0, 0
for j in JOINTS:
    E = out[f"{j}_elems"]
    if not len(E):
        continue
    Nn = out[f"{j}_nodes"]

    def _ix(x):
        x = int(x)
        return (-x - 1) if x < 0 else (x - 1)

    P6 = np.array([[Nn[_ix(x)] for x in row] for row in E])
    dt = wedge_dets(P6)
    mins = dt.min(axis=1)
    nneg = int((mins <= 0).sum())
    nsm = int((mins <= 1e-2).sum())
    _tot_n += nneg
    _tot_s += nsm
    _tot += len(E)
    print(f"  {j:22s} 单元 {len(E):5d}  worst-det(J) min {mins.min():+9.4f}  "
          f"≤0 的 {nneg:4d}  0<x≤1e-2 的 {nsm:4d}")
print(f"  ⇒ 合计 {_tot} 个楔形：≤0 的 {_tot_n} 个（目标 0）、≤1e-2 的 {_tot_s} 个（目标 0）")

tc = sum(r["ncell"] for r in rep.values())
tn = sum(r["nneg"] for r in rep.values())
print(f"① 关节面三角总数 : {tc:6d}   （基线底面 1119；v1 细分 m=4 应 = 4²×三角化数）")
print(f"② 负体积（翻转前）: {tn:6d} / {tc} ({100*tn/max(tc,1):.2f}%)  ← 目标 0")
print(f"③ 最小内角 全局最低: {min(r['amin3'] for r in rep.values()):.1f}°  （基线底 quad 最低 20°）")
