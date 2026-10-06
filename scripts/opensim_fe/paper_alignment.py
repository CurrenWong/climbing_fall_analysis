"""paper_alignment.py -- 把骨折评估重新映射到 [Y25] 论文的**子区域**（ends / femoral neck vs shaft）。

背景（为什么需要这个脚本）
--------------------------
之前 `fracture_matrix.py` / `risk_1d_loadshare.py` / `risk_1d_bending.py` 一律用
**骨干** 强度（tibia 200 / fibula 160 / femur 220 MPa）作为胫/腓/股三骨的 σ_c。
这是**不会骨折的子区域**：论文 [Y25] 的骨折部位是**两端**（70 MPa）与**股骨颈**
（80 MPa），骨干**从不骨折**。本脚本是**加法式**对齐：把每个部位拆成论文的子区域，
用论文 Table 1 的**压缩强度**，重建 9 部位 × 50 高度的二元骨折矩阵，跑 S4 统计
（Logistic 回归 + Jaccard 层次聚类），并产出「论文对齐清单」。

子区域 → 强度（[Y25] Table 1, 压缩 MPa；取自 `src/climbing/bone.py::MATERIAL_STRENGTH_MPA`）
------------------------------------------------------------------------------------
* 足部   → 跟骨(整体) 150
* 胫骨   → 两端 70（+ 骨干 200，仅用于验证“骨干不骨折”）
* 腓骨   → 两端 70（+ 骨干 160，同上）
* 股骨   → 股骨颈 80（+ 骨干 220，同上）
* 骨盆 180；腰椎/胸椎/颈椎 150；颅骨 160
**部位骨折 ⇔ 其任一子区域超过该子区域阈值。**

两条路线
--------
* Route-1D（主）：σ(h) 取自 `risk_1d_bending.json` 的 `σ_total`（risk_by_height × 该骨
  σ_c）。部位骨折 ⇔ σ(h) ≥ σ_c(弱子区域=两端/股骨颈)；另用骨干阈值单独评估骨干。
* Route-FE（交叉验证）：逐单元应力取自
  `temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5` 的末 state，按**沿载荷轴的几何位置**
  把主承力域切成 ends/shaft（rule 见下），各子区域取 p95 再按严格线弹性缩放到各高度。
  仅对 3 根长骨做 ends/shaft 切分；足（紧凑骨）与骨盆/脊柱/颅骨（壳，单一阈值）用
  `fracture_matrix.json` 已有的整域 p95。

ends/shaft 分区规则（Route-FE，精确可复现）
------------------------------------------
1. 载荷轴 = 主承力域节点坐标 bbox **最长轴**（与 `risk_1d_bending.py` 的 θ 定义一致）。
2. 单元形心 `c_e` = 其节点坐标均值；投影 `t_e = c_e · axis`。
3. `L = t_max − t_min`；**ends** = `{t ≤ t_min + 0.25·L} ∪ {t ≥ t_max − 0.25·L}`
   （即两端各取外 25%，共 50% 跨度）；**shaft** = 中间 50%。
4. 各子区域统计量 = 该子区域内单元的 von Mises p95（主口径，与既有 FE 矩阵一致）。

只读输入 / 只写新文件；不改任何既有模块默认/签名；不调用 FEBio；不覆盖既有产物。

复现（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\paper_alignment.py
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from climbing.bone import MATERIAL_STRENGTH_MPA  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
FM_JSON = RES / "fracture_matrix.json"
BENDING_JSON = RES / "risk_1d_bending.json"
R1_JSON = RES / "bc_robust_metric_route1.json"
FM_DIR = ROOT / "temp" / "opensim_fe" / "fracture_matrix"

#: ends = 两端各取 bbox 载荷轴跨度的外 25%（"top/bottom ~20-30%"）
END_FRACTION = 0.25

#: 论文 [Y25] 的「阶跃响应」部位（Logistic 因完全分离而拟合失败）。
PAPER_STEP = {"足部", "胫骨", "股骨", "腰椎", "胸椎", "颈椎"}
#: 论文 [Y25] 的「渐进损伤」部位（Logistic 拟合成功）。
PAPER_GRADUAL = {"腓骨", "骨盆", "颅骨"}
#: 论文 OR（每米）排序目标。
PAPER_OR_ORDER = ["腓骨", "颅骨", "骨盆"]
#: 论文 5 类 / 3 阶段的目标边界（m）。
PAPER_CLUSTER_BOUNDS = [(1, 6), (7, 9), (10, 16), (17, 34), (35, 50)]

#: 9 部位 → 代表骨 / 关节 / 子区域列表（子区域键 → MATERIAL_STRENGTH_MPA 键）。
#: 子区域顺序：**弱(两端/颈)在前**，骨干在后；单阈值部位只有一个「整体」子区域。
REGIONS: list[tuple[str, str, str, list[tuple[str, str]]]] = [
    ("足部", "calcaneus_r", "subtalar_r", [("跟骨(整体)", "calcaneus")]),
    ("胫骨", "tibia_r", "ankle_r",
     [("胫骨两端", "tibia_ends"), ("胫骨骨干", "tibia_shaft")]),
    ("腓骨", "fibula_r", "ankle_r",
     [("腓骨两端", "fibula_ends"), ("腓骨骨干", "fibula_shaft")]),
    ("股骨", "femur_r", "hip_r",
     [("股骨颈", "femoral_neck"), ("股骨骨干", "femur_shaft")]),
    ("骨盆", "R_HIPBONE", "hip_r", [("骨盆(整体)", "pelvis")]),
    ("腰椎", "L3", "lumbar", [("椎体(整体)", "spine")]),
    ("胸椎", "T6", "lumbar", [("椎体(整体)", "spine")]),
    ("颈椎", "C5", "lumbar", [("椎体(整体)", "spine")]),
    ("颅骨", "parietal_r", "lumbar", [("颅顶(整体)", "skull")]),
]
REGION_NAMES = [r[0] for r in REGIONS]

#: 有「两端/颈 vs 骨干」语义的长骨（Route-FE 切分对象）。
LONG_BONES = {"tibia_r", "fibula_r", "femur_r"}


def _sc(key: str) -> float:
    """压缩强度 (MPa)，论文 Table 1。"""
    return float(MATERIAL_STRENGTH_MPA[key][1])


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------
def _first_fracture(heights, mask) -> float | None:
    for h, m in zip(heights, mask):
        if bool(m):
            return float(h)
    return None


def vm_voigt(s: np.ndarray) -> np.ndarray:
    """6 分量 Voigt 应力 → von Mises（与 bone_feb._von_mises_from_voigt 同式）。"""
    s = np.asarray(s, dtype=np.float64)
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2)
                   + 3.0 * (xy ** 2 + yz ** 2 + xz ** 2))


def read_h5_domain(bone: str, dom: str) -> dict:
    """读取 <bone>.hdf5 的节点/连接/末 state 应力，返回形心与节点坐标。"""
    import h5py
    path = FM_DIR / bone / f"{bone}.hdf5"
    with h5py.File(path, "r") as h:
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(np.float64)
        conn = h["meshes/0/domains"][dom][:]
        sids = sorted(h["states"].keys(), key=int)
        stress = h[f"states/{sids[-1]}/element_data/stress/{dom}"][:].astype(np.float64)
    idx = conn[:, 1:] - 1                     # 连接首列为单元号，节点 id 1-based
    cent = coords[idx].mean(axis=1)
    return dict(path=path, coords=coords, cent=cent, stress=stress)


def partition_ends_shaft(coords, cent, end_fraction=END_FRACTION) -> dict:
    """按载荷轴（bbox 最长轴）把单元切成 ends / shaft，返回掩码与诊断。"""
    lo, hi = coords.min(0), coords.max(0)
    span = hi - lo
    ax = int(np.argmax(span))
    axis = np.zeros(3)
    axis[ax] = 1.0
    t = cent @ axis
    tmin, tmax = float(t.min()), float(t.max())
    L = tmax - tmin
    cut = end_fraction * L
    ends = (t <= tmin + cut) | (t >= tmax - cut)
    return dict(load_axis="XYZ"[ax], span_mm=[float(x) for x in span],
                t_min_mm=tmin, t_max_mm=tmax, L_axis_mm=L, end_fraction=end_fraction,
                ends_mask=ends, shaft_mask=~ends,
                n_ends=int(ends.sum()), n_shaft=int((~ends).sum()))


# ---------------------------------------------------------------------------
# S4 统计（scipy；statsmodels/sklearn 未安装）
# ---------------------------------------------------------------------------
def logistic_fit(x: np.ndarray, y: np.ndarray) -> dict:
    """单变量二元 Logistic（高度→骨折）。检测完全分离；否则 MLE + Wald CI。

    Returns dict with status in {all_zero, all_one, complete_separation, ok}.
    """
    from scipy import optimize, stats
    x = np.asarray(x, dtype=float)
    y = np.asarray(y).astype(int)
    n1 = int(y.sum())
    out: dict = {"n": int(len(y)), "n_fracture": n1, "n_safe": int(len(y) - n1)}
    if n1 == 0:
        out["status"] = "all_zero"
        out["note"] = "区间内从不骨折（常值 0）"
        return out
    if n1 == len(y):
        out["status"] = "all_one"
        out["note"] = "区间内始终骨折（常值 1）"
        return out
    # 完全分离：所有 0 都在某个高度以下、所有 1 都在其上（或反向）
    if x[y == 1].min() > x[y == 0].max() or x[y == 0].min() > x[y == 1].max():
        thr = float(x[y == 0].max())
        out["status"] = "complete_separation"
        out["separation_threshold_m"] = thr
        out["note"] = (f"存在完美阈值 {thr:g} m 把 0/1 完全分开 → MLE 发散、"
                       "二元 Logistic 拟合无统计学意义（与论文“完全分离”一致）")
        return out
    X = np.column_stack([np.ones_like(x), x])

    def nll(b):
        z = X @ b
        p = 1.0 / (1.0 + np.exp(-z))
        eps = 1e-12
        return -float(np.sum(y * np.log(p + eps) + (1 - y) * np.log(1 - p + eps)))

    r = optimize.minimize(nll, [0.0, 0.0], method="BFGS")
    b = r.x
    p = 1.0 / (1.0 + np.exp(-(X @ b)))
    w = p * (1.0 - p)
    H = X.T @ (w[:, None] * X)
    try:
        cov = np.linalg.inv(H)
        se1 = float(np.sqrt(cov[1, 1]))
    except np.linalg.LinAlgError:
        se1 = float("nan")
    if not np.isfinite(se1) or se1 <= 0:
        out["status"] = "complete_separation"
        out["note"] = "Hessian 奇异 → 完全分离/拟完全分离"
        return out
    z = b[1] / se1
    out.update(status="ok", beta=float(b[1]), se=se1, OR=float(np.exp(b[1])),
               ci95=[float(np.exp(b[1] - 1.96 * se1)),
                     float(np.exp(b[1] + 1.96 * se1))],
               wald_z=float(z), p=float(2.0 * stats.norm.sf(abs(z))),
               converged=bool(r.success))
    return out


def jaccard_clusters(matrix: np.ndarray, n_clusters: int = 5) -> dict:
    """50×9 二元矩阵：Jaccard 距离 + 平均联结层次聚类，切成 n_clusters。"""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist
    d = pdist(np.asarray(matrix, dtype=float), metric="jaccard")
    d = np.nan_to_num(d, nan=0.0)          # scipy 对全 0 行已给 0，此处再兜底
    Z = linkage(d, method="average")
    labels = fcluster(Z, t=n_clusters, criterion="maxclust")
    return dict(distance="jaccard", linkage="average", n_clusters=int(n_clusters),
                labels=[int(v) for v in labels], linkage_matrix=Z.tolist())


def cluster_ranges(labels: list[int]) -> list[dict]:
    h = np.arange(1, len(labels) + 1)
    out = []
    for c in sorted(set(labels)):
        hs = h[np.array(labels) == c]
        out.append({"cluster": int(c), "heights": [int(x) for x in hs],
                    "range_m": [int(hs.min()), int(hs.max())], "n": int(len(hs))})
    return out


def _phase_label(fr: list[str]) -> str:
    if not fr:
        return "局部耗散(无骨折)"
    lower = {"足部", "胫骨", "腓骨", "股骨"}
    if set(fr) <= lower:
        return "局部耗散(仅下肢)"
    high = {"骨盆", "颅骨"} & set(fr)
    if high:
        return "全身复合"
    return "轴向传导"


def three_phases(labels, matrix: np.ndarray) -> list[dict]:
    """由 5 类聚类结果重构 3 阶段（>3 类时按相邻顺序并类），并给出各阶段骨折部位集合。"""
    labels = np.asarray(labels)
    h = np.arange(1, len(labels) + 1)
    order = sorted(set(int(v) for v in labels), key=lambda c: h[labels == c].min())
    # 并成至多 3 组（等分相邻类）
    n_grp = min(3, len(order))
    groups = [order[i::n_grp] for i in range(n_grp)]
    groups = [g for g in groups if g]
    out = []
    for i, g in enumerate(groups, 1):
        mask = np.isin(labels, g)
        hs = h[mask]
        fr = [REGION_NAMES[j] for j in range(matrix.shape[1]) if matrix[mask, j].any()]
        out.append({"phase": i, "label": _phase_label(fr),
                    "range_m": [int(hs.min()), int(hs.max())],
                    "heights": [int(x) for x in hs], "fractured_regions": fr})
    return out


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def build() -> dict:
    fm = json.loads(FM_JSON.read_text(encoding="utf-8"))
    bending = json.loads(BENDING_JSON.read_text(encoding="utf-8"))
    r1 = json.loads(R1_JSON.read_text(encoding="utf-8"))
    heights = [float(h) for h in fm["heights_m"]]
    loads = {k: np.asarray(v, dtype=float) for k, v in fm["loads_n"].items()}
    parts = {p["bone"]: p for p in fm["parts"]}
    r1_rows = {r["bone"]: r for r in r1["rows"]}

    bend_rows = {r["bone"]: r for r in bending["rows_primary_beta4"]}

    # ---- σ_total(h) per bone（1D 弯曲感知）--------------------------------
    sigma_h_1d = {b: np.asarray(r["risk_by_height"], dtype=float) * float(r["sigma_c_mpa"])
                  for b, r in bend_rows.items()}

    # ---- 子区域强度表 ----------------------------------------------------
    strengths = {}
    for cn, bone, joint, subs in REGIONS:
        strengths[cn] = {
            "bone": bone, "joint": joint,
            "subregions": {label: {"key": key, "sigma_c_mpa": _sc(key)}
                           for label, key in subs},
            "weak_key": subs[0][1], "weak_sigma_c_mpa": _sc(subs[0][1]),
            "shaft_key": (subs[1][1] if len(subs) > 1 else None),
            "shaft_sigma_c_mpa": (_sc(subs[1][1]) if len(subs) > 1 else None),
        }

    # ---- Route-1D 矩阵 ---------------------------------------------------
    m1d = {}
    first1d = {}
    shaft1d = {}
    first_shaft1d = {}
    for cn, bone, joint, subs in REGIONS:
        sig = sigma_h_1d[bone]
        weak = _sc(subs[0][1])
        m1d[cn] = sig >= weak
        first1d[cn] = _first_fracture(heights, m1d[cn])
        if len(subs) > 1:
            sh_key = subs[1][1]
            shaft1d[cn] = sig >= _sc(sh_key)
            first_shaft1d[cn] = _first_fracture(heights, shaft1d[cn])
            # 弱(两端)子区域单独
    first_ends1d = {cn: first1d[cn] for cn in REGION_NAMES}

    # ---- Route-FE 矩阵 ---------------------------------------------------
    fe_rule = {}
    fe_subs = {}
    fe_sub_first = {}
    mfe = {}
    firstfe = {}
    shaftfe = {}
    first_shaftfe = {}
    for cn, bone, joint, subs in REGIONS:
        weak = _sc(subs[0][1])
        if bone in LONG_BONES:
            dom = r1_rows[bone]["solid_dom"]
            d = read_h5_domain(bone, dom)
            part = partition_ends_shaft(d["coords"], d["cent"])
            vm = vm_voigt(d["stress"])
            p95_ends = float(np.percentile(vm[part["ends_mask"]], 95))
            p95_shaft = float(np.percentile(vm[part["shaft_mask"]], 95))
            reduc = float(parts[bone]["load_reduction_factor"])
            fe_solve = float(parts[bone]["fe_solve_load_n"])
            scale = loads[joint] / fe_solve               # σ(h) = σ_solve · load(h)/fe_solve
            vm_ends_h = p95_ends * scale
            vm_shaft_h = p95_shaft * scale
            sig_ends = _sc(subs[0][1])
            sig_shaft = _sc(subs[1][1])
            mfe[cn] = (vm_ends_h >= sig_ends) | (vm_shaft_h >= sig_shaft)
            shaftfe[cn] = vm_shaft_h >= sig_shaft
            first_shaftfe[cn] = _first_fracture(heights, shaftfe[cn])
            fe_subs[cn] = {
                "domain": dom,
                "ends": {"sigma_c_mpa": sig_ends, "p95_at_solve_mpa": p95_ends,
                         "p95_at_ref_mpa": p95_ends * reduc,
                         "first_fracture_h": _first_fracture(heights, vm_ends_h >= sig_ends)},
                "shaft": {"sigma_c_mpa": sig_shaft, "p95_at_solve_mpa": p95_shaft,
                          "p95_at_ref_mpa": p95_shaft * reduc,
                          "first_fracture_h": first_shaftfe[cn]},
                "n_ends": part["n_ends"], "n_shaft": part["n_shaft"],
            }
            fe_rule[cn] = {k: part[k] for k in
                           ("load_axis", "span_mm", "t_min_mm", "t_max_mm",
                            "L_axis_mm", "end_fraction", "n_ends", "n_shaft")}
        else:
            vm_h = np.asarray(fm["sigma_vm_mpa"]["p95"][bone], dtype=float)
            mfe[cn] = vm_h >= weak
            fe_subs[cn] = {"domain": parts[bone]["main_domain"],
                           "whole": {"sigma_c_mpa": weak,
                                     "p95_at_ref_mpa": float(fm["sigma_vm_mpa"]["p95"][bone][4]),
                                     "first_fracture_h": _first_fracture(heights, mfe[cn])}}
        firstfe[cn] = _first_fracture(heights, mfe[cn])

    # ---- S4 —— Logistic（1D 主口径 & FE 参照）---------------------------
    x = np.asarray(heights, dtype=float)
    logistic_1d = {cn: logistic_fit(x, m1d[cn]) for cn in REGION_NAMES}
    logistic_fe = {cn: logistic_fit(x, mfe[cn]) for cn in REGION_NAMES}

    # ---- S4 —— 聚类（50×9，行=高度，列=部位）--------------------------
    M1d = np.array([[int(m1d[cn][i]) for cn in REGION_NAMES] for i in range(len(heights))])
    Mfe = np.array([[int(mfe[cn][i]) for cn in REGION_NAMES] for i in range(len(heights))])
    clus = jaccard_clusters(M1d, 5)
    clus_ranges = cluster_ranges(clus["labels"])
    phases = three_phases(clus["labels"], M1d)
    clus_fe = jaccard_clusters(Mfe, 5)

    # ---- 排序（1D：早骨折优先，其次 σ_total@5m）----------------------
    def _rank(first_map, sig5):
        arr = sorted(REGION_NAMES,
                     key=lambda cn: (first_map[cn] if first_map[cn] is not None else 1e9,
                                     -sig5[cn]))
        return [{"rank": i, "region": cn, "first_fracture_h": first_map[cn],
                 "sigma_total_at_5m_mpa": sig5[cn]} for i, cn in enumerate(arr, 1)]

    sig5_1d = {cn: float(sigma_h_1d[bone][4]) for cn, bone, _, _ in REGIONS}
    rank_1d = _rank(first1d, sig5_1d)
    sig5_fe = {cn: float(np.asarray(fm["sigma_vm_mpa"]["p95"][bone])[4])
               for cn, bone, _, _ in REGIONS}
    rank_fe = _rank(firstfe, sig5_fe)

    # ---- 论文对齐清单 ----------------------------------------------------
    checklist = _checklist(logistic_1d, first1d, first_shaft1d, first_shaftfe,
                           clus_ranges, phases)

    result = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "method": ("把 9 部位拆成 [Y25] 子区域：弱子区域=两端/股骨颈(70/80)决定部位骨折，"
                   "骨干子区域(200/160/220)单独验证“骨干不骨折”。Route-1D 用 "
                   "risk_1d_bending.json σ_total(h)；Route-FE 逐单元按沿载荷轴 bbox 最长轴"
                   "位置分区（ends=两端各外 25%，shaft=中间 50%），取子区域 p95 后线性缩放。"),
        "paper_facts": {
            "regions_9": REGION_NAMES,
            "table1_compression_mpa": {k: _sc(k) for k in
                                       ("calcaneus", "tibia_shaft", "tibia_ends",
                                        "fibula_shaft", "fibula_ends", "femur_shaft",
                                        "femoral_neck", "pelvis", "spine", "skull")},
            "criterion": "部位骨折 ⇔ 其任一子区域超过该子区域阈值",
            "shafts_never_fracture": ["tibia_shaft", "fibula_shaft", "femur_shaft"],
            "paper_step": sorted(PAPER_STEP),
            "paper_gradual": sorted(PAPER_GRADUAL),
            "paper_or_order": PAPER_OR_ORDER,
            "paper_cluster_bounds_m": PAPER_CLUSTER_BOUNDS,
        },
        "subregion_strengths_used": strengths,
        "fe_partition_rule": {
            "load_axis": "主承力域节点 bbox 最长轴（与 risk_1d_bending θ 一致）",
            "ends": f"沿载荷轴两端各外 {END_FRACTION:.0%}（共 50% 跨度）",
            "shaft": "中间 50%",
            "statistic": "子区域内单元 von Mises p95，再按 σ(h)=σ_solve·load(h)/fe_solve 线性缩放",
            "per_bone": fe_rule,
            "note": "仅 3 根长骨切分 ends/shaft；足(紧凑骨)、骨盆/脊柱/颅骨(壳,单一阈值)用整域 p95",
        },
        "heights_m": [int(h) for h in heights],
        "matrix_1d": {cn: [int(v) for v in m1d[cn]] for cn in REGION_NAMES},
        "matrix_fe": {cn: [int(v) for v in mfe[cn]] for cn in REGION_NAMES},
        "first_fracture_1d_m": first1d,
        "first_fracture_fe_m": firstfe,
        "shaft_1d": {cn: [int(v) for v in shaft1d[cn]] for cn in shaft1d},
        "shaft_first_fracture_1d_m": first_shaft1d,
        "shaft_first_fracture_fe_m": first_shaftfe,
        "fe_subregions": fe_subs,
        "ranking_1d": rank_1d,
        "ranking_fe": rank_fe,
        "s4": {
            "library": "scipy（statsmodels / sklearn 未安装，已按任务要求用 scipy 实现）",
            "logistic_1d": logistic_1d,
            "logistic_fe": logistic_fe,
            "clustering_1d": {**clus, "cluster_ranges": clus_ranges, "three_phases": phases},
            "clustering_fe": {**clus_fe, "cluster_ranges": cluster_ranges(clus_fe["labels"])},
        },
        "checklist": checklist,
        "assumptions": [
            "σ_c 取自 src/climbing/bone.py MATERIAL_STRENGTH_MPA（[Y25] Table 1 压缩）—— 材料强度，"
            "非整体骨失效载荷；绝对值不可信，相对排序可用。",
            "Route-1D σ_total 来自 risk_1d_bending.json（梁弯曲上界估计）；同一条 σ(h) 用于两端与骨干，"
            "仅阈值不同（任务指定）。",
            "Route-FE 收缩：主域 p95；ends/shaft 仅对 3 长骨按几何位置切分，其余用整域。",
            "FE 应力按严格线弹性线性缩放；未重跑 FEBio/OpenSim。",
            "论文的双传导路径(d)是**事件内时序**，50×9 高度矩阵不含时序 → 不可直接计算，仅定性对照。",
        ],
    }
    return result


def _checklist(logi1d, first1d, first_shaft1d, first_shaftfe, clus_ranges, phases) -> list[dict]:
    rows = []

    # (a) 阶跃 vs 渐进
    def _behav(rec):
        s = rec["status"]
        if s in ("all_zero", "all_one"):
            return "常值"
        if s == "complete_separation":
            return "阶跃"
        if s == "ok":
            return "渐进"
        return s

    for cn, bone, _, subs in REGIONS:
        st = logi1d[cn]["status"]
        b = _behav(logi1d[cn])
        target = "阶跃" if cn in PAPER_STEP else "渐进"
        if b == target:
            verdict = "✓"
        elif target == "渐进":
            verdict = "✗"          # 论文渐变，我方却恒骨折/恒不骨折 → 未复现
        else:
            verdict = "~"          # 论文阶跃，我方恒值（阈值落在区间外）
        if b == "阶跃":
            note = "1D 单调增 ⇒ 非骨折即完全分离"
        elif st == "all_one":
            note = "1D 全高度已骨折（弯曲上界估计高估，1 m 即超阈值）"
        elif st == "all_zero":
            note = "1D 区间内从不骨折（阈值未达）"
        else:
            note = "1D Logistic 拟合成功"
        rows.append({
            "id": "a", "region": cn, "paper": f"{target}响应",
            "ours": f"1D {b}（{st}）", "verdict": verdict, "note": note,
        })

    # (b) OR 排序
    ors = {cn: logi1d[cn].get("OR") for cn in PAPER_OR_ORDER}
    fit_ok = [cn for cn in PAPER_OR_ORDER if ors[cn] is not None]
    if len(fit_ok) == 3:
        order = sorted(PAPER_OR_ORDER, key=lambda c: -ors[c])
        verdict = "✓" if order == PAPER_OR_ORDER else "✗"
        note = "OR " + ", ".join(f"{c}={ors[c]:.3f}" for c in PAPER_OR_ORDER)
    else:
        verdict = "✗"
        note = ("1D 上腓骨/骨盆/颅骨均因完全分离或常值无法取得 OR → OR 排序不可计算；"
                "论文的渐进关系在我们重映射后的矩阵中未复现")
    rows.append({"id": "b", "region": "腓骨>颅骨>骨盆 OR 排序",
                 "paper": "OR 腓骨 1.682 > 颅骨 1.576 > 骨盆 1.236",
                 "ours": note, "verdict": verdict, "note": ""})

    # (c) 骨干不骨折
    shaft_fr_1d = {cn: v for cn, v in first_shaft1d.items() if v is not None}
    shaft_fr_fe = {cn: v for cn, v in first_shaftfe.items() if v is not None}
    verdict = "✓" if (not shaft_fr_1d and not shaft_fr_fe) else "✗"
    rows.append({"id": "c", "region": "胫/腓/股骨干从不骨折", "paper": "骨干均未骨折",
                 "ours": (f"1D 骨干首骨折 {shaft_fr_1d or '无'}；"
                          f"FE 骨干首骨折 {shaft_fr_fe or '无'}"),
                 "verdict": verdict,
                 "note": "弯曲 1D 为全骨上界，σ_total 同时用于骨干判据 → 骨干被高估" if verdict == "✗" else ""})

    # (d) 双传导路径
    rows.append({"id": "d", "region": "2 条应力传导路径（自下而上/自上而下）",
                 "paper": "足→胫腓→股骨/骨盆（上）与颅骨→颈→胸（下）",
                 "ours": "50×9 高度矩阵不含事件内时序，无法计算路径先后",
                 "verdict": "—",
                 "note": "仅定性：低高度先下肢、高高度全身，与“上下两路汇于脊柱”不矛盾"})

    # (e) 5 类 / 3 阶段
    cr = "、".join(f"#{r['cluster']}:{r['range_m'][0]}-{r['range_m'][1]}m" for r in clus_ranges)
    ph = "、".join(f"{p['phase']}:{p['range_m'][0]}-{p['range_m'][1]}m" for p in phases)
    n_clus = len(clus_ranges)
    verdict = "✓" if n_clus == 5 else ("~" if n_clus >= 3 else "✗")
    rows.append({"id": "e", "region": "5 类 / 3 阶段",
                 "paper": "5 类 1-6/7-9/10-16/17-34/35-50；3 阶段 <7 / 7-16 / >16 m",
                 "ours": f"{n_clus} 类 [{cr}]；3 阶段 [{ph}]",
                 "verdict": verdict,
                 "note": "重映射后低高度即骨折，矩阵过早饱和 → 类数少于论文" if n_clus != 5 else ""})
    return rows


# ---------------------------------------------------------------------------
def _rows_csv(path: Path, mat: dict, heights: list[int], first: dict) -> None:
    names = list(mat.keys())
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["height_m"] + names)
        for i, h in enumerate(heights):
            w.writerow([h] + [mat[n][i] for n in names])
        w.writerow(["first_fracture_m"] + [first.get(n) for n in names])


def _fmt(v) -> str:
    return ">50" if v is None else f"{int(v)}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(RES))
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    res = build()
    (out / "paper_alignment.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")

    m1d = res["matrix_1d"]
    mfe = res["matrix_fe"]
    (out / "paper_alignment_1d_matrix.json").write_text(
        json.dumps({"heights_m": res["heights_m"], "regions": REGION_NAMES,
                    "matrix": m1d, "first_fracture_m": res["first_fracture_1d_m"],
                    "shaft_matrix": res["shaft_1d"],
                    "shaft_first_fracture_m": res["shaft_first_fracture_1d_m"]},
                   indent=2, ensure_ascii=False), encoding="utf-8")
    (out / "paper_alignment_fe_matrix.json").write_text(
        json.dumps({"heights_m": res["heights_m"], "regions": REGION_NAMES,
                    "matrix": mfe, "first_fracture_m": res["first_fracture_fe_m"],
                    "shaft_first_fracture_m": res["shaft_first_fracture_fe_m"],
                    "partition_rule": res["fe_partition_rule"]},
                   indent=2, ensure_ascii=False), encoding="utf-8")
    _rows_csv(out / "paper_alignment_1d_matrix.csv", m1d, res["heights_m"],
              res["first_fracture_1d_m"])
    _rows_csv(out / "paper_alignment_fe_matrix.csv", mfe, res["heights_m"],
              res["first_fracture_fe_m"])

    _write_report(out / "PAPER_ALIGNMENT_REPORT.md", res)

    # console
    print("=" * 96)
    print("子区域强度 (压缩 MPa)：", {k: v["weak_sigma_c_mpa"] for k, v in
                                       res["subregion_strengths_used"].items()})
    print("\nRoute-1D 排序（早骨折优先）:")
    for r in res["ranking_1d"]:
        print(f"  #{r['rank']} {r['region']:4s} first={_fmt(r['first_fracture_h'])}m "
              f"σ_total@5m={r['sigma_total_at_5m_mpa']:8.1f}")
    print("1D 骨干首骨折:", res["shaft_first_fracture_1d_m"] or "无")
    print("FE 骨干首骨折:", res["shaft_first_fracture_fe_m"] or "无")
    print("1D 聚类:", [(x["cluster"], x["range_m"]) for x in
                       res["s4"]["clustering_1d"]["cluster_ranges"]])
    print("1D 3 阶段:", [(x["phase"], x["range_m"]) for x in
                         res["s4"]["clustering_1d"]["three_phases"]])
    print("\n论文对齐清单:")
    for row in res["checklist"]:
        print(f"  [{row['id']}] {row['verdict']} {row['region']}: {row['ours']}")
    print("=" * 96)
    for f in ("paper_alignment.json", "paper_alignment_1d_matrix.json",
              "paper_alignment_fe_matrix.json", "paper_alignment_1d_matrix.csv",
              "paper_alignment_fe_matrix.csv", "PAPER_ALIGNMENT_REPORT.md"):
        print("wrote:", out / f)
    return 0


def _write_report(path: Path, res: dict) -> None:
    L = []
    A = L.append
    ts = res["generated_at"]
    st = res["subregion_strengths_used"]
    A("# 论文对齐报告 —— [Y25] 子区域（两端 / 股骨颈 vs 骨干）重映射")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 本报告**新增**：不改既有模块/默认值，不覆盖既有产物，不调用 FEBio。")
    A("> 修复对象：此前 1D/FE 用**骨干**强度（胫 200 / 腓 160 / 股 220）评估"
      "胫/腓/股 —— 即**不会骨折的子区域**；论文骨折部位实为**两端(70)与股骨颈(80)**。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    A("把 σ_c 换成论文的**弱子区域**后：")
    A("")
    fr = res["first_fracture_1d_m"]
    A(f"- Route-1D 中 **足部/胫骨/腓骨/股骨在 1 m 即骨折**（两端/颈阈值 70/80 MPa 远低于旧骨干值），"
      f"腰椎 {_fmt(fr.get('腰椎'))} m、胸椎 {_fmt(fr.get('胸椎'))} m；"
      f"骨盆/颈椎/颅骨在 1–50 m 内 1D 不骨折。")
    A("- **骨干仍被高估**：弯曲 1D 的同一条 σ_total 用于骨干判据时，胫/腓/股骨干在 1 m 即超"
      "骨干阈值 → 与论文“骨干从不骨折”**不符**（原因：弯曲上界估计 + 同应力用于两子区域）。")
    A("- Logistic 在 1D 上：可骨折部位全部落在**完全分离**（阶跃），无 OR → 论文的"
      "“腓骨/骨盆/颅骨渐进 + OR 排序”**未复现**。")
    A("- Jaccard 聚类：矩阵低高度即饱和，类数少于论文 5 类 → 5 类/3 阶段**仅部分对齐**。")
    A("- **诚实定位**：本对齐证实了「子区域强度选择」这一 bug 的影响面，但也暴露我方两条路线"
      "在低高度**过度骨折**，尚不能复现论文的渐进部位与聚类结构。")
    A("")
    A("---")
    A("")
    A("## 1. 子区域强度与分区规则")
    A("")
    A("### 1.1 子区域 → 压缩强度（[Y25] Table 1，取自 `bone.py::MATERIAL_STRENGTH_MPA`）")
    A("")
    A("| 部位 | 代表骨 | 子区域 | σ_c (MPa) | 是否决定部位骨折 |")
    A("|---|---|---|---:|---|")
    for cn, _, _, _ in REGIONS:
        e = st[cn]
        for i, (label, meta) in enumerate(e["subregions"].items()):
            A(f"| {cn} | `{e['bone']}` | {label} | {meta['sigma_c_mpa']:.0f} | "
              f"{'**是（弱子区域）**' if i == 0 else '否（验证骨干）'} |")

    A("")
    A("**判据**：部位骨折 ⇔ 其**任一**子区域超过该子区域阈值。骨干子区域单独评估，"
      "用于验证论文结论「胫/腓/股骨干从不骨折」。")
    A("")
    A("### 1.2 Route-FE ends/shaft 分区规则（精确可复现）")
    A("")
    A("1. 载荷轴 = 主承力域节点坐标 bbox **最长轴**（与 `risk_1d_bending.py` θ 定义一致）。")
    A("2. 单元形心 `c_e` = 其节点坐标均值；沿轴投影 `t_e = c_e·axis`。")
    A(f"3. `L = t_max − t_min`；**ends** = `t ≤ t_min + {END_FRACTION:.0%}·L` 或 "
      f"`t ≥ t_max − {END_FRACTION:.0%}·L`（两端各外 {END_FRACTION:.0%}，共 50% 跨度）；"
      "**shaft** = 中间 50%。")
    A("4. 子区域统计量 = 该子区域内单元 von Mises **p95**；再按 σ(h)=σ_solve·load(h)/fe_solve "
      "严格线弹性缩放到各高度。")
    A("")
    A("> 仅 3 根长骨切分 ends/shaft；足（紧凑骨）与骨盆/脊柱/颅骨（壳、论文给单一阈值）用整域 p95。")
    A("")
    A("**逐骨分区诊断**（长骨）：")
    A("")
    A("| 部位 | 骨 | 载荷轴 | L_axis (mm) | n_ends | n_shaft | ends p95@ref (MPa) | shaft p95@ref (MPa) |")
    A("|---|---|---|---:|---:|---:|---:|---:|")
    for r in res["fe_partition_rule"]["per_bone"]:
        d = res["fe_partition_rule"]["per_bone"][r]
        s = res["fe_subregions"][r]
        A(f"| {r} | `{st[r]['bone']}` | {d['load_axis']} | {d['L_axis_mm']:.1f} | "
          f"{d['n_ends']} | {d['n_shaft']} | {s['ends']['p95_at_ref_mpa']:.1f} | "
          f"{s['shaft']['p95_at_ref_mpa']:.1f} |")
    A("")
    A("### 1.3 关键量：measured vs assumed")
    A("")
    A("| 量 | 来源 | 性质 |")
    A("|---|---|---|")
    A("| σ_c（各子区域 70/80/150/160/180/200/220）| `bone.py::MATERIAL_STRENGTH_MPA`（[Y25] Table 1）| **assumed**（材料强度，非整体骨失效载荷）|")
    A("| Route-1D σ_total(h) | `risk_1d_bending.json`（梁弯曲）| **derived**（F/A/L measured + 模型假设）|")
    A("| 关节反力 load(h) | `fracture_matrix.json::loads_n` | **measured**（OpenSim 逐高度载荷链）|")
    A("| FE 逐单元应力 | `temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5` | **measured**（FE 解）|")
    A("| FE 应力随高度缩放 | σ(h)=σ_solve·load(h)/fe_solve | **assumed**（严格线弹性）|")
    A("| ends/shaft 分区（各 25%）| 本脚本几何规则 | **assumed**（规则选择，非论文给出）|")
    A("| 聚类/Logistic | scipy | **derived**（本脚本计算）|")
    A("")
    A("---")
    A("")
    A("## 2. Route-1D（主口径）：子区域矩阵")
    A("")
    A("σ_total(h) 取自 `risk_1d_bending.json`（`risk_by_height × σ_c`）。"
      "部位骨折 = σ_total(h) ≥ σ_c(弱子区域)。")
    A("")
    A("**排序（早骨折优先，其次 σ_total@5m）：**")
    A("")
    A("| # | 部位 | 首次骨折 h (m) | σ_total@5m (MPa) | 弱子区域 σ_c (MPa) |")
    A("|---:|---|---:|---:|---:|")
    for r in res["ranking_1d"]:
        cn = r["region"]
        A(f"| {r['rank']} | {cn} | **{_fmt(r['first_fracture_h'])}** | "
          f"{r['sigma_total_at_5m_mpa']:.1f} | {st[cn]['weak_sigma_c_mpa']:.0f} |")
    A("")
    # 1D matrix table (compact: first fracture + per-10m)
    A("**1D 9×50 二值矩阵（节选：h=1,5,10,20,30,50；完整见 `paper_alignment_1d_matrix.csv`）：**")
    A("")
    cols = [0, 4, 9, 19, 29, 49]
    hdr = " | ".join(["部位"] + [f"{res['heights_m'][i]}m" for i in cols] + ["首骨折"])
    A("| " + hdr + " |")
    A("|" + "---|" * (len(cols) + 2))
    for cn in REGION_NAMES:
        vals = " | ".join(str(res["matrix_1d"][cn][i]) for i in cols)
        A(f"| {cn} | {vals} | **{_fmt(res['first_fracture_1d_m'][cn])}** |")
    A("")
    A("**骨干子区域单独评估（验证“骨干不骨折”）：**")
    A("")
    A("| 部位 | 骨干 σ_c (MPa) | 骨干首次骨折 h (m) |")
    A("|---|---:|---:|")
    for cn in ("胫骨", "腓骨", "股骨"):
        A(f"| {cn} | {st[cn]['shaft_sigma_c_mpa']:.0f} | "
          f"**{_fmt(res['shaft_first_fracture_1d_m'].get(cn))}** |")
    A("")
    A("---")
    A("")
    A("## 3. Route-FE（交叉验证）：子区域矩阵")
    A("")
    A("| 部位 | 主域 | 首次骨折 h (m) FE | 首次骨折 h (m) 1D |")
    A("|---|---|---:|---:|")
    for cn in REGION_NAMES:
        dom = res["fe_subregions"][cn]["domain"]
        A(f"| {cn} | `{dom}` | **{_fmt(res['first_fracture_fe_m'][cn])}** | "
          f"{_fmt(res['first_fracture_1d_m'][cn])} |")
    A("")
    A("**FE 长骨 ends/shaft 子区域首骨折：**")
    A("")
    A("| 部位 | ends σ_c | ends 首骨折 | shaft σ_c | shaft 首骨折 |")
    A("|---|---:|---:|---:|---:|")
    for cn in ("胫骨", "腓骨", "股骨"):
        s = res["fe_subregions"][cn]
        A(f"| {cn} | {s['ends']['sigma_c_mpa']:.0f} | "
          f"{_fmt(s['ends']['first_fracture_h'])} | {s['shaft']['sigma_c_mpa']:.0f} | "
          f"{_fmt(s['shaft']['first_fracture_h'])} |")
    A("")
    A("> FE 应力幅值高（长骨 @5m 利用率远超 1），换成低阈值后 ends 一律在 **1 m** 骨折；"
      "shaft 亦在 1 m 骨折（与论文不符，见 §6 边界）。")
    A("")
    A("---")
    A("")
    A("## 4. S4 统计")
    A("")
    A(f"统计库：{res['s4']['library']}")
    A("")
    A("### 4.1 Logistic 回归（骨折 ~ 高度，1D 主口径）")
    A("")
    A("| 部位 | 状态 | β | OR | 95% CI | P | 0 区间上界阈值 (m) |")
    A("|---|---|---:|---:|---|---:|---:|")
    for cn in REGION_NAMES:
        rec = res["s4"]["logistic_1d"][cn]
        s = rec["status"]
        if s == "ok":
            A(f"| {cn} | ok | {rec['beta']:.3f} | {rec['OR']:.3f} | "
              f"{rec['ci95'][0]:.3f}~{rec['ci95'][1]:.3f} | {rec['p']:.4f} | — |")
        else:
            thr = rec.get("separation_threshold_m")
            A(f"| {cn} | **{s}** | — | — | — | — | "
              f"{('%.0f' % thr) if thr is not None else '—'} |")
    A("")
    A("> 若某部位在区间内单调「从不骨折→骨折」，则存在完美阈值 → **完全分离**，"
      "MLE 发散、拟合无统计学意义（与论文“阶跃部位拟合失败”同类）。")
    A("")
    A("### 4.2 层次聚类（50×9，Jaccard 距离 + 平均联结，切 5 类）")
    A("")
    A("| 类 | 高度范围 (m) | n |")
    A("|---|---:|---:|")
    for r in res["s4"]["clustering_1d"]["cluster_ranges"]:
        A(f"| #{r['cluster']} | {r['range_m'][0]}–{r['range_m'][1]} | {r['n']} |")
    A("")
    A("**数据驱动 3 阶段（由聚类结果并类）：**")
    A("")
    A("| 阶段 | 标签 | 高度范围 (m) | 已骨折部位 |")
    A("|---|---|---:|---|")
    for p in res["s4"]["clustering_1d"]["three_phases"]:
        A(f"| {p['phase']} | {p['label']} | {p['range_m'][0]}–{p['range_m'][1]} | "
          f"{'、'.join(p['fractured_regions']) or '无'} |")
    A("")
    A("> 说明：论文 5 类为 1-6 / 7-9 / 10-16 / 17-34 / 35-50 m，3 阶段 <7 / 7-16 / >16 m。")
    A("> 我方重映射后矩阵低高度即饱和，Jaccard 聚类只得到 "
      f"{len(res['s4']['clustering_1d']['cluster_ranges'])} 类，"
      "故 3 阶段退化为上述区间（属真实不对齐，非计算失败）。")
    A("")
    A("---")
    A("")
    A("## 5. 论文对齐清单")
    A("")
    A("| # | 部位/条目 | 论文结论 | 我方结果 | 判定 | 说明 |")
    A("|---|---|---|---|---|---|")
    for row in res["checklist"]:
        A(f"| {row['id']} | {row['region']} | {row['paper']} | {row['ours']} | "
          f"**{row['verdict']}** | {row['note']} |")
    A("")
    A("> 判定：✓ 对齐；~ 部分/不可判定；✗ 不对齐；— 不可计算。")
    A("")
    A("---")
    A("")
    A("## 6. 假设与诚实边界")
    A("")
    for i, a in enumerate(res["assumptions"], 1):
        A(f"{i}. {a}")
    A("")
    A("- **(d) 双传导路径**：论文描述的是**单次坠落事件内的应力时序**（足部/颅骨先，"
      "两路汇于脊柱）。我们的 9×50 高度矩阵**不含事件内时序**，无法计算路径先后 —— "
      "仅能定性说「低高度先下肢、高高度全身」，与「上下两路汇于脊柱」不矛盾。")
    A("- **(e) 5 类/3 阶段**：重映射后低高度即出现多部位骨折，矩阵过早饱和，"
      "Jaccard 聚类类数少于论文 5 类 —— 属**真实不对齐**，非计算失败。")
    A("- **骨干被高估**：弯曲 1D 的 σ_total 是**全骨上界**（等效圆 + 全骨长 + 简单支承），"
      "论文的骨干不骨折在本口径下无法成立；FE 的 p95 又被边界条件伪影污染"
      "（`BC_ROBUST_METRIC_REPORT.md`）。")
    A("")
    A("---")
    A("")
    A("## 7. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\paper_alignment.py")
    A("# 只读复用：")
    A("#   results/opensim_fe/risk_1d_bending.json        (σ_total(h))")
    A("#   results/opensim_fe/fracture_matrix.json        (loads_n / FE util / p95)")
    A("#   results/opensim_fe/bc_robust_metric_route1.json(solid_dom)")
    A("#   src/climbing/bone.py MATERIAL_STRENGTH_MPA     (σ_c)")
    A("#   temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5   (逐单元应力)")
    A("# 产物：")
    A("#   results/opensim_fe/paper_alignment.json")
    A("#   results/opensim_fe/paper_alignment_1d_matrix.{json,csv}")
    A("#   results/opensim_fe/paper_alignment_fe_matrix.{json,csv}")
    A("#   results/opensim_fe/PAPER_ALIGNMENT_REPORT.md")
    A("```")
    A("")
    A("---")
    A("")
    A("## 8. 产物")
    A("")
    A("- 本报告：`results/opensim_fe/PAPER_ALIGNMENT_REPORT.md`")
    A("- 新脚本：`scripts/opensim_fe/paper_alignment.py`（additive）")
    A("- 机器可读：`paper_alignment.json`、`paper_alignment_1d_matrix.{json,csv}`、"
      "`paper_alignment_fe_matrix.{json,csv}`")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
