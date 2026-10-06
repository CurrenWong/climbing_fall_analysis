"""pad_stiffness.py -- S1 场景：地面/垫子刚度 k_pad -> 首折高度 曲线。

背景
----
场景 S1（``docs/论文对齐基线与场景探索计划.md`` §4）问：抱石垫的刚度如何改变
首次骨折高度？物理直觉是**软垫把冲击摊到更长的减速行程 -> 峰值地面/关节反力更低
-> 骨应力更低 -> 首折高度更高**；刚性地面是 ``k -> inf`` 的极限。

冻结基线（``results/opensim_fe/AXIAL_SUBREGION_REPORT.md`` / ``risk_1d_loadshare.json``）
只给了**刚性地面**下的 σ(h)；它**没有**垫子刚度这一维。本脚本做**加法式**后处理，
把两个可辩护的垫子模型贴到冻结基线上，得到 k -> 首折高度：

**路线 M1（主口径，实测转移函数）**
    复用已有 FE 跖面 BC 证据（``PLANTAR_BC_REPORT.md`` / ``plantar_bc_summary.json``）：
    弹性地基（Winkler k）相对 fixed 的 **gauge 应力折减比**
    ``R_FE(k) = gauge_max(spring_k) / gauge_max(fixed)``。
    假设垫子按同一比例折减整条载荷链的峰值应力：
    ``sigma_pad(h; k) = R_FE(k) · sigma_rigid(h)``，首折 = 最小 h 使 sigma_pad >= σ_c(弱子区域)。

**路线 M2（交叉校验，物理模型）**
    把垫子与人体组织视为**串联线性弹簧**：峰值力由能量/冲量-动量给出
    ``F_peak = v0·sqrt(m·k_series)``，``k_series = k_pad·k_body/(k_pad+k_body)``，
    故 ``R_phys(k_pad) = F_pad/F_rigid = sqrt(k_pad/(k_pad+k_body))``。
    用 M2 反推“FE 实测折减比对应多大的物理 k_pad”，并给出量级边界。

只读输入 / 只写新文件；不改既有模块默认/签名；不覆盖既有产物；不调用 FEBio；不重跑 OpenSim。

输入（只读）
-----------
* ``results/opensim_fe/risk_1d_loadshare.json``    -> 纯轴向 risk_by_height + stored σ_c（冻结基线）
* ``results/opensim_fe/fracture_matrix.json``       -> heights_m（1..50）
* ``results/opensim_fe/plantar_bc_summary.json``    -> 实测 spring gauge_max(k) / fixed gauge / FE 外推高度
* ``src/climbing/bone.py::MATERIAL_STRENGTH_MPA``   -> 弱子区域压缩强度
* ``src/climbing/pad.py::hard_surface``             -> k_body 默认（人体等效轴向刚度，1e7 N/m）

产物（新文件）
-------------
* ``results/opensim_fe/pad_stiffness.json``
* ``results/opensim_fe/PAD_STIFFNESS_REPORT.md``
* ``results/opensim_fe/pad_stiffness_curve.png``（曲线，可选）

复现（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\pad_stiffness.py
"""
from __future__ import annotations

import argparse
import json
import math
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
LOADSHARE_JSON = RES / "risk_1d_loadshare.json"
FM_JSON = RES / "fracture_matrix.json"
PLANTAR_JSON = RES / "plantar_bc_summary.json"

#: 承重骨（S1 关注）：(中文名, 代表骨, 弱子区域键)。
#: 弱子区域 = 该部位最早失效的子区域（[Y25] Table 1）。
LOAD_BONES: list[tuple[str, str, str]] = [
    ("足部", "calcaneus_r", "calcaneus"),
    ("腓骨", "fibula_r", "fibula_ends"),
    ("胫骨", "tibia_r", "tibia_ends"),
    ("股骨", "femur_r", "femoral_neck"),
]

#: 人体等效轴向刚度（N/m），取自 ``pad.hard_surface`` 的默认标定值 1.0e7。
K_BODY_N_PER_M = 1.0e7

#: M1 曲线用的 k 网格（N/mm，FE Winkler 刚度）。inf 用 None 表示刚性极限。
K_GRID_N_PER_MM: list[float | None] = [
    10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0,
    10000.0, 30000.0, 100000.0, 1000000.0, None,
]

#: M2 物理路线的 k_pad 网格（N/m），围绕 k_body 量级。
K_PAD_PHYSICS_N_PER_M: list[float] = [1.0e5, 3.0e5, 1.0e6, 3.0e6, 1.0e7, 3.0e7]


def _sc(key: str) -> float:
    """弱子区域压缩强度 (MPa)，[Y25] Table 1（bone.py MATERIAL_STRENGTH_MPA）。"""
    return float(MATERIAL_STRENGTH_MPA[key][1])


def first_fracture(heights: np.ndarray, risk: np.ndarray, risk_threshold: float) -> float | None:
    """最小的高 h 使 ``risk(h) >= risk_threshold``（否则 None 代表 > 扫描上界）。"""
    mask = risk >= risk_threshold
    idx = np.nonzero(mask)[0]
    if idx.size == 0:
        return None
    return float(heights[int(idx[0])])


# ---------------------------------------------------------------------------
# 路线 M1：FE 实测转移函数 R_FE(k)
# ---------------------------------------------------------------------------
def build_fe_transfer(plantar: dict) -> dict:
    """从 plantar_bc_summary.json 构造 R_FE(k) = gauge(spring k)/gauge(fixed)。"""
    fixed_gauge = float(plantar["spring_scan"]["fixed_gauge_max_mpa"])
    ks = [float(k) for k in plantar["spring_scan"]["k_N_per_mm"]]
    gs = [float(g) for g in plantar["spring_scan"]["gauge_max_mpa"]]
    ratios = [g / fixed_gauge for g in gs]
    r_min = min(ratios)
    k_min = ks[int(np.argmin(ratios))]
    return {
        "source": "results/opensim_fe/plantar_bc_summary.json (FE plantar-BC probe)",
        "fixed_gauge_max_mpa": fixed_gauge,
        "k_N_per_mm": ks,
        "gauge_max_mpa": gs,
        "ratio_vs_fixed": ratios,
        "monotone_non_decreasing": bool(plantar["spring_scan"]["monotone_non_decreasing"]),
        "r_min": float(r_min),
        "k_at_r_min_N_per_mm": float(k_min),
        "hardest_k_residual_vs_fixed_pct": float(
            plantar["spring_scan"]["hardest_k_residual_vs_fixed_pct"]
        ),
    }


def r_fe(k_n_per_mm: float | None, table: dict) -> float:
    """R_FE(k)：log10(k) 线性插值；k>=1e6 视为刚性极限 1.0；k<10 平推软端。"""
    if k_n_per_mm is None:
        return 1.0
    k = float(k_n_per_mm)
    ks = np.asarray(table["k_N_per_mm"], dtype=float)
    rs = np.asarray(table["ratio_vs_fixed"], dtype=float)
    if k >= ks[-1]:
        return 1.0
    if k <= ks[0]:
        return float(rs[0])
    return float(np.interp(math.log10(k), np.log10(ks), rs))


# ---------------------------------------------------------------------------
# 路线 M2：串联线性弹簧物理模型
# ---------------------------------------------------------------------------
def r_phys(k_pad_n_per_m: float, k_body_n_per_m: float = K_BODY_N_PER_M) -> float:
    """R_phys(k_pad) = sqrt(k_pad / (k_pad + k_body))。

    推导：峰值力 ``F_peak = v0·sqrt(m·k_series)``（能量 ½mv² = ½k_s x²，F=k_s x），
    串联 ``k_s = k_pad·k_body/(k_pad+k_body)``；对刚性地面 ``R -> 1``，软垫 ``R -> 0``。
    """
    k = max(0.0, float(k_pad_n_per_m))
    if k <= 0.0:
        return 0.0
    return math.sqrt(k / (k + float(k_body_n_per_m)))


def k_pad_equiv(r: float, k_body_n_per_m: float = K_BODY_N_PER_M) -> float:
    """把折减比 R 反解成等效物理 k_pad = k_body·R²/(1−R²)。"""
    rr = float(r)
    if not (0.0 <= rr < 1.0):
        return float("inf")
    return float(k_body_n_per_m) * rr * rr / (1.0 - rr * rr)


def fe_absolute_heights(plantar: dict) -> dict:
    """FE 实测外推：fixed vs 各 spring k 的首折高度（anat 尺度，非冻结 1D 尺度）。"""
    out: dict[str, float | None] = {}
    for key, rec in plantar["first_fracture_extrapolation"].items():
        out[key] = (None if rec.get("h_first_m") is None else float(rec["h_first_m"]))
    return out


# ---------------------------------------------------------------------------
def build() -> dict:
    ls = json.loads(LOADSHARE_JSON.read_text(encoding="utf-8"))
    fm = json.loads(FM_JSON.read_text(encoding="utf-8"))
    plantar = json.loads(PLANTAR_JSON.read_text(encoding="utf-8"))

    heights = np.asarray([float(h) for h in fm["heights_m"]], dtype=float)
    rows = {r["bone"]: r for r in ls["rows_with_share_parallel"]}

    fe_table = build_fe_transfer(plantar)

    # ---- 逐骨：rigid 基线 + M1 曲线 + M2 曲线 --------------------------------
    bone_curves: list[dict] = []
    for cn, bone, weak_key in LOAD_BONES:
        r0 = rows[bone]
        stored_sc = float(r0["sigma_c_mpa"])
        risk = np.asarray(r0["risk_by_height"], dtype=float)
        weak_sc = _sc(weak_key)
        # 冻结基线（R=1）：σ_rigid(h) = risk(h)·stored σ_c；折 ⇔ risk(h) ≥ weak/stored
        ratio_rigid = weak_sc / stored_sc
        h_rigid = first_fracture(heights, risk, ratio_rigid)

        # M1：risk(h) ≥ (weak/stored) / R
        m1: list[dict] = []
        for k in K_GRID_N_PER_MM:
            R = r_fe(k, fe_table)
            thr = ratio_rigid / R if R > 0 else float("inf")
            hf = first_fracture(heights, risk, thr)
            m1.append({
                "k_N_per_mm": (None if k is None else float(k)),
                "R_fe": float(R),
                "risk_threshold": float(thr),
                "first_fracture_h": hf,
                "sigma_at_5m_mpa": float(R * risk[4] * stored_sc),
            })

        # M2（物理）：R_phys(k_pad)；同一公式映射到冻结 σ(h)
        m2: list[dict] = []
        for kp in K_PAD_PHYSICS_N_PER_M:
            R = r_phys(kp)
            thr = ratio_rigid / R if R > 0 else float("inf")
            hf = first_fracture(heights, risk, thr)
            m2.append({
                "k_pad_N_per_m": float(kp),
                "R_phys": float(R),
                "risk_threshold": float(thr),
                "first_fracture_h": hf,
            })

        bone_curves.append({
            "region_cn": cn,
            "bone": bone,
            "weak_subregion": weak_key,
            "weak_sigma_c_mpa": float(weak_sc),
            "sigma_stored_mpa": float(stored_sc),
            "risk_threshold_rigid": float(ratio_rigid),
            "sigma_at_5m_rigid_mpa": float(risk[4] * stored_sc),
            "first_fracture_rigid_h": h_rigid,
            "curve_M1_fe_transfer": m1,
            "curve_M2_physics": m2,
        })

    # ---- 关键点：软垫极限（R_min）与 FE k=1000 ---------------------------------
    R_best = fe_table["r_min"]
    k_best = fe_table["k_at_r_min_N_per_mm"]
    soft_limit: list[dict] = []
    for cn, bone, weak_key in LOAD_BONES:
        r0 = rows[bone]
        stored_sc = float(r0["sigma_c_mpa"])
        risk = np.asarray(r0["risk_by_height"], dtype=float)
        ratio_rigid = _sc(weak_key) / stored_sc
        h_r = first_fracture(heights, risk, ratio_rigid)
        h_b = first_fracture(heights, risk, ratio_rigid / R_best)
        soft_limit.append({
            "region_cn": cn, "bone": bone,
            "first_fracture_rigid_h": h_r,
            "first_fracture_best_pad_h": h_b,
            "delta_m": (None if (h_r is None or h_b is None) else float(h_b - h_r)),
        })

    # ---- 物理反推：FE 最佳折减比对应的 k_pad --------------------------------
    k_pad_equiv_for_fe_best = k_pad_equiv(R_best)
    k_pad_equiv_for_fe_1000 = k_pad_equiv(r_fe(1000.0, fe_table))

    # ---- FE 绝对外推（anat 尺度）对照 ----------------------------------------
    fe_abs = fe_absolute_heights(plantar)
    fe_fixed_h = fe_abs.get("fixed:None")
    fe_k1000_h = fe_abs.get("spring:1000.0")

    # ---- 组装 --------------------------------------------------------------
    result = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "scenario": "S1 垫子/地面刚度 -> 首折高度（抱石）",
        "units": "mm-N-MPa-s；k(N/mm)=FE Winkler；k_pad(N/m)=物理等效",
        "method": (
            "冻结纯轴向基线（risk_1d_loadshare rows_with_share_parallel）上叠加两个垫子模型："
            "M1=FE 实测 gauge 折减比 R_FE(k)=gauge(spring)/gauge(fixed)（主口径），"
            "sigma_pad(h;k)=R_FE(k)·sigma_rigid(h)；"
            "M2=串联线性弹簧 R_phys(k_pad)=sqrt(k_pad/(k_pad+k_body))（物理交叉校验）。"
            "首折 = 最小 h 使 sigma_pad >= σ_c(弱子区域)。"
        ),
        "inputs": {
            "risk_1d_loadshare_json": str(LOADSHARE_JSON.relative_to(ROOT)),
            "fracture_matrix_json": str(FM_JSON.relative_to(ROOT)),
            "plantar_bc_summary_json": str(PLANTAR_JSON.relative_to(ROOT)),
            "material_strength_source": "src/climbing/bone.py::MATERIAL_STRENGTH_MPA",
            "share_policy_primary": "rows_with_share_parallel",
            "k_body_N_per_m": K_BODY_N_PER_M,
            "k_body_source": "src/climbing/pad.py::hard_surface 默认 1.0e7 N/m",
        },
        "heights_m": [int(h) for h in heights],
        "fe_transfer_table": fe_table,
        "physics_model": {
            "formula": "R_phys(k_pad) = sqrt(k_pad / (k_pad + k_body))",
            "derivation": (
                "F_peak = v0·sqrt(m·k_series)，k_series = k_pad·k_body/(k_pad+k_body)"
                "（能量 ½mv²=½k x² 与 F=k x）；R = F_pad/F_rigid = sqrt(k_series/k_body)"
            ),
            "k_body_N_per_m": K_BODY_N_PER_M,
            "k_grid_N_per_m": [float(k) for k in K_PAD_PHYSICS_N_PER_M],
            "k_pad_equiv_for_fe_best_N_per_m": float(k_pad_equiv_for_fe_best),
            "k_pad_equiv_for_fe_k1000_N_per_m": float(k_pad_equiv_for_fe_1000),
        },
        "rigid_baseline": {
            cn: next(b["first_fracture_rigid_h"] for b in bone_curves if b["region_cn"] == cn)
            for cn, _, _ in LOAD_BONES
        },
        "soft_pad_limit": soft_limit,
        "bone_curves": bone_curves,
        "fe_absolute_comparison": {
            "scale": "FE anat 单骨跟骨尺度（非冻结 1D 尺度，绝对高度不可直接比较）",
            "fixed_first_fracture_h": fe_fixed_h,
            "spring_k1000_first_fracture_h": fe_k1000_h,
            "height_ratio_k1000_over_fixed": (
                None if (fe_fixed_h is None or fe_k1000_h is None) else float(fe_k1000_h / fe_fixed_h)
            ),
            "note": "FE 外推高度用 anat 扫描真跑 F_subt(h) 反查；本脚本不重跑 FEBio。",
        },
        "assumptions": [
            "冻结基线（纯轴向 σ(h)=risk_by_height×stored σ_c）与子区域强度（足150/腓两端70/胫两端70/股骨颈80）"
            "逐位沿用 AXIAL_SUBREGION_REPORT.md；本脚本只加垫子这一维。",
            "M1 假设垫子按同一折减比 R_FE(k) 折减整条载荷链峰值应力（1D 链线弹性 ⇒ 载荷与应力成比例）。"
            "R_FE 是**跟骨局部 gauge** 实测（单骨 THUMS AM50、h=5m、a=0、线弹性后处理、无单元删除），"
            "外推到腓/胫/股骨是**假设**，不是实测。",
            "M1 的 R_FE(k) 非单调（软端峰值在内部/关节面边缘，硬端在跖面边缘；k=1e3 最低 0.764）；"
            "k>=1e6 视为刚性极限 R=1；k<10 平推软端 0.838。",
            "M2 用线性弹簧；真实泡沫是平台+压实非线性，R_phys 只应视为量级/方向校验，"
            "若要定量应取相关压缩量下的**割线刚度**。",
            "绝对首折高度不可移植（[Y25] 范式）；本报告只给相对变化与排序。",
            "本脚本只读复用既有产物，不调用 FEBio，不重跑 OpenSim 正动力学。",
        ],
    }
    return result


# ---------------------------------------------------------------------------
def _fmt_h(v) -> str:
    return ">50" if v is None else f"{int(round(v))}"


def _fmt_k(k) -> str:
    if k is None:
        return "∞ (刚性)"
    return f"{k:g}"


def _write_plot(res: dict, out_png: Path) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False
    try:
        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        fig, ax = plt.subplots(figsize=(8.6, 5.2))
        for b in res["bone_curves"]:
            xs, ys = [], []
            for rec in b["curve_M1_fe_transfer"]:
                k = rec["k_N_per_mm"]
                h = rec["first_fracture_h"]
                if k is None or h is None:
                    continue
                xs.append(k)
                ys.append(h)
            if xs:
                ax.semilogx(xs, ys, "-o", ms=4, label=f"{b['region_cn']}(M1)")
        ax.set_xlabel("垫子刚度 k (N/mm, FE Winkler, log)")
        ax.set_ylabel("首折高度 h* (m)")
        ax.set_title("S1 垫子刚度 -> 首折高度（冻结纯轴向基线 + FE 实测折减比）")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(out_png, dpi=140)
        plt.close(fig)
        return True
    except Exception:
        return False


def _write_report(path: Path, res: dict) -> None:
    L: list[str] = []
    A = L.append
    t = res["fe_transfer_table"]
    ph = res["physics_model"]
    rb = res["rigid_baseline"]
    sl = {r["region_cn"]: r for r in res["soft_pad_limit"]}
    fe = res["fe_absolute_comparison"]
    fe_fixed = fe["fixed_first_fracture_h"]
    fe_k1000 = fe["spring_k1000_first_fracture_h"]
    fe_ratio = fe["height_ratio_k1000_over_fixed"]
    km1 = [rec["k_N_per_mm"] for rec in res["bone_curves"][0]["curve_M1_fe_transfer"]]

    A("# S1 场景：垫子/地面刚度 → 首次骨折高度（抱石）")
    A("")
    A(f"> 生成时间 {res['generated_at']}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 本报告**新增**：不改既有模块/默认值，不覆盖既有产物，不调用 FEBio，不重跑 OpenSim。")
    A("> 只读复用：冻结纯轴向基线 + 已有 FE 跖面 BC 实测（`plantar_bc_summary.json`）。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    A(f"- **垫子越软，首折越高**：FE 实测折减比最低点 R_min = **{t['r_min']:.4f}**"
      f"（k = {t['k_at_r_min_N_per_mm']:g} N/mm，即应力降 {(1-t['r_min'])*100:.1f}%）。")
    A(f"- 刚性基线（R=1）：腓骨 {_fmt_h(rb['腓骨'])} m、足部 {_fmt_h(rb['足部'])} m、"
      f"胫骨 {_fmt_h(rb['胫骨'])} m、股骨 {_fmt_h(rb['股骨'])} m。")
    A(f"- 最好垫子（R_min，k≈{t['k_at_r_min_N_per_mm']:g} N/mm）：足部 "
      f"{_fmt_h(sl['足部']['first_fracture_rigid_h'])} → **{_fmt_h(sl['足部']['first_fracture_best_pad_h'])} m**、"
      f"腓骨 {_fmt_h(sl['腓骨']['first_fracture_rigid_h'])} → **{_fmt_h(sl['腓骨']['first_fracture_best_pad_h'])} m**、"
      f"胫骨 {_fmt_h(sl['胫骨']['first_fracture_rigid_h'])} → **{_fmt_h(sl['胫骨']['first_fracture_best_pad_h'])} m**、"
      f"股骨 {_fmt_h(sl['股骨']['first_fracture_rigid_h'])} → **{_fmt_h(sl['股骨']['first_fracture_best_pad_h'])} m**。")
    A(f"- **对比已有 FE 结果**：同一 k≈1e3，FE 单骨跟骨外推把首折从 {fe_fixed:.2f} m "
      f"抬到 {fe_k1000:.2f} m（×{fe_ratio:.2f}）；冻结 1D（本脚本 M1）足部从 "
      f"{_fmt_h(rb['足部'])} m 抬到 {_fmt_h(sl['足部']['first_fracture_best_pad_h'])} m。"
      "两者方向一致，但**绝对高度不同源**（见 §4）。")
    A("")
    A("---")
    A("")
    A("## 1. 地面/垫子当前如何建模（rigid vs compliant）")
    A("")
    A("| 层 | 模块/文件 | 现状 | 刚性/柔性 |")
    A("|---|---|---|---|")
    A("| 1D 接触 | `src/climbing/pad.py` | `simulate_boulder_fall` 的渐进接触；"
      "`on_pad=False` 用 `hard_surface(k=1e7 N/m)` 线性弹簧（人体等效刚度，非真刚性墙） | **刚性等价**（基线） |")
    A("| GRF 生成 | `src/climbing/coupling/opensim_grf.py::ground_reaction` | `on_pad` 默认 False；"
      "把 1D 接触力翻译成左右脚的 `(t, F)` | **刚性地面** |")
    A("| 多体注射 | `src/climbing/coupling/opensim_fall.py::run_dead_drop` | `PrescribedForce` 把 GRF 施加到 "
      "`calcn`，腿关节锁死、肌肉关闭；**模型内没有地面接触单元** | 载荷边界（刚性） |")
    A("| FE 跖面 BC | `src/climbing/coupling/plantar_bc.py::apply_plantar_bc` | `fixed`/`roller`/`spring`(Winkler k)；"
      "`spring` 是可参数化的柔性支承 | **可柔性**（仅 FE 探针用） |")
    A("| 载荷链数据 | `temp/opensim_fe/fracture_matrix/load_chain.npz` | 50 高度的 `grf_peak_n` + 逐关节 "
      "`vert/abs` 峰值；全部由**刚性地面**链产生 | **刚性** |")
    A("")
    A("**结论**：冻结基线（`risk_1d_loadshare` / `AXIAL_SUBREGION_REPORT` / `load_chain.npz`）"
      "是**刚性地面**口径，没有垫子刚度这一维。唯一已有的柔性地面证据是**独立的 FE 跖面 BC 探针**"
      "（`plantar_bc.py` + `plantar_bc_summary.json`），它给出 `k -> gauge 折减比`，正好可用作"
      "把刚性基线迁移到柔性地面的**实测转移函数**。")
    A("")
    A("---")
    A("")
    A("## 2. 垫子模型与公式")
    A("")
    A("### 2.1 路线 M1（主口径）：FE 实测折减比转移")
    A("")
    A("```")
    A("R_FE(k) = gauge_max(spring_k) / gauge_max(fixed)      # 来自 plantar_bc_summary.json，实测")
    A("sigma_pad(h; k) = R_FE(k) · sigma_rigid(h)            # sigma_rigid(h)=risk(h)·stored σ_c")
    A("首折 h*(k) = min h 使 sigma_pad(h;k) >= σ_c(弱子区域)")
    A("等价：risk(h) >= (σ_c_weak / σ_c_stored) / R_FE(k)")
    A("```")
    A("")
    A("| k (N/mm) | gauge_max (MPa) | R_FE = gauge/fixed |")
    A("|---:|---:|---:|")
    A(f"| fixed | {t['fixed_gauge_max_mpa']:.2f} | 1.0000 |")
    for k, g, r in zip(t["k_N_per_mm"], t["gauge_max_mpa"], t["ratio_vs_fixed"]):
        A(f"| {k:g} | {g:.2f} | {r:.4f} |")
    A("")
    A(f"> R_FE(k) **非单调**：软端（k≤1e3）峰值移到内部/关节面边缘（跖面固定端奇异被消掉），"
      f"最低 {t['r_min']:.4f} @ k={t['k_at_r_min_N_per_mm']:g}；硬端峰值回到跖面边缘，"
      f"k=1e6 残差相对 fixed 仅 {t['hardest_k_residual_vs_fixed_pct']:.2f}%（收敛回刚性）。")
    A("")
    A("### 2.2 路线 M2（物理交叉校验）：串联线性弹簧")
    A("")
    A("```")
    A("F_peak = v0 · sqrt(m · k_series)        # 能量 ½mv² = ½k_s x²，F = k_s x")
    A("k_series = k_pad · k_body / (k_pad + k_body)   # 垫子与人体组织串联")
    A("R_phys(k_pad) = F_pad / F_rigid = sqrt(k_pad / (k_pad + k_body))")
    A("k_body = 1.0e7 N/m  (pad.hard_surface 默认)")
    A("```")
    A("")
    A(f"- 反推：FE 最佳折减比 R_min={t['r_min']:.4f} 对应物理 **k_pad ≈ "
      f"{ph['k_pad_equiv_for_fe_best_N_per_m']:.3g} N/m**（≈ k_body 量级）；"
      f"k=1000 N/mm 的 R={r_fe(1000.0, t):.4f} 对应 k_pad ≈ "
      f"{ph['k_pad_equiv_for_fe_k1000_N_per_m']:.3g} N/m。")
    A(f"- **选路说明**：主口径用 **M1（实测）**，因为它直接复用已有 FE 证据且与任务锚点"
      f"（FE 跟骨 5.43 m @ k≈1e3）同源；M2 只做量级/方向校验（真实泡沫非线性，R_phys 会高估软垫收益）。")
    A("")
    A("---")
    A("")
    A("## 3. k → 首折高度曲线（M1 主口径，逐骨）")
    A("")
    A("| 部位 | 弱子区域(σ_c) | 刚性 h* | " +
      " | ".join(_fmt_k(k) for k in km1) + " |")
    A("|---|---|---:|" + "|".join([":---:"] * len(km1)) + " |")
    for b in res["bone_curves"]:
        cells = " | ".join(_fmt_h(rec["first_fracture_h"]) for rec in b["curve_M1_fe_transfer"])
        A(f"| {b['region_cn']} | {b['weak_subregion']} ({b['weak_sigma_c_mpa']:.0f}) | "
          f"**{_fmt_h(b['first_fracture_rigid_h'])}** | {cells} |")
    A("")
    A("> 表头 k 列（N/mm）：" + "、".join(_fmt_k(k) for k in km1) + "。")
    A("> 首折高度由离散的 1..50 m 扫描给出，故为整数米（分辨率 1 m）。")
    A("")
    A("### 3.1 5 m 处的应力（MPa，随 k 变化）")
    A("")
    A("| 部位 | 刚性 σ@5m | R_min σ@5m | 降幅 |")
    A("|---|---:|---:|---:|")
    for b in res["bone_curves"]:
        sig0 = b["sigma_at_5m_rigid_mpa"]
        sigb = sig0 * t["r_min"]
        A(f"| {b['region_cn']} | {sig0:.2f} | {sigb:.2f} | {(1-sigb/sig0)*100:.1f}% |")
    A("")
    A("### 3.2 M2 物理曲线（k_pad N/m → 首折，逐骨）")
    A("")
    A("| 部位 | " + " | ".join(f"{k:.0e}" for k in ph["k_grid_N_per_m"]) + " |")
    A("|---|" + "|".join([":---:"] * len(ph["k_grid_N_per_m"])) + "|")
    for b in res["bone_curves"]:
        cells = " | ".join(_fmt_h(rec["first_fracture_h"]) for rec in b["curve_M2_physics"])
        A(f"| {b['region_cn']} | {cells} |")
    A("")
    A("---")
    A("")
    A("## 4. 对比：刚性基线 / 软垫极限 / 已有 FE 结果")
    A("")
    A("| 部位 | 刚性基线 h* (m) | 最好垫子 h* (m) @ R_min | Δ (m) |")
    A("|---|---:|---:|---:|")
    for r in res["soft_pad_limit"]:
        A(f"| {r['region_cn']} | {_fmt_h(r['first_fracture_rigid_h'])} | "
          f"**{_fmt_h(r['first_fracture_best_pad_h'])}** | "
          f"{'—' if r['delta_m'] is None else f'+{r['delta_m']:.0f}'} |")
    A("")
    A("**已有 FE 单骨跟骨结果（anat 尺度，实测外推）**：")
    A("")
    A("| 量 | 值 |")
    A("|---|---:|")
    A(f"| FE fixed 首折高度 | {fe['fixed_first_fracture_h']:.2f} m |")
    A(f"| FE spring k=1000 N/mm 首折高度 | {fe['spring_k1000_first_fracture_h']:.2f} m |")
    A(f"| 高度比 (k=1000 / fixed) | {fe['height_ratio_k1000_over_fixed']:.2f}× |")
    A("")
    A("> **尺度警告**：FE anat 尺度（fixed 2.11 m）与冻结 1D 尺度（足部 8 m）**绝对高度不同源**"
      "（FE 含弯曲/BC 伪影，1D 为纯轴向），二者**不可直接比较绝对高度**。可比较的是**变化方向与量级**：")
    A(f"> - FE：应力降 {(1-t['r_min'])*100:.1f}% → 高度 ×{fe['height_ratio_k1000_over_fixed']:.2f}（2.11 → 5.43 m）；")
    A(f"> - 冻结 1D（本脚本 M1）：同一折减比使足部 {_fmt_h(sl['足部']['first_fracture_rigid_h'])} → "
      f"{_fmt_h(sl['足部']['first_fracture_best_pad_h'])} m、腓骨 {_fmt_h(sl['腓骨']['first_fracture_rigid_h'])} → "
      f"{_fmt_h(sl['腓骨']['first_fracture_best_pad_h'])} m。")
    A("> 两者方向一致（垫子抬高首折），1D 因 σ(h) 增长较缓，绝对抬高幅度更大。")
    A("")
    A("---")
    A("")
    A("## 5. 敏感性与不确定范围（避免单个误导数字）")
    A("")
    A("1. **FE 转移非单调** ⇒ 应力折减比在 k∈[10, 1e6] 内落在 "
      f"**[{t['r_min']:.4f}（k=1e3）, 1.0000（k=1e6 刚性极限）]**；对应足部首折的**有界范围** "
      f"**[{_fmt_h(rb['足部'])} m（R=1）, {_fmt_h(sl['足部']['first_fracture_best_pad_h'])} m（R_min）]**。"
      "不把“k=1e3 → 5.43 m”当作唯一答案。")
    A("2. **M1 折减比是跟骨局部实测**：外推到腓/胫/股骨为假设；若垫子主要消掉的是**跟骨固定端奇异**"
      "而非全局载荷，则其它骨的收益应更小（保守）。")
    A("3. **M2 物理路线**给的是方向/量级：R_phys(k_pad) 单调递减，"
      f"k_pad∈[{ph['k_grid_N_per_m'][0]:.0e}, {ph['k_grid_N_per_m'][-1]:.0e}] N/m → "
      f"R∈[{r_phys(ph['k_grid_N_per_m'][0]):.3f}, {r_phys(ph['k_grid_N_per_m'][-1]):.3f}]；"
      "真实泡沫非线性（平台+压实）会显著削弱软端收益，故 M2 的软端是**乐观上界**。")
    A("4. **名义应力 1D 局限**：σ 是纯轴向名义压缩，忽略弯曲/剪切/偏心/应力集中"
      "（对细长骨低估、对紧凑骨高估）；本报告只在**相对变化**上使用它。")
    A("")
    A("---")
    A("")
    A("## 6. 假设与诚实边界")
    A("")
    for i, a in enumerate(res["assumptions"], 1):
        A(f"{i}. {a}")
    A("")
    A("- **测量 vs 建模**：M1 的 R_FE(k) 是**实测**（FE 探针）；σ(h) 映射与整链外推是**建模**；"
      "M2 全程**建模**。所有首折高度都是**线弹性比例外推**，不是重跑高度扫描。")
    A("- 本脚本不改变任何既有产物；`risk_1d_loadshare.json`、`plantar_bc_summary.json`、"
      "`fracture_matrix.json` 均为只读输入。")
    A("")
    A("---")
    A("")
    A("## 7. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\pad_stiffness.py")
    A("# 只读复用：")
    A("#   results/opensim_fe/risk_1d_loadshare.json    (冻结纯轴向 risk_by_height / stored σ_c)")
    A("#   results/opensim_fe/fracture_matrix.json      (heights_m 1..50)")
    A("#   results/opensim_fe/plantar_bc_summary.json   (FE 实测 gauge(k) / fixed / 外推高度)")
    A("#   src/climbing/bone.py MATERIAL_STRENGTH_MPA  (弱子区域 σ_c)")
    A("# 产物：")
    A("#   results/opensim_fe/pad_stiffness.json")
    A("#   results/opensim_fe/PAD_STIFFNESS_REPORT.md")
    A("#   results/opensim_fe/pad_stiffness_curve.png")
    A("```")
    A("")
    A("---")
    A("")
    A("## 8. 产物")
    A("")
    A("- 本报告：`results/opensim_fe/PAD_STIFFNESS_REPORT.md`")
    A("- 新脚本：`scripts/opensim_fe/pad_stiffness.py`（additive）")
    A("- 机器可读：`results/opensim_fe/pad_stiffness.json`")
    A("- 曲线：`results/opensim_fe/pad_stiffness_curve.png`")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="S1 垫子刚度 -> 首折高度（冻结纯轴向基线 + FE 实测折减比）")
    ap.add_argument("--out-dir", default=str(RES))
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    res = build()
    js_path = out / "pad_stiffness.json"
    js_path.write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    md_path = out / "PAD_STIFFNESS_REPORT.md"
    _write_report(md_path, res)
    png_ok = _write_plot(res, out / "pad_stiffness_curve.png")

    t = res["fe_transfer_table"]
    print("=" * 96)
    print("S1 垫子刚度 -> 首折高度")
    print(f"FE 实测折减比 R_FE(k): min={t['r_min']:.4f} @ k={t['k_at_r_min_N_per_mm']:g} N/mm")
    print(f"物理 k_pad 反推（R_min） = {res['physics_model']['k_pad_equiv_for_fe_best_N_per_m']:.3g} N/m")
    sl_by_cn = {r["region_cn"]: r for r in res["soft_pad_limit"]}
    print("-" * 96)
    for b in res["bone_curves"]:
        cn = b["region_cn"]
        print(f"  {cn:4s} rigid={_fmt_h(b['first_fracture_rigid_h']):>4}m  "
              f"M1(k=1e3)={_fmt_h(b['curve_M1_fe_transfer'][4]['first_fracture_h']):>4}m  "
              f"M1(R_min)={_fmt_h(sl_by_cn[cn]['first_fracture_best_pad_h']):>4}m")
    fe = res["fe_absolute_comparison"]
    print(f"FE 单骨跟骨（anat 尺度）：fixed={fe['fixed_first_fracture_h']:.2f} m -> "
          f"k=1e3 {fe['spring_k1000_first_fracture_h']:.2f} m "
          f"(×{fe['height_ratio_k1000_over_fixed']:.2f})")
    print("=" * 96)
    print("wrote:", js_path)
    print("wrote:", md_path)
    print("wrote:", out / "pad_stiffness_curve.png", "(ok)" if png_ok else "(skipped)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
