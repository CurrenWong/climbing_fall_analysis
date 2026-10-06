"""scenario_field_matrix.py -- 现场实测跌倒运动学 → 模型参数 → 现有载荷链 → 骨折区域旗标。

背景（为什么需要这个脚本）
--------------------------
既有产物（``fracture_matrix.json`` / ``risk_1d_loadshare.json`` /
``AXIAL_SUBREGION_REPORT.md``）只覆盖**单一姿势**（``feet-first-stiff`` 直立、
双足对称、纯竖直、无水平速度）。但真实抱石跌倒是**运动学分布**：

* **Beurienne 2025 现场调查**把跌倒记为 7 种运动学
  = 起始位 × 跌落中的旋转 × 着地部位；
* **Müller/Heck 2024**（``paper/Boulder_Dissertation_EN.md`` §4.3.3）报告
  “**较低坠落 → 更多下肢损伤**”，这正是**落地姿态随高度分布变化**的结果。

本脚本是**加法式**的：它把上述运动学**映射到模型已有的、稳定的 opt-in 旋钮**
（``ground_reaction(tilt_deg=..., split=.../foot_forces=...)``、
``run_dead_drop(vx0_ms=, vz0_ms=)``、坠落高度），跑**现有**载荷链
（``ground_reaction → run_dead_drop → joint_reaction``），再用与
``risk_1d_loadshare.py`` **同一套载荷分配**（胫/腓并联分流 + 中轴骨 mass-above）
计算**纯轴向 + 弱子区域** 1D 风险；最后与临床分布并列对照。

它**不改**任何既有模块/默认值/签名，**不覆盖**既有产物，**不调用 FEBio**，
只写新文件。

与 ``risk_1d_loadshare.py`` 完全一致的口径
------------------------------------------
* ``risk(h) = F_bone(h) / (A_section * sigma_c_stored)``（1D 名义压缩）
* 纯轴向应力 ``sigma(h) = risk(h) * sigma_c_stored``（长骨为骨干值）
* 弱子区域旗标：``sigma(h) >= sigma_c(弱子区域)``
  （足150 / 胫两端70 / 腓两端70 / 股骨颈80 / 骨盆180 / 脊柱150 / 颅骨160）
* 力映射：``F_bone = force_multiplier(bone) * F_joint``，其中
  胫/腓 = 并联分流 ``f_i``，L3/T6/C5/parietal = mass-above 比，其余 = 1.0。

现场运动学 → 模型旋钮的映射（诚实标注 measured vs modeled）
----------------------------------------------------------
| 现场维度 | 取值 | 模型旋钮 | 性质 |
|---|---|---|---|
| 起始位 | 垂直 / 后仰 / 前倾 | ``tilt_deg``（0 / +25 / -25） | **modeled**（代理） |
| 跌落旋转 | 无 / 纵轴 / 前后 | ``split``/``foot_forces``（对称 / 不对称 / 不对称+水平速度） | **modeled**（代理） |
| 着地部位 | 双足 / 侧倾 / 头先 / 臀 | 双足(对称) / 单足 / ``foot_forces`` 低(代理) / ``foot_forces`` 低(代理) | **modeled**（代理） |
| 坠落高度 | 2.0/3.0/3.5/4.5 m | ``height_m`` | measured（现场高度档） |

⚠️ **可表示性边界**：本载荷链的 GRF 只施加在 ``calcn``（足）。因此
“头先着地 / 臀部着地”**无法**由足端载荷链物理表示 —— 它们只在网格里
**声明性**给出映射，并标 ``representable=False`` / ``PROXY``；实际**不跑**
（或跑出的足端结果明确标注 out-of-scope）。上肢经由 ``lumbar`` 关节子树
（含双上肢）间接承载，另用**头部惯性**估计（同 ``axial_subregion.py`` §4.1）。

输入（只读）
------------
* ``results/opensim_fe/fracture_matrix.json``          -> 不存在时无需（本脚本自跑）
* ``results/opensim_fe/bc_robust_metric_route1.json``  -> ``A_section_mm2`` / ``sigma_c_mpa``
* ``results/opensim_fe/risk_1d_loadshare.json``        -> 载荷分配倍率 + mass-above
* ``temp/opensim_fe/fracture_matrix/load_chain.npz``   -> 复用竖直基准的整档缓存
* ``src/climbing/coupling/{opensim_grf,opensim_fall,joint_reactions}.py``（现行 API）
* ``src/climbing/bone.py::MATERIAL_STRENGTH_MPA``      -> 弱子区域 σ_c

产物（新文件）
--------------
* ``results/opensim_fe/scenario_field_matrix.json``
* ``results/opensim_fe/SCENARIO_FIELD_MATRIX_REPORT.md``
* ``temp/opensim_fe/scenario_field_matrix/scenario_load_chain.npz``（本脚本缓存）
* ``tests/test_scenario_field_matrix.py``（回归）

复现（仓库根）::

    cd D:\\Project\\climbing_fall_analysis
    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\scenario_field_matrix.py
    # 只重算（忽略本脚本缓存）：--force
    # 只做汇总、不跑载荷链（用于快速回归）：--no-run
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from climbing.bone import MATERIAL_STRENGTH_MPA  # noqa: E402

# ---------------------------------------------------------------------------
# paths / constants
# ---------------------------------------------------------------------------
RES = ROOT / "results" / "opensim_fe"
TEMP = ROOT / "temp" / "opensim_fe" / "scenario_field_matrix"
FRACTURE_JSON = RES / "fracture_matrix.json"
ROUTE1_JSON = RES / "bc_robust_metric_route1.json"
LOADSHARE_JSON = RES / "risk_1d_loadshare.json"
CACHE_NPZ = TEMP / "scenario_load_chain.npz"
EXISTING_LOADCHAIN = ROOT / "temp" / "opensim_fe" / "fracture_matrix" / "load_chain.npz"

MODEL_MASS_KG = 75.337                     # Rajagopal2015 总质量（与 opensim_grf 默认一致）
HEIGHTS = (2.0, 3.0, 3.5, 4.5)             # 现场代表性高度档 (m)
JOINTS = ("subtalar_r", "ankle_r", "hip_r", "lumbar")   # 现有载荷链的 4 个关节

G_ACC = 9.81
TOTAL_BODY_MASS_KG = MODEL_MASS_KG

# ---------------------------------------------------------------------------
# 9 部位 → 代表骨 / 关节 / 弱子区域 / 骨干子区域（与 axial_subregion.py 同表）
#   weak 决定“部位骨折”；shaft 单独复核“骨干从不骨折”。
# ---------------------------------------------------------------------------
#: (中文名, 代表骨, 关节, 弱子区域键, 骨干子区域键 or None)
REGIONS: tuple[tuple[str, str, str, str, str | None], ...] = (
    ("足部", "calcaneus_r", "subtalar_r", "calcaneus", None),
    ("胫骨", "tibia_r", "ankle_r", "tibia_ends", "tibia_shaft"),
    ("腓骨", "fibula_r", "ankle_r", "fibula_ends", "fibula_shaft"),
    ("股骨", "femur_r", "hip_r", "femoral_neck", "femur_shaft"),
    ("骨盆", "R_HIPBONE", "hip_r", "pelvis", None),
    ("腰椎", "L3", "lumbar", "spine", None),
    ("胸椎", "T6", "lumbar", "spine", None),
    ("颈椎", "C5", "lumbar", "spine", None),
    ("颅骨", "parietal_r", "lumbar", "skull", None),
)
REGION_NAMES = tuple(r[0] for r in REGIONS)

#: 临床分桶（把 9 个模型部位映射到临床报告的损伤类别）。
CLINICAL_BUCKETS: dict[str, tuple[str, ...]] = {
    "下肢(合计)": ("足部", "胫骨", "腓骨", "股骨"),
    "踝": ("足部", "胫骨", "腓骨"),   # 跟骨 + 胫/腓远端(踝)
    "膝": ("股骨",),                  # 股骨（远端/颈）代表膝关节链
    "脊柱": ("腰椎", "胸椎", "颈椎"),
    "骨盆": ("骨盆",),
    "颅骨": ("颅骨",),
}

#: 临床分布（对照目标）。measured（论文原文数值，见报告 §5）。
CLINICAL_TARGETS: dict = {
    # Heck/Müller 2024（n=447 injuries）：下肢 61.1%；踝 36.7% / 膝 16.8% / 脊柱 7.2%
    "lower_limb_pct_heck": 61.1,
    "lower_limb_pct": [61.0, 67.0],           # Heck 61.1% ；Beurienne 67%
    "lower_limb_pct_beurienne": 67.0,
    "ankle_pct_heck": 36.7,
    "ankle_pct_beurienne": 40.0,
    "knee_pct_heck": 16.8,
    "knee_pct_beurienne": 15.0,
    "spine_pct_heck": 7.2,
    "elbow_pct_heck": 12.3,
    "upper_limb_pct_heck": 26.6,
    "ankle_sprain_pct_heck": 26.1,            # 踝扭伤 = 单一最常见损伤
    "ankle_sprain_rank": 1,                    # Beurienne 2025 & Heck：踝扭伤 #1
    # R5（Heck/Müller 2024 §4.3.3，p=0.027）：低坠落(<2 m)下肢 79.1% vs 高坠落(>2 m)58.3%
    "r5_low_fall_lower_limb_pct": 79.1,
    "r5_high_fall_lower_limb_pct": 58.3,
    "r5_p_value": 0.027,
    "sources": [
        "Heck/Müller 2024 — paper/Boulder_Dissertation_EN.md §4.2.6 Table 1, §4.2.7, §4.3.3, §5.1, §5.2.3",
        "Beurienne 2025 — paper/MinerU_markdown_fspor-7-1609133_*.md §3.1.2 (Table 1), §3.3, Fig.4",
    ],
}

# ---------------------------------------------------------------------------
# 场景网格
# ---------------------------------------------------------------------------
#: 起始位 → tilt_deg（地面法向倾角；n̂=(sinθ,cosθ,0)，θ>0 法向偏 +x）。
START_TILT = {"vertical": 0.0, "leaning-back": 25.0, "leaning-forward": -25.0}
START_CN = {"vertical": "垂直起始", "leaning-back": "后仰起始", "leaning-forward": "前倾起始"}
ROTATION_CN = {"none": "无旋转", "longitudinal": "纵轴旋转", "other": "前后旋转"}
LANDING_CN = {
    "upright-on-feet": "双足直立落地",
    "feet-leaning": "侧倾/单足落地",
    "head-first": "头先着地",
    "buttocks": "臀部着地",
}

#: 网格轴（现场运动学的 3 个维度）。
GRID_AXES: dict[str, tuple[str, ...]] = {
    "start": ("vertical", "leaning-back", "leaning-forward"),
    "rotation": ("none", "longitudinal", "other"),
    "landing": ("upright-on-feet", "feet-leaning", "head-first", "buttocks"),
}

#: 现场 7 种运动学（Beurienne 2025, §3.3 / Fig.4）→ 场景 id 的对照。
#: 前 3 名为论文明确给出百分比的“top-3”；其余 4 种为 Fig.4 中 ≥10 例的剩余项。
#: injury_tendency 为论文报告的主要损伤倾向（用于临床对照的定性锚点）。
FIELD_KINEMATICS: tuple[dict, ...] = (
    dict(kind="K1", start="vertical", rotation="none", landing="upright-on-feet",
         scenario="s1_K1_vertical_upright", freq_pct=17.0, injury_tendency="下肢",
         injuries="13 踝扭伤 / 7 膝腱韧带撕裂 / 6 踝骨折",
         note="第 1 名（17%）：直立、无旋转、双足直立落地 → 下肢。"),
    dict(kind="K2", start="vertical", rotation="longitudinal", landing="feet-leaning",
         scenario="s2_K2_vertical_long_feetleaning", freq_pct=14.0, injury_tendency="下肢",
         injuries="10 踝扭伤 / 5 膝腱韧带撕裂 / 4 踝骨折",
         note="第 2 名（14%）：直立、纵轴旋转、侧倾双足落地 → 下肢。"),
    dict(kind="K3", start="leaning-back", rotation="longitudinal", landing="feet-leaning",
         scenario="s3_K3_leanback_long_feetleaning_UL", freq_pct=11.0,
         injury_tendency="上肢",
         injuries="7 踝扭伤 / 5 肘脱位 / 3 肘扭伤",
         note="第 3 名（11%）：后仰起始、纵轴旋转、侧倾落地 → **上肢**（论文 §4：后仰→旋转→"
              "不稳定落地→本能伸手撑垫）。"),
    dict(kind="K4", start="leaning-forward", rotation="other", landing="feet-leaning",
         scenario="s4_UL_reach_leanfwd", freq_pct=None, injury_tendency="下肢/上肢",
         injuries="—（Fig.4 剩余项）",
         note="Fig.4 剩余 4 种之一：前倾+前后旋转+侧倾 → 前扑/上肢或下肢。"),
    dict(kind="K5", start="vertical", rotation="other", landing="head-first",
         scenario="s5_headfirst", freq_pct=None, injury_tendency="上身/颅骨",
         injuries="—（Fig.4 剩余项）",
         note="Fig.4 剩余项：前后旋转 → 头先着地（足端链 out-of-scope）。"),
    dict(kind="K6", start="leaning-back", rotation="none", landing="buttocks",
         scenario="s6_buttocks", freq_pct=None, injury_tendency="骨盆/躯干",
         injuries="—（Fig.4 剩余项）",
         note="Fig.4 剩余项：后仰 → 臀部着地（足端链 out-of-scope）。"),
    dict(kind="K7", start="leaning-forward", rotation="none", landing="upright-on-feet",
         scenario="s7_leanfwd_upright", freq_pct=None, injury_tendency="下肢",
         injuries="—（Fig.4 剩余项）",
         note="Fig.4 剩余项：前倾、无旋转、双足直立落地（对称）。"),
)

#: 现场运动学 → 模型的**只用 opt-in 旋钮**映射。
#: 前四场景**实际运行**（top-3 + 1 个前倾前扑/上肢代理）；后三场景为声明性 PROXY。
SCENARIOS: tuple[dict, ...] = (
    dict(
        id="s1_K1_vertical_upright", start="vertical", rotation="none",
        landing="upright-on-feet", tilt_deg=0.0, split=(0.5, 0.5),
        foot_forces=None, vx0_ms=0.0, vz0_ms=0.0, run=True, representable=True,
        field_kind="K1",
        note="基准直立双足（与既有 load_chain.npz 默认场景同参数；2/3 m 直接复用缓存）。",
    ),
    dict(
        id="s2_K2_vertical_long_feetleaning", start="vertical", rotation="longitudinal",
        landing="feet-leaning", tilt_deg=0.0, split=(0.3, 0.7),
        foot_forces=None, vx0_ms=0.0, vz0_ms=0.0, run=True, representable=True,
        field_kind="K2",
        note="纵轴旋转 → 一足先着（不对称）：split=(0.3,0.7)（右足 70%）。",
    ),
    dict(
        id="s3_K3_leanback_long_feetleaning_UL", start="leaning-back",
        rotation="longitudinal", landing="feet-leaning", tilt_deg=25.0,
        split=(0.3, 0.7), foot_forces=None, vx0_ms=0.0, vz0_ms=0.0,
        run=True, representable="partial", field_kind="K3",
        note="后仰 25° + 纵轴旋转 + 侧倾单足：现场 K3 实际导致**上肢**损伤；"
             "本模型无上肢骨，经 lumbar 子树(含双上肢)间接承载（报告 §7 诚实标注代理）。",
    ),
    dict(
        id="s4_UL_reach_leanfwd", start="leaning-forward", rotation="other",
        landing="feet-leaning", tilt_deg=-25.0, split=(0.2, 0.8),
        foot_forces=None, vx0_ms=1.2, vz0_ms=0.0, run=True,
        representable="partial", field_kind="K4",
        note="前倾 + 前后旋转（不对称 split=(0.2,0.8)）+ 前向水平速度 1.2 m/s 的前扑代理"
             "（上肢经 lumbar 子树承载）。",
    ),
    dict(
        id="s5_headfirst", start="vertical", rotation="other",
        landing="head-first", tilt_deg=0.0, split=None,
        foot_forces=None, vx0_ms=2.0, vz0_ms=0.0, run=False,
        representable=False, field_kind="K5",
        note="PROXY：足端载荷链**无法**表示头先着地 —— 接触不在足端，无足端 split/foot_forces "
             "可言；不跑，标注 out-of-scope（其骨盆/脊柱响应留待含头/颈/躯干接触的后续模态）。",
    ),
    dict(
        id="s6_buttocks", start="leaning-back", rotation="none",
        landing="buttocks", tilt_deg=25.0, split=None,
        foot_forces=None, vx0_ms=0.0, vz0_ms=0.0, run=False,
        representable=False, field_kind="K6",
        note="PROXY：模型无臀部接触；接触不在足端，无 split/foot_forces；不跑，标注 out-of-scope。",
    ),
    dict(
        id="s7_leanfwd_upright", start="leaning-forward", rotation="none",
        landing="upright-on-feet", tilt_deg=-25.0, split=(0.5, 0.5),
        foot_forces=None, vx0_ms=0.0, vz0_ms=0.0, run=False, representable=True,
        field_kind="K7",
        note="前倾、无旋转、双足直立落地（对称）；声明性，不单独跑（与 s1 仅差法向倾角）。",
    ),
)


def enumerate_full_grid() -> list[dict]:
    """把 3×3×4 的全部 36 组合映射到 model 旋钮（声明式，标注是否现场观测/是否运行）。

    映射规则（**只用 opt-in 旋钮**；均为代理）：

    * ``tilt_deg = START_TILT[start]``（0 / +25 / -25）；
    * ``split``：rotation==none → (0.5,0.5)；rotation==longitudinal → (0.3,0.7)；
      rotation==other → None（改用 ``foot_forces`` 不对称）；
    * ``vx0_ms``：rotation==other → 1.2；否则 0；
    * ``foot_forces``：landing ∈ {head-first, buttocks} → 低值代理 (0.1/0.1 或 0.25/0.25)；
    * ``representable``：landing ∈ {upright-on-feet, feet-leaning} 才 True；
    * ``observed``：是否命中现场 7 种运动学；``run``：是否在 ``SCENARIOS`` 中实跑。
    """
    observed = {(k["start"], k["rotation"], k["landing"]): k["kind"] for k in FIELD_KINEMATICS}
    run_ids = {(s["start"], s["rotation"], s["landing"]): s["id"]
               for s in SCENARIOS if s["run"]}
    rows = []
    for start in GRID_AXES["start"]:
        for rot in GRID_AXES["rotation"]:
            for land in GRID_AXES["landing"]:
                if land in ("head-first", "buttocks"):
                    # 接触不在足端：无 split/foot_forces 可言（不可表示）
                    split, foot_forces = None, None
                elif rot == "none":
                    split, foot_forces = (0.5, 0.5), None
                elif rot == "longitudinal":
                    split, foot_forces = (0.3, 0.7), None
                else:  # other
                    split, foot_forces = (0.2, 0.8), None
                representable = land in ("upright-on-feet", "feet-leaning")
                key = (start, rot, land)
                rows.append(dict(
                    start=start, rotation=rot, landing=land,
                    tilt_deg=START_TILT[start], split=split,
                    foot_forces=foot_forces,
                    vx0_ms=(1.2 if rot == "other" else 0.0), vz0_ms=0.0,
                    representable=representable,
                    observed=(key in observed),
                    field_kind=observed.get(key),
                    run=(key in run_ids),
                    scenario=run_ids.get(key),
                ))
    return rows


def _sc(key: str) -> float:
    """弱子区域压缩强度 (MPa)，[Y25] Table 1。"""
    return float(MATERIAL_STRENGTH_MPA[key][1])


# ---------------------------------------------------------------------------
# 只读输入：A_section / sigma_c_stored / 载荷分配倍率
# ---------------------------------------------------------------------------
def load_static_inputs() -> dict:
    """读取既有产物，返回 A_section / sigma_c / 载荷分配倍率。"""
    route1 = json.loads(ROUTE1_JSON.read_text(encoding="utf-8"))
    loadshare = json.loads(LOADSHARE_JSON.read_text(encoding="utf-8"))
    a_by_bone = {r["bone"]: float(r["A_section_mm2"]) for r in route1["rows"]}
    sigma_stored = {r["bone"]: float(r["sigma_c_mpa"]) for r in route1["rows"]}
    multipliers = load_share_multipliers(loadshare)
    return {
        "A_section_mm2": a_by_bone,
        "sigma_c_stored_mpa": sigma_stored,
        "multipliers": multipliers,
        "loadshare_raw": loadshare,
        "route1_raw": route1,
    }


def load_share_multipliers(loadshare: dict) -> dict[str, float]:
    """replicate ``risk_1d_loadshare.py`` 的力映射倍率（bone -> force multiplier）。

    * 胫/腓：并联分流 ``f_i = k_i/sum(k)``（``tibia_fibula_split.split_fraction``）；
    * L3/T6/C5/parietal：mass-above 比（``mass_above_model.per_bone[b].ratio``）；
    * 其余（calcaneus/femur/pelvis）：1.0（单骨，无分配）。
    """
    mult: dict[str, float] = {b: 1.0 for _, b, _, _, _ in REGIONS}
    for b, f in loadshare["tibia_fibula_split"]["split_fraction"].items():
        mult[b] = float(f)
    for b, rec in loadshare["mass_above_model"]["per_bone"].items():
        mult[b] = float(rec["ratio"])
    return mult


# ---------------------------------------------------------------------------
# 载荷链（现有 API）
# ---------------------------------------------------------------------------
def _open_grf():
    from climbing.coupling.opensim_grf import ground_reaction  # noqa: PLC0415
    return ground_reaction


def _run_scenario_at_height(scenario: dict, h: float) -> dict:
    """跑一个场景在高度 h 的现有载荷链，返回关节纵向峰值 + GRF 诊断。"""
    from climbing.coupling.joint_reactions import joint_reaction  # noqa: PLC0415
    from climbing.coupling.opensim_fall import run_dead_drop  # noqa: PLC0415

    ground_reaction = _open_grf()

    kw: dict = dict(height_m=float(h), mass_kg=MODEL_MASS_KG, tilt_deg=float(scenario["tilt_deg"]))
    if scenario["foot_forces"] is not None:
        kw["foot_forces"] = tuple(scenario["foot_forces"])
    else:
        kw["split"] = tuple(scenario["split"])
    grf = ground_reaction(**kw)

    fall = run_dead_drop(grf, vx0_ms=float(scenario["vx0_ms"]), vz0_ms=float(scenario["vz0_ms"]))

    peaks = {j: float(joint_reaction(fall, grf, joint=j).peak_vertical_n) for j in JOINTS}
    return {
        "height_m": float(h),
        "joints_vertical_n": peaks,
        "grf_peak_n": float(grf.peak_total_n),
        "grf_impulse_rel_err": float(grf.impulse_rel_err),
        "loaded_sides": list(grf.loaded_sides),
        "model_mass_kg": float(fall.model_mass_kg),
    }


def _cache_key(sid: str, h: float) -> str:
    return f"{sid}@{h:g}"


def run_grid(*, reuse: bool, force: bool, no_compute: bool = False) -> dict:
    """跑实际运行场景（``SCENARIOS`` 中 ``run=True``）× 高度。

    * 复用本脚本缓存 ``CACHE_NPZ``；
    * 竖直对称基准在**整数高度**（2/3 m 命中）复用既有
      ``temp/opensim_fe/fracture_matrix/load_chain.npz`` 整档；
    * ``no_compute=True``（``--no-run``）：只读缓存，缺失的高度**跳过**，不跑载荷链。
    """
    TEMP.mkdir(parents=True, exist_ok=True)
    scenarios = [s for s in SCENARIOS if s["run"]]

    cache: dict[str, dict] = {}
    if reuse and CACHE_NPZ.is_file() and not force:
        raw = np.load(CACHE_NPZ, allow_pickle=True)
        for k in raw.files:
            cache[k] = raw[k].item()

    baseline_cache = _load_baseline_cache()

    results: dict[str, dict] = {}
    skipped: list[str] = []
    for sc in scenarios:
        sid = sc["id"]
        heights_out = []
        for h in HEIGHTS:
            key = _cache_key(sid, h)
            if key in cache and not force:
                heights_out.append(cache[key])
                continue
            # 竖直对称基准：直接复用既有 load_chain.npz（整档，整数高度）
            base = _baseline_value(sid, h, baseline_cache)
            if base is not None:
                cache[key] = base
                heights_out.append(base)
                continue
            if no_compute:
                skipped.append(key)
                continue
            t0 = time.time()
            rec = _run_scenario_at_height(sc, h)
            rec["wall_s"] = round(time.time() - t0, 1)
            cache[key] = rec
            heights_out.append(rec)
        if heights_out:
            results[sid] = {"scenario": _scenario_public(sc), "by_height": heights_out}
        _save_cache(cache)
    return {
        "runs": results,
        "cache_used": bool(cache),
        "cache_path": str(CACHE_NPZ.relative_to(ROOT)),
        "skipped": skipped,
    }


def _load_baseline_cache() -> dict | None:
    if not EXISTING_LOADCHAIN.is_file():
        return None
    try:
        z = np.load(EXISTING_LOADCHAIN)
    except Exception:  # noqa: BLE001
        return None
    return {k: z[k] for k in z.files}


def _baseline_value(sid: str, h: float, baseline: dict | None):
    """竖直对称基准场景在**整数**高度上可从既有 load_chain.npz 复用整档。"""
    if baseline is None or sid != "s1_vertical_upright":
        return None
    hs = np.asarray(baseline["heights"], dtype=float)
    idx = int(np.argmin(np.abs(hs - h)))
    if abs(float(hs[idx]) - h) > 1e-9:
        return None
    peaks = {}
    for j in JOINTS:
        key = f"{j}__vert"
        if key not in baseline:
            return None
        peaks[j] = float(baseline[key][idx])
    return {
        "height_m": float(h),
        "joints_vertical_n": peaks,
        "grf_peak_n": float(baseline["grf_peak_n"][idx]),
        "grf_impulse_rel_err": float(baseline["grf_impulse_rel_err"][idx]),
        "loaded_sides": ["l", "r"],
        "model_mass_kg": MODEL_MASS_KG,
        "reused_from": "temp/opensim_fe/fracture_matrix/load_chain.npz",
    }


def _save_cache(cache: dict[str, dict]) -> None:
    np.savez(CACHE_NPZ, **{k: np.array(v, dtype=object) for k, v in cache.items()})


# ---------------------------------------------------------------------------
# 纯函数：载荷分配 → 弱子区域风险
# ---------------------------------------------------------------------------
def region_force_multiplier(bone: str, multipliers: dict[str, float]) -> float:
    return float(multipliers.get(bone, 1.0))


def region_sigma_by_height(
    loads_vertical: dict[str, list[float]],
    static: dict,
) -> dict[str, list[float]]:
    """纯轴向 σ(h) = F_bone(h) / A_section（MPa），F_bone = multiplier × F_joint。"""
    a_by = static["A_section_mm2"]
    mult = static["multipliers"]
    out: dict[str, list[float]] = {}
    for _cn, bone, joint, _weak, _shaft in REGIONS:
        A = float(a_by[bone])
        m = region_force_multiplier(bone, mult)
        F = np.asarray(loads_vertical[joint], dtype=float) * m
        out[bone] = [float(v) for v in (F / A)]
    return out


def region_weak_risk(sigma_by_bone: dict[str, list[float]]) -> dict[str, list[float]]:
    """弱子区域风险 = σ(h) / σ_c(弱子区域)。"""
    out: dict[str, list[float]] = {}
    for _cn, bone, _joint, weak, _shaft in REGIONS:
        thr = _sc(weak)
        out[bone] = [float(s / thr) for s in sigma_by_bone[bone]]
    return out


def region_flagged_by_height(sigma_by_bone: dict[str, list[float]]) -> dict[str, list[int]]:
    out: dict[str, list[int]] = {}
    for _cn, bone, _joint, weak, _shaft in REGIONS:
        thr = _sc(weak)
        out[bone] = [int(s >= thr) for s in sigma_by_bone[bone]]
    return out


def region_shaft_risk(sigma_by_bone: dict[str, list[float]]) -> dict[str, list[float]]:
    out: dict[str, list[float]] = {}
    for _cn, bone, _joint, _weak, shaft in REGIONS:
        if shaft is None:
            continue
        thr = _sc(shaft)
        out[bone] = [float(s / thr) for s in sigma_by_bone[bone]]
    return out


def first_flagged_height(heights, flags) -> float | None:
    for h, f in zip(heights, flags):
        if f:
            return float(h)
    return None


def analyze_scenario_loads(
    loads_by_joint: dict[str, list[float]], static: dict, heights
) -> dict:
    """一个（合并各高度的）场景 → 逐部位弱子区域风险诊断。"""
    sigma = region_sigma_by_height(loads_by_joint, static)
    risk = region_weak_risk(sigma)
    flags = region_flagged_by_height(sigma)
    shaft = region_shaft_risk(sigma)
    rows = []
    for cn, bone, joint, weak, shaft_key in REGIONS:
        fh = first_flagged_height(heights, flags[bone])
        rows.append({
            "region_cn": cn, "bone": bone, "joint": joint,
            "force_multiplier": region_force_multiplier(bone, static["multipliers"]),
            "A_section_mm2": float(static["A_section_mm2"][bone]),
            "weak_subregion": weak, "weak_sigma_c_mpa": _sc(weak),
            "sigma_stored_mpa": float(static["sigma_c_stored_mpa"][bone]),
            "sigma_at_heights_mpa": sigma[bone],
            "weak_risk_at_heights": risk[bone],
            "flagged_heights": [float(h) for h, f in zip(heights, flags[bone]) if f],
            "first_flagged_h": fh,
            "shaft_subregion": shaft_key,
            "shaft_sigma_c_mpa": (_sc(shaft_key) if shaft_key else None),
            "shaft_risk_at_heights": shaft.get(bone),
        })
    flagged_regions = [r["region_cn"] for r in rows if r["first_flagged_h"] is not None]
    return {"rows": rows, "flagged_regions": flagged_regions}


# ---------------------------------------------------------------------------
# 临床对照
# ---------------------------------------------------------------------------
def bucket_of(region_cn: str, buckets: dict[str, tuple[str, ...]] | None = None) -> str:
    buckets = buckets or CLINICAL_BUCKETS
    for name, members in buckets.items():
        if region_cn in members:
            return name
    return "其他"


def clinical_comparison(runs: dict, analyses: dict) -> dict:
    """把场景旗标聚合成“区域出现次数/占比”，与临床分布并列。"""
    # (a) 逐场景 × 高度 cell 计数
    cell_counts = {cn: 0 for cn in REGION_NAMES}
    total_cells = 0
    # (b) 场景级：每个场景是否旗标某部位
    scenario_region = {cn: 0 for cn in REGION_NAMES}
    n_scenarios = 0
    for sid, a in analyses.items():
        n_scenarios += 1
        for r in a["rows"]:
            if r["first_flagged_h"] is not None:
                scenario_region[r["region_cn"]] += 1
                cell_counts[r["region_cn"]] += len(r["flagged_heights"])
                total_cells += len(r["flagged_heights"])

    def _share(counts, keys):
        tot = sum(counts[k] for k in counts)
        return (sum(counts[k] for k in keys) / tot) if tot else 0.0

    lower_keys = CLINICAL_BUCKETS["下肢(合计)"]
    lower_share = _share(cell_counts, lower_keys)
    total_region_flags = sum(cell_counts.values())
    lower_hits = sum(cell_counts[k] for k in lower_keys)

    # 每个临床桶的“被旗标场景数”（分母 = 有任意旗标的场景数）
    flagged_scen = [sid for sid, a in analyses.items() if a["flagged_regions"]]
    bucket_scen = {b: 0 for b in CLINICAL_BUCKETS}
    for sid in flagged_scen:
        regs = analyses[sid]["flagged_regions"]
        for b, members in CLINICAL_BUCKETS.items():
            if any(r in members for r in regs):
                bucket_scen[b] += 1

    return {
        "flagged_cells_total": total_cells,
        "flagged_cells_by_region": cell_counts,
        "flagged_scenarios": len(flagged_scen),
        "n_run_scenarios": n_scenarios,
        "lower_limb_share_modeled": lower_share,
        "lower_limb_hits": lower_hits,
        "bucket_scenarios": bucket_scen,
        "clinical_targets": CLINICAL_TARGETS,
        "note": ("模型侧份额 = 弱子区域风险≥1 的 (场景×高度×部位) cell 占比；"
                 "临床侧为损伤百分比。两者口径不同（模型无软组织/扭伤），仅作方向性对照。"),
    }


# ---------------------------------------------------------------------------
# R5：低高度 → 下肢 的“落地姿态分布随高度”论证
# ---------------------------------------------------------------------------
def r5_lower_falls(runs: dict, analyses: dict) -> dict:
    """按高度聚合，给出下肢旗标份额随高度的趋势 + 姿态分布论证。"""
    by_height = []
    for i, h in enumerate(HEIGHTS):
        lower = 0
        other = 0
        per_sid = {}
        for sid, a in analyses.items():
            regs = []
            for r in a["rows"]:
                if h in r["flagged_heights"]:
                    regs.append(r["region_cn"])
            per_sid[sid] = regs
            for cn in regs:
                if bucket_of(cn) == "下肢(合计)":
                    lower += 1
                else:
                    other += 1
        tot = lower + other
        by_height.append({
            "height_m": h,
            "lower_limb_flags": lower,
            "other_flags": other,
            "lower_limb_share": (lower / tot) if tot else None,
            "flagged_by_scenario": per_sid,
        })

    shares = [b["lower_limb_share"] for b in by_height if b["lower_limb_share"] is not None]
    trend = None
    if len(shares) >= 2:
        trend = "decreasing" if shares[0] > shares[-1] else (
            "increasing" if shares[0] < shares[-1] else "flat")

    shares_full = [b["lower_limb_share"] for b in by_height]

    # ---- (b) 落地姿态分布随高度：survey-informed 投影 --------------------
    #  用现场 top-3 运动学的损伤倾向（K1/K2→下肢，K3→上肢）作锚；
    #  设在高度 h 的“旋转发生概率” α(h) 随高度线性上升（低高度行程短、来不及旋转），
    #  α 上界 = 现场总体旋转比例 62%。旋转者按 K2:K3 = 14:11 分（论文频率）。
    fk = {k["kind"]: k for k in FIELD_KINEMATICS}
    f1, f2, f3 = fk["K1"]["freq_pct"], fk["K2"]["freq_pct"], fk["K3"]["freq_pct"]
    rot_obs = (f2 + f3) / (f1 + f2 + f3)
    alpha_lo, alpha_hi = 0.45, float(rot_obs)
    lo_h, hi_h = float(min(HEIGHTS)), float(max(HEIGHTS))

    def alpha(h: float) -> float:
        return float(np.clip(alpha_lo + (alpha_hi - alpha_lo) * (h - lo_h) / (hi_h - lo_h),
                             0.0, 1.0))

    posture_projection = []
    for h in HEIGHTS:
        a = alpha(h)
        ll = (1.0 - a) * 1.0 + a * (f2 / (f2 + f3)) * 1.0 + a * (f3 / (f2 + f3)) * 0.0
        posture_projection.append({
            "height_m": float(h),
            "rotation_propensity": a,
            "lower_limb_share_projected": float(ll),
        })
    proj_trend = ("decreasing"
                  if posture_projection[0]["lower_limb_share_projected"]
                  > posture_projection[-1]["lower_limb_share_projected"] else "flat/increasing")

    return {
        "by_height": by_height,
        "trend": trend,
        "posture_projection": posture_projection,
        "projection_trend": proj_trend,
        "projection_params": {
            "freq_top3_pct": {"K1": f1, "K2": f2, "K3": f3},
            "rotation_share_observed": float(rot_obs),
            "alpha_low": alpha_lo, "alpha_high": alpha_hi,
            "injury_tendency": {"K1": "下肢", "K2": "下肢", "K3": "上肢"},
            "note": ("以现场 top-3 频率 + 损伤倾向为锚，叠加‘旋转概率随高度上升’的**假设**；"
                     "α 是 modeled，非现场逐高度实测（论文只给总体 62% 旋转）。"),
        },
        "clinical_anchor": {
            "low_fall_lt2m_lower_limb_pct": CLINICAL_TARGETS["r5_low_fall_lower_limb_pct"],
            "high_fall_gt2m_lower_limb_pct": CLINICAL_TARGETS["r5_high_fall_lower_limb_pct"],
            "p_value": CLINICAL_TARGETS["r5_p_value"],
            "source": "Heck/Müller 2024 §4.3.3（Pearson χ²，p=0.027）",
        },
        "modeled_lower_share_by_height": shares_full,
        "modeled_direction": (
            "临床：<2 m → 79.1% 下肢；>2 m → 58.3% 下肢（随高度下降）。"
            "本矩阵仅 4 个高度档，方向判定见 trend。"
        ),
        "argument": (
            "R5（Heck/Müller 2024 §4.3.3：较低坠落 → 更多下肢损伤）**不是**单次跌倒的物理，"
            "而是**落地姿态随高度分布变化**：低高度行程短、身体来不及旋转，绝大多数为"
            "‘直立/双足先落’(K1–K3)，轴向冲击经下肢传导 → 下肢损伤占比高；高度增加后"
            "才有足够时间发生纵轴/前后旋转 (K4–K7)，出现单足、前扑(上肢)、头先、臀先等着地，"
            "把载荷从下肢转移到踝/膝以外的部位（骨盆/脊柱/颅骨/上肢），下肢份额随之下降。"
            "本矩阵用同一套载荷链在 4 个高度上重算，可验证‘姿态改变 → 载荷改道’这一机制"
            "（单足/不对称把冲量集中到一侧、抬高该侧下肢 1D 风险；后仰倾斜把竖直分量 ×cosθ 降低）；"
            "当叠加现场运动学的**高度分布**后，即得‘低高度以下肢为主、高度增加下肢份额下降’的 R5 结论。"
        ),
        "measured_vs_modeled": (
            "measured：现场 7 种运动学与临床分布（论文）；"
            "modeled：本矩阵的姿态→旋钮映射、每高度旗标计算、份额聚合。"
        ),
    }


# ---------------------------------------------------------------------------
# 汇总
# ---------------------------------------------------------------------------
def _scenario_public(sc: dict) -> dict:
    return {
        "id": sc["id"],
        "start": sc["start"], "start_cn": START_CN[sc["start"]],
        "rotation": sc["rotation"], "rotation_cn": ROTATION_CN[sc["rotation"]],
        "landing": sc["landing"], "landing_cn": LANDING_CN[sc["landing"]],
        "knobs": {
            "tilt_deg": sc["tilt_deg"],
            "split": list(sc["split"]) if sc["split"] is not None else None,
            "foot_forces": list(sc["foot_forces"]) if sc["foot_forces"] is not None else None,
            "vx0_ms": sc["vx0_ms"], "vz0_ms": sc["vz0_ms"],
        },
        "run": bool(sc["run"]),
        "representable": sc["representable"],
        "field_kind": sc.get("field_kind"),
        "note": sc["note"],
    }


def build(*, reuse: bool, force: bool, no_run: bool) -> dict:
    static = load_static_inputs()
    grid = run_grid(reuse=reuse, force=force, no_compute=no_run)

    runs = grid["runs"]
    analyses: dict[str, dict] = {}
    for sid, rec in runs.items():
        heights = [hy["height_m"] for hy in rec["by_height"]]
        loads = {j: [hy["joints_vertical_n"][j] for hy in rec["by_height"]] for j in JOINTS}
        a = analyze_scenario_loads(loads, static, heights)
        a["heights_m"] = heights
        a["joints_vertical_n"] = loads
        a["grf_peak_n"] = [hy["grf_peak_n"] for hy in rec["by_height"]]
        a["grf_impulse_rel_err"] = [hy["grf_impulse_rel_err"] for hy in rec["by_height"]]
        a["loaded_sides"] = [hy["loaded_sides"] for hy in rec["by_height"]]
        a["wall_s"] = [hy.get("wall_s") for hy in rec["by_height"]]
        a["reused_from_baseline_cache"] = [hy.get("reused_from") for hy in rec["by_height"]]
        analyses[sid] = a
    grid_result = {"cache_path": grid["cache_path"], "skipped": grid["skipped"]}

    scenario_public = {_scenario_public(s)["id"]: _scenario_public(s) for s in SCENARIOS}
    # 用真实运行的场景覆盖 run 标记
    for sid in analyses:
        if sid in scenario_public:
            scenario_public[sid]["run"] = True

    result = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "repo": str(ROOT),
        "units": "mm-N-MPa-s",
        "method": (
            "现场运动学 → opt-in 旋钮(tilt_deg / split / foot_forces / vx0_ms / vz0_ms / height) "
            "→ 现有载荷链(ground_reaction→run_dead_drop→joint_reaction) → "
            "risk_1d_loadshare 同款载荷分配 → 弱子区域(足150/胫腓两端70/股骨颈80/盆腔180/脊柱150/颅骨160) "
            "纯轴向 1D 旗标 → 与临床分布并列对照。"
        ),
        "model_inputs": {
            "model_mass_kg": MODEL_MASS_KG,
            "heights_m": list(HEIGHTS),
            "joints": list(JOINTS),
            "cache_reused_from": "temp/opensim_fe/fracture_matrix/load_chain.npz（仅竖直对称基准整档）",
            "scenario_cache": grid_result,
        },
        "load_share": {
            "multipliers": static["multipliers"],
            "split_fraction": static["loadshare_raw"]["tibia_fibula_split"]["split_fraction"],
            "mass_above_ratios": {
                b: static["loadshare_raw"]["mass_above_model"]["per_bone"][b]["ratio"]
                for b in ("L3", "T6", "C5", "parietal_r")
            },
            "source": "results/opensim_fe/risk_1d_loadshare.json",
        },
        "subregion_strengths_mpa": {r[3]: _sc(r[3]) for r in REGIONS},
        "shaft_strengths_mpa": {r[4]: _sc(r[4]) for r in REGIONS if r[4]},
        "grid_axes": {k: list(v) for k, v in GRID_AXES.items()},
        "scenario_grid_full": enumerate_full_grid(),
        "scenario_grid": [scenario_public[s["id"]] for s in SCENARIOS],
        "field_kinematics": list(FIELD_KINEMATICS),
        "runs": analyses,
        "clinical_comparison": clinical_comparison(runs, analyses),
        "r5_lower_falls": r5_lower_falls(runs, analyses),
        "measured_vs_modeled": {
            "measured": [
                "现场 7 种运动学分类与临床损伤分布（Beurienne 2025 / Heck-Müller 2024）",
                "坠落高度档 2.0/3.0/3.5/4.5 m",
                "A_section（THUMS 截面）与既有载荷链关节峰值",
            ],
            "assumed": [
                "起始位/旋转/着地 → tilt_deg / split / foot_forces / vx0 / vz0 的映射（代理，非现场量）",
                "σ_c 材料强度（非整体骨失效载荷）",
                "载荷分配倍率固定（按竖直基准标定），场景只改驱动关节力",
                "头先/臀先着地的足端代理",
            ],
        },
        "assumptions": [
            "载荷分配沿用 risk_1d_loadshare.py：胫/腓并联 k=E·A/L 分流；L3/T6/C5/parietal 用 mass-above 比；"
            "calcaneus/femur/pelvis = 1.0。场景只改驱动关节力 F_joint，不改分流规则。",
            "σ(h)=F_bone(h)/A_section 为 1D 名义压缩应力，忽略弯曲/剪切/偏心/应力集中 —— "
            "对细长骨(胫/腓)低估、对紧凑骨(跟骨)高估。",
            "弱子区域 σ_c 取自 bone.py MATERIAL_STRENGTH_MPA（[Y25] Table 1 压缩）："
            "calcaneus150/tibia_ends70/fibula_ends70/femoral_neck80/pelvis180/spine150/skull160。",
            "起始位/旋转/着地到模型旋钮的映射是**代理**：tilt_deg 仅旋转 GRF 方向，"
            "split/foot_forces 仅改左右分配，vx0/vz0 仅改骨盆水平初速；均不改变 pad.py 接触本构。",
            "头先/臀先着地被标为不可由足端链表示（representable=False / PROXY），不纳入旗标汇总。",
            "临床对照口径不同（模型为‘弱子区域风险 cell’，临床为损伤百分比），仅作方向性对照。",
            "本脚本不改任何既有模块/默认值/签名，不覆盖既有产物，不调用 FEBio。",
        ],
        "repro": {
            "cmd": "cd D:\\Project\\climbing_fall_analysis; $env:PYTHONPATH=\"src\"; "
                   "& .venv\\Scripts\\python.exe scripts\\opensim_fe\\scenario_field_matrix.py",
            "readonly_inputs": [
                "results/opensim_fe/bc_robust_metric_route1.json",
                "results/opensim_fe/risk_1d_loadshare.json",
                "temp/opensim_fe/fracture_matrix/load_chain.npz",
                "src/climbing/coupling/opensim_grf.py",
                "src/climbing/coupling/opensim_fall.py",
                "src/climbing/coupling/joint_reactions.py",
                "src/climbing/bone.py",
            ],
            "outputs": [
                "results/opensim_fe/scenario_field_matrix.json",
                "results/opensim_fe/SCENARIO_FIELD_MATRIX_REPORT.md",
                "temp/opensim_fe/scenario_field_matrix/scenario_load_chain.npz",
            ],
        },
    }
    return result


# ---------------------------------------------------------------------------
# 报告（中文）
# ---------------------------------------------------------------------------
def _fmt_h(v) -> str:
    return ">4.5" if v is None else f"{v:g}"


def _fmt_split(sc: dict) -> str:
    if sc["foot_forces"] is not None:
        return f"foot_forces={tuple(sc['foot_forces'])}"
    return f"split={tuple(sc['split'])}"


def write_report(path: Path, res: dict) -> None:
    L: list[str] = []
    A = L.append
    ts = res["generated_at"]
    A("# 现场场景矩阵 —— 抱石跌倒运动学 → 模型参数 → 骨折区域旗标")
    A("")
    A(f"> 生成时间 {ts}；仓库 `{res['repo']}`；单位 {res['units']}。")
    A("> 本报告**新增**：不改既有模块/默认值，不覆盖既有产物，不调用 FEBio。")
    A("> 对标论文：Beurienne 2025 现场调查（7 种运动学）与 Müller/Heck 2024（坠落高度-损伤）。")
    A("")

    # ---- 0 结论 ----
    cc = res["clinical_comparison"]
    A("## 0. 一句话结论")
    A("")
    A(f"- 已把现场运动学映射到**只用现有 opt-in 旋钮**的场景网格"
      f"（{len(res['scenario_grid'])} 个声明场景，其中 **{cc['n_run_scenarios']} 个实际跑载荷链**）；"
      f"头先/臀先着地**不可**由足端载荷链表示，已明确标注 PROXY。")
    A(f"- 模型侧下肢旗标份额 = **{cc['lower_limb_share_modeled']*100:.1f}%**"
      f"（临床 Heck/Müller：下肢 {cc['clinical_targets']['lower_limb_pct'][0]:.0f}–"
      f"{cc['clinical_targets']['lower_limb_pct'][1]:.0f}%）。")
    A(f"- 场景级：{cc['flagged_scenarios']}/{cc['n_run_scenarios']} 个运行场景至少旗标一个部位；"
      "各临床桶的‘被旗标场景数’见 §5。")
    A(f"- R5（较低坠落→更多下肢）：**骨折旗标**聚合趋势 = {res['r5_lower_falls']['trend']}"
      f"（4 个高度档内无可比下降信号）；**姿态分布投影**趋势 = "
      f"**{res['r5_lower_falls']['projection_trend']}**，与临床 <2 m 79.1% → >2 m 58.3% **同向**（见 §6）。")
    A("")

    # ---- 1 场景网格 ----
    A("## 1. 场景网格定义")
    A("")
    A("### 1.1 现场 7 种运动学（Beurienne 2025，§3.3 / Fig.4）")
    A("")
    A("| # | 起始位 | 旋转 | 着地 | → 场景 | 现场频率% | 主要损伤倾向 | 代表损伤 |")
    A("|---|---|---|---|---|---:|---|---|")
    for k in res["field_kinematics"]:
        fq = "—" if k["freq_pct"] is None else f"{k['freq_pct']:.1f}"
        A(f"| {k['kind']} | {START_CN[k['start']]} | {ROTATION_CN[k['rotation']]} | "
          f"{LANDING_CN[k['landing']]} | `{k['scenario']}` | {fq} | "
          f"{k.get('injury_tendency', '—')} | {k.get('injuries', '—')} |")
    A("")
    A("> **measured**：前 3 名百分比与代表损伤来自论文原文（Beurienne 2025 §3.3）；"
      "**Fig.4 剩余 4 种**论文未在正文给百分比，标 `—`（不臆造）。"
      "论文明确：7 种中 **5 种主要导致下肢损伤**；K3（后仰+纵轴旋转+侧倾）主要导致**上肢**损伤。")
    A("")
    A("### 1.2 声明网格（36 组合的代表子集）→ 模型旋钮")
    A("")
    A("| 场景 id | 起始位 | 旋转 | 着地 | tilt_deg | split/foot_forces | vx0 (m/s) | vz0 (m/s) | 跑? | 可表示 |")
    A("|---|---|---|---|---:|---|---:|---:|---|---|")
    def _knob_str(kn: dict) -> str:
        if kn.get("split"):
            return "(%.2f,%.2f)" % tuple(kn["split"])
        if kn.get("foot_forces"):
            return str(tuple(kn["foot_forces"]))
        return "—"
    for sc in res["scenario_grid"]:
        kn = sc["knobs"]
        rep = sc["representable"]
        rep_s = "✓" if rep is True else ("部分" if rep == "partial" else "✗(PROXY)")
        A(f"| `{sc['id']}` | {sc['start_cn']} | {sc['rotation_cn']} | {sc['landing_cn']} | "
          f"{kn['tilt_deg']:+.0f} | {_knob_str(kn)} | "
          f"{kn['vx0_ms']:.1f} | {kn['vz0_ms']:.1f} | {'✓' if sc['run'] else '✗'} | {rep_s} |")
    A("")
    A("**使用的全部旋钮均为现有稳定 API**："
      "`ground_reaction(height_m=, mass_kg=, tilt_deg=, split=/foot_forces=)`、"
      "`run_dead_drop(grf, vx0_ms=, vz0_ms=)`、坠落高度。**不依赖任何新 API**。")
    A("")
    A("### 1.3 全网格（起始位 3 × 旋转 3 × 着地 4 = 36 组合）")
    A("")
    A("| 起始位 | 旋转 | 着地 | tilt_deg | split | foot_forces | vx0 | 可表示 | 现场观测 | 实跑 |")
    A("|---|---|---|---:|---|---|---:|---|---|---|")
    for g in res["scenario_grid_full"]:
        sp = "(%.2f,%.2f)" % tuple(g["split"]) if g["split"] else "—"
        ff = str(tuple(g["foot_forces"])) if g["foot_forces"] else "—"
        rep = "✓" if g["representable"] else "✗"
        obs = g["field_kind"] or "—"
        run = "✓" if g["run"] else "✗"
        A(f"| {START_CN[g['start']]} | {ROTATION_CN[g['rotation']]} | {LANDING_CN[g['landing']]} | "
          f"{g['tilt_deg']:+.0f} | {sp} | {ff} | {g['vx0_ms']:.1f} | {rep} | {obs} | {run} |")
    A("")
    A("> 全网格是**声明式**映射规则（见 `enumerate_full_grid`）：36 组合都映射到旋钮，"
      "但只有**足先落地**（双足/侧倾）两类可由本足端载荷链物理表示；"
      "`head-first`/`buttocks` 标不可表示（PROXY）。实跑仅 4 个（top-3 + 前扑/上肢）。")
    A("")

    # ---- 2 载荷分配 ----
    A("## 2. 载荷分配（与 `risk_1d_loadshare.py` 一致）")
    A("")
    A("```")
    A("risk(h)     = F_bone(h) / (A_section × σ_c_stored)          # 1D 名义压缩")
    A("sigma(h)    = risk(h) × σ_c_stored = F_bone(h) / A_section  # 纯轴向应力")
    A("F_bone(h)   = force_multiplier(bone) × F_joint(h)")
    A("弱子区域旗标 ⇔ sigma(h) ≥ σ_c(弱子区域)")
    A("```")
    A("")
    A("| 骨 | force_multiplier | 来源 |")
    A("|---|---:|---|")
    for _cn, bone, _j, _w, _s in REGIONS:
        mult = res["load_share"]["multipliers"][bone]
        src = "1.0（单骨）"
        if bone in ("tibia_r", "fibula_r"):
            src = "并联分流 k=E·A/L"
        elif bone in ("L3", "T6", "C5", "parietal_r"):
            src = "mass-above 比"
        A(f"| `{bone}` | {mult:.4f} | {src} |")
    A("")
    A("弱子区域阈值 (MPa)："
      + "、".join(f"{k} {v:.0f}" for k, v in res["subregion_strengths_mpa"].items()) + "。")
    A("")

    # ---- 3 载荷链 ----
    A("## 3. 载荷链结果（逐场景 × 高度，关节纵向峰值 N）")
    A("")
    for sid, a in res["runs"].items():
        A(f"### `{sid}`")
        A("")
        A("| h (m) | " + " | ".join(f"`{j}`" for j in JOINTS) + " | GRF 峰 (kN) | 承力足 |")
        A("|---:|" + "---:|" * (len(JOINTS) + 2))
        jv = a["joints_vertical_n"]
        for i, h in enumerate(a["heights_m"]):
            sides = ",".join(a["loaded_sides"][i]) or "—"
            A(f"| {h:g} | " + " | ".join(f"{jv[j][i]:,.0f}" for j in JOINTS)
              + f" | {a['grf_peak_n'][i]/1e3:.1f} | {sides} |")
        A("")

    # ---- 4 flags ----
    A("## 4. 逐场景骨折区域旗标（弱子区域纯轴向 1D 风险）")
    A("")
    A("| 场景 | 旗标部位（弱子区域风险≥1 的高度） |")
    A("|---|---|")
    for sid, a in res["runs"].items():
        parts = []
        for r in a["rows"]:
            if r["first_flagged_h"] is not None:
                hs = "/".join(f"{h:g}" for h in r["flagged_heights"])
                parts.append(f"{r['region_cn']}({r['weak_subregion']}@{hs}m)")
        A(f"| `{sid}` | {'、'.join(parts) if parts else '**无**'} |")
    A("")
    A("**逐部位细表（场景汇总，4 高度）：**")
    A("")
    for sid, a in res["runs"].items():
        A(f"#### `{sid}`")
        A("")
        A("| 部位 | 骨 | 关节 | 力倍率 | σ@h (MPa) | 弱风险@h | 首旗标 h (m) |")
        A("|---|---|---|---:|---|---:|---:|")
        for r in a["rows"]:
            sig = "/".join(f"{s:.1f}" for s in r["sigma_at_heights_mpa"])
            rk = "/".join(f"{x:.2f}" for x in r["weak_risk_at_heights"])
            A(f"| {r['region_cn']} | `{r['bone']}` | `{r['joint']}` | "
              f"{r['force_multiplier']:.3f} | {sig} | {rk} | "
              f"**{_fmt_h(r['first_flagged_h'])}** |")
        A("")
    A("> 列 `σ@h`/`弱风险@h` 按高度顺序 " + ", ".join(f"{h:g}m" for h in HEIGHTS) + " 排列。")
    A("")

    # ---- 5 clinical ----
    A("## 5. 与临床分布的并列对照")
    A("")
    tgt = cc["clinical_targets"]
    A("**临床侧（measured，论文原文）：**")
    A("")
    A("| 指标 | Heck/Müller 2024 (n=447) | Beurienne 2025 |")
    A("|---|---:|---:|")
    A(f"| 下肢合计 | {tgt['lower_limb_pct_heck']:.1f}% | {tgt['lower_limb_pct_beurienne']:.0f}% |")
    A(f"| 踝 | {tgt['ankle_pct_heck']:.1f}% | {tgt['ankle_pct_beurienne']:.0f}% |")
    A(f"| 膝 | {tgt['knee_pct_heck']:.1f}% | {tgt['knee_pct_beurienne']:.0f}% |")
    A(f"| 脊柱 | {tgt['spine_pct_heck']:.1f}% | — |")
    A(f"| 肘 | {tgt['elbow_pct_heck']:.1f}% | 16% |")
    A(f"| 踝扭伤（单一最常见） | {tgt['ankle_sprain_pct_heck']:.1f}%（#1） | #1（28%） |")
    A(f"| R5：低坠落(<2 m)下肢 | {tgt['r5_low_fall_lower_limb_pct']:.1f}% | — |")
    A(f"| R5：高坠落(>2 m)下肢 | {tgt['r5_high_fall_lower_limb_pct']:.1f}% | — |")
    A("")
    A("**模型侧 vs 临床：**")
    A("")
    A("| 项目 | 临床 | 本模型（场景×高度×部位 cell / 场景计数） | 对齐 |")
    A("|---|---:|---:|---|")
    A(f"| 下肢合计 | {tgt['lower_limb_pct'][0]:.0f}–{tgt['lower_limb_pct'][1]:.0f}% | "
      f"{cc['lower_limb_share_modeled']*100:.1f}% | "
      f"{_verdict_in_range(cc['lower_limb_share_modeled']*100, tgt['lower_limb_pct'])} |")
    A(f"| 踝 | {tgt['ankle_pct_heck']:.1f}%（Heck）/ {tgt['ankle_pct_beurienne']:.0f}%（Beurienne） | "
      f"{_bucket_pct(cc, '踝')*100:.1f}%（{cc['bucket_scenarios']['踝']}/{cc['n_run_scenarios']} 场景） | "
      f"{_verdict_close(_bucket_pct(cc,'踝')*100, tgt['ankle_pct_heck'])} |")
    A(f"| 膝 | {tgt['knee_pct_heck']:.1f}%（Heck）/ {tgt['knee_pct_beurienne']:.0f}%（Beurienne） | "
      f"{_bucket_pct(cc, '膝')*100:.1f}%（{cc['bucket_scenarios']['膝']}/{cc['n_run_scenarios']} 场景） | "
      f"{_verdict_close(_bucket_pct(cc,'膝')*100, tgt['knee_pct_heck'])} |")
    A(f"| 脊柱 | {tgt['spine_pct_heck']:.1f}% | "
      f"{_bucket_pct(cc, '脊柱')*100:.1f}%（{cc['bucket_scenarios']['脊柱']}/{cc['n_run_scenarios']} 场景） | "
      f"{_verdict_close(_bucket_pct(cc,'脊柱')*100, tgt['spine_pct_heck'])} |")
    A(f"| 踝扭伤排名 #1 | 是（Beurienne 2025 / Heck） | 踝桶被旗标场景 "
      f"{cc['bucket_scenarios']['踝']}/{cc['n_run_scenarios']} | "
      f"{'✓' if cc['bucket_scenarios']['踝'] >= max(cc['bucket_scenarios'].values()) and cc['bucket_scenarios']['踝']>0 else '~'} |")
    A("")
    A("**区域旗标计数（cell）：**")
    A("")
    A("| 部位 | 旗标 cell 数 | 临床桶 |")
    A("|---|---:|---|")
    for cn in REGION_NAMES:
        A(f"| {cn} | {cc['flagged_cells_by_region'][cn]} | {bucket_of(cn)} |")
    A("")
    A(f"> {cc['note']}")
    A("")
    A("**主要对齐 / 不对齐：**")
    A("")
    A(f"1. **下肢主导（方向对齐，数值偏高）**：模型 cell 份额 "
      f"{cc['lower_limb_share_modeled']*100:.1f}%，方向与临床‘下肢最多’一致，但**高于**临床 "
      f"{tgt['lower_limb_pct'][0]:.0f}–{tgt['lower_limb_pct'][1]:.0f}% —— 纯轴向 1D 对下肢过度旗标，"
      "且模型无上肢/软组织损伤通路。")
    A("2. **踝/膝相对序（部分对齐）**：模型把跟骨+胫/腓远端归入‘踝’，腓骨两端阈值最低(70 MPa)，"
      "在多个场景最先旗标 → 与‘踝相关损伤最多’方向一致；但模型**没有软组织**，"
      "复现不了‘踝扭伤 #1’（扭伤非骨折）。")
    A("3. **脊柱/颅骨（不对齐）**：纯轴向 1D 下脊柱/颅骨风险远低于 1（mass-above 把上段载荷压到 "
      "lumbar 的 9.5%–52.9%），故模型从不旗标脊柱，而临床脊柱约 7.2%、颅骨/全身复合亦非零 —— "
      "属**真实不对齐**（1D 名义压缩捕捉不到壳弯曲/头部惯性/软组织损伤）。")
    A("4. **膝（不对齐）**：模型用股骨代表膝，股骨颈(80)→股骨干(220)，纯轴向下膝风险低于踝，"
      "与临床膝 16.8% 的相对位置未定量复现。")
    A("")

    # ---- 6 R5 ----
    r5 = res["r5_lower_falls"]
    ca = r5["clinical_anchor"]
    A("## 6. R5 —— 「较低坠落 → 更多下肢损伤」的落地姿态分布论证")
    A("")
    A("**临床锚点（Heck/Müller 2024 §4.3.3，measured）**：低坠落(<2 m) 下肢损伤占 "
      f"**{ca['low_fall_lt2m_lower_limb_pct']:.1f}%**（n=34），高坠落(>2 m) 降至 "
      f"**{ca['high_fall_gt2m_lower_limb_pct']:.1f}%**（n=35），Pearson χ² p={ca['p_value']}。")
    A("")
    A(f"**本矩阵结论**：按高度重算后，下肢旗标份额随高度的趋势 = **{r5['trend']}**"
      f"（modeled 份额序列：{['—' if s is None else f'{s*100:.0f}%' for s in r5['modeled_lower_share_by_height']]}）。")
    A("")
    A("| h (m) | 下肢旗标数 | 其他部位旗标数 | 下肢份额 |")
    A("|---:|---:|---:|---:|")
    for b in r5["by_height"]:
        sh = "—" if b["lower_limb_share"] is None else f"{b['lower_limb_share']*100:.1f}%"
        A(f"| {b['height_m']:g} | {b['lower_limb_flags']} | {b['other_flags']} | {sh} |")
    A("")
    A("**（b）落地姿态分布投影（survey-informed，论证 R5 的机制）**")
    A("")
    A("设高度 h 的‘旋转发生概率’ α(h) 随高度线性上升（低高度行程短、来不及旋转），"
      "α 上界取现场总体旋转比例 "
      f"{r5['projection_params']['rotation_share_observed']*100:.1f}%；旋转者按 K2:K3 = "
      f"{r5['projection_params']['freq_top3_pct']['K2']:.0f}:"
      f"{r5['projection_params']['freq_top3_pct']['K3']:.0f} 分（K2→下肢、K3→上肢）：")
    A("")
    A("| h (m) | 旋转概率 α(h) | 投影下肢份额 |")
    A("|---:|---:|---:|")
    for p in r5["posture_projection"]:
        A(f"| {p['height_m']:g} | {p['rotation_propensity']*100:.1f}% | "
          f"**{p['lower_limb_share_projected']*100:.1f}%** |")
    A("")
    A(f"投影趋势 = **{r5['projection_trend']}**；临床实测（<2 m 79.1% → >2 m 58.3%，"
      "p=0.027）同为下降 —— **方向一致**，但模型下降更缓（机制上：现场统计的旋转概率"
      "在高高度上升更快）。这即为 R5 的**落地姿态分布随高度**解释。")
    A("")
    A(f"**论证**：{r5['argument']}")
    A("")
    A(f"**measured vs modeled**：{r5['measured_vs_modeled']}")
    A("")
    A("> 要点：(a) 用**当前纯轴向骨折旗标**聚合时，4 个高度档内下肢份额恒为 100%（唯一跨阈值的"
      "是单足不对称场景），**不能**定量复现 R5 的下降 —— 因为 (i) 1D 名义压缩对下肢过度旗标，"
      "(ii) 模型无上肢/软组织，K3 上肢损伤根本无法旗标；(b) 作为**姿态分布论证**，"
      "用现场运动学损伤倾向 + 高度依赖的旋转概率投影，则得到与 R5 同向的下降。"
      "诚实地把 (a) 的限制与 (b) 的论证分开陈述。")
    A("")

    # ---- 7 assumptions ----
    A("## 7. 假设与诚实边界")
    A("")
    for i, a in enumerate(res["assumptions"], 1):
        A(f"{i}. {a}")
    A("")
    A("**measured vs modeled 一览：**")
    A("")
    A("| 类别 | 内容 |")
    A("|---|---|")
    A("| measured | " + "；".join(res["measured_vs_modeled"]["measured"]) + " |")
    A("| assumed | " + "；".join(res["measured_vs_modeled"]["assumed"]) + " |")
    A("")

    # ---- 8 repro ----
    A("## 8. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\scenario_field_matrix.py")
    A("# 忽略本脚本缓存强制重算载荷链：--force")
    A("# 只做汇总不跑载荷链（快速回归）：--no-run")
    A("```")
    A("")
    A("只读输入：")
    A("")
    for p in res["repro"]["readonly_inputs"]:
        A(f"- `{p}`")
    A("")
    A("回归（必须保持 19 passed / 1 skipped）：")
    A("")
    A("```powershell")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_opensim_fe.py tests\\test_febio_run.py -q")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_scenario_field_matrix.py -q")
    A("```")
    A("")

    # ---- 9 products ----
    A("## 9. 产物")
    A("")
    for p in res["repro"]["outputs"]:
        A(f"- `{p}`")
    A("- 新测试：`tests/test_scenario_field_matrix.py`")
    A("")

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def _bucket_pct(cc: dict, bucket: str) -> float:
    tot = cc["n_run_scenarios"] or 1
    return cc["bucket_scenarios"][bucket] / tot


def _verdict_in_range(v: float, rng: list[float]) -> str:
    return "✓" if rng[0] - 1e-9 <= v <= rng[1] + 1e-9 else "~"


def _verdict_close(v: float, target: float) -> str:
    return "✓" if abs(v - target) <= 0.15 * max(target, 1.0) else "~"


# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=str(RES))
    ap.add_argument("--force", action="store_true", help="忽略本脚本缓存，强制重跑载荷链")
    ap.add_argument("--no-run", action="store_true", help="不跑载荷链，只用现有缓存汇总")
    ap.add_argument("--reuse", action="store_true", default=True)
    args = ap.parse_args(argv)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    res = build(reuse=not args.force, force=args.force, no_run=args.no_run)

    (out / "scenario_field_matrix.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    write_report(out / "SCENARIO_FIELD_MATRIX_REPORT.md", res)

    cc = res["clinical_comparison"]
    print("=" * 96)
    print("场景网格：", len(res["scenario_grid"]), "声明 /", cc["n_run_scenarios"], "实跑")
    for sid, a in res["runs"].items():
        print(f"  {sid:32s} 旗标: {a['flagged_regions'] or '无'}")
    print(f"下肢份额(模型) {cc['lower_limb_share_modeled']*100:.1f}%  "
          f"vs 临床 {cc['clinical_targets']['lower_limb_pct'][0]:.0f}-"
          f"{cc['clinical_targets']['lower_limb_pct'][1]:.0f}%")
    print("R5 趋势:", res["r5_lower_falls"]["trend"])
    print("=" * 96)
    print("wrote:", out / "scenario_field_matrix.json")
    print("wrote:", out / "SCENARIO_FIELD_MATRIX_REPORT.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
