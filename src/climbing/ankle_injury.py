"""Phase P2 —— 踝（P2 2D 模型的损伤判据）。

为什么需要这个模块
------------------
``injury.assess_boulder_fall`` 评估"全身压力 + HIC + 脊柱/下肢轴向力"，
``bone.assess_sites`` 评估"7 个骨部位骨折（material-strength × area）"，
二者都**不**识别"踝"这条 [B25] 第一大损伤通路（踝 40%，含踝扭伤 28% +
踝骨折 8%；下肢合计 67%）。

本模块把踝从结果里读出来，按与 ``injury.VERDICTS`` 同源的分档规则输出
``AnkleVerdict``。这是 P1 → P2 验证的关键钩子：把"踝"维度加入
``phase_p1_montecarlo.py`` 的损伤部位分布，让原本为 0 的踝列变得非空。

与上下游的契约
-------------
* **上游**（``pad2d.simulate_boulder_fall_2d``，并行 agent 正在写）：在
  ``BoulderFallResult2D`` 上输出 ``ankle_mode`` / ``peak_ankle_beta_rad``
  / ``peak_inversion_moment_nmm`` / ``peak_ligament_strain`` /
  ``ankle_slide_frac`` / ``ankle_bone_utilization`` / ``ankle_fracture`` /
  ``ankle_sprain``，并在 ``meta`` 里镜像同名标量。
  ``Posture2D`` 上输出 ``ankle_beta0_rad`` / ``beta_from_theta`` /
  ``mu_slide`` / ``ankle_d_lateral_m`` / ``ankle_d_medial_m`` /
  ``ligament_r_mm``（默认除 ``ligament_r_mm=22.0`` 外都为 0.0 ⇒ 默认路径
  下踝不参与，与 1D 模型 bit-identical —— 不变量 1，硬门槛）。
* **本模块**：仅读属性，不导入 ``pad2d``。若 ``result`` 上没有
  ``peak_ligament_strain``（1D ``BoulderFallResult``），抛 ``ValueError``，
  行为与 ``bone.py:301-307``（缺 ``impact_window_s``）同源。
* **下游**（``scripts/phase_p1_montecarlo.py``）：把踝列加入
  ``trial_rows`` / ``site_rows``，重新计算 ``site_distribution``，
  写入 ``results/p1_montecarlo.{txt,csv}`` 与 ``results/figures/p1_coverage.png``。

分档规则（与 ``injury.VERDICTS`` 同源 ``["低","中","高","极高"]``）
---------------------------------------------------------------
按最严重先吞的顺序：

* 骨折 ∧ 利用率 ≥ 2.0                              → 极高
* 骨折                                              → 高
* 扭伤 ∧ 峰值韧带应变 ≥ 0.2                        → 高（高应变扭伤）
* 扭伤                                              → 中
* 其余                                              → 低

0.2 的应变阈值选取：ATFL 失效应变 0.14（[R6]），CFL 失效 0.13，
0.2 ≈ 1.4× ATFL 失效 —— 越阈值后还有 ~40% 的过伸空间；
用于区分"刚好扭伤"（中）与"严重扭伤"（高）。

单位 mm–N·mm–rad–无量纲。
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "AnkleVerdict",
    "assess_ankle",
    "ankle_band",
    "ANKLE_MODE_VALUES",
]


# 踝模式的允许值（与 pad2d 合约中 BoulderFallResult2D.ankle_mode 同步）
# 参考 P2 输出合约："none" / "stuck" / "sliding" / "mixed"。
ANKLE_MODE_VALUES: frozenset[str] = frozenset({"none", "stuck", "sliding", "mixed"})


# --------------------------------------------------------------------------
# 数据类
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class AnkleVerdict:
    """一次 P2 仿真的踝评估结果。

    Fields
    ------
    posture : str
        姿势键名（与 ``pad.POSTURES`` 对齐）。
    height_m : float
        坠落高度（m）。
    beta0_rad : float
        初始踝旋后 / 外翻角（rad）。正 = 内翻（inversion），负 = 外翻（eversion）。
        取自 ``Posture2D.ankle_beta0_rad``，优先 ``result.meta``，缺省 0.0。
    mode : str
        踝运动学模式，集合为 :data:`ANKLE_MODE_VALUES`。
    peak_beta_rad : float
        ``|beta(t)|`` 的最大值（rad），取自 ``peak_ankle_beta_rad``。
    peak_moment_nmm : float
        峰值内翻力矩（N·mm），取自 ``peak_inversion_moment_nmm``。
    peak_ligament_strain : float
        峰值韧带应变（无量纲），取自 ``peak_ligament_strain``。
        ATFL 失效 ≈ 0.14，[R6]；CFL 失效 ≈ 0.13。
    bone_utilization : float
        踝骨利用率（峰值 / 阈值），取自 ``ankle_bone_utilization``。
        ≥ 1.0 ⇒ 骨折（``ankle_fracture=True``）。
    fracture : bool
        是否踝骨骨折（取自 ``ankle_fracture``）。
    sprain : bool
        是否韧带扭伤（取自 ``ankle_sprain``）。
    band : str
        损伤档位："低" / "中" / "高" / "极高"。由 :func:`ankle_band` 派生。
    note : str
        自由文本注释（机制、来源、契约异常）。
    """

    posture: str
    height_m: float
    beta0_rad: float
    mode: str
    peak_beta_rad: float
    peak_moment_nmm: float
    peak_ligament_strain: float
    bone_utilization: float
    fracture: bool
    sprain: bool
    band: str
    note: str = ""


# --------------------------------------------------------------------------
# 分档
# --------------------------------------------------------------------------
def ankle_band(
    sprain: bool,
    fracture: bool,
    bone_utilization: float,
    peak_ligament_strain: float,
) -> str:
    """按踝损伤严重度生成档位（与 ``injury.VERDICTS`` 同源）。

    规则（顺序敏感，最严重者先吞）：

    * 骨折 ∧ 利用率 ≥ 2.0                       → 极高
    * 骨折                                       → 高
    * 扭伤 ∧ 峰值韧带应变 ≥ 0.2                  → 高（高应变扭伤）
    * 扭伤                                       → 中
    * 其余                                       → 低

    Parameters
    ----------
    sprain, fracture
        同 :attr:`AnkleVerdict.fracture` / :attr:`AnkleVerdict.sprain`。
    bone_utilization
        踝骨折利用率（peak / threshold）。骨折 ⇒ 通常 ≥ 1.0。
    peak_ligament_strain
        峰值韧带应变。扭伤 ⇒ 通常 ≥ 0.14（ATFL 失效）。
    """
    if fracture and bone_utilization >= 2.0:
        return "极高"
    if fracture:
        return "高"
    if sprain and peak_ligament_strain >= 0.2:
        return "高"
    if sprain:
        return "中"
    return "低"


# --------------------------------------------------------------------------
# 鸭式探测
# --------------------------------------------------------------------------
# 探针字段：``peak_ligament_strain`` 是 P2 合约中**无重名**的浮点字段，
# 1D ``BoulderFallResult`` 没有，且合约明确要求它被镜像到 ``meta``。
# 用它做存在性判定，避免误判 1D 残留字段（``peak_force_n`` 等共用名）。
_PROBE_ATTR = "peak_ligament_strain"


def _peek_ankle_data(result) -> bool:
    """检测 ``result`` 是否携带 P2 踝数据（duck-typed，不导入 ``pad2d``）。"""
    if hasattr(result, _PROBE_ATTR):
        return True
    meta = getattr(result, "meta", None)
    if isinstance(meta, dict) and _PROBE_ATTR in meta:
        return True
    return False


def _read_field(result, key: str, default):
    """先在 ``result`` 上找属性，再到 ``meta``；都没有则返回 ``default``。"""
    if hasattr(result, key):
        return getattr(result, key)
    meta = getattr(result, "meta", None)
    if isinstance(meta, dict) and key in meta:
        return meta[key]
    return default


# --------------------------------------------------------------------------
# 主评估
# --------------------------------------------------------------------------
def assess_ankle(result, posture: str = "") -> AnkleVerdict:
    """对一次 P2 仿真做踝评估（duck-typed）。

    Parameters
    ----------
    result
        ``pad2d.simulate_boulder_fall_2d`` 的返回值（``BoulderFallResult2D``）。
        若结果是 1D ``BoulderFallResult``（无 P2 踝自由度）⇒ 抛 ``ValueError``。
    posture
        姿势键名（仅用于报告 / 诊断；缺省取 ``result.posture``）。

    Returns
    -------
    AnkleVerdict

    Raises
    ------
    ValueError
        结果对象无踝数据（``peak_ligament_strain`` 既不在属性里也不在
        ``meta`` 里）。镜像 ``src/climbing/bone.py:301-307`` 的行为
        （缺 ``impact_window_s`` 时抛 ``ValueError``）。
    """
    if not _peek_ankle_data(result):
        raise ValueError(
            "踝数据缺失：peak_ligament_strain 不在 result 上也不在 "
            "result.meta 中 —— 本评估需要 P2 踝自由度（Posture2D 的 "
            "ankle_beta0_rad / beta_from_theta 等 + BoulderFallResult2D 的 "
            "ankle_* 字段）。当前 result 是 1D BoulderFallResult，没有这些；"
            "请确认调用的是 pad2d.simulate_boulder_fall_2d（带 Posture2D 参数）。"
        )

    # ---- 读字段（含合约异常处理） ------------------------------------
    raw_mode = str(_read_field(result, "ankle_mode", "none"))
    if raw_mode not in ANKLE_MODE_VALUES:
        # 合约变更：只警示不断流，把未知模式回退为 'none' 并在 note 里记录
        mode_note = f"（未知模式 {raw_mode!r}，已回退到 'none'）"
        mode = "none"
    else:
        mode = raw_mode
        mode_note = ""

    fracture = bool(_read_field(result, "ankle_fracture", False))
    sprain = bool(_read_field(result, "ankle_sprain", False))
    bone_util = float(_read_field(result, "ankle_bone_utilization", 0.0))
    peak_strain = float(_read_field(result, "peak_ligament_strain", 0.0))
    peak_beta = float(_read_field(result, "peak_ankle_beta_rad", 0.0))
    peak_moment = float(_read_field(result, "peak_inversion_moment_nmm", 0.0))
    beta0 = float(_read_field(result, "ankle_beta0_rad", 0.0))
    height = float(getattr(result, "height_m", 0.0))

    # ---- 档位 ---------------------------------------------------------
    band = ankle_band(sprain, fracture, bone_util, peak_strain)

    # ---- 注释 ---------------------------------------------------------
    bits: list[str] = []
    if mode != "none":
        bits.append(f"模式={mode}")
    if sprain:
        bits.append(f"扭伤（peak_ligament_strain={peak_strain:.3f}）")
    if fracture:
        bits.append(f"骨折（utilization={bone_util:.2f}）")
    if not bits:
        bits.append("踝未越阈")
    base_note = "；".join(bits) + mode_note

    return AnkleVerdict(
        posture=posture or getattr(result, "posture", ""),
        height_m=height,
        beta0_rad=beta0,
        mode=mode,
        peak_beta_rad=peak_beta,
        peak_moment_nmm=peak_moment,
        peak_ligament_strain=peak_strain,
        bone_utilization=bone_util,
        fracture=fracture,
        sprain=sprain,
        band=band,
        note=base_note,
    )