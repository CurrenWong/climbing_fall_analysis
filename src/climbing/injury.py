"""Phase B-2 —— 抱石落地**分部位损伤评估**。

为什么需要这个模块
------------------
Phase B-1 只报"峰值力 / 峰值加速度 / HIC"，漏掉了抱石摔伤里最要命的一类：

    **力不大，但压强大。**

头朝下落地时峰值力只有 20.7 kN（六种姿势里最低），看上去"最安全"；
但接触面积只有 0.078 m²，峰值接触压力 **264 kPa**，而那 6.4 kg 的头
独自扛下整个冲击 —— 加速���高达 85 g。反过来，平拍落地峰值力 61 kN
（全场最高），但接触面积 0.60 m²，压力只有 **102 kPa**，是全场最低。

单看力，排序会完全反过来。所以本模块引入三个正交的判据：

1. **峰值接触压力** ``p = F(t) / A(t)`` —— 局部软组织/骨损伤的直接驱动量
2. **HIC** —— 脑损伤（沿用 metrics 的标准分档）
3. **脊柱轴向力** —— 椎体压缩/爆裂骨折

数据性质警告
------------
本模块的阈值来自公开生物力学文献的**量级**估计，用于区分
"这个量级值得关注"和"这个量级很危险"，**不是临床判据**。
个体差异、体位、既往损伤、疲劳都会大幅移动实际阈值。
任何单个数字都应理解为"值得注意的量级"，而非"超过就必然骨折"。

用法
----
    from climbing.injury import assess_boulder_fall
    rep = assess_boulder_fall(result, posture_name)
    print(rep.verdict)          # 总体判读
    for part in rep.parts:      # 分部位
        print(part.name, part.metric, part.band)

⚠️ 绝对判读 vs 相对排序（重要）
---------------------------------
本模块的**绝对分档目前偏保守**，3 m 落差下几乎所有姿势都会顶到"极高"，
跨高度/质量的区分度不足。原因不是阈值随手定的，而是双质点模型的结构限制：

* 躯干除了屈曲弹簧外**没有第二条通路**，现实中"腿屈完后躯干直接落在垫上"
  的过程在这里无法表达，躯干屈曲行程被系统性低估（受控屈膝实测只走到
  8.7 cm，而姿势允许 40 cm），于是 ``a_torso ≈ v²/(2x)`` 被高估约 1.5-2 倍。
* 由此脊柱/下肢惯性载荷 ``m·a`` 整体偏高，普遍跨过"极高"档（16 kN）；
  低落差端也一样失真（0.2 m 掉落就判"中"）。

**结论**：
- ✅ **相对排序**（哪种姿势更危险、哪个部位先受伤）—— 可用且稳定
- ✅ **机制解释**（接触面积 / 屈曲行程 / 质量集中）—— 可用
- ❌ **绝对判读分档**（"3 m = 极高"、"0.2 m 就开始受伤"）—— 暂不可用

要修好绝对判读，需要 Phase B-3 的多刚体接触（躯干直接接触垫面 +
Hill 型肌肉力-速关系），而不是继续调阈值。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .metrics import HIC_BANDS, RISK_G, hic_from_result

__all__ = [
    "Band",
    "PartAssessment",
    "InjuryReport",
    "assess_boulder_fall",
    "PRESSURE_BANDS",
    "BONE_BANDS",
    "VERDICTS",
]


# --------------------------------------------------------------------------
# 分档
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Band:
    """一条分档规则。"""

    lo: float
    hi: float
    label: str

    def contains(self, v: float) -> bool:
        return self.lo <= v < self.hi


def _band_of(v: float, table) -> str:
    for b in table:
        if b.contains(v):
            return b.label
    return table[-1].label


def _as_bands(table) -> list:
    """把 ``metrics`` 里的旧式 ``(lo, hi, label)`` 元组转成本模块的 ``Band``。

    ``metrics.HIC_BANDS`` / ``metrics.RISK_G`` 是公开 API（``smoke.py`` 和
    早期脚本都在用），改成 ``Band`` 会破坏它们，所以在这里做适配而不是
    回头改 metrics。
    """
    out = []
    for item in table:
        if isinstance(item, Band):
            out.append(item)
        else:
            lo, hi, label = item
            out.append(Band(lo, hi, label))
    return out


HIC_BANDS_B = _as_bands(HIC_BANDS)
RISK_G_B = _as_bands(RISK_G)


def _rank(table, label: str) -> int:
    for i, b in enumerate(table):
        if b.label == label:
            return i
    return 0


#: 峰值接触压力分档 (kPa)。瞬态冲击下软组织/骨的局部损伤量级。
#: 注意：压疮文献里的 4-32 kPa 是**持续**受压，与毫秒级冲击不是一回事，
#: 这里用的是冲击载荷下"软组织挫伤/撕裂开始明显"的量级。
PRESSURE_BANDS = [
    Band(0.0, 150.0, "低"),
    Band(150.0, 400.0, "中"),
    Band(400.0, 800.0, "高"),
    Band(800.0, float("inf"), "极高"),
]

#: 骨骼轴向峰值载荷分档 (kN)。用于脊柱（椎体爆裂骨折）与下肢（胫骨骨折）。
#: 依据：活体腰椎/胫骨在轴向压缩下的失效量级约 8-16 kN
#: （离体试验更高，活体更接近 8-12 kN）。
BONE_BANDS = [
    Band(0.0, 6.0, "低"),
    Band(6.0, 10.0, "中"),
    Band(10.0, 16.0, "高"),
    Band(16.0, float("inf"), "极高"),
]

#: 总体判读（取所有部位里最严重的一档）
VERDICTS = ["低", "中", "高", "极高"]


# --------------------------------------------------------------------------
# 分部位评估
# --------------------------------------------------------------------------
@dataclass
class PartAssessment:
    """一个身体部位的损伤评估结果。"""

    name: str
    metric: float
    metric_name: str
    unit: str
    band: str
    mechanism: str = ""      # 为什么是这个部位 / 这个机制
    detail: dict = field(default_factory=dict)

    def __str__(self) -> str:  # pragma: no cover - 仅供打印
        return (f"{self.name:<10s} {self.metric_name}={self.metric:8.1f} {self.unit:<4s} "
                f"[{self.band}]")


@dataclass
class InjuryReport:
    posture: str
    surface: str
    height_m: float
    mass_kg: float
    parts: list[PartAssessment]
    verdict: str
    peak_pressure_kpa: float
    peak_force_kn: float
    peak_accel_g: float
    hic: float

    def worst_parts(self, n: int = 3) -> list[PartAssessment]:
        return sorted(self.parts, key=lambda p: -_rank(
            {"低": 0, "中": 1, "高": 2, "极高": 3}.get(p.band, 0)))[:n]


# --------------------------------------------------------------------------
# 主评估
# --------------------------------------------------------------------------
def assess_boulder_fall(result, posture_name: str = "") -> InjuryReport:
    """对一次 ``BoulderFallResult`` 做分部位损伤评估。

    Parameters
    ----------
    result
        ``pad.simulate_boulder_fall`` 的返回值。
    posture_name
        姿势键名（``POSTURES`` 的键），用于查接触部位表。若为空则用
        ``result.meta['posture_obj']``。
    """
    post = result.meta.get("posture_obj")
    name = posture_name or (post.name if post is not None else result.posture)
    m = float(result.mass_kg)

    # ---- 峰值接触压力：局部损伤的直接驱动量 -----------------------------
    # 压力 = 力 / **该时刻的接触面积**。接触面积是时变的（脚 -> 臀 -> 背
    # 逐渐展开），所以峰值压力出现在**面积还小、力已经上来了**的早期，
    # 常常不是峰值力那一瞬间。
    area = np.maximum(np.asarray(result.contact_area_m2, dtype=float), 1e-9)
    force = np.maximum(np.asarray(result.pad_force_n, dtype=float), 0.0)
    pressure = force / area
    p_peak_kpa = float(np.max(pressure)) / 1e3
    i_ppeak = int(np.argmax(pressure))

    # ---- 骨骼轴向载荷 ---------------------------------------------------
    # 关键：这里用**惯性载荷** ``F = m·a``，而不是屈曲弹簧力 f_flex。
    # f_flex 里混着 ``c_flex·rate`` 这个纯阻尼项（3e3 N·s/m，屈曲速率
    # 15 m/s 时就有 45 kN），把它当成脊柱轴向力会让"腿绷直"的姿势脊柱
    # 载荷虚高到 44-48 kN，判读全部顶到"极高"，分档完全失去区分度。
    # 真正压在脊柱/胫骨上的是"该段身体质量 × 它自己的减速加速度"。
    m_up = float(result.meta.get("m_up", result.mass_kg * 0.7))
    m_low = float(result.meta.get("m_low", result.mass_kg * 0.3))
    a_up = float(np.max(result.accel_torso_g))
    a_low = float(np.max(result.accel_leg_g))
    spine_kN = m_up * 9.80665 * a_up / 1e3
    leg_kN = m_low * 9.80665 * a_low / 1e3

    # ---- HIC -------------------------------------------------------------
    h = hic_from_result(result)

    # ---- 各部位 ----------------------------------------------------------
    parts: list[PartAssessment] = []

    first = post.first_contact if post is not None else "?"
    is_head = bool(post.first_contact_is_head) if post is not None else False

    # (1) 第一接触部位 —— 压力直接作用在这里
    parts.append(PartAssessment(
        name=f"接触部位({first})",
        metric=p_peak_kpa,
        metric_name="峰值接触压力",
        unit="kPa",
        band=_band_of(p_peak_kpa, PRESSURE_BANDS),
        mechanism=(
            f"第一接触部位是{first}，峰值时刻接触面积仅 "
            f"{1e4 * float(np.max(area[np.argmax(pressure)])):.0f} cm²，"
            f"{result.peak_force_kn:.1f} kN 的力压在上面 = {p_peak_kpa:.0f} kPa。"
            f"面积小是这里最容易忽略的伤害来源。"
        ),
        detail={"peak_force_kn": result.peak_force_kn,
                "area_cm2_at_ppeak": 1e4 * float(area[np.argmax(pressure)])},
    ))

    # (2) 头/颈
    if is_head:
        parts.append(PartAssessment(
            name="头/颈",
            metric=h,
            metric_name="HIC",
            unit="",
            band=_band_of(h, HIC_BANDS_B),
            mechanism=("头是唯一的第一接触点，m_low 就是头本身（约 8% 体重），"
                       "全部冲击由它独自承担；HIC 算在头部的加速度上。"),
            detail={"peak_accel_g": float(result.peak_accel_leg_g)},
        ))
    else:
        # 头不是第一接触点，但会跟着躯干一起减速。头部在屈曲链末端，
        # 实际加速度略高于躯干 CoM，用一个温和的放大系数。
        a_head = float(result.peak_accel_torso_g) * 1.15
        parts.append(PartAssessment(
            name="头/颈",
            metric=h,
            metric_name="HIC",
            unit="",
            band=_band_of(h, HIC_BANDS_B),
            mechanism=("头不是第一接触点，HIC 取躯干代表点加速度（略低估，"
                       "真实头加速度在屈曲链末端更高）。"),
            detail={"peak_accel_g": a_head},
        ))

    # (3) 脊柱 —— 躯干惯性载荷
    parts.append(PartAssessment(
        name="脊柱",
        metric=spine_kN,
        metric_name="轴向峰值载荷",
        unit="kN",
        band=_band_of(spine_kN, BONE_BANDS),
        mechanism=(f"脊柱承受躯干（{m_up:.0f} kg）的减速惯性载荷，"
                   f"峰值 = m·a = {spine_kN:.1f} kN（躯干 {a_up:.1f} g）。"
                   + ("屈曲行程不足，躯干的速度必须靠刚性止挡吃��，"
                      "载荷直接沿脊柱传下去。"
                      if post is not None and post.x_flex_max < 0.10 else
                      "有可用的屈曲行程，载荷被肌肉/关节分散吸收。")),
        detail={"peak_accel_g": a_up, "m_up_kg": m_up},
    ))

    # (4) 下肢 —— 第一接触段（足/踝/胫）的轴向载荷
    parts.append(PartAssessment(
        name="下肢(足/踝)",
        metric=leg_kN,
        metric_name="轴向峰值载荷",
        unit="kN",
        band=_band_of(leg_kN, BONE_BANDS),
        mechanism=(f"第一接触段（{first}，{m_low:.0f} kg）的减速惯性载荷 "
                   f"= {leg_kN:.1f} kN（{a_low:.1f} g）。"
                   + ("接触面积小 + 腿基本不屈曲时，这个载荷全部由"
                      "跖骨/踝的极小截面承担，是扭伤与骨折的典型机制。"
                      if post is not None and post.m_low_frac < 0.15 else
                      "接触面积较大且有屈曲行程，载荷被分散。")),
        detail={"peak_accel_g": a_low, "m_low_kg": m_low,
                "x_flex_max_m": post.x_flex_max if post is not None else None},
    ))

    # ---- 总体判读 --------------------------------------------------------
    band_rank = {"低": 0, "中": 1, "高": 2, "极高": 3}
    verdict = max((p.band for p in parts), key=lambda b: band_rank.get(b, 0))

    return InjuryReport(
        posture=name,
        surface=result.pad_label,
        height_m=result.height_m,
        mass_kg=m,
        parts=parts,
        verdict=verdict,
        peak_pressure_kpa=p_peak_kpa,
        peak_force_kn=result.peak_force_kn,
        peak_accel_g=float(result.peak_primary_g),
        hic=h,
    )
