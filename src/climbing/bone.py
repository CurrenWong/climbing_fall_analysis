"""分部位骨骼损伤评估 —— Phase B-2 补充（方案 v2 · P0）。

为什么需要这个模块
------------------
v1 的 ``injury.py`` 只对"脊柱/下肢"给一个笼统的力阈值（6/10/16 kN）。
两篇文献表明这抹掉了两个关键事实：

1. **[Y25] Table 1**：骨的**两端**（松质骨为主）比**骨干**（皮质骨）弱约 3 倍。
   腓骨两端 30/70 MPa vs 腓骨骨干 90/160 MPa；股骨颈 40/80 vs 股骨骨干 140/220。
   骨折恰恰最爱发生在两端（踝、膝、髋）。
2. **骨折是可定位的**：不同部位有各自的阈值与主导失效模式（压缩/拉伸）。

本模块把"一根骨一个阈值"细化成**部位级**评估，并明确每条通路的载荷来源。

载荷来源映射（双质点 B-2 模型）
------------------------------
在 ``pad.py`` 的双质点模型里，只有两条内在信号可用：

* ``f_pad``：地面/垫面反力 —— 直接作用在**第一接触部位**
* ``f_flex``：经屈曲段上传的轴向力 —— 由牛顿第三定律
  ``f_flex = m_up · (g − a_torso)``，即**脊柱/骨盆/股骨/胫腓骨的轴向载荷**

因此：

| 部位 | 载荷 |
|---|---|
| 跟骨、足弓 | ``f_pad`` |
| 胫骨（远端/骨干）、腓骨、股骨颈、骨盆、脊柱椎体 | ``f_flex`` |

⚠️ 颅骨**不在**本模块 —— 它的失效由接触面积与加速度主导，
标准做法是 HIC（见 ``metrics.py``），不在轴向力通路上。

关于阈值来源的重要说明
----------------------
[Y25] Table 1 给的是**材料强度（MPa）**，用于和有限元 von Mises 应力对比。
把它 × 截面面积换算成整体骨失效载荷是**近似**，且对部分部位会显著高估
（足、骨盆、脊柱尤甚）——因为整体骨的失效还涉及足弓、小梁骨、弯曲与应力集中。

因此本模块：
* 保留 [Y25] 的**材料强度**作为 ``MATERIAL_STRENGTH_MPA`` 参考表
* 用**显式记载的承载截面**换算成力阈值，并在 ``note`` 里标注置信度
* **相对排序（两端 vs 骨干）比绝对值可靠** —— 这是可以放心用的部分

独立佐证：跟骨骨折的**力判据**已被验证可用 —— Barnes 等（IRCOBI 2019）
比较全身 PMHS 与足-腿构件试验，确认"足跟接触力"是严重跟骨骨折的有效
预测量（Voo 等的 HIPC）。本模块用峰值力而非应力，正是这个理由。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "BoneSite",
    "BONE_SITES",
    "MATERIAL_STRENGTH_MPA",
    "SiteLoad",
    "SiteVerdict",
    "assess_sites",
    "knockdown_report",
]


# --------------------------------------------------------------------------
# [Y25] Table 1 —— 材料强度参考（MPa）
# --------------------------------------------------------------------------
def _mpa(sigma_t: float, sigma_c: float) -> tuple[float, float]:
    return (float(sigma_t), float(sigma_c))


#: [Y25] 表 1 原始值：材料级拉伸/压缩强度 (MPa)。
#: **这是材料强度，不是整体骨失效载荷。** 见模块 docstring。
MATERIAL_STRENGTH_MPA: dict[str, tuple[float, float]] = {
    "calcaneus": _mpa(100, 150),      # 足部
    "tibia_shaft": _mpa(120, 200),    # 胫骨骨干
    "tibia_ends": _mpa(30, 70),       # 胫骨两端
    "fibula_shaft": _mpa(90, 160),    # 腓骨骨干
    "fibula_ends": _mpa(30, 70),      # 腓骨两端
    "femur_shaft": _mpa(140, 220),    # 股骨骨干
    "femoral_neck": _mpa(40, 80),     # 股骨颈
    "pelvis": _mpa(100, 180),         # 骨盆
    "spine": _mpa(70, 150),           # 脊柱椎体
    "skull": _mpa(80, 160),           # 颅骨（本模块不评估，见 docstring）
}


# --------------------------------------------------------------------------
# 部位定义
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class BoneSite:
    """一个可评估的骨骼部位。

    Attributes
    ----------
    key
        ``MATERIAL_STRENGTH_MPA`` 的键（材料强度参考）。
    name_cn
        中文名。
    area_mm2
        **承载截面面积**（mm²），用于把材料强度换算成力阈值。
        这是本模块最不确定的参数（个体差异可达 ±30%）。
    mode
        ``"compression"`` / ``"tension"`` —— 该部位在实际坠落中的主导失效模式。
    load_source
        ``"pad"``（用 f_pad）或 ``"flex"``（用 f_flex）。
    source
        面积取值的依据。
    confidence
        ``"low"`` / ``"medium"`` —— 换算置信度。
    """

    key: str
    name_cn: str
    area_mm2: float
    mode: str
    load_source: str
    source: str
    load_fraction: float = 1.0
    confidence: str = "medium"
    note: str = ""

    @property
    def sigma_mpa(self) -> float:
        """该部位按主导失效模式选取的材料强度 (MPa)。"""
        st, sc = MATERIAL_STRENGTH_MPA[self.key]
        return sc if self.mode == "compression" else st

    def load_from(self, peak_pad_n: float, peak_flex_n: float) -> float:
        """该部位（单根骨）的实际峰值载荷 (N)。

        ``f_pad`` / ``f_flex`` 是**整条通路**的总力，而阈值是**单根骨**的。
        倍数分配由 ``load_fraction`` 给出 —— 见各部位的 source 说明。
        """
        gross = peak_pad_n if self.load_source == "pad" else peak_flex_n
        return gross * self.load_fraction

    @property
    def force_threshold_n(self) -> float:
        """材料强度 × 承载截面 → 名义失效载荷 (N)。

        ``1 MPa = 1 N/mm²``，故 ``MPa × mm² = N``。
        """
        return self.sigma_mpa * self.area_mm2

    @property
    def threshold_kn(self) -> float:
        return self.force_threshold_n / 1e3


#: 轴向载荷通路上的部位（足先落地场景）。
#:
#: 截面面积取自人体测量学/生物力学文献的**量级**（皮质骨为主的骨干取较小值，
#: 含松质骨的骨端取较大值）。**置信度整体偏低** —— 这是本模块最大的不确定源，
#: 也是 ``scripts/validate_vs_fem.py`` 要量化的对象。
BONE_SITES: dict[str, BoneSite] = {
    "calcaneus": BoneSite(
        key="calcaneus", name_cn="跟骨（足部）",
        area_mm2=900.0, mode="compression", load_source="pad",
        load_fraction=0.50,       # 双足着地，单侧跟骨承担约一半地面反力
        source="跟骨体横截面量级；单侧承重按双足均分取 0.5", confidence="low",
        note="荷载集中在跟骨、非均布；材料强度×截面积会显著高估。"
             "独立佐证：足跟接触力已是经验证的跟骨骨折预测量（Voo HIPC）。",
    ),
    "tibia_distal": BoneSite(
        key="tibia_ends", name_cn="胫骨远端",
        area_mm2=450.0, mode="compression", load_source="flex",
        load_fraction=0.45,       # 单侧胫骨 ≈ 单腿载荷的 90%，双足均分后 0.45
        source="胫骨远端横截面；单侧胫骨取整条通路约 45%", confidence="medium",
        note="胫骨远端骨折是足先落地的经典损伤。",
    ),
    "tibia_mid": BoneSite(
        key="tibia_shaft", name_cn="胫骨骨干",
        area_mm2=350.0, mode="compression", load_source="flex",
        load_fraction=0.45,
        source="胫骨中段皮质骨环截面；同上 45%", confidence="medium",
    ),
    "fibula_ends": BoneSite(
        key="fibula_ends", name_cn="腓骨两端",
        area_mm2=80.0, mode="compression", load_source="flex",
        load_fraction=0.05,       # 关键：腓骨在完整小腿中只承担 ~10% 轴向载荷
        source="腓骨远端横截面（细）；单侧腓骨按并联分流取 0.05", confidence="medium",
        note="[Y25] 中骨折风险随高度上升最快（OR=1.682）—— 相对脆弱的非主承重骨。"
             "⚠️ 载荷分配是重点：整条小腿的轴向载荷主要由胫骨承担，"
             "腓骨只分担约 10%（单侧 5%）。若两骨吃同样的力，"
             "腓骨会被严重高估（实测 1 m 落垫就误判腓骨骨折 2.25×）。",
    ),
    "femoral_neck": BoneSite(
        key="femoral_neck", name_cn="股骨颈",
        area_mm2=500.0, mode="compression", load_source="flex",
        load_fraction=0.50,       # 单侧股骨承担约一半
        source="股骨颈横截面；单侧按双足均分取 0.5", confidence="low",
        note="强度低（40/80 MPa）+ 弯曲/剪切复合受力，实际失效载荷远低于纯压缩估算。",
    ),
    "pelvis": BoneSite(
        key="pelvis", name_cn="骨盆",
        area_mm2=1200.0, mode="compression", load_source="flex",
        load_fraction=1.0,        # 骨盆是单一结构，承担两腿汇合后的全部载荷
        source="骨盆环承载截面量级；单一致结构取 1.0", confidence="low",
        note="[Y25] 中风险随高度上升最慢（OR=1.236）—— 承重能力最强。",
    ),
    "lumbar_spine": BoneSite(
        key="spine", name_cn="腰椎椎体",
        area_mm2=900.0, mode="compression", load_source="flex",
        load_fraction=1.0,        # 脊柱是单一轴向柱
        source="L1 椎体横截面；单一柱取 1.0", confidence="low",
        note="[Y25] 属'阶跃响应'部位：低于临界高度绝对安全，高于则骤增。",
    ),
}


# --------------------------------------------------------------------------
# 评估
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class SiteLoad:
    """一个部位的实际载荷（N）与阈值（N）。"""

    key: str
    name_cn: str
    peak_force_n: float
    threshold_n: float
    utilization: float          # peak / threshold
    confidence: str
    note: str = ""

    @property
    def peak_force_kn(self) -> float:
        return self.peak_force_n / 1e3

    @property
    def threshold_kn(self) -> float:
        return self.threshold_n / 1e3

    @property
    def fractured(self) -> bool:
        return self.utilization >= 1.0

    @property
    def band(self) -> str:
        u = self.utilization
        if u >= 1.0:
            return "超过阈值"
        if u >= 0.7:
            return "接近阈值"
        if u >= 0.4:
            return "中等"
        return "低"


@dataclass
class SiteVerdict:
    """一次仿真的分部位骨骼判读。"""

    posture: str
    height_m: float
    on_pad: bool
    peak_pad_n: float = 0.0
    peak_flex_n: float = 0.0
    sites: list[SiteLoad] = field(default_factory=list)

    @property
    def fractured_sites(self) -> list[SiteLoad]:
        return [s for s in self.sites if s.fractured]

    @property
    def worst(self) -> SiteLoad | None:
        return max(self.sites, key=lambda s: s.utilization) if self.sites else None

    def summary(self) -> str:
        head = (f"{self.posture} @ {self.height_m:.1f} m "
                f"({'软垫' if self.on_pad else '刚性地面'})")
        if not self.fractured_sites:
            w = self.worst
            return (f"{head}: 无部位超过阈值"
                    + (f"（最高 {w.name_cn} {w.utilization:.2f}）" if w else ""))
        names = "、".join(f"{s.name_cn}({s.utilization:.2f})"
                          for s in sorted(self.fractured_sites,
                                          key=lambda s: -s.utilization))
        return f"{head}: 超过阈值 —— {names}"


def assess_sites(result, posture: str = "") -> SiteVerdict:
    """对一次 ``BoulderFallResult`` 做分部位骨骼评估。

    Parameters
    ----------
    result
        ``pad.simulate_boulder_fall`` 的返回值。需要 ``impact_window`` 切片内的
        ``pad_force_n`` 与 ``flex_force_n``。若结果里没有 ``flex_force_n``，
        则用 ``m_up·(g − a_torso)`` 由牛顿第三定律反算。
    posture
        姿势名（仅用于报告）。

    Returns
    -------
    SiteVerdict
    """
    meta = getattr(result, "meta", {}) or {}
    win = meta.get("impact_window_s")
    if win is None:
        raise ValueError(
            "result.meta['impact_window_s'] 缺失：分部位评估必须在冲击窗口内取峰值，"
            "否则会被回弹重绷污染（见 阶段总结.md §2.4.2 约定①）。"
        )
    t = np.asarray(result.t_s)
    t0, t1 = float(win[0]), float(win[1])
    sl = slice(int(np.searchsorted(t, t0)), int(np.searchsorted(t, t1)) + 1)

    f_pad = np.asarray(result.pad_force_n)[sl]
    peak_pad = float(np.max(np.abs(f_pad)))

    # f_flex：结果里已由 pad.py 在冲击窗口内算好（meta['f_flex_peak_n']）。
    # 若缺失才由躯干方程反算：f_flex = m_up·(g − a_torso)（牛顿第三定律）。
    peak_flex = meta.get("f_flex_peak_n")
    if peak_flex is None:
        a_torso = np.asarray(result.accel_torso_g)[sl] * 9.81
        m_up = meta.get("m_up")
        if m_up is None:
            m = float(result.mass_kg)
            frac = meta.get("m_low_frac")
            if frac is None:
                raise ValueError("需要 meta['m_up'] 或 meta['m_low_frac'] 才能反算 f_flex")
            m_up = m * (1.0 - frac)
        peak_flex = float(np.max(np.abs(m_up * (9.81 - a_torso))))
    peak_flex = float(peak_flex)

    verdict = SiteVerdict(posture=posture or result.posture,
                          height_m=float(result.height_m),
                          on_pad=bool(meta.get("on_pad", True)),
                          peak_pad_n=peak_pad, peak_flex_n=peak_flex)
    for key, site in BONE_SITES.items():
        peak = site.load_from(peak_pad, peak_flex)
        thr = site.force_threshold_n
        verdict.sites.append(SiteLoad(
            key=key, name_cn=site.name_cn,
            peak_force_n=peak, threshold_n=thr,
            utilization=peak / thr if thr > 0 else float("inf"),
            confidence=site.confidence, note=site.note,
        ))
    return verdict


# --------------------------------------------------------------------------
# 换算可靠性诊断
# --------------------------------------------------------------------------
def knockdown_report(f_pad_peak_n: float) -> str:
    """给出"材料强度换算"相对"已知整体骨失效载荷"的偏差诊断。

    用跟骨做例子（它是唯一有经验证的力判据的部位）：
    把本模块算出的名义阈值与文献中跟骨骨折的量级（约 4–8 kN）对比，
    输出所需的结构性折减系数（knockdown factor）。

    这个系数**不是**要拿去硬乘，而是用来量化"材料强度路线"的偏差幅度。
    """
    site = BONE_SITES["calcaneus"]
    nominal = site.force_threshold_n
    lo, hi = 4.0e3, 8.0e3
    k_hi, k_lo = nominal / hi, nominal / lo
    lines = [
        "—— 换算可靠性诊断（跟骨）——",
        f"名义阈值 = σ_c × A = {site.sigma_mpa:.0f} MPa × {site.area_mm2:.0f} mm² "
        f"= {nominal/1e3:.1f} kN",
        f"文献跟骨骨折量级 = {lo/1e3:.0f}–{hi/1e3:.0f} kN",
        f"所需结构性折减 = {k_hi:.1f}× … {k_lo:.1f}×",
        "",
        f"结论：材料强度 × 承载截面 对跟骨**高估约 {k_hi:.0f}–{k_lo:.0f} 倍**。",
        "原因是跟骨失效不是均匀柱体压缩 —— 涉及足弓、小梁骨、以及",
        "荷载在距下关节的非均布。故：",
        "  ✅ 部位间**相对排序**（两端 vs 骨干、腓骨 vs 胫骨）可用",
        "  ❌ 绝对力阈值需整体骨标定后才可用",
        f"（仿真当前峰值 f_pad = {f_pad_peak_n/1e3:.1f} kN）",
    ]
    return "\n".join(lines)
