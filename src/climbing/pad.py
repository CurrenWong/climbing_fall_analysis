"""Phase B-1 —— 抱石坠落 1D 模型（软垫 + 渐进接触 + 屈膝缓冲）。

为什么不是"质点 + 双线性弹簧"
------------------------------
方案初稿把攀岩者当作质点压在软垫上。这会漏掉抱石伤害的**核心机制**：

1. **接触面积是渐进的**。脚先着地时只有 ~0.06 m² 接触，软垫把力集中在
   踝/足；人继续塌下去，接触面积才逐渐展开到臀、背。力-位移曲线因此
   是*时变*的。
2. **肌肉/关节屈曲是第二个吸能器**。屈膝落地的缓冲主要来自股四头肌和
   髋的离心收缩，不是软垫。软垫 + 关节屈曲的**串联行程**才是决定
   峰值力的东西。
3. **头部/躯干的减速晚于双脚**。同样的垫反力，作用在小质量上给出的
   加速度远大于作用在全身质量上。头朝下落地之所以致命，正因为头是
   第一接触点且质量占比小。

所以本模块用**双质点串联模型**：

    腿/足 (m_low)  ──[ 屈曲弹簧 k_flex, 行程 x_flex_max ]──  躯干+头 (m_up)
        │
        └── 软垫（渐进接触面积 A(κ)，非线��� σ(ξ)）

输出：峰值垫反力、双/躯干峰值加速度（→ HIC）、最大压缩、是否压穿、
以及**未被弹性吸收的残余能量**（= 必然转化为组织损伤的那部分）。
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import solve_ivp

from . import G

__all__ = [
    "CrashPad",
    "Posture",
    "BoulderFallResult",
    "POSTURES",
    "simulate_boulder_fall",
    "hard_surface",
]


# --------------------------------------------------------------------------
# 软垫泡沫
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class CrashPad:
    """Crash pad 泡沫的压缩应力-应变本构 + 厚度。

    应力-应变用**平台 + 压实**两段（开放孔 PU 泡沫的典型形态）::

        σ(ξ) = s_pl·(1 - e^(-ξ/ξ_pl))  +  s_d·(ξ/ξ_d)^n  +  E_bot·max(ξ-ξ_d, 0)

    - 平台段（0 → ~0.15 应变）：开孔结构逐层压溃，应力上升很慢。
    - 压实段（→ ``xi_dense_max``）：孔壁贴死，曲线急剧上扬。
    - 超过 ``xi_dense_max``：泡沫压成实心层，只剩很陡的线性支撑。

    **标定依据（2026-10 重标）**
    ------------------------
    30 kg/m³ 级开放孔 PU 攀岩垫的压缩应力-应变量级：

    ==========  ================
    压缩应变 ξ   应力
    ==========  ================
    10 %         10-20 kPa
    30 %         30-45 kPa
    50 %         60-120 kPa
    70 %         150-300 kPa
    85 %         300-500 kPa
    ==========  ================

    旧参数（``xi_dense_max=0.60``、``n_dense=5``、``s_dense=350 kPa``、
    ``bottoming_stiffness=50``）算出的曲线是 50%→174 kPa、60%→384 kPa、
    75%→2853 kPa —— **在 60% 就进入近乎刚性的壁**。后果是 3 m 落差对
    几乎所有姿势都会压穿，峰值接触压力被抬到 470-520 kPa，
    于是"落在 20 cm 软垫"和"直接摔在混凝土"的压力差只剩 4%
    （472 vs 494 kPa），软垫等于没起作用。这是**标定缺陷，不是物理结论**。

    Attributes
    ----------
    thickness_m
        总厚度。运动馆常见 10-30 cm，多层堆叠取总和。
    s_plateau_pa
        平台应力 (Pa)。30 kg/m³ 的开孔 PU 约 25-45 kPa。
    s_dense_pa
        压实段起点应力 (Pa)，即 ``xi_dense_max`` 处的应力。
    xi_plateau
        平台段特征应变。
    n_dense
        压实段指数，越大越"突然"。
    xi_dense_max
        压实极限应变（压穿点）。PU 泡沫可压缩到 70-85% 才成为实心层。
    bottoming_modulus_pa
        压穿后的线性支撑模量 (Pa)。泡沫压实后仍能继续压缩一小段，
        不应变成刚性墙。

    ⚠️ 这些是**文献量级**，不是本项目的实测标定。相对比较可用，
    绝对力值应视为量级估计。
    """

    thickness_m: float = 0.20
    s_plateau_pa: float = 35_000.0
    s_dense_pa: float = 365_000.0
    xi_plateau: float = 0.15
    n_dense: float = 3.56
    xi_dense_max: float = 0.85
    bottoming_modulus_pa: float = 5.0e6  # 压穿后的线性支撑模量
    linear_modulus_pa: float = 0.0     # 额外的线性 σ = E·ξ 项 (Pa)，泡沫用 0

    # -- 应力-应变 ---------------------------------------------------------
    def stress_pa(self, xi: np.ndarray | float) -> np.ndarray:
        xi = np.maximum(np.asarray(xi, dtype=float), 0.0)
        plateau = self.s_plateau_pa * (1.0 - np.exp(-xi / self.xi_plateau))
        frac = np.clip(xi / self.xi_dense_max, 0.0, None)
        dense = self.s_dense_pa * np.power(frac, self.n_dense)
        out = plateau + dense
        # 线性项：给"纯弹性"的接触面用（见 hard_surface）
        out = out + self.linear_modulus_pa * xi
        # 压穿之后：实心泡沫层的陡峭但有限的支撑
        out = out + self.bottoming_modulus_pa * np.maximum(xi - self.xi_dense_max, 0.0)
        return out

    def force_n(self, compression_m: np.ndarray | float, area_m2: float) -> np.ndarray:
        """给定压缩量和**当前接触面积**返回反力 (N)。"""
        xi = np.maximum(np.asarray(compression_m, dtype=float), 0.0) / self.thickness_m
        return self.stress_pa(xi) * area_m2

    def capacity_j(self, area_m2: float, xi_limit: float | None = None) -> float:
        """整块垫在给定面积下、压穿前能吸收的能量 (J)。"""
        xi_max = self.xi_dense_max if xi_limit is None else min(xi_limit, self.xi_dense_max)
        xi = np.linspace(0.0, xi_max, 2000)
        return float(np.trapezoid(self.stress_pa(xi), xi) * self.thickness_m * area_m2)

    @property
    def max_compression_m(self) -> float:
        return self.xi_dense_max * self.thickness_m


def hard_surface(k_n_per_m: float = 1.0e7) -> CrashPad:
    """把硬地面（混凝土/木地板）折算成同接口的"垫"。

    地板本身近乎刚性，真正变形的是人体组织，所以这里用的是
    **人体在高速冲击下的等效轴向刚度** ``k_n_per_m``。

    取值说明（重要，见 README「已知局限」）
    ------------------------------------
    公开文献里"人体等效刚度"有两个量级不同的数：

    * **准静态 / 慢加载**：1-8 MN/m。坐位、立位缓慢下压的实验值。
    * **动态冲击（7-8 m/s）**：10-30 MN/m。80 kg 的人要在 ~5 cm 内
      从 7.7 m/s 停住，等效刚度就是 ``m·a/x = 80×593/0.05 ≈ 1e7 N/m``。

    早期版本默认取 4e6（准静态量级），结果"离垫摔到混凝土"算出来的
    峰值力**低于落在软垫上**，安全排序整个颠倒 —— 因为一个 3 m 落差的
    动态冲击被当成了慢加载。默认改为 1.0e7 (N/m)，落在动态区间下沿。

    ⚠️ 这是全模型**标定最弱**的一个参数。``tests/`` 里只有"硬地面必须
    比软垫更硬/峰值更高"这种定性断言，没有对标具体实验数据。
    有实测跌落数据后应优先重标这个值。

    线性项走 ``linear_modulus_pa``：``σ = E·ξ = k·x``（``thickness=1``
    时 ξ 就是压缩量 m）。早期版本试图用 ``s_plateau_pa`` 表达同一个刚度，
    但那会走 ``s_pl·(1-e^(-ξ/ξ_pl))`` 的**饱和指数**曲线 ——
    在实际压缩量（几厘米）范围内它远低于 k·x，进一步把硬地面算软。
    """
    return CrashPad(
        thickness_m=1.0,
        s_plateau_pa=0.0,
        s_dense_pa=0.0,
        linear_modulus_pa=k_n_per_m,   # σ = k·x，thickness=1 时 ξ=x
        xi_plateau=1.0,
        n_dense=1.0,
        xi_dense_max=1.0,
        bottoming_modulus_pa=0.0,
    )


# --------------------------------------------------------------------------
# 落地姿势
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Posture:
    """一种落地姿势。

    Attributes
    ----------
    m_low_frac
        腿/足占总质量的比例。头朝下落地时头部是第一接触点，
        此时"第一接触质量"很小，峰值加速度会被放大。
    k_flex
        屈曲弹簧刚度 (N/m)：踝+膝+髋的串联等效刚度。
    x_flex_max
        屈曲可用行程 (m)。
    contact_schedule
        ``(塌陷进度 κ, 接触面积 m²)`` 的单调递增折线。κ 是
        "相对释放点的下落量 / 总塌陷行程"。
    collapse_ref_m
        总塌陷行程（κ=1 时的下落量），仅用于把 κ 定义成无量纲量。
    first_contact
        首次接触部位。
    """

    name: str
    m_low_frac: float
    k_flex: float
    x_flex_max: float
    contact_schedule: tuple[tuple[float, float], ...]
    collapse_ref_m: float
    first_contact: str
    note: str = ""
    first_contact_is_head: bool = False

    @property
    def m_low_frac_check(self) -> float:
        return min(1.0, max(0.0, self.m_low_frac))

    def split_mass(self, mass_kg: float) -> tuple[float, float]:
        """返回 ``(m_low, m_up)``，两者都严格为正。

        ``m_low_frac`` 落在 ``[0, 1]`` 闭区间端点会让其中一个质量为 0，
        躯干/腿的方程随即退化成除以 0（实现里被 ``max(m, 1e-6)`` 掩盖成
        一个巨大的、但不炸的假加速度）。这比直接抛异常危险得多 ——
        早期版本就因为 ``m_low_frac=1.00`` 让躯干一路自由落体 20 m，
        而峰值力看起来"正常"，能量账悄悄全错。因此这里硬性拦截。
        """
        f = self.m_low_frac
        if not (0.0 < f < 1.0):
            raise ValueError(
                f"姿势 {self.name!r} 的 m_low_frac={f} 越界："
                f"必须严格落在 (0, 1) 开区间，否则 m_low 或 m_up 会为 0，"
                f"屈曲弹簧无法建立、躯干/腿的加速度方程退化。"
            )
        m_low = mass_kg * f
        m_up = mass_kg * (1.0 - f)
        if m_low <= 0.0 or m_up <= 0.0:
            raise ValueError(
                f"姿势 {self.name!r} 质量分配退化: m_low={m_low}, m_up={m_up} (m={mass_kg})"
            )
        return m_low, m_up

    def area_at(self, kappa: float) -> float:
        """按塌陷进度插值接触面积 (m²)。"""
        pts = self.contact_schedule
        ks = [p[0] for p in pts]
        areas = [p[1] for p in pts]
        kappa = float(np.clip(kappa, ks[0], ks[-1]))
        return float(np.interp(kappa, ks, areas))


#: 典型落地姿势库。面积取人体测量学量级（双足 ~0.06 m²，臀 ~0.25，全背 ~0.55）。
#:
#: ``m_low_frac`` 用人体测量学约束（80 kg 成年男性量级）：
#: 双腿合计约 18-19%，头颈约 8%，骨盆+大腿约 50%，躯干约 43%。
#: **必须严格落在 (0, 1) 开区间** —— 取 1.0 会让 ``m_up = 0``，
#: 躯干方程退化成 ``g - f_flex/1e-6``，屈曲弹簧永远不建立，
#: 躯干会一路自由落体（实测 2 s 内掉 20 m），能量账全错。
POSTURES: dict[str, Posture] = {
    "flat-flop": Posture(
        name="平拍（背/臀同时着垫）",
        m_low_frac=0.85,           # 背臀整体接触、有效减速质量≈全身
        k_flex=4.0e5,             # 几乎不屈曲
        x_flex_max=0.02,
        contact_schedule=((0.0, 0.30), (0.15, 0.50), (0.6, 0.62), (1.0, 0.70)),
        collapse_ref_m=0.20,
        first_contact="背/臀",
        note="最差：接触面积立刻最大、行程最短，峰值力最高。",
    ),
    "butt-impact": Posture(
        name="坐落（臀先着垫）",
        m_low_frac=0.50,           # 骨盆+大腿先着垫，躯干+头在其上
        k_flex=1.2e5,
        x_flex_max=0.08,
        contact_schedule=((0.0, 0.14), (0.2, 0.24), (0.7, 0.36), (1.0, 0.42)),
        collapse_ref_m=0.28,
        first_contact="臀/骶",
        note="脊柱压缩风险；接触面积比平拍小，峰值力略低。",
    ),
    "feet-first-stiff": Posture(
        name="脚先落但腿绷直",
        m_low_frac=0.18,           # 双足+小腿 = 腿部总质量
        k_flex=1.1e6,             # 腿几乎不弯
        x_flex_max=0.03,
        contact_schedule=((0.0, 0.045), (0.1, 0.07), (0.5, 0.18), (1.0, 0.30)),
        collapse_ref_m=0.22,
        first_contact="双足外缘",
        note="踝扭伤/跟骨骨折经典姿势：足部接触面积极小。",
    ),
    "controlled-drop": Posture(
        name="屈膝缓冲（受控下跳）",
        m_low_frac=0.18,
        k_flex=2.2e5,
        x_flex_max=0.40,
        contact_schedule=((0.0, 0.050), (0.25, 0.09), (0.5, 0.20), (0.8, 0.45), (1.0, 0.58)),
        collapse_ref_m=0.52,
        first_contact="双足（全掌）",
        note="最优：吸能行程长、接触面积渐进展开，峰值力最低。",
    ),
    "tuck-roll": Posture(
        name="团身翻滚（脚→臀→背）",
        m_low_frac=0.20,           # 双足先接触，屈膝后大腿加入
        k_flex=1.6e5,
        x_flex_max=0.30,
        contact_schedule=((0.0, 0.055), (0.3, 0.13), (0.6, 0.32), (0.85, 0.52), (1.0, 0.60)),
        collapse_ref_m=0.46,
        first_contact="双足",
        note="接近最优：能量分散到多个部位，滚转动能进一步耗散。",
    ),
    "head-first": Posture(
        name="头朝下（失控）",
        m_low_frac=0.08,           # 头颈 = 唯一第一接触点
        k_flex=8.0e4,
        x_flex_max=0.25,
        contact_schedule=((0.0, 0.055), (0.3, 0.07), (1.0, 0.10)),
        collapse_ref_m=0.30,
        first_contact="头顶",
        note="最危险：极小接触面积 + 头部无缓冲结构。",
        first_contact_is_head=True,   # HIC 必须取头部（= m_low）的加速度
    ),
    # ---- Phase B-2 新增：软垫也救不了的姿势 ------------------------------
    "toe-point": Posture(
        name="前脚掌点地（垫偏/垫边）",
        m_low_frac=0.06,           # 前脚掌，是人体上最小的承重接触面之一
        k_flex=9.0e5,             # 踝关节几乎不屈
        x_flex_max=0.02,
        # 跖骨前缘接触，面积比全掌小一个量级
        contact_schedule=((0.0, 0.018), (0.2, 0.030), (0.6, 0.075), (1.0, 0.13)),
        collapse_ref_m=0.16,
        first_contact="跖骨前缘",
        note="常见于落点没对准垫子。接触面积极小 -> 压力极高，"
             "典型后果是跖骨骨折/前脚掌挫伤。",
    ),
    "shoulder-first": Posture(
        name="肩先着垫（横向失衡）",
        m_low_frac=0.16,           # 单侧肩 + 上臂
        k_flex=3.0e5,
        x_flex_max=0.10,
        contact_schedule=((0.0, 0.075), (0.25, 0.13), (0.7, 0.26), (1.0, 0.34)),
        collapse_ref_m=0.24,
        first_contact="单侧肩",
        note="接触面积极小且不居中，肩关节/锁骨承压；"
             "横向落地时身体会绕接触点翻滚，力矩放大。",
    ),
    "one-leg-awkward": Posture(
        name="单脚先落（踝外翻）",
        m_low_frac=0.10,           # 单腿
        k_flex=1.3e6,             # 支撑腿绷直
        x_flex_max=0.025,
        # 单脚，且落地姿态歪斜，接触集中在足外缘
        contact_schedule=((0.0, 0.030), (0.15, 0.045), (0.5, 0.10), (1.0, 0.18)),
        collapse_ref_m=0.20,
        first_contact="单足外缘",
        note="单腿 + 踝外翻是扭伤/腓骨骨折的经典机制；"
             "另一条腿无法分担，屈曲行程也用不上。",
    ),
    "superman-late": Posture(
        name="伸直腿滞后落地（够垫失败）",
        m_low_frac=0.20,
        k_flex=1.6e6,             # 双腿绷直
        x_flex_max=0.02,
        contact_schedule=((0.0, 0.060), (0.12, 0.10), (0.55, 0.30), (1.0, 0.52)),
        collapse_ref_m=0.26,
        first_contact="双足弓",
        note="想屈膝但来不及/够不到垫子。接触面积不小，"
             "但**没有可用的屈曲行程**，力全部由踝和脊柱承担。",
    ),
}


# --------------------------------------------------------------------------
# 结果容器
# --------------------------------------------------------------------------
@dataclass
class BoulderFallResult:
    posture: str
    pad_label: str
    height_m: float
    mass_kg: float

    t_s: np.ndarray
    z_foot_m: np.ndarray         # 首次接触点相对释放点的下落 (m)
    z_torso_m: np.ndarray        # 躯干 CoM 相对释放点的下落 (m)
    pad_compression_m: np.ndarray
    contact_area_m2: np.ndarray
    pad_force_n: np.ndarray
    accel_leg_g: np.ndarray      # 正 = 向上减速
    accel_torso_g: np.ndarray    # 正 = 向上减速

    impact_speed_ms: float
    impact_energy_j: float
    peak_force_n: float
    peak_accel_torso_g: float
    peak_accel_leg_g: float
    max_compression_m: float
    bottomed_out: bool
    energy_into_pad_j: float
    energy_into_flex_j: float
    energy_residual_j: float     # 未被弹性吸收 = 必须由组织损伤承担
    hic: float = 0.0
    meta: dict = field(default_factory=dict)

    @property
    def peak_force_kn(self) -> float:
        return self.peak_force_n / 1e3

    @property
    def compression_ratio(self) -> float:
        """压缩量 / 垫厚，用于判断是否压穿。"""
        return self.max_compression_m / max(1e-9, self.meta.get("pad_thickness_m", 1.0))

    @property
    def absorbed_fraction(self) -> float:
        tot = self.energy_into_pad_j + self.energy_into_flex_j
        return (tot / self.impact_energy_j) if self.impact_energy_j > 0 else 0.0

    @property
    def primary_accel_g(self) -> np.ndarray:
        """"受伤部位"的加速度时程 —— HIC 应该算在**第一接触部位**上。

        头朝下落地时 ``m_low`` 就是头（6.4 kg），而躯干有 73.6 kg。
        早期版本让 ``hic_from_result`` 一律取 ``accel_torso_g``，结果
        头朝下落地算出 HIC≈82、看起来是全场最安全的姿势 —— 完全反了。
        真实风险恰恰是那 6.4 kg 的头要独自扛下整个冲击。
        """
        if self.meta.get("head_is_contact"):
            return self.accel_leg_g
        return self.accel_torso_g

    @property
    def peak_primary_g(self) -> float:
        return float(np.max(self.primary_accel_g))


# --------------------------------------------------------------------------
# 求解
# --------------------------------------------------------------------------
def simulate_boulder_fall(
    *,
    height_m: float = 3.0,
    mass_kg: float = 80.0,
    posture: str | Posture = "controlled-drop",
    pad: CrashPad | None = None,
    on_pad: bool = True,
    g: float = G,
    dt_sample: float = 2e-4,
    t_max: float = 2.0,
    max_ode_step: float = 1e-4,
    rtol: float = 1e-8,
    atol: float = 1e-11,
) -> BoulderFallResult:
    """模拟一次抱石坠落。

    Parameters
    ----------
    height_m
        释放点距垫面的高度 (m)。自由落体 v = √(2gh)。
    on_pad
        True = 落在软垫上；False = **离垫坠落**，改用 ``hard_surface()``。
    posture
        姿势名（``POSTURES`` 的键）或自定义 ``Posture``。
    """
    post = POSTURES[posture] if isinstance(posture, str) else posture
    pad = (pad or CrashPad()) if on_pad else hard_surface()

    m = float(mass_kg)
    m_low, m_up = post.split_mass(m)
    # 躯干 CoM 在双足之上的初始高度差（站立时约 0.85 m）
    leg_ext0 = 0.85
    k_flex = post.k_flex
    x_flex_max = post.x_flex_max
    # 屈曲行程用尽后的止挡刚度。
    #
    # 这里的建模困境值得记下来：真实世界里"腿屈完了"之后，躯干是**直接
    # 落在软垫上**继续被吸收的。但本模型是双质点串联，躯干除了屈曲弹簧
    # 以外**没有任何其他通路** —— 它只能压在这个止挡上。于是：
    #   * 止挡太硬 (5e7) -> 12 kg 躯干撞钢板，报 430 g
    #   * 止挡 2e6      -> 仍然产生尖峰，平拍/坐落报到 150-160 g
    # 止挡刚度大到能吃掉全部能量，会让结果对 x_flex_max 极度敏感；太小
    # 又让"腿"形同虚设。这里取 3e5（与屈曲刚度同量级），把它理解为
    # "腿屈完之后的**渐进**塌陷"而不是刚性止挡，剩余能量继续由垫承担。
    # 真正正确的解法是 Phase B-3 的多刚体接触（躯干直接接触垫面）。
    k_stop = 3.0e5

    # 屈曲阻尼 (N·s/m)，代表肌腱/肌肉的粘性损耗。
    # 早期版本取 3.0e3，这是个**虚高的人为参数**：冲击时屈曲速率可达
    # 5-8 m/s，单这一项就产生 15-24 kN，占了"屈膝缓冲"总峰值力的
    # 大部分，把躯干加速度顶到 40 g、脊柱惯性载荷顶到 26 kN（判读全
    # 顶"极高"，分档失去区分度）。
    # 真实的肌肉吸能主要是**离心收缩做功**（力-位移-功），不是粘性
    # 阻尼；整条腿的等效粘性系数在 100-500 N·s/m 量级。这里取 400。
    # Phase B-2 的遗留项：用 Hill 型力-速关系 + 肌肉做功上限替代它。
    c_flex = 400.0

    # ---- 几何约定 ------------------------------------------------------
    # z 是"自释放点起的向下位移"，向下为正。躯干在双足**之上**，
    # 所以躯干的 z 比脚小 leg_ext0。腿长 = z_foot - z_torso（脚越深腿越长）。
    # 屈曲量 xf = 初始腿长 - 当前腿长，xf 增大 = 塌陷。
    # 注意初值必须是 z_torso = -leg_ext0；早期版本写成 +leg_ext0，
    # 等于把躯干放到了脚下方 0.85 m（埋在垫里），几何整个反了。
    def flex_disp(z_foot, z_torso):
        """腿的屈曲量 = 初始腿长 - 当前腿长 (m)，>0 表示在塌陷。"""
        return leg_ext0 - (z_foot - z_torso)

    def flex_force(z_foot, z_torso, vz_foot, vz_torso):
        """返回 ``(力的大小, 屈曲量)``。力的大小恒为非负，方向由调用方按牛顿第三定律分配。"""
        xf = flex_disp(z_foot, z_torso)
        if xf <= 0.0:
            return 0.0, 0.0
        if xf <= x_flex_max:
            f = k_flex * xf
        else:
            f = k_flex * x_flex_max + k_stop * (xf - x_flex_max)
        # dxf/dt = -(vz_foot - vz_torso)：躯干朝静止的脚砸下去时为正
        rate = vz_torso - vz_foot
        f_damp = c_flex * max(0.0, rate)     # 只在加载（塌陷）时耗能
        return f + f_damp, xf

    def rhs(_t, y):
        zf, vz_f, zt, vz_t = y
        # 垫面接触：首次接触点的下落超过 height 才压缩垫子
        comp = zf - height_m
        if comp > 0.0:
            kappa = comp / post.collapse_ref_m
            area = post.area_at(kappa)
            f_pad = float(pad.force_n(comp, area))
        else:
            area = 0.0
            f_pad = 0.0

        f_flex, _ = flex_force(zf, zt, vz_f, vz_t)

        # ---- 牛顿第三定律 ----------------------------------------------
        # f_flex 是**内力**，必须一上一下相互抵消：
        #   屈曲弹簧把躯干往上顶（+），同时把下段质量往下压（-）。
        # 早期版本把 f_flex 也算进下段质量的"向上"方向，
        # 合外力变成 f_pad + 2·f_flex —— 系统凭空获得能量，
        # 总动量不再守恒。这也正是加严 m_up 守卫后 ODE 直接发散的原因：
        # 之前 f_flex 恒为 0，这个错误被掩盖着。
        #
        #   m_low: 垫力向上(+)，屈曲反力向下(-)
        #   m_up : 屈曲力向上(+)，重力向下
        # 两者相加 -> 合力 = f_pad，符合 CoM 方程 m·a_com = f_pad - m·g。
        # m_low / m_up 已由 Posture.split_mass() 保证严格为正，
        # 这里不做 max(..., 1e-6) 兜底 —— 那种兜底会把除零伪装成一个
        # 量级看似正常的假加速度，让错误静默通过。
        az_f = g - f_pad / m_low + f_flex / m_low
        az_t = g - f_flex / m_up
        return [vz_f, az_f, vz_t, az_t]

    # ---- 终止事件：躯干停止且垫子已完全卸载 -------------------------------
    def settled(_t, y):
        comp = y[0] - height_m
        return 1.0 if (comp > 0.0 and abs(y[3]) < 0.05 and abs(y[1]) < 0.05) else -1.0

    settled.terminal = True
    settled.direction = 1.0

    sol = solve_ivp(
        rhs,
        (0.0, t_max),
        # y0 = [z_foot, vz_foot, z_torso, vz_torso]；躯干在脚上方 -> z_torso = -leg_ext0
        [0.0, 0.0, -leg_ext0, 0.0],
        events=[settled],
        max_step=max_ode_step,
        rtol=rtol,
        atol=atol,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(f"积分失败: {sol.message}")

    t_out = np.arange(0.0, sol.t[-1] + dt_sample, dt_sample)
    t_out = t_out[t_out <= sol.t[-1] + 1e-12]
    if t_out.size < 4 or t_out[-1] < sol.t[-1] - 1e-9:
        t_out = np.append(t_out, sol.t[-1])
    zf, zt = sol.sol(t_out)[0], sol.sol(t_out)[2]

    comp = np.maximum(0.0, zf - height_m)
    kappa = comp / post.collapse_ref_m
    area = np.where(comp > 0.0, np.array([post.area_at(k) for k in kappa]), 0.0)
    f_pad = np.where(comp > 0.0, pad.force_n(comp, np.maximum(area, 1e-9)), 0.0)
    f_pad = np.maximum(f_pad, 0.0)

    # 与 rhs() 用同一套几何/牛顿第三定律约定
    xf = np.maximum(0.0, leg_ext0 - (zf - zt))
    rate_flex = np.maximum(0.0, np.gradient(xf, t_out))
    vz_f = np.gradient(zf, t_out)
    vz_t = np.gradient(zt, t_out)
    f_flex = np.where(
        xf > 0.0,
        np.where(xf <= x_flex_max, k_flex * xf, k_flex * x_flex_max + k_stop * (xf - x_flex_max))
        + c_flex * rate_flex,
        0.0,
    )

    # 悬挂/减速方向朝上为正。屈曲力对下段质量是**向下**的（牛顿第三定律），
    # 早期版本写成 +，导致 a_leg 被严重高估（实测 1000+ g）。
    # 注意：这里必须 /g 归一到"g"——字段名以 _g 结尾，且 rope.py 的
    # accel_g 同样以 g 为单位。早期版本漏掉了这步除法，pad 侧报的是 m/s^2，
    # 与 rope 侧差 9.8 倍，两张表直接对比会得出完全错误的结论。
    a_leg = ((f_pad - f_flex) / m_low - g) / g
    a_torso = (f_flex / m_up - g) / g

    # ---- 冲击窗口：只取"第一次接触 -> 第一次分离" ------------------------
    # 峰值力 / 峰值加速度 / 能量账都必须在**单次冲击**内统计。3 m 落到
    # 20 cm 垫上会压穿并弹起，系统随后进入第二轮自由落体+再冲击；
    # 若对全程取 max，第二次甚至第三次的振荡会盖掉首次冲击的真实峰值，
    # 而"能量吸收"取全程 running-integral 的最大值则会算出
    # E_pad + E_flex > m·g·h（吸能超过输入能量）这种自相矛盾的结果。
    in_contact = comp > 0.0
    if in_contact.any():
        i0 = int(np.argmax(in_contact))
        gap = np.where(~in_contact[i0:])[0]
        i1 = (i0 + int(gap[0])) if gap.size else t_out.size
    else:
        i0, i1 = 0, 0
    i1 = max(i1, i0 + 2)
    i1 = min(i1, t_out.size)
    sl = slice(i0, i1)

    peak_force = float(np.max(f_pad[sl])) if f_pad[sl].size else 0.0
    max_comp = float(np.max(comp[sl]))
    bottomed = bool(max_comp > pad.max_compression_m * 0.999)

    # 能量账：单次冲击窗口内对**加载段**做 ∫F·dx（功率的时间累积）。
    # 必须取同一个时刻的快照：早期版本让 E_pad 和 E_flex 各自取
    # 全程 running-integral 的最大值，两个峰值发生在不同时刻，
    # 相加就超过了真实同时吸收量。np.trapezoid(f, x) 更糟 ——
    # 回弹段 x 减小会把能量又减回来，弹性体净积分为 0，
    # 于是 e_pad 算出来是 -0.098 J（负数），"absorbed = -0.0%"。
    t_win = t_out[sl]
    if t_win.size < 2:
        e_pad = e_flex = 0.0
    else:
        dt_w = np.diff(t_win)
        p_pad = f_pad[sl] * vz_f[sl]
        p_flex = f_flex[sl] * np.gradient(xf, t_out)[sl]
        e_pad_t = np.concatenate([[0.0], np.cumsum(0.5 * (p_pad[:-1] + p_pad[1:]) * dt_w)])
        e_flex_t = np.concatenate([[0.0], np.cumsum(0.5 * (p_flex[:-1] + p_flex[1:]) * dt_w)])
        tot_t = e_pad_t + e_flex_t
        k = int(np.argmax(tot_t))          # 总吸能最大的那一刻，两项取同一快照
        e_pad = float(max(0.0, e_pad_t[k]))
        e_flex = float(max(0.0, e_flex_t[k]))
    e_impact = m * g * height_m          # 释放到静止释放的总势能

    v_impact = float(np.sqrt(2.0 * g * height_m))

    return BoulderFallResult(
        posture=post.name,
        pad_label=("Crash pad" if on_pad else "硬地面（离垫）"),
        height_m=height_m,
        mass_kg=m,
        t_s=t_out,
        z_foot_m=zf,
        z_torso_m=zt,
        pad_compression_m=comp,
        contact_area_m2=area,
        pad_force_n=f_pad,
        accel_leg_g=a_leg,
        accel_torso_g=a_torso,
        impact_speed_ms=v_impact,
        impact_energy_j=e_impact,
        peak_force_n=peak_force,
        peak_accel_torso_g=float(np.max(a_torso[sl])),
        peak_accel_leg_g=float(np.max(a_leg[sl])),
        max_compression_m=max_comp,
        bottomed_out=bottomed,
        energy_into_pad_j=e_pad,
        energy_into_flex_j=e_flex,
        energy_residual_j=max(0.0, e_impact - e_pad - e_flex),
        meta={"pad_thickness_m": pad.thickness_m, "k_flex": k_flex,
              "x_flex_max": x_flex_max, "on_pad": on_pad, "posture_obj": post,
              "head_is_contact": post.first_contact_is_head,
              # 供 injury.py 估算脊柱轴向载荷
              "f_flex_peak_n": float(np.max(f_flex[sl])) if f_flex[sl].size else 0.0,
              "m_low": m_low, "m_up": m_up,
              "impact_window_s": (float(t_out[i0]), float(t_out[i1 - 1]))},
    )
