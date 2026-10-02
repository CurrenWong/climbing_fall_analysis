"""Phase B-3 —— 多段人体链落地模型（躯干可直接接触垫面 + Hill 型肌肉）。

B-2 留下的两个结构性问题
------------------------
Phase B-2 用的是**双质点串联**：下肢（m_low）压软垫 -> 屈曲弹簧 -> 躯干（m_up）。
后果是：

1. **躯干除屈曲弹簧外没有第二条通路。** 现实里"腿屈完了，躯干直接落在
   垫上继续被吸收"这个过程在双质点模型里无法表达 —— 躯干只能压在一个
   止挡上。受控屈膝在 B-2 里屈曲行程只走到 **8.7 cm**（姿势允许 40 cm），
   于是 ``a_torso ≈ v²/2x`` 被高估 1.5-2 倍。
2. **损伤绝对判读全部顶"极高"。** 由上条的加速度高估传导到脊柱/下肢的
   惯性载荷 ``m·a``，3 m 落差下无一姿势落在"低"档。

本模块把人体换成**沿竖直方向的 6 节点链**::

    5 头顶  ── 颈/头 (link4)
    4 肩/臀 ── 躯干   (link3)
    3 髋    ── 大腿   (link2)
    2 膝    ── 小腿   (link1)
    1 踝    ── 足     (link0)
    0 足底  ── 接触面

关键改进：

* **任意节点都能接触垫面。** 腿屈完后 3(髋/臀)/4(肩) 会直接落到垫上，
  由软垫继续吸能 —— 这正是 B-2 表达不了的那一段。
* **关节有真实行程上限。** 膝 ~0.5 m、髋 ~0.45 m、踝 ~0.15 m，腰椎很硬。
  B-2 只有"屈曲弹簧 + 止挡"两段。
* **Hill 型力-速关系替代粘性阻尼。** 肌肉吸能主要靠离心收缩做功，
  不是 ``c·ė`` 那种粘性阻尼。B-2 的 ``c_flex=400 N·s/m`` 是个凑出来的
  等效值，而它曾经（取 3000 时）单独贡献了峰值力的绝大部分。

数据性质
--------
各段质量取 80 kg 成年男性的公开人体测量学量级（腿 ~18%、躯干+骨盆 ~47%、
头颈 ~8%、上肢 ~16%）。关节刚度、行程、肌肉参数是**量级估计**。

⚠️ **当前状态：结构正确，参数未标定，数值不可用**
--------------------------------------------------
`verify_b3.py` 的两条验收判据**都未通过**，不要引用本模块输出的绝对数值：

* 判据 A 失败：屈膝缓冲的总腿屈曲行程只有 **6.8 cm**（要求 20-50 cm）
* 判据 B 失败：躯干峰值加速度 68.6 g，**高于** B-2 的 37.2 g
  （本该至少降低 25%）

已确认**正确**的部分（可以放心继承）：
* 牛顿第三定律：内力一上一下配平
* 状态布局：交错 [z0,v0,z1,v1,...]，rhs 返回同布局
* 能量守恒：实测 2292 J vs 起始 2354 J，误差 **2.6%**
* 关节**双向**：受压有力、过伸有限（否则人体会散架）
* 多节点接触：髋/肩在腿屈完后能直接落到垫上

未解决的**物理问题**（诊断结论）：
要让 80 kg 的人从 7.67 m/s 停下来，需要约 28 kN 的减速力；而当前模型里
垫峰值只有 3 kN、膝关节峰值 2.5 kN，合计约 7 g 减速 —— **差约 10 倍**。
把肌肉等长力 f_iso 放大 12 倍几乎不改变结果（终态动能 560→523 J、
屈曲 6.8→3.8 cm），说明瓶颈不在肌肉参数。

更可能的根因是**接触时序**：足节点只有 1.2 kg，撞上软垫后在约 10 ms 内
就弹回，其余 79 kg 还没来得及"接上"，链条被自己的抗拉止挡拉回去弹跳。
真实落地是质量逐渐接入的（脚→小腿→大腿→臀），需要**渐进接触**或
真正的多刚体接触求解（MuJoCo），不是调参数能解决的。

下一步应做的是：把足/踝段的质量与接触面积按人体测量学重新分配，
并给链加入**接触面随塌陷展开**的机制（复刻 B-2 里证明有效的那部分），
再重新标定。在此之前 B-2 的 `pad.py` 仍是主用模型。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from . import G
from .pad import CrashPad, hard_surface, POSTURES

__all__ = [
    "Link",
    "ChainPosture",
    "CHAIN_POSTURES",
    "ChainFallResult",
    "simulate_chain_fall",
    "default_chain",
]


# ==========================================================================
# 结构定义
# ==========================================================================
@dataclass(frozen=True)
class Link:
    """链上的一段（关节）。

    Attributes
    ----------
    name
    length_m
        静止长度（m）。
    k_elast_pa
        弹性刚度 (N/m)。代表韧带/被动软组织。
    f_iso_n
        等长收缩力 (N)。Hill 模型的基准张力。
    v_max_ms
        肌纤维最大缩短速度 (m/s)，决定力-速关系的陡度。
    x_max_m
        可用屈曲行程 (m)。超过后由止挡承担。
    k_stop_n_per_m
        行程用尽后的止挡刚度 (N/m)。
    m_frac
        该段占总质量的比例。
    node_index
        该段**上端**节点编号（链从足底 node 0 往上数）。
    """

    name: str
    length_m: float
    k_elast_pa: float
    f_iso_n: float
    v_max_ms: float
    x_max_m: float
    k_stop_n_per_m: float
    m_frac: float
    node_index: int
    k_tension_n_per_m: float = 2.0e4   # 抗拉刚度（韧带/肌腱）
    x_ext_max_m: float = 0.02          # 向后（过伸）极限 (m)

    def force_n(self, x: float, rate: float, f_iso_scale: float = 1.0,
                x_travel: float | None = None) -> float:
        """关节力 (N)。

        **返回值带符号**：>0 = 压缩（把两端推开），<0 = 拉伸（把两端拉拢）。

        抗拉这一支不能省。最初写成 ``if x <= 0: return 0``，等于关节
        可以无限伸长 —— 脚撞上垫后只有 1.18 kg 被减速，其余 79 kg 继续
        下落把踝关节拉到 **-54 cm**（人体散架：末态头节点在 z=-4.6 m、
        还在以 3.1 m/s 下落），终态仍有 429 J 动能，能量缺口 2001 J。
        真实关节由韧带/肌腱限制过伸，踝/膝的反向行程只有几毫米到 2 cm。
        """
        # ---- 拉伸侧（x < 0）：韧带抗拉 ----------------------------------
        if x < 0.0:
            over = x + self.x_ext_max_m          # 超过过伸极限的量
            f = self.k_tension_n_per_m * x       # 基础抗拉（可为负 = 拉回）
            if over < 0.0:                       # 超过极限，止挡
                f += 5.0e6 * over
            return min(f, 0.0)

        xt = self.x_max_m if x_travel is None else x_travel
        xt = max(xt, 1e-3)
        # 被动弹性：韧带/软组织
        f = self.k_elast_pa * x

        # ---- Hill 型力-速关系 ------------------------------------------
        # 缩短（rate<0）时张力下降：F = F_iso (1 - v/v_max)
        # 拉长（rate>0，即离心收缩）时张力**升高**且不降为零：
        # 真实肌肉离心收缩能产生约 1.5-2 倍等长力，这正是落地吸能的主力。
        # 这也是 B-2 的粘性阻尼项 c_flex 所无法表达的。
        v = -rate                     # v > 0 = 缩短
        if v >= 0.0:
            hill = 1.0 - min(v / self.v_max_ms, 0.9)
        else:
            hill = 1.0 + min(-v / self.v_max_ms, 1.0) * 0.8
        # 被动张力随肌肉长度增加（力-长度关系）
        passive = min(x / max(0.5 * xt, 1e-6), 1.0)
        f += self.f_iso_n * f_iso_scale * hill * passive

        # ---- 行程用尽后的止挡 -------------------------------------------
        if x > xt:
            f += self.k_stop_n_per_m * (x - xt)
        return f


def joint_forces(links: tuple[Link, ...], x: np.ndarray, rate: np.ndarray,
                 f_iso_scale: float, x_avail: np.ndarray) -> np.ndarray:
    """一次性算出全部关节的力 (T, n_links)。

    与 ``rhs`` 里用的是**同一个** ``force_n``，保证后处理与积分一致。
    """
    out = np.zeros_like(x)
    for i, lk in enumerate(links):
        f = np.zeros(x.shape[0])
        for k in range(x.shape[0]):
            f[k] = lk.force_n(x[k, i], rate[k, i], f_iso_scale, x_avail[i])
        out[:, i] = f
    return out


#: 6 节点链。node 0 = 足底，node 5 = 头顶。
DEFAULT_LINKS: tuple[Link, ...] = (
    Link("足/踝", 0.07, 3.0e5, 250.0, 1.2, 0.15, 2.0e6, 0.015, 1),
    Link("小腿/膝", 0.42, 6.0e4, 900.0, 0.9, 0.50, 3.0e5, 0.065, 2),
    Link("大腿/髋", 0.42, 8.0e4, 1100.0, 0.9, 0.45, 3.0e5, 0.100, 3),
    Link("躯干/腰椎", 0.55, 2.5e5, 700.0, 0.5, 0.20, 5.0e5, 0.470, 4),
    Link("颈/头", 0.22, 2.0e5, 150.0, 0.8, 0.10, 2.0e5, 0.350, 5),
)


def default_chain() -> tuple[Link, ...]:
    return DEFAULT_LINKS


# ==========================================================================
# 姿势
# ==========================================================================
@dataclass(frozen=True)
class ChainPosture:
    """一种落地姿势在链模型上的参数。

    与 B-2 的 ``Posture`` 一一对应（名字/接触部位），额外给出链特有的
    参数：哪些节点能接触垫面、初始关节预屈曲。

    Attributes
    ----------
    contact_nodes
        能接触垫面的节点编号。``(0,)`` = 只有脚；``(0, 3, 4)`` = 脚 + 髋/臀 + 肩。
    preflex_m
        各关节的初始预压缩量 (m)。"平拍"落地时腿已经是弯的。
    scale_f_iso
        肌肉激活系数（相对默认）。腿绷直/不屈曲的姿势接近 0（没有主动缓冲）。
    contact_area_m2
        每个接触节点的接触面积 (m²)。
    """

    name: str
    base: str                      # 对应 B-2 POSTURES 的键
    contact_nodes: tuple[int, ...]
    preflex_m: tuple[float, ...]  # 与 LINKS 等长
    scale_f_iso: float
    contact_area_m2: dict[int, float]
    note: str = ""

    def area_at(self, node: int) -> float:
        return self.contact_area_m2.get(node, 0.0)


#: 与 ``pad.POSTURES`` 同名的链姿势。
CHAIN_POSTURES: dict[str, ChainPosture] = {
    "flat-flop": ChainPosture(
        name="平拍（背/臀同时着垫）", base="flat-flop",
        contact_nodes=(0, 3, 4), preflex_m=(0.02, 0.10, 0.08, 0.02, 0.0),
        scale_f_iso=0.25,
        contact_area_m2={0: 0.30, 3: 0.32, 4: 0.20},
        note="腿几乎不屈曲，接触面立刻铺到臀/背 —— 压力大但行程短。"),
    "butt-impact": ChainPosture(
        name="坐落（臀先着垫）", base="butt-impact",
        contact_nodes=(0, 3, 4), preflex_m=(0.02, 0.06, 0.06, 0.02, 0.0),
        scale_f_iso=0.45, contact_area_m2={0: 0.16, 3: 0.26, 4: 0.18},
        note="臀先着垫，脊柱压缩风险。"),
    "feet-first-stiff": ChainPosture(
        name="脚先落但腿绷直", base="feet-first-stiff",
        contact_nodes=(0, 3), preflex_m=(0.0, 0.0, 0.0, 0.0, 0.0),
        scale_f_iso=0.12, contact_area_m2={0: 0.045, 3: 0.24},
        note="腿几乎不弯，接触面积极小。"),
    "controlled-drop": ChainPosture(
        name="屈膝缓冲（受控下跳）", base="controlled-drop",
        contact_nodes=(0, 3, 4), preflex_m=(0.01, 0.06, 0.06, 0.01, 0.0),
        scale_f_iso=1.0, contact_area_m2={0: 0.050, 3: 0.30, 4: 0.26},
        note="最优：肌肉充分激活，膝/髋有 0.5/0.45 m 可用行程。"),
    "tuck-roll": ChainPosture(
        name="团身翻滚（脚→臀→背）", base="tuck-roll",
        contact_nodes=(0, 3, 4), preflex_m=(0.02, 0.14, 0.16, 0.03, 0.02),
        scale_f_iso=0.85, contact_area_m2={0: 0.055, 3: 0.22, 4: 0.30},
        note="接近最优：已预屈曲，接触面积随塌陷逐段展开。"),
    "head-first": ChainPosture(
        name="头朝下（失控）", base="head-first",
        contact_nodes=(5, 0, 3, 4), preflex_m=(0.0, 0.0, 0.0, 0.0, 0.0),
        scale_f_iso=0.20, contact_area_m2={0: 0.05, 3: 0.20, 4: 0.16, 5: 0.055},
        note="头是唯一第一接触点（node 5），质量占比小。"),
    "toe-point": ChainPosture(
        name="前脚掌点地（垫偏/垫边）", base="toe-point",
        contact_nodes=(0, 3, 4), preflex_m=(0.01, 0.0, 0.0, 0.0, 0.0),
        scale_f_iso=0.15, contact_area_m2={0: 0.018, 3: 0.26, 4: 0.24},
        note="跖骨前缘接触，面积极小。"),
    "shoulder-first": ChainPosture(
        name="肩先着垫（横向失衡）", base="shoulder-first",
        contact_nodes=(4, 0, 3), preflex_m=(0.01, 0.05, 0.05, 0.02, 0.01),
        scale_f_iso=0.40, contact_area_m2={0: 0.06, 3: 0.24, 4: 0.075},
        note="肩是第一接触点。"),
    "one-leg-awkward": ChainPosture(
        name="单脚先落（踝外翻）", base="one-leg-awkward",
        contact_nodes=(0, 3, 4), preflex_m=(0.01, 0.0, 0.0, 0.0, 0.0),
        scale_f_iso=0.18, contact_area_m2={0: 0.030, 3: 0.24, 4: 0.22},
        note="单脚 + 踝外翻，另一条腿无法分担。"),
    "superman-late": ChainPosture(
        name="伸直腿滞后落地（够垫失败）", base="superman-late",
        contact_nodes=(0, 3, 4), preflex_m=(0.0, 0.0, 0.0, 0.0, 0.0),
        scale_f_iso=0.10, contact_area_m2={0: 0.060, 3: 0.28, 4: 0.26},
        note="想屈膝但来不及/够不到垫子，肌肉没时间激活。"),
}


# ==========================================================================
# 结果容器
# ==========================================================================
@dataclass
class ChainFallResult:
    posture: str
    pad_label: str
    height_m: float
    mass_kg: float

    t_s: np.ndarray
    node_z: np.ndarray            # (T, 6) 各节点向下位移
    node_v: np.ndarray            # (T, 6) 各节点速度
    joint_x: np.ndarray           # (T, 5) 各关节压缩量
    joint_force_n: np.ndarray     # (T, 5) 各关节力
    pad_force_node: np.ndarray    # (T, 6) 各节点的垫反力

    # 关节总屈曲行程
    knee_x: np.ndarray
    hip_x: np.ndarray
    ankle_x: np.ndarray
    total_leg_flex_m: float       # 踝+膝+髋串联行程

    # 峰值量
    peak_force_n: float
    peak_pressure_kpa: float
    peak_accel_head_g: float
    peak_accel_trunk_g: float
    peak_accel_foot_g: float
    max_spine_load_kn: float
    max_knee_load_kn: float
    max_hip_load_kn: float
    impact_speed_ms: float
    impact_energy_j: float
    bottomed_out: bool
    first_contact_node: int
    meta: dict = field(default_factory=dict)

    @property
    def peak_force_kn(self) -> float:
        return self.peak_force_n / 1e3

    def node_accel_g(self, node: int) -> np.ndarray:
        return np.gradient(self.node_v[:, node], self.t_s) / G

    def primary_accel_g(self) -> np.ndarray:
        """HIC 应该算在哪一段上 —— 第一接触的那个节点。"""
        return np.gradient(self.node_v[:, self.first_contact_node], self.t_s) / G


# ==========================================================================
# 求解
# ==========================================================================
def simulate_chain_fall(
    *,
    height_m: float = 3.0,
    mass_kg: float = 80.0,
    posture: str | ChainPosture = "controlled-drop",
    pad: CrashPad | None = None,
    on_pad: bool = True,
    links: tuple[Link, ...] = DEFAULT_LINKS,
    g: float = G,
    dt_sample: float = 5e-4,
    t_max: float = 2.0,
    max_ode_step: float = 2e-5,
    rtol: float = 1e-8,
    atol: float = 1e-11,
    radau: bool = True,
) -> ChainFallResult:
    """模拟一次抱石坠落（多段链模型）。

    Parameters
    ----------
    height_m
        释放点距**足底**的高度 (m)。足底触地时下落量 = height_m。
    posture
        ``CHAIN_POSTURES`` 的键或自定义 ``ChainPosture``。
    on_pad
        False = 离垫坠落，改用 ``hard_surface()``。
    """
    post = CHAIN_POSTURES[posture] if isinstance(posture, str) else posture
    pad = (pad or CrashPad()) if on_pad else hard_surface()
    m = float(mass_kg)
    n_links = len(links)
    n_nodes = n_links + 1

    # ---- 节点质量 ------------------------------------------------------
    m_nodes = np.zeros(n_nodes)
    for lk in links:
        m_nodes[lk.node_index] += m * lk.m_frac
    m_nodes[0] = m * links[0].m_frac          # 足底节点带一段足的质量
    m_nodes[-1] = m * links[-1].m_frac       # 头顶节点带头/上肢质量
    # 归一化到总质量
    m_nodes *= m / m_nodes.sum()

    # ---- 初始几何 + 关节的参考静止长度 --------------------------------
    # **预屈曲必须体现在"参考静止长度"上，不能只是压缩量偏置。**
    # 若把 preflex 当作压缩量偏置（最初的写法），t=0 时关节就带着 kN 级
    # 预应力：1.18 kg 的足节点拿到 **3043 m/s²** 的初始加速度，身体在
    # 释放瞬间被弹开；而只改 z0 又没用 —— gap 相对**伸展**长度 L 算，
    # 压缩量还是 pf_i。
    # 正确做法：蹲姿就是"静止长度更短"，压缩量从这个姿态的参考构型算起。
    L_rest = np.array([lk.length_m for lk in links])
    for i in range(n_links):
        pf = post.preflex_m[i] if i < len(post.preflex_m) else 0.0
        L_rest[i] -= pf
    # 该姿态下**剩余**可用行程（预屈曲已经用掉一部分）
    x_avail = np.array([lk.x_max_m for lk in links]) - np.array(
        [post.preflex_m[i] if i < len(post.preflex_m) else 0.0
         for i in range(n_links)])

    z0 = np.zeros(n_nodes)
    z0[0] = -height_m
    for i in range(n_links):
        z0[i + 1] = z0[i] - L_rest[i]      # 向上为负 z

    contact_nodes = post.contact_nodes

    def joint_state(y):
        """返回 (压缩量, 屈曲速率)，形状都是 (n_links,) 或 (n_links, T)。

        节点 i 在下、i+1 在上，z 向下为正，所以 gap = z_i - z_{i+1}。
        压缩量相对**该姿态的参考静止长度** ``L_rest`` 计算（见上），
        所以 t=0 时所有关节 x=0，身体是纯粹的自由落体。
        """
        z = y[0::2]
        v = y[1::2]
        gap = z[:-1] - z[1:]
        L = L_rest
        if z.ndim > 1:                       # (n_links, T)
            L = L.reshape(-1, 1)
        x = L - gap
        rate = -(v[:-1] - v[1:])             # 压缩 = 上方节点相对下方下移
        return x, rate

    def rhs(_t, y):
        z = y[0::2]
        v = y[1::2]
        x, rate = joint_state(y)

        f_joint = np.zeros(n_links)
        for i, lk in enumerate(links):
            f_joint[i] = lk.force_n(x[i], rate[i], post.scale_f_iso,
                                    x_travel=x_avail[i])

        a = np.zeros(n_nodes)
        for i in range(n_nodes):
            # 关节力：被压缩的 link 把下端节点往下推、上端往上拉
            if i < n_links:
                a[i] += f_joint[i]             # link i 压 node i 向下
            if i > 0:
                a[i] -= f_joint[i - 1]         # link i-1 压 node i 向上
        a += g

        # ---- 垫接触：任意接触节点都能被垫托住 ------------------------
        for nd in contact_nodes:
            comp = z[nd] - 0.0                 # 垫面在 z=0
            if comp > 0.0:
                area = post.area_at(nd)
                f = float(pad.force_n(comp, area))
                a[nd] -= f / m_nodes[nd]        # 垫力向上

        # 状态是**交错**布局 [z0, v0, z1, v1, ...]，返回必须同布局。
        # 最初写成 np.concatenate([v, a])（分块布局），与交错的 y0 不匹配，
        # 第一步之后状态就错位成 z/v 混排，然后在第 ~440 步炸掉。
        out = np.zeros(2 * n_nodes)
        out[0::2] = v
        out[1::2] = a
        return out

    # ---- 终止：所有节点静止 --------------------------------------------
    # 必须同时要求**已经接触过垫面**。最初只判速度全零，而初始速度就是零，
    # 事件在 t=0 立刻触发，积分直接结束（t_out 只剩 2 个点）。
    def settled(_t, y):
        z = y[0::2]
        v = y[1::2]
        touched = any(z[nd] > 0.0 for nd in contact_nodes)
        return 1.0 if (touched and np.all(np.abs(v) < 0.05)) else -1.0

    settled.terminal = True
    settled.direction = 1.0

    y0 = np.zeros(2 * n_nodes)
    y0[0::2] = z0
    y0[1::2] = 0.0

    # ---- 积分器 --------------------------------------------------------
    # 链模型是**刚性**的：抗拉止挡 5e6 N/m 作用在 1.2 kg 的足节点上
    # （ω ≈ 2000 rad/s），加上软垫的压实段和 Hill 项，显式 RK45 在
    # max_step=2e-5 下仍然失稳 —— 实测末态动能 705 J、能量缺口 3130 J，
    # 也就是系统在**凭空获得能量**，整条链弹跳不止。
    # 所以默认用隐式的 Radau（自带自适应刚性处理），显式只留给调试。
    method = "Radau" if radau else "RK45"
    sol = solve_ivp(rhs, (0.0, t_max), y0, events=[settled], method=method,
                    max_step=max_ode_step, rtol=rtol, atol=atol,
                    dense_output=True)
    if not sol.success:
        raise RuntimeError(f"积分失败: {sol.message}")

    t_out = np.arange(0.0, sol.t[-1] + dt_sample, dt_sample)
    t_out = t_out[t_out <= sol.t[-1] + 1e-12]
    if t_out.size < 4 or t_out[-1] < sol.t[-1] - 1e-9:
        t_out = np.append(t_out, sol.t[-1])

    Y = sol.sol(t_out)              # (2*n_nodes, T) —— **节点优先**
    x, rate = joint_state(Y)        # (n_links, T)
    # 转成时间优先，与 f_pad_node / t_out 对齐。
    # 保持节点优先再写 `z[:, nd]` 会切到**时间**轴上（返回长度 n 的数组），
    # 而不是取某个节点。
    z = Y[0::2].T                   # (T, n_nodes)
    v = Y[1::2].T                   # (T, n_nodes)
    x = x.T                         # (T, n_links)
    rate = rate.T                   # (T, n_links)

    # ---- 各节点垫反力与压力 --------------------------------------------
    f_pad_node = np.zeros((t_out.size, n_nodes))
    for nd in contact_nodes:
        comp = z[:, nd]
        f_pad_node[:, nd] = np.where(comp > 0.0,
                                    pad.force_n(np.maximum(comp, 0.0),
                                                post.area_at(nd)), 0.0)
    f_total = f_pad_node.sum(axis=1)
    press = np.zeros_like(f_total)
    for nd in contact_nodes:
        area = max(post.area_at(nd), 1e-9)
        press = np.maximum(press, f_pad_node[:, nd] / area)
    p_peak_kpa = float(np.max(press)) / 1e3

    # ---- 冲击窗口：第一次接触 -> 第一次分离 ----------------------------
    any_contact = f_total > 0
    if any_contact.any():
        i0 = int(np.argmax(any_contact))
        gap = np.where(~any_contact[i0:])[0]
        i1 = (i0 + int(gap[0])) if gap.size else t_out.size
    else:
        i0, i1 = 0, 0
    i1 = min(max(i1, i0 + 2), t_out.size)
    sl = slice(i0, i1)

    # ---- 峰值与载荷 ----------------------------------------------------
    a_head = np.gradient(v[:, n_nodes - 1], t_out) / g
    a_trunk = np.gradient(v[:, 4], t_out) / g
    a_foot = np.gradient(v[:, 0], t_out) / g

    peak_force = float(np.max(f_total[sl])) if f_total[sl].size else 0.0
    # 脊柱轴向载荷 = 头部节点（含上肢质量）的减速惯性力峰值
    m_head = float(m_nodes[n_nodes - 1])
    a_head_abs = np.abs(np.gradient(v[:, n_nodes - 1], t_out))
    spine_kn = float(np.max(a_head_abs[sl])) * m_head / 1e3
    # 关节力（与 rhs 用同一个 force_n，保证前后处理一致）
    rate_all = np.gradient(x, t_out, axis=0)      # x 是 (T, n_links)
    fj = joint_forces(links, x, rate_all, post.scale_f_iso, x_avail)
    knee_kn = float(np.max(fj[sl, 1])) / 1e3 if fj[sl].size else 0.0
    hip_kn = float(np.max(fj[sl, 2])) / 1e3 if fj[sl].size else 0.0

    max_comp = 0.0
    for nd in contact_nodes:
        max_comp = max(max_comp, float(np.max(np.maximum(z[sl, nd], 0.0))))
    bottomed = max_comp > pad.max_compression_m * 0.999

    # 首次接触节点 = 垫力最先起来的那一段
    first = int(np.argmax(np.any(f_pad_node[sl] > 0, axis=0))) if f_total[sl].size else 0

    return ChainFallResult(
        posture=post.name,
        pad_label=("Crash pad" if on_pad else "硬地面（离垫）"),
        height_m=height_m, mass_kg=m,
        t_s=t_out, node_z=z, node_v=v, joint_x=x, pad_force_node=f_pad_node,
        knee_x=x[:, 1], hip_x=x[:, 2], ankle_x=x[:, 0],
        joint_force_n=fj,
        total_leg_flex_m=float(np.max(x[sl, 0] + x[sl, 1] + x[sl, 2]))
        if x[sl].size else 0.0,
        peak_force_n=peak_force,
        peak_pressure_kpa=p_peak_kpa,
        peak_accel_head_g=float(np.max(np.abs(a_head[sl]))),
        peak_accel_trunk_g=float(np.max(np.abs(a_trunk[sl]))),
        peak_accel_foot_g=float(np.max(np.abs(a_foot[sl]))),
        max_spine_load_kn=spine_kn,
        max_knee_load_kn=knee_kn,
        max_hip_load_kn=hip_kn,
        impact_speed_ms=float(np.sqrt(2.0 * g * height_m)),
        impact_energy_j=m * g * height_m,
        bottomed_out=bottomed,
        first_contact_node=first,
        meta={"posture_obj": post, "links": links,
              "head_is_contact": (first == n_nodes - 1)},
    )
