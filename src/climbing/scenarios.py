"""数据驱动的坠落场景库 —— 方案 v2 · P1。

为什么需要这个模块
------------------
v1 的"场景"是我们**手写**的：10 种姿势 + 拍脑袋的高度/权重。
[B25]（Beurienne 2025, *Front Sports Act Living* 7:1609133，n=245）
给出了真实频次，本模块把它变成**可追溯到文献的权重 + 可复现的采样器**。

它解决一个具体问题：
> 把"哪种姿势最危险"升级成
> **"按真实频次加权，最常导致损伤的是哪条通路"**。

⚠️ 三条必须写在前面的诚实声明
------------------------------
1. **采样假设各维度独立。** [B25] 只给了边际分布（Table 3）与 3 个联合场景
   （Figure 4），**没给完整的联合分布**。原文明确指出旋转与落地姿势相关
   （"a rotation during the fall more frequently occurred…"），
   所以本采样器**高估了不常见组合**。联合场景另有 ``NAMED_SCENARIOS`` 表，
   用于覆盖度计算（那里用实测联合概率，不用独立性假设）。

2. **"高度带 → 米"是映射假设，不是实测值。** [B25] 只给
   墙底 20% / 墙中 31% / 墙顶 47%，没有米数。见 ``HEIGHT_BANDS``。

3. **"落地姿势 → POSTURES"是映射假设。** [B25] 的 5 个落地类别与
   ``pad.POSTURES`` 的 10 种姿势不是一一对应，见 ``LANDING_TO_POSTURE``
   逐条的 rationale。这是本模块**最需要被质疑**的部分。

来源标记
--------
每条 Category 都带 ``source``。凡 ``source`` 里含 ``"假设"`` 的，
都是**本项目的假设**而非文献实测值 —— 换更细的人群数据即可替换。

数据陷阱（我踩过）
------------------
[B25] 有两张容易混的表：

* **Table 3** —— 全体 n=245 的边际分布  ← **权重取这张**
* **Table 4** —— 按受伤部位拆成 Upper limbs(76) / Lower limbs(201) 两列

初稿误把 **Table 4 的 Lower-limbs 列**当全体分布，于是得到"脚先落地 87%"
（真值 73%）、"无旋转 39%"（真值 30%）。差 9 个百分点足以改变结论。
详见 [`.learnings/LEARNINGS.md` LRN-20261003-008](../.learnings/LEARNINGS.md)。
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

__all__ = [
    "Category",
    "Dimension",
    "DIMENSIONS",
    "HEIGHT_BANDS",
    "WALL_HEIGHT_M",
    "LANDING_TO_POSTURE",
    "PostureMapping",
    "NAMED_SCENARIOS",
    "OBSERVED_INJURY_LOCATION",
    "OBSERVED_INJURY_TYPE",
    "FallScenario",
    "SimulationSpec",
    "ScenarioSampler",
    "MASS_MODEL",
    "coverage_report",
]

_B25 = "[B25] Beurienne 2025 Front Sports Act Living 7:1609133 (n=245)"


# --------------------------------------------------------------------------
# 基本容器
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Category:
    """一个取值及其实测占比。"""

    label: str
    prob: float
    ci: tuple[float, float] | None = None      # 实测 95% CI（用于验证采样器）
    source: str = _B25

    @property
    def is_unknown(self) -> bool:
        return self.label.startswith("unknown")


@dataclass(frozen=True)
class Dimension:
    """一个场景维度（一组互斥取值的分布）。"""

    key: str
    name_cn: str
    cats: tuple[Category, ...]
    source: str = _B25
    note: str = ""

    def _kept(self, drop_unknown: bool) -> list[Category]:
        ks = [c for c in self.cats if not (drop_unknown and c.is_unknown)]
        if not ks:
            raise ValueError(f"维度 {self.key!r} 去掉 unknown 后为空")
        return ks

    def probs(self, drop_unknown: bool = True) -> tuple[list[str], np.ndarray]:
        """返回 ``(labels, probs)``，已按需要剔掉 unknown 并重新归一。"""
        ks = self._kept(drop_unknown)
        p = np.array([c.prob for c in ks], dtype=float)
        return [c.label for c in ks], p / p.sum()

    def raw_total(self) -> float:
        """原始占比之和（含 unknown），应约为 1 —— 用于校验抄录无误。"""
        return float(sum(c.prob for c in self.cats))

    def measured_ci(self, label: str) -> tuple[float, float] | None:
        for c in self.cats:
            if c.label == label:
                return c.ci
        return None


# --------------------------------------------------------------------------
# [B25] Table 3 —— 全体 n=245，边际分布
# --------------------------------------------------------------------------
DIMENSIONS: dict[str, Dimension] = {
    "wall_type": Dimension(
        key="wall_type", name_cn="墙型",
        cats=(
            Category("vertical", 0.45, (0.388, 0.512)),
            Category("steep", 0.29, (0.240, 0.354)),
            Category("slab", 0.12, (0.084, 0.165)),
            Category("roof_overhang", 0.07, (0.044, 0.108)),
            Category("corner", 0.05, (0.028, 0.084)),
            Category("unknown", 0.02, (0.008, 0.047)),
        ),
        note="垂直 45% + 倾斜 29% = 74% 的事故发生在最常见的两种墙上。",
    ),
    "height_band": Dimension(
        key="height_band", name_cn="坠落高度带",
        cats=(
            Category("top", 0.47, (0.412, 0.536)),
            Category("middle", 0.31, (0.260, 0.363)),   # 原文 CI 误印为与 bottom 相同
            Category("bottom", 0.20, (0.158, 0.259)),
            Category("unknown", 0.02, (0.002, 0.030)),
        ),
        note="墙顶坠落占比最高，且 [B25] 回归显示墙顶坠落的严重伤风险"
             "比墙底高约 5 倍（墙底 OR 0.18）。",
    ),
    "start_position": Dimension(
        key="start_position", name_cn="起始姿势",
        cats=(
            Category("vertical", 0.62, (0.554, 0.675)),
            Category("leaning_backward", 0.26, (0.210, 0.320)),
            Category("leaning_forward", 0.09, (0.063, 0.137)),
            Category("unknown", 0.03, (0.014, 0.058)),
        ),
        note="「上肢损伤者 41% 是后仰起始」出自 Table 4 的上肢列（31/76），"
             "**不能**与全体 26% 混用。",
    ),
    "rotation": Dimension(
        key="rotation", name_cn="坠落中的旋转",
        cats=(
            Category("longitudinal", 0.31, (0.259, 0.375)),
            Category("antero_posterior", 0.12, (0.086, 0.168)),
            Category("multiple", 0.09, (0.057, 0.127)),
            Category("back_transverse", 0.07, (0.044, 0.108)),
            Category("front_transverse", 0.03, (0.014, 0.058)),
            Category("without", 0.30, (0.248, 0.362)),
            Category("unknown", 0.08, (0.053, 0.123)),
        ),
        note="有旋转合计 62%。**无旋转只有 30%** —— 这个数是 1D 模型"
             "「可表达比例」的分子。",
    ),
    "landing_position": Dimension(
        key="landing_position", name_cn="落地姿势",
        cats=(
            Category("on_feet_leaning", 0.39, (0.333, 0.454)),
            Category("standing_on_feet", 0.34, (0.286, 0.404)),
            Category("head_first", 0.12, (0.084, 0.165)),
            Category("back_or_front_to_pad", 0.06, (0.037, 0.099)),
            Category("on_buttocks", 0.05, (0.031, 0.089)),
            Category("unknown", 0.04, (0.017, 0.063)),
        ),
        note="脚先落地合计 **73%**（原文摘要作 74%）。"
             "注意 Table 4 的下肢列把这两个数写成 44%/43%（=87%），是子群体。",
    ),
    "cause": Dimension(
        key="cause", name_cn="坠落诱因",
        cats=(
            Category("foot_slip", 0.29, (0.237, 0.350)),
            Category("dyno", 0.22, (0.173, 0.276)),
            Category("missed_hold", 0.21, (0.162, 0.263)),
            Category("route_finish", 0.06, (0.037, 0.099)),
            Category("balance_loss", 0.04, (0.022, 0.073)),
            Category("hands_let_go", 0.04, (0.019, 0.068)),
            Category("downclimbing_trouble", 0.03, (0.017, 0.063)),
            Category("other", 0.06, (0.037, 0.099)),
            Category("unknown", 0.05, (0.006, 0.041)),
        ),
        note="脚滑 + dyno + 抓空 = 72%，与原文一致。",
    ),
    "landing_surface": Dimension(
        key="landing_surface", name_cn="落点表面",
        cats=(
            Category("pads", 0.94, (0.911, 0.969)),
            Category("between_pads_and_wall", 0.03, (0.014, 0.058)),
            Category("between_2_pads", 0.02, (0.008, 0.047)),
            Category("bump_into_someone", 0.01, (0.001, 0.023)),
        ),
        note="94% 落在护垫上 → 那 6% 才走 ``on_pad=False``（硬地面）。",
    ),
}


# --------------------------------------------------------------------------
# 高度带 → 米（**假设**）
# --------------------------------------------------------------------------
#: 抱石墙典型高度。[B25] 原文：「short walls around 4.5 meters high」。
WALL_HEIGHT_M: float = 4.5

#: 高度带 → (下限, 上限) 米。**三档的米数是本项目假设**，
#: [B25] 只给了三档的占比（20/31/47%），没有给米数。
#: 划档依据：把 4.5 m 墙按三等分对应"墙底/墙中/墙顶"，
#: 下界 0.5 m 表示"离地很近但仍是坠落"。
HEIGHT_BANDS: dict[str, tuple[float, float]] = {
    "bottom": (0.5, 1.5),
    "middle": (1.5, 3.0),
    "top": (3.0, 4.5),
}


# --------------------------------------------------------------------------
# 落地姿势 → POSTURES（**假设**）
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class PostureMapping:
    """[B25] 的一个落地类别映射到 ``pad.POSTURES`` 的哪几种姿势。"""

    postures: tuple[tuple[str, float], ...]
    rationale: str
    ligament_risk: bool
    source: str = "本项目假设（[B25] 的 5 类落地 ≠ POSTURES 的 10 种姿势）"

    def __post_init__(self):
        tot = sum(w for _, w in self.postures)
        if abs(tot - 1.0) > 1e-9:
            raise ValueError(f"PostureMapping 权重和须为 1，实为 {tot}")


LANDING_TO_POSTURE: dict[str, PostureMapping] = {
    "on_feet_leaning": PostureMapping(
        postures=(("one-leg-awkward", 0.45), ("toe-point", 0.30),
                  ("feet-first-stiff", 0.25)),
        rationale="【B25】「在脚上且倾斜」= 身体倾斜着落在脚上 → 足部以外缘/"
                  "局部先着垫，是**踝旋后**的经典机制（原文明确写：护垫下沉"
                  "→ 踝旋后 → 韧带负荷）。映射成小接触面积姿势以**近似力的后果**，"
                  "但注意：本模型没有踝关节，**算不出旋后本身**。",
        ligament_risk=True,
    ),
    "standing_on_feet": PostureMapping(
        postures=(("feet-first-stiff", 0.60), ("controlled-drop", 0.40)),
        rationale="【B25】「站在脚上」= 双脚对称且直立着地。**是否屈膝缓冲"
                  "文献没给** —— 假设偏绷直，因为 85% 是**非自主**坠落"
                  "（来不及预判，难以主动屈曲）。这一条最需要实测数据替换。",
        ligament_risk=False,
    ),
    "head_first": PostureMapping(
        postures=(("head-first", 1.0),),
        rationale="1:1 对应。",
        ligament_risk=False,
    ),
    "back_or_front_to_pad": PostureMapping(
        postures=(("flat-flop", 1.0),),
        rationale="1:1 对应（背/面朝垫 = 平拍）。",
        ligament_risk=False,
    ),
    "on_buttocks": PostureMapping(
        postures=(("butt-impact", 1.0),),
        rationale="1:1 对应。",
        ligament_risk=False,
    ),
}

#: 落点表面 → 是否落在软垫上。垫与墙之间 / 两垫之间都是**坚硬缝隙**。
SURFACE_TO_ON_PAD: dict[str, bool] = {
    "pads": True,
    "between_pads_and_wall": False,
    "between_2_pads": False,
    "bump_into_someone": True,
}


# --------------------------------------------------------------------------
# 质量模型（**假设**）
# --------------------------------------------------------------------------
#: [B25] 人群：62% 女性、平均年龄 29.2±7.7。[B25] **没有报告体重**，
#: 故下面用一般人群人体测量学量级。**这是假设，不是 [B25] 实测。**
MASS_MODEL: dict[str, tuple[float, float, float]] = {
    # label: (prob, mean_kg, sd_kg)
    "female": (0.62, 62.0, 9.0),
    "male": (0.36, 78.0, 10.0),
    "unknown": (0.02, 70.0, 10.0),
}
MASS_CLIP_KG: tuple[float, float] = (40.0, 130.0)


# --------------------------------------------------------------------------
# [B25] Figure 4 —— 实测的联合场景（前 3 个有数字）
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class NamedScenario:
    """[B25] Figure 4 中一个实测的**联合**场景。"""

    rank: int
    share: float                 # 占全部场景的比例（实测）
    start_position: str
    rotation: str
    landing_position: str
    injury_location: str
    detail: str

    @property
    def rotation_free(self) -> bool:
        return self.rotation == "without"


#: 只有前 3 名在正文里有明确数字；其余 4 个在 Figure 4 图内，未能提取。
NAMED_SCENARIOS: tuple[NamedScenario, ...] = (
    NamedScenario(1, 0.17, "vertical", "without", "standing_on_feet",
                  "lower_limb",
                  "13 踝扭伤 / 7 膝肌腱韧带断裂 / 6 踝骨折"),
    NamedScenario(2, 0.14, "vertical", "longitudinal", "on_feet_leaning",
                  "lower_limb",
                  "10 踝扭伤 / 5 膝肌腱韧带断裂 / 4 踝骨折"),
    NamedScenario(3, 0.11, "leaning_backward", "longitudinal", "on_feet_leaning",
                  "upper_limb",
                  "7 踝扭伤 / 5 肘脱位 / 3 肘扭伤"),
)


# --------------------------------------------------------------------------
# [B25] Table 1 —— 实测的损伤结局（P1 的对标靶子）
# --------------------------------------------------------------------------
OBSERVED_INJURY_LOCATION: dict[str, tuple[float, tuple[float, float], int]] = {
    # label: (share, CI, n)
    "lower_limb": (0.67, (0.613, 0.719), 201),
    "upper_limb": (0.25, (0.207, 0.304), 76),
    "trunk": (0.04, (0.023, 0.068), 12),
    "head_neck": (0.03, (0.018, 0.060), 10),
    "unknown": (0.01, (0.002, 0.024), 2),
}

OBSERVED_INJURY_TYPE: dict[str, tuple[float, tuple[float, float], int]] = {
    "sprain": (0.36, (0.304, 0.411), 107),
    "fracture": (0.23, (0.182, 0.276), 68),
    "dislocation": (0.11, (0.082, 0.154), 34),
    "tendon_ligament_rupture": (0.11, (0.079, 0.150), 33),
    "tendonitis": (0.02, (0.011, 0.047), 7),
    "tbi": (0.01, (0.005, 0.034), 4),
    "bruises_soft_tissue": (0.03, (0.016, 0.056), 9),
    "cutaneous": (0.01, (0.005, 0.034), 4),
    "other": (0.07, (0.046, 0.104), 21),
    "unknown": (0.05, (0.028, 0.077), 14),
}

#: 最常受伤的**具体**部位（[B25] Table 1 正文）。
OBSERVED_BODY_PART_TOP: dict[str, float] = {
    "ankle": 0.40,
    "elbow": 0.16,
    "knee": 0.15,
}
#: 具体损伤（占比最大者）。
OBSERVED_SPECIFIC_INJURY: dict[str, float] = {
    "ankle_sprain": 0.28,
    "ankle_fracture": 0.08,
    "elbow_dislocation": 0.07,
}


# --------------------------------------------------------------------------
# 场景与仿真规格
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class FallScenario:
    """一次坠落场景（[B25] 的 5 个运动学维度 + 环境维度）。"""

    wall_type: str
    height_band: str
    start_position: str
    rotation: str
    landing_position: str
    cause: str
    landing_surface: str

    @property
    def on_pad(self) -> bool:
        return SURFACE_TO_ON_PAD[self.landing_surface]

    @property
    def rotation_free(self) -> bool:
        """本模型**唯一**有资格讨论的场景类型：坠落中无旋转。"""
        return self.rotation == "without"

    @property
    def foot_first(self) -> bool:
        return self.landing_position in ("standing_on_feet", "on_feet_leaning")

    @property
    def upper_limb_mechanism(self) -> bool:
        """后仰起始 → 旋转 → 落地不稳 → 本能用上肢撑垫（[B25] 讨论）。

        代理变量只用**起始姿势**：后仰占全体 26%，但在上肢损伤者中占 41%
        （Table 4 上肢列 31/76）→ 粗 OR ≈ 2。**这是单变量代理，不是回归结果。**
        """
        return self.start_position == "leaning_backward"

    def to_dict(self) -> dict[str, str]:
        return {
            "wall_type": self.wall_type,
            "height_band": self.height_band,
            "start_position": self.start_position,
            "rotation": self.rotation,
            "landing_position": self.landing_position,
            "cause": self.cause,
            "landing_surface": self.landing_surface,
        }


@dataclass(frozen=True)
class SimulationSpec:
    """一个场景落到"1D 模型能吃的输入"上的结果。"""

    scenario: FallScenario
    posture: str
    height_m: float
    mass_kg: float
    on_pad: bool
    #: 模型能否正确表达"旋转"。False ⇒ 载荷是**忽略旋转**的近似。
    rotation_representable: bool
    #: 该场景的损伤是否以**韧带/扭伤**为主 —— 模型产生不了这类损伤。
    ligament_dominant: bool
    #: 该场景是否属于"上肢撑垫"机制 —— 模型里没有上肢。
    upper_limb_mechanism: bool

    @property
    def weight_exact(self) -> bool:
        """模型给出的结论**可信**的充要条件（本模块的口径）。

        必须同时满足：无旋转可表达、且损伤不是韧带主导、且不涉及上肢。
        """
        return (self.rotation_representable
                and not self.ligament_dominant
                and not self.upper_limb_mechanism)


# --------------------------------------------------------------------------
# 采样器
# --------------------------------------------------------------------------
class ScenarioSampler:
    """按 [B25] Table 3 的实测占比采样坠落场景。

    Examples
    --------
    >>> s = ScenarioSampler(seed=0)
    >>> sc = s.sample()
    >>> spec = s.sample_simulation()
    >>> spec.posture in __import__("climbing.pad", fromlist=["POSTURES"]).POSTURES
    True
    """

    def __init__(self, seed: int | None = None, drop_unknown: bool = True):
        self.rng = np.random.default_rng(seed)
        self.drop_unknown = drop_unknown
        self._prepared: dict[str, tuple[list[str], np.ndarray]] = {}
        for k, d in DIMENSIONS.items():
            self._prepared[k] = d.probs(drop_unknown)
        # MASS_MODEL: label -> (prob, mean, sd)
        # 早期版本写成 zip(*[(v[0], v[1]) for v in ...])，把 **prob 当成了 label、
        # mean 当成了 prob** → rng.choice 拿到浮点数当标签，KeyError: '0.62'。
        self._mass_labels = list(MASS_MODEL.keys())
        self._mass_p = np.array([v[0] for v in MASS_MODEL.values()], dtype=float)
        self._mass_p = self._mass_p / self._mass_p.sum()

    # -- 内部 ---------------------------------------------------------
    def _pick(self, key: str) -> str:
        labels, p = self._prepared[key]
        return str(self.rng.choice(labels, p=p))

    def _mass(self) -> float:
        lab = str(self.rng.choice(self._mass_labels, p=self._mass_p))
        _, mu, sd = MASS_MODEL[lab]
        return float(np.clip(self.rng.normal(mu, sd), *MASS_CLIP_KG))

    # -- 对外 ---------------------------------------------------------
    def sample(self) -> FallScenario:
        return FallScenario(
            wall_type=self._pick("wall_type"),
            height_band=self._pick("height_band"),
            start_position=self._pick("start_position"),
            rotation=self._pick("rotation"),
            landing_position=self._pick("landing_position"),
            cause=self._pick("cause"),
            landing_surface=self._pick("landing_surface"),
        )

    def sample_many(self, n: int) -> list[FallScenario]:
        return [self.sample() for _ in range(n)]

    def height_m(self, band: str) -> float:
        lo, hi = HEIGHT_BANDS[band]
        return float(self.rng.uniform(lo, hi))

    def to_simulation(self, sc: FallScenario, mass_kg: float | None = None) -> SimulationSpec:
        """把一个场景映射成 ``simulate_boulder_fall`` 的参数。"""
        m = LANDING_TO_POSTURE[sc.landing_position]
        names = [n for n, _ in m.postures]
        weights = np.array([w for _, w in m.postures], dtype=float)
        posture = str(self.rng.choice(names, p=weights / weights.sum()))
        return SimulationSpec(
            scenario=sc,
            posture=posture,
            height_m=self.height_m(sc.height_band),
            mass_kg=self._mass() if mass_kg is None else float(mass_kg),
            on_pad=sc.on_pad,
            rotation_representable=sc.rotation_free,
            ligament_dominant=m.ligament_risk,
            upper_limb_mechanism=sc.upper_limb_mechanism,
        )

    def sample_simulation(self, mass_kg: float | None = None) -> SimulationSpec:
        return self.to_simulation(self.sample(), mass_kg=mass_kg)


# --------------------------------------------------------------------------
# 覆盖度：模型能代表真实世界的多少
# --------------------------------------------------------------------------
def _p(key: str, label: str) -> float:
    """某维度某取值的**归一化后**概率（剔掉 unknown）。"""
    labels, probs = DIMENSIONS[key].probs(drop_unknown=True)
    return float(probs[labels.index(label)])


def coverage_report() -> str:
    """给出"1D 模型覆盖真实场景/损伤的比例"，并逐条说明缺口。

    三种口径，从宽到严：

    1. **场景覆盖（宽）** —— 模型有对应落地姿势的比例。近乎全覆盖，
       因为 ``POSTURES`` 里脚/头/背/臀都有。
    2. **场景覆盖（严）** —— **无旋转**的比例。旋转是模型完全没有的自由度。
    3. **损伤覆盖（口径最硬）** —— 模型能产生**损伤类型 × 部位**的比例。
       模型只有骨折判据、只有 7 个骨部位。
    """
    p_no_rot = _p("rotation", "without")
    p_feet = _p("landing_position", "standing_on_feet") + \
        _p("landing_position", "on_feet_leaning")
    p_feet_no_rot = _p("landing_position", "standing_on_feet") * p_no_rot + \
        _p("landing_position", "on_feet_leaning") * p_no_rot   # 独立性假设

    # ---- 损伤类型覆盖 ----
    #  模型能产出：骨折（在下肢 7 个骨部位）。
    #  [B25]：骨折占 23%，其中 71% 在下肢（Table 2）→ 模型覆盖 23%×71%。
    frac = OBSERVED_INJURY_TYPE["fracture"][0]
    frac_ll = 0.71                                  # Table 2 骨折的部位分布
    covered_type = frac * frac_ll
    #  头/颈 3% 由 HIC 覆盖（另一条判据，不在 bone.py 里）
    covered_head = OBSERVED_INJURY_LOCATION["head_neck"][0]

    lines = [
        "=" * 74,
        "  覆盖度报告 —— 1D 模型 vs [B25] 真实世界",
        "=" * 74,
        "",
        "【口径 1】场景覆盖（宽）：模型有对应落地姿势的场景",
        f"  落地姿势映射覆盖 100% 的 [B25] 落地类别",
        f"  （脚上直立 / 脚上倾斜 / 头先 / 背 / 臀 → POSTURES 均存在）",
        "  ⚠️ 但「有姿势」不等于「算得对」 —— 见口径 2、3。",
        "",
        "【口径 2】场景覆盖（严）：坠落中**无旋转**的比例",
        f"  P(无旋转) = {p_no_rot:.3f}   ← 这是模型唯一没有丢掉自由度的场景",
        f"  P(脚先) = {p_feet:.3f}",
        f"  P(无旋转 ∧ 脚先) ≈ {p_feet_no_rot:.3f}  （独立性假设；"
        f"原文指出旋转与落地相关，故这是**上界**）",
        "",
        "【口径 3】损伤覆盖（最硬）：模型能产出的「损伤类型 × 部位」",
        f"  [B25] 骨折占比 {frac:.2f}，其中下肢占 {frac_ll:.0%}（Table 2）",
        f"  → 模型覆盖 {frac:.2f} × {frac_ll:.2f} = {covered_type:.3f}",
        f"  + 头/颈 {covered_head:.2f}（由 HIC 覆盖，另一条判据）",
        f"  → 合计约 **{covered_type + covered_head:.3f}**",
        "",
        "  模型**无法产出**的损伤类型（合计约 58%）：",
    ]
    for k in ("sprain", "dislocation", "tendon_ligament_rupture", "tendonitis"):
        sh = OBSERVED_INJURY_TYPE[k][0]
        lines.append(f"    ❌ {k:<26} {sh:.2f}")
    lines += [
        "",
        "  模型**无法定位**的常见部位：",
        "    ❌ 踝（40%，含扭伤 28% + 骨折 8%）—— 无踝关节自由度",
        "    ❌ 肘（16%）—— 无上肢",
        "    ❌ 膝（15%）—— 无膝关节韧带",
        "    ⚠️ 踝骨折（8%）**部分**可由 tibia_distal / fibula_ends 承担，",
        "       但模型无法把它辨识成「踝」（会报成胫骨远端/腓骨骨折）。",
        "",
        "  模型**能**定位的部位（都有骨折判据）：",
        "    ✅ 跟骨、胫骨远端、胫骨骨干、腓骨两端、股骨颈、骨盆、腰椎椎体",
        "",
        "【[B25] 自己给的 7 个场景】只有 1 个无旋转：",
    ]
    for ns in NAMED_SCENARIOS:
        mark = "✅ 可仿真" if ns.rotation_free else "❌ 有旋转"
        lines.append(f"    #{ns.rank}  {ns.share:.0%}  {mark}  "
                     f"{ns.start_position}/{ns.rotation}/{ns.landing_position}"
                     f"  → {ns.injury_location}")
    lines += [
        "",
        "  ⚠️ 那个唯一可仿真的 #1（17%）里，最主要的损伤是 **13 例踝扭伤**",
        "     —— 仍然落在模型算不出的类别里。",
        "",
        "=" * 74,
        f"  结论：模型能给出**可信**结论的坠落比例 ≲ {p_no_rot:.0%}（无旋转），",
        f"        能产出的损伤比例 ≈ {covered_type + covered_head:.0%}。",
        "        其余部分是**结构性缺口**，不是标定误差 —— 除非加自由度，",
        "        调参永远补不上。",
        "=" * 74,
    ]
    return "\n".join(lines)
