"""
boulder_fall_1d.py — 抱石坠落摔伤分析 · Phase B-1：1D 简化模型

物理模型
--------
把攀岩者当成**质点**，软垫当成**双线性非线性弹簧**（先软后硬，模拟泡沫压实）。
坐标系：x 从垫表面开始，向下为正（自由落体向下加速，x 增大表示陷入软垫）。

  - 自由落体段 (x <= 0)：m * dv/dt = -m*g   （重力向下，取向下为正 → dv/dt = +g）
  - 软垫压缩段 (x >  0)：m * dv/dt = m*g - F_pad(x)
      其中 F_pad 是软垫向上反力（阻碍陷落），故取负号。

动力学方程（向下为正）：
    dx/dt = v
    若 x <= 0:  dv/dt = g
    若 x  > 0:  dv/dt = g - F_pad(x)/m

双线性软垫力-位移：
    F_pad(x) = k1 * x                       ,  x <= x_switch
             = k1*x_switch + k2*(x-x_switch),  x >  x_switch
动能全部由软垫吸收时，峰值压缩 x_max 满足：
    m*g*x_max + 0.5*m*v_impact^2 = E_pad(x_max)
（上式用于能量守恒自检）

安全阈值（来自仿真方案文档，B 类抱石场景）：
    软垫峰值反力   < 12 kN
    峰值加速度     < 20 g
    最大压缩位移   < 0.15 m
    软垫吸能占比   > 90%
"""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp


# ---------------------------------------------------------------------------
# 软垫模型
# ---------------------------------------------------------------------------
@dataclass
class CrashPad:
    """双线性软垫（开放孔 PU + 封闭孔 PE 两层泡沫的等效模型）。

    Parameters
    ----------
    k1 : 初始（软）刚度 N/m，对应上层软泡沫
    k2 : 压实（硬）刚度 N/m，对应下层硬泡沫
    x_switch : 软硬转换位移 m（约 5 cm，泡沫开始压实）
    max_compression : 软垫总可用压缩行程 m（超过即压穿到底→地面）
    """

    k1: float = 50_000.0
    k2: float = 300_000.0
    x_switch: float = 0.05
    max_compression: float = 0.15
    k_ground: float = 5_000_000.0  # 压穿后的硬地面等效刚度 N/m

    def force(self, x: float) -> float:
        """软垫反力 (N)，向上为正（阻碍陷落）。x<=0 时返回 0。

        当 x > max_compression 时叠加硬地面刚度（模拟压穿到底）。
        """
        if x <= 0.0:
            return 0.0
        if x <= self.x_switch:
            f = self.k1 * x
        else:
            f = self.k1 * self.x_switch + self.k2 * (x - self.x_switch)
        if x > self.max_compression:
            f += self.k_ground * (x - self.max_compression)
        return f

    def energy(self, x) -> float:
        """软垫从 0 压缩到 x 所吸收的能量 (J) = ∫ F dx（含压穿后地面项）。

        支持标量或 numpy 数组输入。
        """
        x = np.asarray(x, dtype=float)
        e = np.zeros_like(x)
        soft = x <= self.x_switch
        e[soft] = 0.5 * self.k1 * x[soft] ** 2
        hard = (~soft) & (x > 0.0)
        xh = x[hard]
        e[hard] = (0.5 * self.k1 * self.x_switch**2
                   + self.k1 * self.x_switch * (xh - self.x_switch)
                   + 0.5 * self.k2 * (xh - self.x_switch) ** 2)
        bottom = x > self.max_compression
        if np.any(bottom):
            xb = x[bottom]
            e[bottom] += self.k_ground * (xb - self.max_compression) ** 2 / 2
        return e if e.shape != () else float(e)

    def force_at(self, x: np.ndarray) -> np.ndarray:
        return np.array([self.force(xi) for xi in np.atleast_1d(x)])


# 常见软垫预设（用于对比演示）
PAD_PRESETS = {
    # 薄垫：单层、压实快 → 高风险
    "thin": CrashPad(k1=60_000.0, k2=400_000.0, x_switch=0.04, max_compression=0.10),
    # 标准垫：开放孔+封闭孔双层，常见室内抱石馆配置
    "standard": CrashPad(k1=50_000.0, k2=300_000.0, x_switch=0.05, max_compression=0.15),
    # 厚垫：加厚多层，行程充足 → 显著降低峰值
    "thick": CrashPad(k1=25_000.0, k2=120_000.0, x_switch=0.08, max_compression=0.25),
}


# ---------------------------------------------------------------------------
# 仿真参数与结果
# ---------------------------------------------------------------------------
@dataclass
class SimConfig:
    mass: float = 80.0           # 攀岩者质量 kg（含装备）
    height: float = 3.0          # 坠落高度 m（从垫表面到起跳/脱手点）
    g: float = 9.81              # 重力加速度 m/s^2
    pad: CrashPad = field(default_factory=CrashPad)
    t_max: float = 2.0           # 仿真时长 s（自适应：足够覆盖自由落体+回弹）
    max_step: float = 0.0005     # 积分最大步长 s（软垫急停需要小步长）


@dataclass
class SimResult:
    t: np.ndarray
    x: np.ndarray                # 位移（垫面为 0，向下正）
    v: np.ndarray                # 速度（向下正）
    a: np.ndarray                # 加速度（向下正，含重力；向下为 +）
    pad_force: np.ndarray        # 软垫反力（向上为正 → 此处存绝对值）
    impact_velocity: float       # 撞击垫瞬间速度 m/s
    peak_force: float            # 软垫峰值反力 N
    peak_accel_g: float          # 峰值"净"加速度（= 力/m，不含重力）g 数
    peak_total_accel_g: float    # 峰值总加速度（含重力分量）g 数
    max_compression: float       # 最大压缩 m
    kinetic_energy: float        # 撞击动能 J
    energy_absorbed: float       # 软垫实际吸收能量 J
    energy_ratio: float          # 吸能占比（= 吸收/总机械能）
    bottomed_out: bool           # 是否压穿到底
    cfg: SimConfig = None


# ---------------------------------------------------------------------------
# 动力学
# ---------------------------------------------------------------------------
def _deriv(t, y, cfg: SimConfig):
    x, v = y
    pad = cfg.pad
    # 软垫反力（向上，阻碍陷落）；向下为正，故加速度减去反力项
    if x <= 0.0:
        dvdt = cfg.g
    else:
        f = pad.force(x)
        dvdt = cfg.g - f / cfg.mass
    return [v, dvdt]


def _pad_force_from_state(x, v, cfg: SimConfig):
    """根据状态反算软垫反力（用于时程记录）。"""
    if x <= 0.0:
        return 0.0
    return cfg.pad.force(x)


def simulate(cfg: SimConfig) -> SimResult:
    """运行一次 1D 抱石坠落仿真。"""
    # 初值：站在 height 处（相对于垫面向上为负），速度为 0
    # 但我们用 x 从垫面起算向下为正，所以起跳时 x = -height, v = 0
    y0 = [-cfg.height, 0.0]

    # 事件：到达垫面（x 由负穿 0）
    def hit_pad(t, y, *args):
        return y[0]
    hit_pad.terminal = False
    hit_pad.direction = 1

    sol = solve_ivp(
        _deriv,
        (0.0, cfg.t_max),
        y0,
        args=(cfg,),
        method="RK45",
        max_step=cfg.max_step,
        rtol=1e-9,
        atol=1e-11,
        events=hit_pad,
        dense_output=True,
    )

    t = sol.t
    x = sol.y[0]
    v = sol.y[1]

    # 加速度时程（直接由动力学重算，数值更稳）
    a = np.array([_deriv(0.0, [xi, vi], cfg)[1] for xi, vi in zip(x, v)])
    pad_force = np.array([_pad_force_from_state(xi, vi, cfg) for xi, vi in zip(x, v)])

    # 撞击速度：x 第一次过 0 时的速度
    cross = np.where(np.diff(np.sign(x)) > 0)[0]
    if len(cross) > 0:
        i = cross[0]
        impact_velocity = abs(v[i])
    else:
        # 没碰到垫（高度太小或数值问题）——取最低点速度
        impact_velocity = abs(v[np.argmin(x)])

    peak_force = float(np.max(pad_force)) if np.any(pad_force > 0) else 0.0
    peak_accel_net = peak_force / cfg.mass / cfg.g if peak_force > 0 else 0.0  # 净（力/m）
    # 总加速度（含重力）：向下为正，落垫时 a = g - F/m，但 F 向上抵消，
    # 人体实际承受的"减速度"大小 = F/m；加上重力后总 = g - F/m（向下），
    # 取绝对值更大的那个作为冲击 g 数（峰值冲击阶段 F>>mg，故 ≈ F/m）
    a_net = pad_force / cfg.mass / cfg.g
    peak_total_accel_g = float(np.max(np.abs(a_net)))

    max_compression = float(np.max(x)) if np.any(x > 0) else 0.0
    bottomed_out = max_compression >= cfg.pad.max_compression - 1e-6

    kinetic_energy = 0.5 * cfg.mass * impact_velocity**2
    # 软垫吸收能量 = 压缩到底/最低点时的软垫储能（含重力做功）
    # 总机械能守恒：E_total = KE_impact + mg*x_max；软垫吸收 = pad.energy(x_max)
    energy_absorbed = cfg.pad.energy(max_compression)
    # 总输入能量（撞击动能 + 压缩段重力做功）
    total_mech = kinetic_energy + cfg.mass * cfg.g * max_compression
    energy_ratio = energy_absorbed / total_mech if total_mech > 0 else 0.0

    return SimResult(
        t=t, x=x, v=v, a=a, pad_force=pad_force,
        impact_velocity=impact_velocity,
        peak_force=peak_force,
        peak_accel_g=peak_accel_net,
        peak_total_accel_g=peak_total_accel_g,
        max_compression=max_compression,
        kinetic_energy=kinetic_energy,
        energy_absorbed=energy_absorbed,
        energy_ratio=energy_ratio,
        bottomed_out=bottomed_out,
        cfg=cfg,
    )


# ---------------------------------------------------------------------------
# 能量守恒自检
# ---------------------------------------------------------------------------
def energy_check(cfg: SimConfig, res: SimResult) -> dict:
    """验证能量守恒：撞击动能 + 压缩段重力势能 = 软垫吸收 + 残余动能。

    最深处 x_max 速度≈0，故：
        KE_impact + m*g*x_max  ≈  pad.energy(x_max)
    返回误差占比与是否压穿。
    """
    ke_impact = res.kinetic_energy
    grav_work = cfg.mass * cfg.g * res.max_compression
    lhs = ke_impact + grav_work
    rhs = res.energy_absorbed
    residual = lhs - rhs
    rel_err = abs(residual) / lhs if lhs > 0 else 0.0
    return {
        "ke_impact_J": ke_impact,
        "grav_work_J": grav_work,
        "pad_absorbed_J": rhs,
        "residual_J": residual,
        "rel_err": rel_err,
        "bottomed_out": res.bottomed_out,
        "pass": (rel_err < 0.02) and (not res.bottomed_out),
    }


# ---------------------------------------------------------------------------
# 安全阈值与评估
# ---------------------------------------------------------------------------
SAFE_THRESHOLDS = {
    "peak_force_kN": 12.0,       # 软垫峰值反力安全上限 (kN)
    "peak_accel_g": 20.0,        # 峰值加速度安全上限 (g)
    "max_compression_cm": 15.0,  # 最大压缩位移安全上限 (cm)
    "energy_ratio_min": 0.90,    # 软垫吸能占比下限
}


def assess(res: SimResult) -> dict:
    """对照安全阈值给出各项 PASS/FAIL 与整体风险等级。"""
    peak_force_kN = res.peak_force / 1000.0
    peak_accel_g = res.peak_accel_g
    max_comp_cm = res.max_compression * 100.0
    er = res.energy_ratio

    checks = {
        "peak_force": peak_force_kN <= SAFE_THRESHOLDS["peak_force_kN"],
        "peak_accel": peak_accel_g <= SAFE_THRESHOLDS["peak_accel_g"],
        "max_compression": max_comp_cm <= SAFE_THRESHOLDS["max_compression_cm"],
        "energy_ratio": er >= SAFE_THRESHOLDS["energy_ratio_min"],
        "not_bottomed": not res.bottomed_out,
    }
    n_fail = sum(1 for v in checks.values() if not v)

    if res.bottomed_out:
        risk = "CRITICAL"   # 压穿到底 → 直接撞地面
    elif n_fail == 0:
        risk = "SAFE"
    elif n_fail <= 2:
        risk = "MODERATE"
    else:
        risk = "HIGH"
    return {"checks": checks, "n_fail": n_fail, "risk": risk,
            "peak_force_kN": peak_force_kN, "peak_accel_g": peak_accel_g,
            "max_compression_cm": max_comp_cm, "energy_ratio": er}


# ---------------------------------------------------------------------------
# 参数扫描
# ---------------------------------------------------------------------------
def scan_heights(heights=(2.0, 3.0, 4.0, 5.0), mass=80.0, pad=None, **kw) -> list[SimResult]:
    pad = pad or CrashPad()
    results = []
    for h in heights:
        cfg = SimConfig(mass=mass, height=h, pad=pad, **kw)
        results.append(simulate(cfg))
    return results


def scan_matrix(heights=(2.0, 3.0, 4.0, 5.0), masses=(60.0, 80.0, 100.0),
                pad_presets=("standard", "thick"), **kw) -> list[SimResult]:
    """高度 × 质量 × 软垫预设 的批量扫描。"""
    results = []
    for pad_name in pad_presets:
        pad = PAD_PRESETS[pad_name]
        for m in masses:
            for h in heights:
                cfg = SimConfig(mass=m, height=h, pad=pad, **kw)
                res = simulate(cfg)
                res.cfg = cfg
                setattr(res, "pad_name", pad_name)
                setattr(res, "mass", m)
                results.append(res)
    return results


if __name__ == "__main__":
    # 快速自检
    cfg = SimConfig()
    res = simulate(cfg)
    chk = energy_check(cfg, res)
    a = assess(res)
    print(f"撞击速度: {res.impact_velocity:.2f} m/s")
    print(f"撞击动能: {res.kinetic_energy:.0f} J")
    print(f"峰值力:   {a['peak_force_kN']:.1f} kN")
    print(f"峰值加速度(净): {a['peak_accel_g']:.1f} g")
    print(f"最大压缩: {a['max_compression_cm']:.1f} cm")
    print(f"吸能占比: {a['energy_ratio']*100:.1f}%")
    print(f"压穿?: {res.bottomed_out}")
    print(f"风险等级: {a['risk']}")
    print(f"能量自检: rel_err={chk['rel_err']*100:.3f}%, pass={chk['pass']}")
    print("\n--- 软垫预设对比 (80kg) ---")
    for name in ("thin", "standard", "thick"):
        pad = PAD_PRESETS[name]
        r = simulate(SimConfig(mass=80, height=3.0, pad=pad))
        aa = assess(r)
        print(f"{name:9s}: Fpk={aa['peak_force_kN']:.1f}kN  apk={aa['peak_accel_g']:.1f}g  "
              f"xmax={aa['max_compression_cm']:.1f}cm  risk={aa['risk']}")
