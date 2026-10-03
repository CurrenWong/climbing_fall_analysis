"""生成仿真流程图 -> results/figures/flowchart.png

用 matplotlib 画，避免依赖 mermaid 渲染器（离线可看）。
"""
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

OUT = pathlib.Path(__file__).resolve().parents[1] / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

C_IN = "#e3f2fd"; C_BOX = "#fff8e1"; C_POST = "#e8f5e9"
C_METRIC = "#f3e5f5"; C_VERDICT = "#fce4ec"; C_CAL = "#fff3e0"
C_EDGE = "#455a64"

fig, ax = plt.subplots(figsize=(15, 21))
ax.set_xlim(0, 10); ax.set_ylim(0, 21); ax.axis("off")


def box(x, y, w, h, text, color, fs=9.5, bold=False):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.12,rounding_size=0.18",
        facecolor=color, edgecolor=C_EDGE, linewidth=1.3, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, zorder=3, linespacing=1.5,
            fontweight="bold" if bold else "normal")


def arrow(x1, y1, x2, y2, label="", ls="-", color=C_EDGE, lw=1.4,
          rad=0.0, lx=None, ly=None):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=15,
        linestyle=ls, color=color, linewidth=lw, zorder=1,
        connectionstyle=f"arc3,rad={rad}"))
    if label:
        ax.text(lx if lx is not None else (x1 + x2) / 2,
                ly if ly is not None else (y1 + y2) / 2, label,
                ha="center", va="center", fontsize=8.5, color=color,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none"))


# ---- 1. 输入 ------------------------------------------------------------
box(2.0, 19.6, 6.0, 0.95,
    "① 输入参数\n"
    "绳索: 质量 / 绳长 L0 / 坠落系数 FF / 保护器库\n"
    "抱石: 落差 h / 质量 / 垫参数 / 姿势库 POSTURES",
    C_IN, bold=True)

# ---- 2. 建模：两条通道 --------------------------------------------------
box(0.3, 16.6, 4.5, 2.5,
    "②a 绳索建模  (rope.py)\n\n"
    "状态  y = [x, v, L_free]\n"
    "几何约束  e = x − h_free + (L_free − L0)\n"
    "绳本构  T(ε) = T_ref·(ε/ε_ref)ⁿ\n"
    "保护器两态: 锁死 / 滑动(限幅 t_hold)\n"
    "阻尼  c·ė（仅拉伸时耗能）\n"
    "dv/dt = g − T/m",
    C_BOX)

box(5.2, 16.6, 4.5, 2.5,
    "②b 抱石建模  (pad.py)\n\n"
    "状态  y = [z_足, v_足, z_躯, v_躯]\n"
    "双质点  m_low / m_up（split_mass 守卫）\n"
    "屈曲  xf = L0 − (z_足 − z_躯)\n"
    "垫本构  σ(ξ) 平台+压实+有限压穿支撑\n"
    "渐进接触  A(κ),  κ = comp / collapse_ref\n"
    "牛顿第三定律（见右框）",
    C_BOX)

box(0.35, 14.35, 4.4, 1.9,
    "垫反力与力的配平\n\n"
    "az_f = g − f_pad/m_low + f_flex/m_low\n"
    "az_t = g − f_flex/m_up\n\n"
    "两式相加 = f_pad/m − g   质心方程",
    "#fffde7", fs=8.5)

box(5.25, 14.35, 4.4, 1.9,
    "冲击窗口定义（两个模型共用）\n\n"
    "首次接触  t0 → 首次分离  t1\n\n"
    "峰值只在窗口内取 max，\n"
    "避免回弹重绷 / 二次冲击污染",
    "#fffde7", fs=8.5)

arrow(5.0, 19.6, 2.6, 19.15)
arrow(5.0, 19.6, 7.4, 19.15)
arrow(2.55, 16.6, 2.55, 16.25)
arrow(7.45, 16.6, 7.45, 16.25)

# ---- 3. 求解 ------------------------------------------------------------
box(0.3, 12.3, 9.4, 1.5,
    "③ 数值求解  solve_ivp\n\n"
    "绳索: RK45,  max_step 2e−4,  rtol 1e−8\n"
    "抱石: RK45,  max_step 1e−4,  rtol 1e−8,  atol 1e−11\n"
    "（B-3 链模型改用隐式 Radau —— 刚性接触）",
    C_POST, bold=True)
arrow(5.0, 19.6, 5.0, 13.8, label="组装 ODE", ly=13.95)

# ---- 4. 后处理 ----------------------------------------------------------
box(0.3, 9.0, 9.4, 2.9,
    "④ 后处理  窗口内统一处理\n\n"
    "峰值    max(力) / max(加速度) / max(接触压力 p = F/A)  —— 窗口内，非全程\n"
    "能量账  ∫F·dx 用**功率累积的峰值**（回弹段会把能量减回来，\n"
    "        np.trapezoid(f, x) 对弹性体净积分为 0 → 算出负数）\n"
    "        E_pad 与 E_flex 取**同一时刻**快照，相加才不超过输入\n"
    "单位    accel 统一除 g 归一（字段名以 _g 结尾就必须真是 g）",
    C_POST)
arrow(5.0, 12.3, 5.0, 11.9)

# ---- 5. 指标 ------------------------------------------------------------
box(0.3, 6.6, 4.5, 1.9,
    "⑤a 绳索指标\n\n"
    "峰值张力 / 加速度 / 伸长率\n"
    "放绳量 / 绳总滑动 / 滑动态时长\n"
    "保护者侧张力 = e^(μ·θ) 放大",
    C_METRIC)

box(5.2, 6.6, 4.5, 1.9,
    "⑤b 抱石指标\n\n"
    "峰值垫反力 / **接触压力 F/A**\n"
    "各段加速度（第一接触部位）\n"
    "屈曲行程 / 能量分配 / 残余",
    C_METRIC)
arrow(2.55, 9.0, 2.55, 8.5)
arrow(7.45, 9.0, 7.45, 8.5)

# ---- 6. 判读 ------------------------------------------------------------
box(0.3, 4.2, 4.5, 1.9,
    "⑥a 绳索判读\n\n"
    "UIAA 101: 峰值 < 12 kN\n"
    "峰值对 FF(>1) 严格单调\n"
    "保护器峰值 = min(绳峰值, t_hold)",
    C_VERDICT)

box(5.2, 4.2, 4.5, 1.9,
    "⑥b 抱石判读（分部位）\n\n"
    "接触部位  压力 150/400/800 kPa\n"
    "头/颈     HIC 240/700/1000\n"
    "脊柱/下肢 骨载荷 6/10/16 kN = m·a",
    C_VERDICT)
arrow(2.55, 6.6, 2.55, 6.1)
arrow(7.45, 6.6, 7.45, 6.1)

# ---- 7. 标定 & 验证 -----------------------------------------------------
box(0.3, 1.5, 4.5, 2.2,
    "⑦ 标定  (calibrate_rope.py)\n\n"
    "UIAA 条件网格搜索\n"
    "目标 = 加权 RSS(峰值 0.7 / 伸长 0.3)\n"
    "※ 相对误差**相加**会正负抵消成假优解\n"
    "※ 标定完必须**写回默认值**并检查\n"
    "  下游功能是否“活过来”",
    C_CAL, fs=8.5)

box(5.2, 1.5, 4.5, 2.2,
    "⑧ 验证  (tests/ 60 用例)\n\n"
    "正交不变量：牛顿第三定律 / 能量守恒 /\n"
    "单调性 / 已知解(HIC 闭式) / 单位 /\n"
    "反函数自洽 / 质量分配守卫\n"
    "※ 约一半用例是锁死已修的静默 bug",
    C_CAL, fs=8.5)

arrow(2.55, 4.2, 2.55, 3.7, ls="--", label="标定")
arrow(7.45, 4.2, 7.45, 3.7, ls="--", label="验证")
arrow(4.8, 2.6, 5.2, 2.6, ls=":", rad=0, lw=1.2)
ax.text(5.0, 2.15, "闭环", ha="center", fontsize=8, color=C_EDGE)

ax.text(5.0, 20.75, "攀岩坠落仿真 · 完整流程", ha="center", fontsize=15,
        fontweight="bold")

p = OUT / "flowchart.png"
fig.savefig(p, dpi=115, bbox_inches="tight")
print("saved:", p)
