"""当前模型（`climbing.pad` / `climbing.rope`）的图表输出。

产出到 `results/figures/`：

    b2_timeseries_<posture>.png   单场景时程（力/压力/加速度/屈曲/能量）
    b2_posture_compare.png        姿势对比 —— 峰值力 vs 峰值压力
    b2_pressure_matrix.png        姿势 x 高度 的峰值接触压力热力图
    b2_pad_constitutive.png       软垫本构曲线 + 标定依据
    rope_constitutive.png         绳本构曲线 + UIAA 标定点
    rope_ff_scan.png              峰值力/伸长率 对 FF 的扫描
    rope_device_compare.png       保护器对比

⚠️ 所有阈值/参考线都是**文献量级**，不是临床判据，图上已标注。
旧的质点模型图表在 `legacy/output/`，已废弃。
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from climbing.rope import RopeModel, simulate_rope_fall, default_rope
from climbing.belay import DEVICES
from climbing.pad import (
    POSTURES, CrashPad, hard_surface, simulate_boulder_fall,
)
from climbing.metrics import hic_from_result

OUT = pathlib.Path(__file__).resolve().parents[1] / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# 中文字体（Windows 常见）
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 110
plt.rcParams["savefig.bbox"] = "tight"

GREEN, RED, BLUE, ORANGE, GREY = "#2e7d32", "#c62828", "#1565c0", "#ef6c00", "#616161"


def sim(posture, height=3.0, mass=80.0, pad=None, on_pad=True, fast=False):
    kw = {}
    if fast:
        # 扫描/热力图只需要峰值压力的大致量级，放宽积分精度换速度。
        # ��要 ~4x 加速，峰值压力的偏差在 1-2% 以内，不影响分档排序。
        kw = dict(rtol=1e-6, atol=1e-9, max_ode_step=4e-4, dt_sample=1e-3)
    b = simulate_boulder_fall(height_m=height, mass_kg=mass, posture=posture,
                              pad=pad, on_pad=on_pad, **kw)
    b.hic = hic_from_result(b)
    return b


def peak_pressure(b):
    return float(np.max(b.pad_force_n / np.maximum(b.contact_area_m2, 1e-9))) / 1e3


# ==========================================================================
# 1) 单场景时程
# ==========================================================================
def fig_timeseries(posture="controlled-drop", height=3.0, mass=80.0):
    b = sim(posture, height, mass)
    post = b.meta["posture_obj"]
    # 只截"第一次接触 -> 第一次分离"这一段。冲击过程只有几十毫秒，
    # 铺满 2 s 全程的话有用的细节全被压成一根竖线。
    comp = b.pad_compression_m > 0
    idx = np.where(comp)[0]
    if idx.size:
        i0 = idx[0]
        gap = np.where(~comp[i0:])[0]
        i1 = (i0 + gap[0]) if gap.size else min(b.t_s.size, i0 + 2000)
        pad_t = 0.004
        m = (b.t_s >= b.t_s[i0] - pad_t) & (b.t_s <= b.t_s[min(i1, b.t_s.size - 1)] + pad_t)
    else:
        m = np.ones_like(b.t_s, dtype=bool)
    t = b.t_s[m]
    press = (b.pad_force_n / np.maximum(b.contact_area_m2, 1e-9) / 1e3)[m]

    fig, ax = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    a = ax[0]
    a.plot(t, b.z_foot_m[m], label="足 (第一接触点)", color=BLUE)
    a.plot(t, b.z_torso_m[m], label="躯干 CoM", color=ORANGE)
    a.set_ylabel("下落位移 (m)")
    a.legend(fontsize=8)
    a.set_title(f"{post.name} · {height} m / {mass:.0f} kg / 20 cm 软垫")

    a = ax[1]
    a.plot(t, b.pad_force_n[m] / 1e3, label="垫反力 (kN)", color=GREEN)
    a.plot(t, (b.accel_leg_g[m] * mass * 9.80665 / 1e3), label="下段 m·a (kN)",
           color=BLUE, ls="--", lw=1)
    a.set_ylabel("力 (kN)")
    a.legend(fontsize=8)

    a = ax[2]
    a.plot(t, press, color=RED, label="峰值接触压力")
    a.axhline(400, ls=":", c=RED, lw=1)
    a.text(t[0], 405, "400 kPa（分档上界，量级）", fontsize=7, color=RED)
    a.set_ylabel("压力 (kPa)")
    a.legend(fontsize=8)

    a = ax[3]
    a.plot(t, b.accel_leg_g[m], label="下段/接触部位 (g)", color=BLUE)
    a.plot(t, b.accel_torso_g[m], label="躯干 (g)", color=ORANGE)
    a.plot(t, b.primary_accel_g[m], color=RED, lw=2, label="受伤部位 (HIC 用)")
    a.set_ylabel("加速度 (g)")
    a.set_xlabel("时间 (s)")
    a.legend(fontsize=8)

    fig.suptitle(
        f"峰值 {b.peak_force_kn:.1f} kN · 压力 {peak_pressure(b):.0f} kPa · "
        f"接触 {b.peak_primary_g:.0f} g · HIC {b.hic:.0f}",
        fontsize=11, y=0.995)
    fig.text(0.5, 0.005,
             "⚠️ 接触段力/加速度曲线有明显数值振荡（刚性接触 + RK45），"
             "峰值可用，波形细节待 Phase B-3 改善",
             ha="center", fontsize=7.5, color=RED)
    p = OUT / f"b2_timeseries_{posture}.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  {p.name}")


# ==========================================================================
# 2) 姿势对比 —— 峰值力 vs 峰值压力
# ==========================================================================
def fig_posture_compare(height=3.0, mass=80.0):
    names = list(POSTURES)
    F, P, A, H, CT = [], [], [], [], []
    for n in names:
        b = sim(n, height, mass)
        F.append(b.peak_force_kn)
        P.append(peak_pressure(b))
        A.append(b.peak_primary_g)
        H.append(b.hic)
        CT.append(POSTURES[n].first_contact)
    order = np.argsort(P)[::-1]          # 按压力降序，最危险的在上
    names = [names[i] for i in order]
    F, P, A, H, CT = [np.array(v)[order] for v in (F, P, A, H, CT)]
    y = np.arange(len(names))
    labels = [f"{POSTURES[n].name}\n({CT[i]})" for i, n in enumerate(names)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7), sharey=True)

    ax1.barh(y, F, color=BLUE, alpha=.85)
    ax1.set_yticks(y, labels, fontsize=8)
    ax1.invert_yaxis()
    ax1.set_xlabel("峰值垫反力 (kN)")
    ax1.set_title("按峰值力排序", fontsize=11)
    for i, v in enumerate(F):
        ax1.text(v + 0.5, i, f"{v:.1f}", va="center", fontsize=8)

    cols = [RED if p >= 400 else (ORANGE if p >= 150 else GREEN) for p in P]
    ax2.barh(y, P, color=cols, alpha=.85)
    ax2.set_xlabel("峰值接触压力 (kPa)")
    ax2.set_title("按峰值压力排序（伤害的真实驱动量）", fontsize=11)
    for i, v in enumerate(P):
        ax2.text(v + 6, i, f"{v:.0f}", va="center", fontsize=8)
    for x, lab, c in [(150, "中", ORANGE), (400, "高", RED)]:
        ax2.axvline(x, ls=":", c=c, lw=1)
        ax2.text(x + 4, len(names) - 0.3, f"{lab}档 {x}", fontsize=7, color=c)

    fig.suptitle(
        f"{height} m / {mass:.0f} kg / 20 cm 软垫 —— 注意两个排序完全不同：\n"
        "峰值力最小的姿势未必压力最小（接触面积才是关键）",
        fontsize=12)
    p = OUT / "b2_posture_compare.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  {p.name}")


# ==========================================================================
# 3) 姿势 x 高度 的峰值压力热力图
# ==========================================================================
def fig_pressure_matrix(mass=80.0, on_pad=True):
    names = list(POSTURES)
    heights = np.arange(0.5, 5.01, 0.5)
    M = np.zeros((len(names), len(heights)))
    t0 = __import__("time").time()
    for i, n in enumerate(names):
        for j, h in enumerate(heights):
            M[i, j] = peak_pressure(sim(n, float(h), mass, on_pad=on_pad, fast=True))
        done = (i + 1) * len(heights)
        total = len(names) * len(heights)
        print(f"    [{done:3d}/{total}] {POSTURES[n].name:<20s} "
              f"({__import__('time').time()-t0:5.1f}s)", flush=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.pcolormesh(heights, np.arange(len(names)), M, cmap="YlOrRd",
                       norm=LogNorm(vmin=max(M.min(), 10), vmax=M.max()))
    ax.set_yticks(np.arange(len(names)),
                  [POSTURES[n].name for n in names], fontsize=8)
    ax.set_xlabel("落差 (m)")
    ax.set_ylabel("落地姿势")
    ax.set_title(f"峰值接触压力 (kPa) · {mass:.0f} kg · "
                 f"{'20 cm 软垫' if on_pad else '硬地面'}", fontsize=12)
    fig.colorbar(im, ax=ax, label="峰值接触压力 (kPa, 对数色标)")
    for i in range(len(names)):
        for j in range(0, len(heights), 2):
            ax.text(heights[j], i, f"{M[i,j]:.0f}", ha="center", va="center",
                    fontsize=6, color="black")
    p = OUT / ("b2_pressure_matrix.png" if on_pad else "b2_pressure_matrix_hard.png")
    fig.savefig(p)
    plt.close(fig)
    print(f"  {p.name}")


# ==========================================================================
# 4) 软垫本构
# ==========================================================================
def fig_pad_constitutive():
    xi = np.linspace(0, 0.95, 400)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    pads = {"soft/10cm": CrashPad(thickness_m=0.10),
            "standard/20cm": CrashPad(),
            "thick/30cm": CrashPad(thickness_m=0.30)}
    for lbl, pd in pads.items():
        ax1.plot(100 * xi, pd.stress_pa(xi) / 1e3, label=lbl)
    # 30 kg/m3 开放孔 PU 的公开量级区间
    ref = [(10, 10, 20), (30, 30, 45), (50, 60, 120), (70, 150, 300)]
    ax1.fill_between([r[0] for r in ref], [r[1] for r in ref], [r[2] for r in ref],
                     color=GREY, alpha=.2, label="PU 泡沫文献量级")
    ax1.axvline(100 * CrashPad().xi_dense_max, ls="--", c=RED, lw=1)
    ax1.text(85, 250, "压穿点 85%", fontsize=8, color=RED)
    ax1.set_xlabel("压缩应变 ξ (%)")
    ax1.set_ylabel("应力 (kPa)")
    ax1.set_title("软垫本构（重标后对上开放孔 PU 量级）", fontsize=11)
    ax1.legend(fontsize=8)
    ax1.grid(alpha=.3)

    # 吸能能力
    for lbl, pd in pads.items():
        ax2.plot(100 * xi, pd.capacity_j(1.0) * 0 + np.array(
            [pd.stress_pa(x) * pd.thickness_m for x in xi]) / 1e3, label=lbl)
    ax2.axvline(100 * CrashPad().xi_dense_max, ls="--", c=RED, lw=1)
    ax2.set_xlabel("压缩应变 ξ (%)")
    ax2.set_ylabel("单位面积吸能 (kJ/m²)")
    ax2.set_title("压实前的能量吸收能力（1 m² 接触）", fontsize=11)
    ax2.legend(fontsize=8)
    ax2.grid(alpha=.3)

    fig.suptitle("⚠️ 泡沫参数为文献量级，非本项目实测标定", fontsize=10, color=RED)
    p = OUT / "b2_pad_constitutive.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  {p.name}")


# ==========================================================================
# 5) 绳索
# ==========================================================================
def fig_rope_constitutive():
    rope = default_rope()
    eps = np.linspace(0, 0.5, 300)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    ax1.plot(100 * eps, rope.tension(eps) / 1e3, color=BLUE, lw=2)
    ax1.axhline(12.0, ls="--", c=RED, lw=1.2)
    ax1.text(1, 12.3, "UIAA 101 判据 12 kN", color=RED, fontsize=8)
    r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.77, t_hold_n=1e9)
    ax1.plot([r.max_elongation_pct], [r.peak_tension_kn], "o", color=GREEN, ms=9,
             zorder=5, label=f"UIAA 条件实测点\n{r.peak_tension_kn:.2f} kN @ {r.max_elongation_pct:.1f}%")
    ax1.fill_between([0, 30], 0, 9, color=GREEN, alpha=.08)
    ax1.set_xlabel("绳伸长率 ε (%)")
    ax1.set_ylabel("张力 (kN)")
    ax1.set_title("绳本构 + UIAA 101 标定", fontsize=11)
    ax1.legend(fontsize=8)
    ax1.grid(alpha=.3)

    scan = []
    for ff in np.arange(0.25, 2.01, 0.05):
        rr = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=float(ff),
                                t_hold_n=1e9, t_impact_max=1.2)
        scan.append((ff, rr.peak_tension_kn, rr.max_elongation_pct))
    ff = np.array([s[0] for s in scan])
    kN = np.array([s[1] for s in scan])
    pc = np.array([s[2] for s in scan])

    ax2.plot(ff, kN, color=BLUE, lw=2, label="峰值力")
    ax2.set_xlabel("坠落系数 FF")
    ax2.set_ylabel("峰值力 (kN)", color=BLUE)
    ax2.tick_params(axis="y", labelcolor=BLUE)
    ax2.axhline(12.0, ls="--", c=RED, lw=1.2)
    ax2.axvline(1.77, ls=":", c=GREY, lw=1)
    ax2.text(1.79, 1.0, "UIAA\nFF=1.77", fontsize=8, color=GREY)
    ax3 = ax2.twinx()
    ax3.plot(ff, pc, color=ORANGE, lw=2, label="伸长率")
    ax3.set_ylabel("伸长率 (%)", color=ORANGE)
    ax3.tick_params(axis="y", labelcolor=ORANGE)
    ax2.set_title("对 FF 的扫描（峰值力严格单调；FF≤1 无自由落体段故恒定）", fontsize=11)
    ax2.grid(alpha=.3)

    fig.suptitle("⚠️ 幂律本构无法同时拟合峰值力与伸长率，见 README「已知局限」",
                 fontsize=10, color=RED)
    p = OUT / "rope_constitutive.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  {p.name}")


def fig_rope_devices():
    names, peaks, pay, slide = [], [], [], []
    for n, d in DEVICES.items():
        r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.77,
                               device=n, t_hold_n=d.t_hold_n)
        names.append(d.name)
        peaks.append(r.peak_tension_kn)
        pay.append(r.max_payout_m)
        slide.append(r.total_rope_slide_m)
    x = np.arange(len(names))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    cols = [GREEN if p == 8.17 else ORANGE for p in peaks]
    ax1.bar(x, peaks, color=cols, alpha=.85)
    ax1.axhline(12.0, ls="--", c=RED, lw=1.2)
    ax1.text(-0.4, 12.4, "UIAA 12 kN", color=RED, fontsize=8)
    ax1.set_xticks(x, names, rotation=20, ha="right", fontsize=8)
    ax1.set_ylabel("坠落者侧峰值力 (kN)")
    ax1.set_title("保护器峰值力（= min(绳峰值, 保持力)）", fontsize=11)
    for i, v in enumerate(peaks):
        ax1.text(i, v + 0.3, f"{v:.2f}", ha="center", fontsize=8)

    w = 0.38
    ax2.bar(x - w/2, pay, w, label="放绳量 (m)", color=ORANGE, alpha=.85)
    ax2.bar(x + w/2, slide, w, label="绳总滑动 (m)", color=BLUE, alpha=.85)
    ax2.set_xticks(x, names, rotation=20, ha="right", fontsize=8)
    ax2.set_ylabel("长度 (m)")
    ax2.set_title("放绳量与绳磨损（Grigri 峰值低但放绳最多）", fontsize=11)
    ax2.legend(fontsize=8)

    p = OUT / "rope_device_compare.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  {p.name}")


def main():
    print(f"输出目录: {OUT}")
    fig_pad_constitutive()
    fig_timeseries("controlled-drop")
    fig_timeseries("toe-point")
    fig_timeseries("flat-flop")
    fig_posture_compare()
    fig_pressure_matrix(on_pad=True)
    fig_pressure_matrix(on_pad=False)
    fig_rope_constitutive()
    fig_rope_devices()
    print("完成")


if __name__ == "__main__":
    main()
