"""
visualize_b1.py — 抱石坠落 Phase B-1 可视化

提供三类图：
  1. plot_single_timeseries  ：单场景时程（速度 / 加速度 / 压缩 / 软垫反力）
  2. plot_height_comparison   ：多高度峰值力 / 加速度 / 压缩对比柱状图
  3. plot_pad_force_curve     ：软垫力-位移曲线 + 安全阈值标注
  4. plot_matrix_heatmap      ：高度×质量 峰值力热力图（可选 pad 预设）

所有图保存为 PNG 到 output/，并可选弹窗。
"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")  # 无界面保存
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

# 中文字体支持（避免 CJK 字形缺失警告）
for _cand in ("Microsoft YaHei", "SimHei", "Noto Sans SC", "Arial Unicode MS"):
    if any(_cand in f.name for f in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = [_cand, "DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        break

# 允许以脚本方式直接运行（把项目根加入 sys.path）
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.boulder_fall_1d import (
    CrashPad, PAD_PRESETS, SAFE_THRESHOLDS, SimConfig, SimResult,
    assess, energy_check, scan_heights, simulate,
)

OUT = _ROOT / "output"
OUT.mkdir(exist_ok=True)

# 风险配色
RISK_COLOR = {"SAFE": "#2ca02c", "MODERATE": "#ff7f0e",
              "HIGH": "#d62728", "CRITICAL": "#8b0000"}


# ---------------------------------------------------------------------------
# 1. 单场景时程图
# ---------------------------------------------------------------------------
def plot_single_timeseries(res: SimResult, show: bool = False, save: bool = True):
    cfg = res.cfg
    fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
    t = res.t * 1000.0  # ms
    x_cm = res.x * 100.0

    # 速度
    axes[0].plot(t, res.v, color="#1f77b4", lw=1.6)
    axes[0].axhline(0, color="grey", lw=0.6)
    axes[0].set_ylabel("速度 (m/s)")
    axes[0].set_title(f"抱石坠落时程  |  {cfg.mass:.0f}kg, 高度 {cfg.height:.1f}m,  "
                      f"撞击 {res.impact_velocity:.1f} m/s", fontsize=12)
    axes[0].grid(alpha=0.3)

    # 加速度（净减速度 = 力/m，向下为正，取反得向上减速为正冲击）
    a_net_g = res.pad_force / cfg.mass / 9.81
    axes[1].plot(t, a_net_g, color="#d62728", lw=1.6)
    axes[1].axhline(SAFE_THRESHOLDS["peak_accel_g"], color="orange", ls="--", lw=1,
                    label=f"阈值 {SAFE_THRESHOLDS['peak_accel_g']:.0f} g")
    axes[1].set_ylabel("减速度 (g)")
    axes[1].legend(loc="upper right", fontsize=9)
    axes[1].grid(alpha=0.3)

    # 压缩位移
    axes[2].plot(t, x_cm, color="#2ca02c", lw=1.6)
    axes[2].axhline(SAFE_THRESHOLDS["max_compression_cm"], color="orange", ls="--", lw=1,
                    label=f"行程上限 {SAFE_THRESHOLDS['max_compression_cm']:.0f} cm")
    axes[2].axhline(cfg.pad.max_compression * 100, color="red", ls=":", lw=1,
                    label=f"垫行程 {cfg.pad.max_compression*100:.0f} cm")
    axes[2].set_ylabel("压缩 (cm)")
    axes[2].legend(loc="upper right", fontsize=9)
    axes[2].grid(alpha=0.3)

    # 软垫反力
    axes[3].plot(t, res.pad_force / 1000.0, color="#9467bd", lw=1.6)
    axes[3].axhline(SAFE_THRESHOLDS["peak_force_kN"], color="orange", ls="--", lw=1,
                    label=f"安全上限 {SAFE_THRESHOLDS['peak_force_kN']:.0f} kN")
    axes[3].set_ylabel("软垫反力 (kN)")
    axes[3].set_xlabel("时间 (ms)")
    axes[3].legend(loc="upper right", fontsize=9)
    axes[3].grid(alpha=0.3)

    a = assess(res)
    fig.suptitle(f"风险等级: {a['risk']}  |  峰值力 {a['peak_force_kN']:.1f} kN,  "
                 f"峰值减速 {a['peak_accel_g']:.1f} g,  吸能 {a['energy_ratio']*100:.0f}%",
                 fontsize=13, color=RISK_COLOR[a["risk"]], fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    if save:
        fname = OUT / f"timeseries_m{cfg.mass:.0f}_h{cfg.height:.1f}.png"
        fig.savefig(fname, dpi=130)
        plt.close(fig)
        return fname
    if show:
        plt.show()
    return None


# ---------------------------------------------------------------------------
# 2. 多高度峰值对比柱状图
# ---------------------------------------------------------------------------
def plot_height_comparison(results, pad_name="standard", show=False, save=True):
    # 按高度排序
    results = sorted(results, key=lambda r: r.cfg.height)
    heights = [r.cfg.height for r in results]
    fpk = [assess(r)["peak_force_kN"] for r in results]
    apk = [assess(r)["peak_accel_g"] for r in results]
    xmax = [assess(r)["max_compression_cm"] for r in results]
    risks = [assess(r)["risk"] for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    bar_colors = [RISK_COLOR[r] for r in risks]

    axes[0].bar([str(h) for h in heights], fpk, color=bar_colors)
    axes[0].axhline(SAFE_THRESHOLDS["peak_force_kN"], color="orange", ls="--", lw=1.5,
                    label=f"安全 {SAFE_THRESHOLDS['peak_force_kN']:.0f} kN")
    axes[0].set_title("软垫峰值反力 vs 高度")
    axes[0].set_xlabel("坠落高度 (m)")
    axes[0].set_ylabel("峰值反力 (kN)")
    axes[0].legend()
    for i, v in enumerate(fpk):
        axes[0].text(i, v + 1, f"{v:.0f}", ha="center", fontsize=9)

    axes[1].bar([str(h) for h in heights], apk, color=bar_colors)
    axes[1].axhline(SAFE_THRESHOLDS["peak_accel_g"], color="orange", ls="--", lw=1.5,
                    label=f"安全 {SAFE_THRESHOLDS['peak_accel_g']:.0f} g")
    axes[1].set_title("峰值减速度 vs 高度")
    axes[1].set_xlabel("坠落高度 (m)")
    axes[1].set_ylabel("减速度 (g)")
    axes[1].legend()
    for i, v in enumerate(apk):
        axes[1].text(i, v + 1, f"{v:.0f}", ha="center", fontsize=9)

    axes[2].bar([str(h) for h in heights], xmax, color=bar_colors)
    axes[2].axhline(SAFE_THRESHOLDS["max_compression_cm"], color="orange", ls="--", lw=1.5,
                    label=f"上限 {SAFE_THRESHOLDS['max_compression_cm']:.0f} cm")
    axes[2].set_title("最大压缩位移 vs 高度")
    axes[2].set_xlabel("坠落高度 (m)")
    axes[2].set_ylabel("压缩 (cm)")
    axes[2].legend()
    for i, v in enumerate(xmax):
        axes[2].text(i, v + 0.5, f"{v:.0f}", ha="center", fontsize=9)

    fig.suptitle(f"抱石坠落峰值指标对比  ·  软垫预设: {pad_name}  ·  "
                 f"{results[0].cfg.mass:.0f} kg", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    if save:
        fname = OUT / f"height_comparison_{pad_name}_m{results[0].cfg.mass:.0f}.png"
        fig.savefig(fname, dpi=130)
        plt.close(fig)
        return fname
    if show:
        plt.show()
    return None


# ---------------------------------------------------------------------------
# 3. 软垫力-位移曲线
# ---------------------------------------------------------------------------
def plot_pad_force_curve(pad: CrashPad, impact_energy_J: float = None,
                         show=False, save=True, fname_tag="pad"):
    xs = np.linspace(0, pad.max_compression * 1.2, 400)
    fs = pad.force_at(xs)
    ks = pad.energy(xs)

    fig, ax1 = plt.subplots(figsize=(9, 6))
    ax1.plot(xs * 100, fs / 1000, color="#9467bd", lw=2, label="软垫反力 F(x)")
    ax1.axvline(pad.x_switch * 100, color="grey", ls=":", lw=1,
                label=f"压实点 {pad.x_switch*100:.0f} cm")
    ax1.axvline(pad.max_compression * 100, color="red", ls="--", lw=1.2,
                label=f"行程上限 {pad.max_compression*100:.0f} cm")
    ax1.set_xlabel("压缩位移 (cm)")
    ax1.set_ylabel("反力 (kN)", color="#9467bd")
    ax1.tick_params(axis="y", labelcolor="#9467bd")
    ax1.grid(alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(xs * 100, ks, color="#2ca02c", lw=1.6, ls="-", alpha=0.7,
             label="累计吸能 E(x)")
    ax2.set_ylabel("累计吸收能量 (J)", color="#2ca02c")
    ax2.tick_params(axis="y", labelcolor="#2ca02c")

    if impact_energy_J is not None:
        ax2.axhline(impact_energy_J, color="orange", ls="--", lw=1.2,
                    label=f"撞击动能 {impact_energy_J:.0f} J")
        # 求交点（能量=动能）
        idx = np.where(ks >= impact_energy_J)[0]
        if len(idx):
            x_at = xs[idx[0]] * 100
            ax2.axvline(x_at, color="orange", ls=":", lw=1)

    # 合并图例
    lines1, lab1 = ax1.get_legend_handles_labels()
    lines2, lab2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, lab1 + lab2, loc="upper left", fontsize=9)

    ax1.set_title("软垫力-位移与能量曲线（双线性泡沫模型）", fontsize=12, fontweight="bold")
    fig.tight_layout()
    if save:
        fname = OUT / f"pad_curve_{fname_tag}.png"
        fig.savefig(fname, dpi=130)
        plt.close(fig)
        return fname
    if show:
        plt.show()
    return None


# ---------------------------------------------------------------------------
# 4. 高度×质量 峰值力热力图
# ---------------------------------------------------------------------------
def plot_matrix_heatmap(results, pad_name="standard", show=False, save=True):
    masses = sorted({getattr(r, "mass", r.cfg.mass) for r in results})
    heights = sorted({r.cfg.height for r in results})
    grid = np.full((len(masses), len(heights)), np.nan)
    risk_grid = np.empty((len(masses), len(heights)), dtype=object)
    for r in results:
        mi = masses.index(getattr(r, "mass", r.cfg.mass))
        hi = heights.index(r.cfg.height)
        grid[mi, hi] = assess(r)["peak_force_kN"]
        risk_grid[mi, hi] = assess(r)["risk"]

    fig, ax = plt.subplots(figsize=(9, 6))
    im = ax.imshow(grid, cmap="RdYlGn_r", aspect="auto", vmin=0, vmax=60)
    ax.set_xticks(range(len(heights)))
    ax.set_xticklabels([f"{h:.0f}m" for h in heights])
    ax.set_yticks(range(len(masses)))
    ax.set_yticklabels([f"{m:.0f}kg" for m in masses])
    ax.set_xlabel("坠落高度")
    ax.set_ylabel("质量")
    for i in range(len(masses)):
        for j in range(len(heights)):
            v = grid[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.0f}\n{risk_grid[i,j][:4]}", ha="center", va="center",
                        fontsize=8, color="black")
    fig.colorbar(im, ax=ax, label="峰值反力 (kN)")
    ax.set_title(f"峰值反力热力图 · 软垫: {pad_name}\n(绿=安全<12kN, 红=高危险)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    if save:
        fname = OUT / f"heatmap_{pad_name}.png"
        fig.savefig(fname, dpi=130)
        plt.close(fig)
        return fname
    if show:
        plt.show()
    return None


if __name__ == "__main__":
    # 演示生成所有图
    pad = PAD_PRESETS["standard"]
    res = simulate(SimConfig(mass=80, height=3.0, pad=pad))
    f1 = plot_single_timeseries(res)
    f2 = plot_height_comparison(scan_heights(pad=pad))
    f3 = plot_pad_force_curve(pad, impact_energy_J=res.kinetic_energy, fname_tag="standard")
    print("saved:", f1, f2, f3)
