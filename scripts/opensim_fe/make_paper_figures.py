# -*- coding: utf-8 -*-
r"""make_paper_figures.py — 论文图 1–9 出版级绘图（additive，只读复用既有 JSON 产物).

本脚本**只读** results/opensim_fe/ 下的既有机器可读产物，生成 results/opensim_fe/figures/
下的 8 张 PNG（>=200 dpi）+ FIGURE_CAPTIONS.md。不调用 FEBio、不重跑 OpenSim、不修改任何既有
文件/模块。

图 1（管线总览）+ 图 9（距骨骨折机制示意）为纯文献示意图，无 JSON 数据依赖；
图 2–7 走 JSON 数据。
  fig1  本脚本内置（管线骨架, 见 §2.1 / §2.5）
  fig2  results/opensim_fe/axial_subregion.json        (route_a_axial_ends + paper_facts)
  fig3  results/opensim_fe/nonvertical_s2.json         (joint_load_factors_5m + per_bone_risk_5m)
  fig4  results/opensim_fe/nonvertical_s5.json         (scenario_hard)
  fig5  results/opensim_fe/pad_supination_tradeoff.json(sweep.cells + sweep.threshold_by_height)
  fig6  results/opensim_fe/ankle_ligament_criterion.json(scenario.rows + critical_thresholds)
  fig7  results/opensim_fe/ankle_external_rotation.json(scenario_hard；图 4 的外旋镜像 / SER)
  fig9  本脚本内置（距骨骨折机制矢状面示意，论文 §3.6）
  fig8  本脚本内置（三踝骨折链式机制示意图, 见 §3.6 / §4.2 — SER 四阶段链式）

中文字体: 优先 Microsoft YaHei / SimHei；若均不可用则回退英文轴标签（中文仍保留在图注 md）。
运行:
  cd D:\Project\climbing_fall_analysis
  $env:PYTHONPATH="src"
  & .venv\Scripts\python.exe scripts\opensim_fe\make_paper_figures.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 无界面后端，便于批处理

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager as fm

# --------------------------------------------------------------------------- #
# 路径
# --------------------------------------------------------------------------- #
REPO = Path(__file__).resolve().parents[2]
RES = REPO / "results" / "opensim_fe"
FIGDIR = RES / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

DPI = 300
FIGSIZE_WIDE = (9.2, 4.6)
FIGSIZE_SQ = (7.6, 5.2)

# --------------------------------------------------------------------------- #
# 中文字体探测
# --------------------------------------------------------------------------- #
_CJK_CANDIDATES = [
    "Microsoft YaHei",
    "Microsoft YaHei UI",
    "SimHei",
    "SimSun",
    "Noto Sans CJK SC",
    "Source Han Sans SC",
    "PingFang SC",
]


def _pick_cjk_font() -> str | None:
    available = {f.name for f in fm.fontManager.ttflist}
    for name in _CJK_CANDIDATES:
        if name in available:
            return name
    return None


_CJK = _pick_cjk_font()
if _CJK:
    plt.rcParams["font.sans-serif"] = [_CJK, "DejaVu Sans"]
    plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = DPI
plt.rcParams["savefig.dpi"] = DPI
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.30
plt.rcParams["grid.linestyle"] = ":"

# 字体路由标记（写入图注文件）
FONT_ROUTE = (
    f"CJK (Microsoft YaHei 等可用: {_CJK})" if _CJK else "English fallback (无可用 CJK 字体)"
)


def L(cn: str, en: str) -> str:
    """按字体路由返回中文或英文标签。"""
    return cn if _CJK else en


# --------------------------------------------------------------------------- #
# 配色（色盲友好）
# --------------------------------------------------------------------------- #
C_FIB = "#d62728"      # 腓骨 高风险
C_FOOT = "#1f77b4"     # 足部
C_TIB = "#2ca02c"      # 胫骨
C_FEM = "#ff7f0e"      # 股骨
C_OTHER = "#7f7f7f"    # 其余
C_SYM = "#4c72b0"      # 对称
C_SINGLE = "#dd8452"   # 单脚
C_AX = "#4878cf"       # 轴向分量
C_LAT = "#c44e52"      # 侧向分量
C_CONS = "#9467bd"     # 保守口径
C_ELEV = "#8c564b"     # 高承载口径


def load_json(name: str) -> dict:
    with open(RES / name, encoding="utf-8") as fh:
        return json.load(fh)


def _annotate(ax, text, xy, xytext, color="black", arrow=True):
    ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        color=color,
        fontsize=9,
        ha="center",
        va="center",
        arrowprops=dict(arrowstyle="->", color=color, lw=1.1) if arrow else None,
    )


# =========================================================================== #
# fig2 — 与 [Y25] 对齐：首次骨折高度
# =========================================================================== #
def fig2_alignment() -> None:
    d = load_json("axial_subregion.json")
    rows = {r["region_cn"]: r for r in d["route_a_axial_ends"]}
    band = d["paper_facts"]["fibula_ends_fracture_m"]  # [7, 9]

    order = ["腓骨", "足部", "胫骨", "股骨", "胸椎", "腰椎", "颈椎", "骨盆", "颅骨"]
    en_map = {
        "腓骨": "Fibula", "足部": "Foot", "胫骨": "Tibia", "股骨": "Femur",
        "胸椎": "T-spine", "腰椎": "L-spine", "颈椎": "C-spine",
        "骨盆": "Pelvis", "颅骨": "Skull",
    }
    names = [L(k, en_map[k]) for k in order]

    CAP = 50.0
    heights = []
    capped = []
    for k in order:
        h = rows[k]["first_fracture_h"]
        if h is None:
            heights.append(CAP)
            capped.append(True)
        else:
            heights.append(float(h))
            capped.append(False)
    heights = np.array(heights)
    capped = np.array(capped)

    colors = []
    color_map = {"腓骨": C_FIB, "足部": C_FOOT, "胫骨": C_TIB, "股骨": C_FEM}
    for k, c in zip(order, capped):
        colors.append(color_map.get(k, C_OTHER) if not c else "#bfbfbf")

    fig, ax = plt.subplots(figsize=FIGSIZE_WIDE)

    # [Y25] 7–9 m 阈值带
    ax.axhspan(band[0], band[1], color="#d62728", alpha=0.12, zorder=0)
    ax.axhline(band[0], color="#d62728", lw=1.0, ls="--", zorder=1)
    ax.axhline(band[1], color="#d62728", lw=1.0, ls="--", zorder=1)

    x = np.arange(len(order))
    bars = ax.bar(x, heights, color=colors, edgecolor="black", linewidth=0.7, zorder=3)
    # >50 的条加斜纹，标识“截断”
    for b, c in zip(bars, capped):
        if c:
            b.set_hatch("//")
            b.set_edgecolor("#555555")

    # 数值标注
    for xi, (h, c, k) in enumerate(zip(heights, capped, order)):
        label = ">50" if c else f"{int(h) if float(h).is_integer() else h:g}"
        ax.text(xi, h + 0.6, label, ha="center", va="bottom", fontsize=10,
                fontweight="bold" if k == "腓骨" else "normal",
                color=C_FIB if k == "腓骨" else "black")

    # 腓骨在带内 — 注解
    fib_i = order.index("腓骨")
    _annotate(
        ax,
        L("腓骨首骨折 7 m\n落在 [Y25] 7–9 m 区间",
          "Fibula first fracture 7 m\nwithin [Y25] 7-9 m band"),
        xy=(fib_i, rows["腓骨"]["first_fracture_h"]),
        xytext=(fib_i + 1.35, 31),
        color=C_FIB,
    )

    ax.text(
        len(order) - 0.4, band[1] + 0.4,
        L("[Y25] 腓骨两端首骨折区间 7–9 m", "[Y25] fibula-ends first-fracture band 7-9 m"),
        ha="right", va="bottom", color="#d62728", fontsize=9,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right")
    ax.set_ylabel(L("首次骨折高度 (m)", "First-fracture height (m)"))
    ax.set_ylim(0, 56)
    ax.set_title(
        L("图 2  纯轴向+弱子区域 1D 首骨折高度与 [Y25] 对齐检验（标 measured/modeled）",
          "Fig.2  Pure-axial + weak-subregion 1D first-fracture height vs [Y25]"),
        fontsize=11,
    )
    ax.text(
        0.995, 0.02,
        L("灰斜纹 = 1–50 m 内未越阈（>50，外推）\n数据源: axial_subregion.json（modeled 1D）",
          "Hatched grey = no fracture within 1-50 m (>50, extrapolated)\nsource: axial_subregion.json (modeled 1D)"),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color="#444444",
    )
    fig.tight_layout()
    out = FIGDIR / "fig2_alignment.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[ok] {out.name}  ({out.stat().st_size} bytes)")


# =========================================================================== #
# fig3 — 单脚 vs 对称：承力侧关节反力 + 因子
# =========================================================================== #
def fig3_single_foot() -> None:
    d = load_json("nonvertical_s2.json")
    fac = d["joint_load_factors_5m"]
    joints = ["subtalar_r", "ankle_r", "hip_r", "lumbar"]
    en = {"subtalar_r": "Subtalar", "ankle_r": "Ankle", "hip_r": "Hip", "lumbar": "Lumbar"}

    sym = np.array([fac[j]["sym_vertical_n"] for j in joints])
    single = np.array([fac[j]["single_vertical_n"] for j in joints])
    factors = np.array([fac[j]["factor"] for j in joints])

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(11.2, 4.6), gridspec_kw={"width_ratios": [2.4, 1.0]}
    )

    x = np.arange(len(joints))
    w = 0.38
    b1 = ax.bar(x - w / 2, sym, w, label=L("对称 (0.5/0.5)", "Symmetric (0.5/0.5)"),
                color=C_SYM, edgecolor="black", linewidth=0.6)
    b2 = ax.bar(x + w / 2, single, w, label=L("单脚 (0/1)", "Single foot (0/1)"),
                color=C_SINGLE, edgecolor="black", linewidth=0.6)
    for bs in (b1, b2):
        ax.bar_label(bs, fmt="%.0f", fontsize=8, padding=2)

    for xi, f in zip(x, factors):
        ax.text(xi, max(sym[xi], single[xi]) * 1.10, f"×{f:.3f}",
                ha="center", va="bottom", fontsize=10, fontweight="bold",
                color=C_SINGLE if f > 1 else "#444444")

    ax.set_xticks(x)
    ax.set_xticklabels([L(j, en[j]) for j in joints])
    ax.set_ylabel(L("关节反力纵向分量 $F_y$ @5 m (N)", "Joint reaction $F_y$ @5 m (N)"))
    ax.set_ylim(0, single.max() * 1.28)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_title(L("(a) 承力侧关节反力：单脚 vs 对称 @5 m", "(a) Loaded-side joint reaction: single vs symmetric @5 m"),
                 fontsize=10)

    # 右：因子条
    cols = [C_SINGLE if f > 1 else "#8c8c8c" for f in factors]
    bars = ax2.barh(np.arange(len(joints)), factors, color=cols,
                    edgecolor="black", linewidth=0.6)
    ax2.axvline(1.0, color="black", ls=":", lw=1.2)
    ax2.bar_label(bars, fmt="%.3f", fontsize=9, padding=3)
    ax2.set_yticks(np.arange(len(joints)))
    ax2.set_yticklabels([L(j, en[j]) for j in joints])
    ax2.invert_yaxis()
    ax2.set_xlim(0, max(factors) * 1.28)
    ax2.set_xlabel(L("载荷因子 (单脚/对称)", "Load factor (single/sym)"))
    ax2.set_title(L("(b) 载荷因子", "(b) Load factor"), fontsize=10)
    ax2.text(0.98, 0.03, L("下肢 ≈×2；中轴 ≈×1", "Lower limb ≈×2; axial ≈×1"),
             transform=ax2.transAxes, ha="right", va="bottom", fontsize=8, color="#444444")

    fig.suptitle(
        L("图 3  单脚落地对承力侧关节反力的放大（数据源 nonvertical_s2.json，measured @5 m）",
          "Fig.3  Single-foot amplification of loaded-side joint reactions"),
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    out = FIGDIR / "fig3_single_foot.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[ok] {out.name}  ({out.stat().st_size} bytes)")


# =========================================================================== #
# fig4 — 踝旋后：腓骨远端 risk vs β
# =========================================================================== #
def fig4_supination() -> None:
    d = load_json("nonvertical_s5.json")
    rows = d["scenario_hard"]  # 4 heights × 3 betas
    heights = sorted({r["height_m"] for r in rows})
    betas = sorted({r["supination_deg"] for r in rows})
    data = {(r["height_m"], r["supination_deg"]): r for r in rows}

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.6))

    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(1, len(heights) - 1)) for i in range(len(heights))]

    # (a) 侧向 risk（log）
    for ci, (h, c) in enumerate(zip(heights, colors)):
        y = [max(data[(h, b)]["fibula_lateral_risk"], 1e-3) for b in betas]
        ax.plot(betas, y, "-o", color=c, lw=1.8, ms=5, label=f"h = {h:g} m")
    ax.axhline(1.0, color="#d62728", ls="--", lw=1.4)
    ax.text(betas[-1], 1.0 * 1.25, L("阈值 risk = 1", "threshold risk = 1"),
            ha="right", va="bottom", color="#d62728", fontsize=9)
    ax.set_yscale("log")
    ax.set_xticks(betas)
    ax.set_xlabel(L("旋后角 β (°)", "Supination angle β (°)"))
    ax.set_ylabel(L("腓骨远端侧向 risk（log）", "Fibula distal lateral risk (log)"))
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title(L("(a) 腓骨远端侧向 risk 随 β（刚性地面，上界式）",
                   "(a) Fibula distal lateral risk vs β (rigid floor, upper bound)"),
                 fontsize=10)
    # 临界角注记
    bc = np.mean([r["critical_inversion_deg_this_h"] for r in rows])
    _annotate(ax, L(f"临界 β_c ≈ {bc:.2f}°\n(上界式，不可直接引用)",
                    f"critical β_c ≈ {bc:.2f}°\n(upper bound)"),
              xy=(0.0, 1.0), xytext=(12, 0.06), color="#7f2704")

    # (b) F_lat vs β（同一 β 下各高度差异很小 → 几乎重合）
    for h, c in zip(heights, colors):
        y = [data[(h, b)]["F_lat_n"] for b in betas]
        ax2.plot(betas, y, "-s", color=c, lw=1.6, ms=5, label=f"h = {h:g} m")
    ax2.set_xticks(betas)
    ax2.set_xlabel(L("旋后角 β (°)", "Supination angle β (°)"))
    ax2.set_ylabel(L("冠状面横向力 $F_{lat}$ (N)", "Frontal lateral force $F_{lat}$ (N)"))
    ax2.legend(fontsize=8, loc="upper left")
    ax2.set_title(L("(b) 横向力 $F_{lat}=F_y\\sinβ$ 随 β", "(b) Lateral force $F_{lat}$ vs β"),
                  fontsize=10)
    ax2.text(0.98, 0.03,
             L("β=0 时 F_lat=0 → 逐位复现纯轴向", "β=0 ⇒ F_lat=0 (recovers pure axial)"),
             transform=ax2.transAxes, ha="right", va="bottom", fontsize=8, color="#444444")

    fig.suptitle(
        L("图 4  踝旋后使腓骨远端在真实高度内越阈（数据源 nonvertical_s5.json，modeled 上界式 1D）",
          "Fig.4  Ankle supination drives fibula distal over-threshold within real heights"),
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    out = FIGDIR / "fig4_supination.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[ok] {out.name}  ({out.stat().st_size} bytes)")


# =========================================================================== #
# fig5 — 垫子刚度 k × 旋后：total_risk(k)，无内部最优
# =========================================================================== #
def fig5_pad_tradeoff() -> None:
    d = load_json("pad_supination_tradeoff.json")
    cells = d["sweep"]["cells"]
    heights = sorted({c["height_m"] for c in cells})
    k_grid = d["sweep"]["k_grid_n_per_mm"]  # 含 None(=刚性)
    k_vals = [k for k in k_grid if k is not None]
    # 刚性极限用最大的有限 k 之后的“哨兵”位置表示，便于在 log 轴右端画 ∞
    K_INF = max(k_vals) * 4.0
    by = defaultdict(list)
    for c in cells:
        by[c["height_m"]].append(c)
    for h in by:
        by[h].sort(key=lambda c: (c["k_n_per_mm"] is None, c["k_n_per_mm"] or 0.0))

    def kx(c):
        return K_INF if c["k_n_per_mm"] is None else c["k_n_per_mm"]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.8))

    cmap = plt.get_cmap("plasma")
    colors = [cmap(0.10 + 0.70 * i / max(1, len(heights) - 1)) for i in range(len(heights))]

    # (a) total_risk(k) 各高度
    for h, c in zip(heights, colors):
        xs = [kx(cell) for cell in by[h]]
        ys = [cell["total_risk"] for cell in by[h]]
        ax.plot(xs, ys, "-o", color=c, lw=1.8, ms=4.5, label=f"h = {h:g} m")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(1.0, color="#d62728", ls="--", lw=1.4)
    ax.text(k_vals[0], 1.0 * 1.2, L("阈值 total_risk = 1", "threshold total_risk = 1"),
            ha="left", va="bottom", color="#d62728", fontsize=9)
    ax.set_xticks([50, 1e3, 1e4, 1e5, 1e6, K_INF])
    ax.set_xticklabels(["50", "1e3", "1e4", "1e5", "1e6", L("∞(刚性)", "∞(rigid)")])
    ax.set_xlabel(L("垫子刚度 k (N/mm)", "Pad stiffness k (N/mm)"))
    ax.set_ylabel(L("总风险 total_risk（log）", "total_risk (log)"))
    ax.legend(fontsize=8, loc="upper right")
    # 无内部最优 + 右向箭头
    _annotate(ax,
              L("无内部最优\nk* → 刚性极限", "No interior optimum\nk* → rigid limit"),
              xy=(K_INF, 0.72), xytext=(2e4, 0.09), color="#2b2b2b")
    ax.set_title(L("(a) total_risk(k)：单调递减、无内部极小",
                   "(a) total_risk(k): monotone decreasing, no interior min"),
                 fontsize=10)

    # (b) 分解：h=4.5 m 轴向/侧向/合计
    h_rep = max(heights)
    cells_rep = by[h_rep]
    xs = [kx(c) for c in cells_rep]
    ax2.plot(xs, [c["axial_risk"] for c in cells_rep], "-o", color=C_AX, lw=1.8, ms=4,
             label=L("轴向 axial", "axial"))
    ax2.plot(xs, [max(c["lateral_risk"], 1e-3) for c in cells_rep], "-s", color=C_LAT, lw=1.8, ms=4,
             label=L("侧向 lateral", "lateral"))
    ax2.plot(xs, [c["total_risk"] for c in cells_rep], "-^", color="#333333", lw=2.0, ms=4.5,
             label=L("合计 total", "total"))
    # k_th（total=1）
    kth = next((r["k_threshold_n_per_mm"] for r in d["sweep"]["threshold_by_height"]
                if r["height_m"] == h_rep), None)
    if kth:
        ax2.axvline(kth, color="#d62728", ls=":", lw=1.3)
        kth_txt = f"k_th≈{kth:.2e} N/mm"
        ax2.text(kth, 0.03, kth_txt, rotation=90, ha="right", va="bottom",
                 color="#d62728", fontsize=8)
    ax2.axhline(1.0, color="#d62728", ls="--", lw=1.1)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xticks([50, 1e3, 1e4, 1e5, 1e6, K_INF])
    ax2.set_xticklabels(["50", "1e3", "1e4", "1e5", "1e6", L("∞(刚性)", "∞(rigid)")])
    ax2.set_xlabel(L("垫子刚度 k (N/mm)", "Pad stiffness k (N/mm)"))
    ax2.set_ylabel(L("risk 分量（log）", "risk components (log)"))
    ax2.legend(fontsize=8, loc="upper right")
    ax2.set_title(L(f"(b) 分量分解 (h = {h_rep:g} m)", f"(b) Component breakdown (h = {h_rep:g} m)"),
                  fontsize=10)
    ax2.text(0.98, 0.03,
             L("软垫省下的轴向 << 它引入的侧向",
               "axial saving << induced lateral risk"),
             transform=ax2.transAxes, ha="right", va="bottom", fontsize=8, color="#444444")

    fig.suptitle(
        L("图 5  垫子刚度 × 旋后两难：无内部最优、k*→刚性极限"
          "（数据源 pad_supination_tradeoff.json；R_FE measured，β(k) modeled）",
          "Fig.5  Pad stiffness × supination tradeoff: no interior optimum, k*→rigid limit"),
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    out = FIGDIR / "fig5_pad_tradeoff.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[ok] {out.name}  ({out.stat().st_size} bytes)")


# =========================================================================== #
# fig6 — 踝韧带判据：ATFL 关节级 risk vs β（保守 / 高承载）
# =========================================================================== #
def fig6_ligament() -> None:
    d = load_json("ankle_ligament_criterion.json")
    rows = d["scenario"]["rows"]
    heights = sorted({r["height_m"] for r in rows})
    betas = sorted({r["supination_deg"] for r in rows})
    thr = d["critical_thresholds"]
    m_cons = thr["joint_moment_thresholds_nm"]["conservative_no_preload"]
    m_elev = thr["joint_moment_thresholds_nm"]["elevated_2kn_preload"]
    crit = thr["our_model_critical_angle_deg_at_21kn"]

    by = {(r["height_m"], r["supination_deg"], r["preload_model"]): r for r in rows}

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), sharey=True)
    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(1, len(heights) - 1)) for i in range(len(heights))]

    panels = [
        ("conservative", axes[0], C_CONS,
         L(f"保守口径（无预载，$M_{{fail}}$={m_cons:.0f} N·m）",
           f"Conservative (no preload, $M_{{fail}}$={m_cons:.0f} N·m)")),
        ("elevated", axes[1], C_ELEV,
         L(f"高承载口径（2 kN 预载，$M_{{fail}}$={m_elev:.0f} N·m）",
           f"Elevated (2 kN preload, $M_{{fail}}$={m_elev:.0f} N·m)")),
    ]

    for model, ax, base, title in panels:
        for h, c in zip(heights, colors):
            y = [max(by[(h, b, model)]["atfl_joint_risk"], 1e-3) for b in betas]
            ax.plot(betas, y, "-o", color=c, lw=1.8, ms=4.5, label=f"h = {h:g} m")
        ax.axhline(1.0, color="#d62728", ls="--", lw=1.4)
        ax.text(betas[-1], 1.05, L("扭伤阈值 risk=1", "sprain threshold risk=1"),
                ha="right", va="bottom", color="#d62728", fontsize=8.5)
        ax.set_yscale("log")
        ax.set_xticks(betas)
        ax.set_xlabel(L("旋后角 β (°)", "Supination angle β (°)"))
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7.5, loc="lower right")

    axes[0].set_ylabel(L("ATFL 关节级 risk（log）", "ATFL joint-level risk (log)"))

    # 临界角注记
    axes[0].annotate(
        L(f"模型临界 β≈{crit['conservative']:.2f}°（21 kN GRF）",
          f"model critical β≈{crit['conservative']:.2f}°"),
        xy=(crit["conservative"], 1.0), xytext=(11, 0.05),
        color="#7f2704", fontsize=8.5, ha="center",
        arrowprops=dict(arrowstyle="->", color="#7f2704", lw=1.0),
    )
    axes[1].annotate(
        L(f"模型临界 β≈{crit['elevated']:.2f}°", f"model critical β≈{crit['elevated']:.2f}°"),
        xy=(crit["elevated"], 1.0), xytext=(15, 0.06),
        color="#4d2010", fontsize=8.5, ha="center",
        arrowprops=dict(arrowstyle="->", color="#4d2010", lw=1.0),
    )

    fig.suptitle(
        L("图 6  踝韧带判据：ATFL 关节级扭伤 risk 随 β（保守 / 高承载两口径）"
          "（数据源 ankle_ligament_criterion.json）",
          "Fig.6  Ankle ligament criterion: ATFL joint-level sprain risk vs β"),
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    out = FIGDIR / "fig6_ligament.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[ok] {out.name}  ({out.stat().st_size} bytes)")


# =========================================================================== #
# fig1 — 管线总览（四段水平流 + 二判据分支）
# =========================================================================== #
def fig1_pipeline() -> None:
    r"""论文图 1：管线总览（论文初稿 §2.1 + §2.5）。

    4 段水平流: GRF(pad) → OpenSim 多体正动力学 → 载荷传递 → FEBio 单骨应力.
    末尾分支到 2 个判据盒（骨骼判据 / 踝韧带判据）。
    顶部横幅: 非垂直落地 6 自由度.
    底部注记: 两端 1D 接口 → 升为 3D（中段本已 3D）.

    视觉约定（与图 2–6 同款）：FancyBboxPatch 圆角矩形 + 浅填色 + 细黑边；
    节点文字黑，注释灰；连接线纯黑实线 1.2pt 配箭头样式；分支用一对角连接器。
    """
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(11.6, 5.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 60)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- 配色（浅 fill + 黑边，色盲友好） ----
    C_STAGE_FILL = "#eaf2fb"     # 阶段浅蓝
    C_STAGE_EDGE = "#0a0a0a"     # 阶段边框
    C_OUT_FILL = "#fdebd3"       # 判据浅橙
    C_OUT_EDGE = "#0a0a0a"
    C_BANNER_FILL = "#f4f4f4"    # 顶部横幅浅灰
    C_BANNER_EDGE = "#7a7a7a"
    C_NOTE_FILL = "#fbfbfb"
    C_NOTE_EDGE = "#bdbdbd"
    C_BRANCH = "#0a0a0a"          # 分支连线
    C_ARROW = "#0a0a0a"

    # ---- 阶段 4 盒（横向排列） ----
    stage_w, stage_h = 22.0, 9.0
    stage_y = 30.0                 # 阶段盒中心 y
    stage_xs = [3.0, 27.0, 51.0, 75.0]  # 阶段盒左下 x

    stages_cn = [
        # (主标题, 副标题)
        (L("地/垫反力 GRF", "Ground/Pad Reaction GRF"),
         L("pad.py：1D 双质点 → 每足 3D 力向量", "pad.py: 1D 2-mass → per-foot 3D")),
        (L("OpenSim 多体正动力学", "OpenSim Multibody FD"),
         L("关节 3D 反力（subtalar_reaction）", "Joint 3D reactions (subtalar_reaction)")),
        (L("载荷传递", "Load Transfer"),
         L("关节 wrench → FE 边界（load_transfer）", "Joint wrench → FE BC (load_transfer)")),
        (L("FEBio", "FEBio"),
         L("单骨 von Mises 应力", "Single-bone von Mises stress")),
    ]

    for (x, (title, sub)) in zip(stage_xs, stages_cn):
        box = FancyBboxPatch(
            (x, stage_y - stage_h / 2), stage_w, stage_h,
            boxstyle="round,pad=0.4,rounding_size=1.2",
            linewidth=1.0, edgecolor=C_STAGE_EDGE, facecolor=C_STAGE_FILL,
            zorder=2,
        )
        ax.add_patch(box)
        ax.text(
            x + stage_w / 2, stage_y + 2.2, title,
            ha="center", va="center", fontsize=11.5, fontweight="bold",
            color="#101010", zorder=3,
        )
        ax.text(
            x + stage_w / 2, stage_y - 2.6, sub,
            ha="center", va="center", fontsize=8.5,
            color="#333333", zorder=3,
        )

    # ---- 主干箭头 1→2→3→4 ----
    arrow_kw = dict(
        arrowstyle="-|>", mutation_scale=18, linewidth=1.4,
        color=C_ARROW, shrinkA=2, shrinkB=2, zorder=4,
    )
    for i in range(3):
        x0 = stage_xs[i] + stage_w
        x1 = stage_xs[i + 1]
        y0 = stage_y
        arr = FancyArrowPatch((x0, y0), (x1, y0), **arrow_kw)
        ax.add_patch(arr)

    # ---- 分支：从 4 右侧 (FEBio) 出来 → 两判据盒 ----
    out_w, out_h = 22.0, 7.5
    out_xs = [3.0, 75.0]    # 左下 x（与阶段盒左右两端对齐：1 ↔ 骨骼判据、4 ↔ 踝韧带判据）
    out_y_top = 51.0          # 判据顶部 y
    out_y_bot = out_y_top - out_h

    outputs_cn = [
        (L("骨骼判据", "Bone Criterion"),
         L("按子区域材料强度", "Subregion material strength")),
        (L("踝韧带判据", "Ankle Ligament Criterion"),
         L("ATFL/CFL 关节级内翻力矩", "ATFL/CFL joint-level inversion moment")),
    ]

    # 分支起点：阶段 4 中心 y = stage_y，分支先向上，再分叉
    branch_start = (stage_xs[3] + stage_w / 2, stage_y)         # FEBio 中点
    branch_junction = (stage_xs[3] + stage_w / 2, stage_y + 6)  # 向上汇合点
    branch_targets = [(out_xs[0] + out_w / 2, out_y_bot),
                        (out_xs[1] + out_w / 2, out_y_bot)]

    # 阶段 4 → 汇合点
    ax.add_line(Line2D([branch_start[0], branch_junction[0]],
                       [branch_start[1], branch_junction[1]],
                       color=C_BRANCH, lw=1.4, zorder=3))
    # 汇合点 → 两判据盒（带箭头）
    for tgt in branch_targets:
        arr = FancyArrowPatch(
            branch_junction, tgt,
            arrowstyle="-|>", mutation_scale=18, linewidth=1.4,
            color=C_ARROW, shrinkA=0, shrinkB=2, zorder=4,
        )
        ax.add_patch(arr)

    # ---- 两判据盒 ----
    for (x, (title, sub)) in zip(out_xs, outputs_cn):
        box = FancyBboxPatch(
            (x, out_y_bot), out_w, out_h,
            boxstyle="round,pad=0.4,rounding_size=1.2",
            linewidth=1.0, edgecolor=C_OUT_EDGE, facecolor=C_OUT_FILL,
            zorder=2,
        )
        ax.add_patch(box)
        ax.text(
            x + out_w / 2, out_y_bot + 4.4, title,
            ha="center", va="center", fontsize=10.5, fontweight="bold",
            color="#101010", zorder=3,
        )
        ax.text(
            x + out_w / 2, out_y_bot + 2.0, sub,
            ha="center", va="center", fontsize=8.5,
            color="#333333", zorder=3,
        )

    # ---- 顶部横幅 ----
    banner_w, banner_h = 96.0, 6.5
    banner_x, banner_y = 2.0, 52.5
    banner_box = FancyBboxPatch(
        (banner_x, banner_y), banner_w, banner_h,
        boxstyle="round,pad=0.3,rounding_size=0.8",
        linewidth=0.8, edgecolor=C_BANNER_EDGE, facecolor=C_BANNER_FILL,
        zorder=2,
    )
    ax.add_patch(banner_box)
    banner_text = L(
        "非垂直落地 6 自由度（S1–S4：①方向 ②地面倾角 ③姿势 ④单脚 ⑤力矩；⑥未做）",
        "Non-vertical landing 6 DoF (S1-S4: 1.direction 2.tilt 3.posture 4.single 5.moment; 6.todo)",
    )
    ax.text(
        banner_x + banner_w / 2, banner_y + banner_h / 2 + 0.8, "非垂直落地 6 自由度",
        ha="center", va="center", fontsize=11, fontweight="bold",
        color="#101010", zorder=3,
    )
    ax.text(
        banner_x + banner_w / 2, banner_y + banner_h / 2 - 1.6,
        L("S1–S4：①方向 ②地面倾角 ③姿势 ④单脚 ⑤力矩；⑥未做",
          "S1-S4: 1.direction 2.tilt 3.posture 4.single 5.moment; 6.todo"),
        ha="center", va="center", fontsize=9, color="#444444", zorder=3,
    )

    # ---- 底部注记 ----
    note_w, note_h = 96.0, 5.2
    note_x, note_y = 2.0, 0.5
    note_box = FancyBboxPatch(
        (note_x, note_y), note_w, note_h,
        boxstyle="round,pad=0.3,rounding_size=0.7",
        linewidth=0.8, edgecolor=C_NOTE_EDGE, facecolor=C_NOTE_FILL,
        zorder=2,
    )
    ax.add_patch(note_box)
    ax.text(
        note_x + note_w / 2, note_y + note_h / 2,
        L("两端 1D 接口 → 升为 3D（中段本已 3D）",
          "Both-end 1D interfaces → promoted to 3D (middle is already 3D)"),
        ha="center", va="center", fontsize=10, fontweight="bold",
        color="#202020", zorder=3,
    )

    # ---- 图题 ----
    ax.text(
        50.0, 12.5,
        L("图 1  管线总览：四段水平流 + 末端分支至骨骼/踝韧带判据",
          "Fig.1  Pipeline overview: 4-stage horizontal flow + branch to bone / ankle-ligament criteria"),
        ha="center", va="center", fontsize=11.5, fontweight="bold",
        color="#101010",
    )

    fig.tight_layout(pad=0.4)
    out = FIGDIR / "fig1_pipeline.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[ok] {out.name}  ({out.stat().st_size} bytes)")


# =========================================================================== #
# fig7 — 踝外旋 / SER：腓骨远端外旋 risk 随高度（图 4 旋后的横截面镜像）
# =========================================================================== #
def fig7_external_rotation() -> None:
    r"""论文图 7：踝外旋（external rotation / 横截面扭转）使腓骨远端 / 外踝在
    真实抱石高度内越阈 —— 图 4（旋后 / 冠状面）的横截面镜像。

    数据源：`results/opensim_fe/ankle_external_rotation.json::scenario_hard`
    (on_pad=False, 24 行 = 4 高度 × 6 外旋角)。

    设计意图：
      * x 轴 = 高度 h (m)，4 个真实抱石高度 ∈ {2, 3, 3.5, 4.5}；
      * y 轴 = 腓骨远端外旋 risk（log 轴，5 个量级跨度）；
      * 5 条曲线，θ ∈ {0°, 5°, 10°, 15°, 30°}；
      * 整个 x 域均落在「真实抱石高度 ≤ 4.5 m」带内（橙色淡填充）—— 暗示全部测试高度都已
        属于真实场景；
      * 红色虚线 `risk = 1` 阈值；
      * 临界外旋角 `θ_c ≈ 0.89–1.06°`（上界式 1D）以小注记给出。
    """
    d = load_json("ankle_external_rotation.json")
    rows = d["scenario_hard"]
    heights = sorted({r["height_m"] for r in rows})
    # 仅绘制 5 个 θ ∈ {0, 5, 10, 15, 30}，与正文 §5 / §5.1 一致
    thetas = [0.0, 5.0, 10.0, 15.0, 30.0]
    by = {(r["height_m"], r["external_rotation_deg"]): r for r in rows}

    fig, ax = plt.subplots(figsize=FIGSIZE_SQ)

    # ---- 「真实抱石高度 ≤4.5 m」淡橙色带 ----
    ax.axvspan(
        min(heights), max(heights),
        color="#ff7f0e", alpha=0.08, zorder=0,
        label=L("真实抱石高度 ≤ 4.5 m 区域", "Real bouldering band h ≤ 4.5 m"),
    )

    # ---- 阈值 risk = 1 ----
    ax.axhline(1.0, color="#d62728", ls="--", lw=1.4, zorder=1)
    ax.text(
        max(heights) + 0.05, 1.0 * 1.18,
        L("阈值 risk = 1", "threshold risk = 1"),
        ha="right", va="bottom", color="#d62728", fontsize=9,
    )

    # ---- 5 条曲线：按 θ 从灰（安全）→ 红（越阈）渐进 ----
    # 注：风险随 θ 单调增大；θ=0 退化为 0（逐位复现纯轴向）。
    theta_palette = [
        "#7f7f7f",   # θ=0° 灰（基线 / 纯轴向）
        "#fdae61",   # θ=5° 浅橙（临界下方）
        "#f16913",   # θ=10° 橙（首次越阈）
        "#d62728",   # θ=15° 红
        "#8b0a50",   # θ=30° 暗紫红（深越阈）
    ]
    markers = ["o", "s", "^", "D", "v"]

    for ci, theta in enumerate(thetas):
        xs, ys = [], []
        for h in heights:
            r = by[(h, theta)]
            # θ=0° 的 risk 严格为 0，log 轴下用 ε 显示，仅示意曲线位置
            v = r["fibula_external_rotation_risk"]
            ys.append(max(v, 1e-3) if v > 0 else 1e-3)
            xs.append(h)
        # θ=0 时所有 y 都被钳到 1e-3 → 画一条灰线压底
        ax.plot(
            xs, ys,
            marker=markers[ci], color=theta_palette[ci],
            lw=1.9 if theta == 0 else 2.1,
            ms=8 if theta != 0 else 6,
            markeredgecolor="black", markeredgewidth=0.6,
            label=L(f"θ = {int(theta)}°", f"θ = {int(theta)}°"),
            zorder=3,
        )
        # 在曲线右端写一个轻量级 "×N" 数值标签，方便读出量级
        if theta > 0:
            y_last = ys[-1]
            ax.text(
                max(heights) + 0.06, y_last,
                f"×{y_last:.0f}",
                color=theta_palette[ci], fontsize=8.5,
                fontweight="bold", va="center", ha="left",
            )

    # ---- 临界角注记（与图 4 同款，色改褐系以区分旋后） ----
    # §8: θ_c ∈ {1.06, 1.01, 0.98, 0.89}° @ {2, 3, 3.5, 4.5} m
    crit_c = [by[(h, 0.0)]["critical_external_rotation_deg_this_h"] for h in heights]
    tc_min, tc_max = min(crit_c), max(crit_c)
    ax.annotate(
        L(
            f"临界 θ_c ≈ {tc_min:.2f}–{tc_max:.2f}°\n(上界式 1D，绝对值不可引用)",
            f"critical θ_c ≈ {tc_min:.2f}–{tc_max:.2f}°\n(upper-bound 1D, abs. not citable)",
        ),
        xy=(min(heights), 1.0),
        xytext=(min(heights) + 0.18, 0.18),
        color="#5a1f02", fontsize=8.8, ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color="#5a1f02", lw=1.1),
    )

    # ---- θ=0 安全 vs θ≥10° 全部越阈 注解 ----
    ax.text(
        min(heights) - 0.05, 0.18,
        L("θ = 0°（纯轴向）\n真实高度内安全",
          "θ = 0° (pure axial)\nsafe inside real heights"),
        ha="left", va="center", fontsize=9, color="#444444",
        bbox=dict(boxstyle="round,pad=0.35", fc="#fbfbfb", ec="#bbbbbb", lw=0.7),
    )
    ax.text(
        max(heights) - 0.10, 18,
        L("θ ≥ 10°\n真实高度内全部越阈",
          "θ ≥ 10°\nall over-threshold"),
        ha="right", va="center", fontsize=10, fontweight="bold", color="#7f0a1c",
        bbox=dict(boxstyle="round,pad=0.40", fc="#fff2f2", ec="#d62728", lw=0.9),
    )

    # 轴 / 标题
    ax.set_xticks(heights)
    ax.set_xticklabels([f"{h:g}" for h in heights])
    ax.set_xlabel(L("抱石落地高度 h (m)", "Bouldering landing height h (m)"))
    ax.set_ylabel(L("腓骨远端 / 外踝 外旋 risk（log）", "Fibula distal / lateral malleolus SER risk (log)"))
    ax.set_yscale("log")
    ax.set_ylim(5e-4, 60)
    ax.set_xlim(min(heights) - 0.25, max(heights) + 1.05)

    # 图例（两列）
    leg = ax.legend(
        loc="upper left", fontsize=9, framealpha=0.92,
        ncol=2, title=L("外旋角 θ", "External rotation θ"),
        title_fontsize=9,
    )
    leg.get_title().set_fontweight("bold")

    # 标题
    ax.set_title(
        L(
            "图 7  踝外旋（external rotation）使腓骨远端 / 外踝在真实抱石高度内越阈"
            "\n（数据源 ankle_external_rotation.json，modeled 上界式 1D，"
            "刚性地面，d = 30 mm）",
            "Fig.7  Ankle external rotation (SER) drives fibula distal over-threshold "
            "within real bouldering heights (rigid floor, upper-bound 1D)",
        ),
        fontsize=11, pad=12,
    )

    # 底部 caveat
    fig.text(
        0.5, 0.01,
        L(
            "诚实边界：上界式 / 1D（材料强度 + 中段截面 + 扭转 + 线性叠加）—— "
            "仅方向性 / 是否越阈，绝对值不可引用；真实外伤角 5–15°，远高于 θ_c ≈ 1°。",
            "Caveat: upper-bound / 1D screening only — directional / threshold-crossing, "
            "absolute values not citable.",
        ),
        ha="center", va="bottom", fontsize=8.0, color="#444444", style="italic",
    )

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = FIGDIR / "fig7_external_rotation.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[ok] {out.name}  ({out.stat().st_size} bytes)")


# =========================================================================== #
# fig8 — 三踝骨折链式机制（Lauge-Hansen SER；§3.6 / §4.2）
# =========================================================================== #
def fig8_trimalleolar_ser() -> None:
    r"""论文图 8：三踝骨折（trimalleolar）的链式机制示意（Lauge-Hansen 旋后-外旋 SER）。

    4 阶段链式（同一外旋下依序发生）：
        SER-I    AITFL（下胫腓前韧带）撕裂             [22,23]
        SER-II   外踝（腓骨远端）斜形骨折              [22,23] — 本文 §3.4 / 图 7 S6 覆盖
        SER-III  后踝骨折（PITFL 撕脱 + 距骨撞击）    [26,27]
        SER-IV   内踝骨折 / 三角韧带撕裂               [26,27]

    视觉约定（与图 1 同款）：
        * FancyBboxPatch 圆角矩形 + 浅填色 + 细黑边；
        * 阶段盒浅蓝填充（#eaf2fb）；本文建模覆盖盒浅橙（#fdebd3，与图 1 输出盒同款）；
          文献支撑盒浅灰（#f4f4f4，与图 1 顶部横幅同款）；
        * 黑色实线连接箭头（1→2→3→4）；外旋箭头红色（#c11b1b）；
        * 中文贯穿（Microsoft YaHei）；底部脚注 italic。
    """
    from matplotlib.patches import (
        FancyBboxPatch, FancyArrowPatch, Rectangle, Polygon, Ellipse,
    )

    fig, ax = plt.subplots(figsize=(12.0, 7.6))
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 76)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- 配色（与图 1 同款） ----
    C_STAGE_FILL = "#eaf2fb"      # 阶段浅蓝
    C_STAGE_EDGE = "#0a0a0a"
    C_COVER_FILL = "#fdebd3"      # 本文建模 = 浅橙（与图 1 输出盒同款）
    C_COVER_EDGE = "#0a0a0a"
    C_LIT_FILL = "#f4f4f4"        # 文献支撑 = 浅灰（与图 1 横幅同款）
    C_LIT_EDGE = "#7a7a7a"
    C_BANNER_FILL = "#f4f4f4"
    C_BANNER_EDGE = "#7a7a7a"
    C_NOTE_FILL = "#fbfbfb"
    C_NOTE_EDGE = "#bdbdbd"
    C_BONE_FILL = "#ecebe5"       # 骨：浅米色
    C_BONE_EDGE = "#2a2a2a"
    C_LIG_OK = "#5a5a5a"          # 完整韧带 = 灰
    C_LIG_BAD = "#c11b1b"         # 撕裂 / 骨折 = 红
    C_FRAG_FILL = "#ffd6d6"       # 撕脱骨片 = 浅粉
    C_FRAG_EDGE = "#9a1414"
    C_ER = "#c11b1b"              # 外旋箭头 = 红
    C_TITLE = "#101010"

    # ---- 顶部横幅（机制说明，非图题） ----
    banner_w, banner_h = 116.0, 7.0
    banner_x, banner_y = 2.0, 67.0
    banner_box = FancyBboxPatch(
        (banner_x, banner_y), banner_w, banner_h,
        boxstyle="round,pad=0.3,rounding_size=0.8",
        linewidth=0.8, edgecolor=C_BANNER_EDGE, facecolor=C_BANNER_FILL,
        zorder=2,
    )
    ax.add_patch(banner_box)
    ax.text(
        banner_x + banner_w / 2, banner_y + banner_h / 2 + 1.1,
        L("三踝骨折的链式机制（Lauge-Hansen 旋后-外旋 SER）",
          "Trimalleolar fracture chain mechanism (Lauge-Hansen SER)"),
        ha="center", va="center", fontsize=11.5, fontweight="bold",
        color=C_TITLE, zorder=3,
    )
    ax.text(
        banner_x + banner_w / 2, banner_y + banner_h / 2 - 1.7,
        L("脚绕胫骨长轴外旋 → 距骨外侧楔形面将腓骨远端沿横截面扭转 + 前后弯曲 → 四阶段依序失效",
          "Foot ER about tibial long axis → talar lateral wedge twists distal fibula "
          "→ four SER stages fail in sequence"),
        ha="center", va="center", fontsize=9.0, color="#444444", zorder=3,
    )

    # ---- 外旋箭头（覆盖 4 阶段） ----
    er_y = 61.5
    er_x0, er_x1 = 8.0, 112.0
    er_arrow = FancyArrowPatch(
        (er_x0, er_y), (er_x1, er_y),
        arrowstyle="-|>", mutation_scale=22, linewidth=1.8,
        color=C_ER, shrinkA=2, shrinkB=2, zorder=3,
    )
    ax.add_patch(er_arrow)
    ax.text(
        (er_x0 + er_x1) / 2, er_y + 1.6,
        L("脚绕胫骨长轴外旋（external rotation of foot about tibial long axis，θ ≈ 5–15°）",
          "External rotation of foot about tibial long axis (θ ≈ 5–15°)"),
        ha="center", va="bottom", fontsize=9.5, fontweight="bold", color=C_ER,
        zorder=4,
    )

    # ---- 踝部解剖标识（背侧观键） ----
    ax.text(
        60.0, 57.6,
        L("踝部示意（背侧观）：胫骨 / 内踝（内侧）· 腓骨 / 外踝（外侧）· 距骨 · 跟骨",
          "Ankle schematic (dorsal view): tibia / medial malleolus · "
          "fibula / lateral malleolus · talus · calcaneus"),
        ha="center", va="center", fontsize=7.6, color="#666666", zorder=3,
    )

    # ---- 4 阶段盒参数 ----
    stage_w, stage_h = 26.5, 33.0
    stage_y_bot = 22.5
    stage_y_top = stage_y_bot + stage_h
    stage_xs = [3.5, 31.5, 59.5, 87.5]

    stages_cn = [
        ("SER-I", L("AITFL（下胫腓前韧带）撕裂",
                    "AITFL (anterior inferior tibiofibular ligament) tear")),
        ("SER-II", L("外踝（腓骨远端）斜形骨折",
                     "Lateral malleolus (distal fibula) spiral fracture")),
        ("SER-III", L("后踝骨折（PITFL 撕脱 + 距骨撞击）",
                      "Posterior malleolus fracture (PITFL avulsion + talar impact)")),
        ("SER-IV", L("内踝骨折 / 三角韧带撕裂",
                     "Medial malleolus fracture / deltoid ligament tear")),
    ]
    source_tags = [
        L("[22,23] Lauge-Hansen 1950 / Yde 1980\n（链式分类）",
          "[22,23] Lauge-Hansen 1950 / Yde 1980\n(chain classification)"),
        L("[22,23]（外踝臂 — 本文 §3.4 / 图 7 S6 覆盖）",
          "[22,23] (lateral arm — covered by §3.4 / Fig.7 S6)"),
        L("[26] Haraguchi & Armiger 2020\n（尸体，PITFL + 距骨撞击）",
          "[26] Haraguchi & Armiger 2020\n(cadaveric, PITFL + axial impact)"),
        L("[27] Zhang 等 2022\n（SER 有限元，旋转应力 > 竖向）",
          "[27] Zhang et al. 2022\n(SER FE, rotational > vertical stress)"),
    ]

    # ---- 踝部 motif 绘制 (dorsal view: anterior=top, posterior=bottom, lateral=left, medial=right) ----
    def _draw_ankle(ax, x0, y0, w, h, stage):
        # Motif 局部坐标：原点 = (x0+w/2, y0+h/2)，范围 ±9
        cx_box = x0 + w / 2
        cy_box = y0 + h / 2
        s = min(w, h) / 18.0

        def P(px, py):
            return (cx_box + px * s, cy_box + py * s)

        bone_kw = dict(facecolor=C_BONE_FILL, edgecolor=C_BONE_EDGE,
                       linewidth=1.1, zorder=4)

        # Tibia: 居中矩形
        x_t, y_t = P(0, 3.5)
        ax.add_patch(Rectangle((x_t - 2.5 * s, y_t - 4.5 * s),
                                5.0 * s, 9.0 * s, **bone_kw))
        # Fibula: 胫骨左侧（外侧）较窄矩形
        x_f, y_f = P(-5.25, 3.5)
        ax.add_patch(Rectangle((x_f - 0.75 * s, y_f - 4.5 * s),
                                1.5 * s, 9.0 * s, **bone_kw))
        # Posterior malleolus (intact in stages 1, 2)：胫骨底缘后突
        if stage < 3:
            x_p, y_p = P(0, -1.85)
            ax.add_patch(Rectangle((x_p - 2.5 * s, y_p - 0.7 * s),
                                    5.0 * s, 1.4 * s, **bone_kw))
        # Medial malleolus：胫骨内侧（右侧）远端突起（stage 4 保留实体，其上叠加骨折线）
        x_m, y_m = P(3.0, -1.0)
        ax.add_patch(Rectangle((x_m - 0.5 * s, y_m - 1.2 * s),
                                1.0 * s, 2.4 * s, **bone_kw))
        # Talus: 胫骨 / 腓骨下方的梯形
        tal_pts = [P(-5.0, -2.5), P(5.0, -2.5),
                   P(4.0, -1.0), P(-4.0, -1.0)]
        ax.add_patch(Polygon(tal_pts, closed=True, **bone_kw))
        # Calcaneus: 距骨下方的椭圆
        x_c, y_c = P(0, -4.8)
        ax.add_patch(Ellipse((x_c, y_c), 9.5 * s, 4.6 * s, **bone_kw))

        # --- AITFL（胫腓骨之间，靠近前侧/顶部） ---
        a1, a2 = P(-4.5, 7.0), P(-2.5, 7.0)
        if stage == 1:
            # 当下撕裂：灰线 + 红色 X
            ax.plot([a1[0], a2[0]], [a1[1], a2[1]],
                    color=C_LIG_OK, lw=2.0, zorder=5)
            ax.plot([a1[0] + 0.6 * s, a2[0] - 0.6 * s],
                    [a1[1] + 0.45 * s, a2[1] - 0.45 * s],
                    color=C_LIG_BAD, lw=3.0, zorder=6)
            ax.plot([a1[0] + 0.6 * s, a2[0] - 0.6 * s],
                    [a1[1] - 0.45 * s, a2[1] + 0.45 * s],
                    color=C_LIG_BAD, lw=3.0, zorder=6)
        else:
            # 此前已撕裂：红色虚线
            ax.plot([a1[0], a2[0]], [a1[1], a2[1]],
                    color=C_LIG_BAD, lw=2.2, linestyle=(0, (3, 2)), zorder=5)
        ax.text(*P(-3.5, 7.9), "AITFL", ha="center", va="bottom",
                fontsize=7.0, color=C_LIG_BAD, zorder=7)

        # --- PITFL（胫腓骨之间，AITFL 稍下） ---
        p1, p2 = P(-4.5, 5.5), P(-2.5, 5.5)
        if stage == 3:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                    color=C_LIG_OK, lw=2.0, zorder=5)
            ax.plot([p1[0] + 0.6 * s, p2[0] - 0.6 * s],
                    [p1[1] + 0.45 * s, p2[1] - 0.45 * s],
                    color=C_LIG_BAD, lw=3.0, zorder=6)
            ax.plot([p1[0] + 0.6 * s, p2[0] - 0.6 * s],
                    [p1[1] - 0.45 * s, p2[1] + 0.45 * s],
                    color=C_LIG_BAD, lw=3.0, zorder=6)
            ax.text(*P(-3.5, 4.9), "PITFL", ha="center", va="bottom",
                    fontsize=7.0, color=C_LIG_BAD, zorder=7)
        elif stage == 4:
            # 此前已撕裂
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                    color=C_LIG_BAD, lw=2.2, linestyle=(0, (3, 2)), zorder=5)
            ax.text(*P(-3.5, 4.9), "PITFL", ha="center", va="bottom",
                    fontsize=7.0, color=C_LIG_BAD, zorder=7)
        else:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                    color=C_LIG_OK, lw=2.0, zorder=5)
            ax.text(*P(-3.5, 4.9), "PITFL", ha="center", va="bottom",
                    fontsize=7.0, color="#555555", zorder=7)

        # --- Stage 2+：外踝（腓骨远端）螺旋 / 斜形骨折 ---
        if stage >= 2:
            fx_x = -5.25
            n_seg = 9
            ys = np.linspace(0.0, 2.8, n_seg)
            xs = fx_x + 0.45 * np.sin(np.linspace(0, 3.2 * np.pi, n_seg))
            pts = [P(x, y) for x, y in zip(xs, ys)]
            ax.plot([p[0] for p in pts], [p[1] for p in pts],
                    color=C_LIG_BAD, lw=2.5, zorder=6)
            ax.text(*P(-6.9, 1.4),
                    L("外踝\n斜折", "Lat.\nspiral"),
                    ha="right", va="center", fontsize=7.5, color=C_LIG_BAD,
                    fontweight="bold", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.25", fc="#fff2f2",
                              ec=C_LIG_BAD, lw=0.6))

        # --- Stage 3+：后踝撕脱骨片（脱离胫骨底缘） ---
        if stage >= 3:
            frag_pts = [P(-2.2, -3.6), P(2.2, -3.6),
                        P(2.2, -2.8), P(-2.2, -2.8)]
            ax.add_patch(Polygon(frag_pts, closed=True,
                                 facecolor=C_FRAG_FILL, edgecolor=C_FRAG_EDGE,
                                 linewidth=1.1, zorder=5))
            # 撕脱裂纹（虚线连接到原位）
            ax.plot([P(0, -2.5)[0], P(0, -2.8)[0]],
                    [P(0, -2.5)[1], P(0, -2.8)[1]],
                    color=C_FRAG_EDGE, lw=1.5, linestyle=(0, (2, 2)),
                    zorder=6)
            ax.text(*P(4.5, -3.2),
                    L("后踝撕脱\n骨片", "Post.\navul. fx"),
                    ha="left", va="center", fontsize=7.5, color=C_LIG_BAD,
                    fontweight="bold", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.25", fc="#fff2f2",
                              ec=C_LIG_BAD, lw=0.6))

        # --- Stage 4+：内踝（胫骨内侧远端突起）斜形骨折 ---
        if stage >= 4:
            crack_pts = [P(2.55, 0.2), P(3.45, -2.2)]
            ax.plot([crack_pts[0][0], crack_pts[1][0]],
                    [crack_pts[0][1], crack_pts[1][1]],
                    color=C_LIG_BAD, lw=2.5, zorder=6)
            ax.text(*P(5.2, -1.0),
                    L("内踝\n骨折", "Med.\nfx"),
                    ha="left", va="center", fontsize=7.5, color=C_LIG_BAD,
                    fontweight="bold", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.25", fc="#fff2f2",
                              ec=C_LIG_BAD, lw=0.6))

    # ---- 绘制各阶段盒 ----
    for i, ((x, (sn, label)), src) in enumerate(zip(
            [(stage_xs[k], stages_cn[k]) for k in range(4)], source_tags)):
        # 阶段盒
        box = FancyBboxPatch(
            (x, stage_y_bot), stage_w, stage_h,
            boxstyle="round,pad=0.4,rounding_size=1.2",
            linewidth=1.0, edgecolor=C_STAGE_EDGE, facecolor=C_STAGE_FILL,
            zorder=2,
        )
        ax.add_patch(box)
        # 阶段编号（顶部居中，粗体）
        ax.text(
            x + stage_w / 2, stage_y_bot + stage_h - 2.8, sn,
            ha="center", va="center", fontsize=12.0, fontweight="bold",
            color=C_TITLE, zorder=3,
        )
        # 阶段损伤标签（编号下方）
        ax.text(
            x + stage_w / 2, stage_y_bot + stage_h - 6.5, label,
            ha="center", va="center", fontsize=8.6, color="#222222",
            zorder=3,
        )
        # 踝部 motif（中部）
        motif_x = x + 1.8
        motif_y = stage_y_bot + 5.5
        motif_w = stage_w - 3.6
        motif_h = stage_h - 13.5
        _draw_ankle(ax, motif_x, motif_y, motif_w, motif_h, stage=i + 1)
        # 文献来源标签（盒底）
        ax.text(
            x + stage_w / 2, stage_y_bot + 2.4, src,
            ha="center", va="center", fontsize=7.4, color="#444444",
            style="italic", zorder=3,
        )

    # ---- 连接箭头 1→2→3→4（与图 1 同款） ----
    arrow_kw = dict(
        arrowstyle="-|>", mutation_scale=20, linewidth=1.5,
        color="#0a0a0a", shrinkA=2, shrinkB=2, zorder=4,
    )
    for i in range(3):
        x0 = stage_xs[i] + stage_w
        x1 = stage_xs[i + 1]
        y0 = stage_y_bot + stage_h / 2
        arr = FancyArrowPatch((x0, y0), (x1, y0), **arrow_kw)
        ax.add_patch(arr)

    # ---- 覆盖度注解 1：本文建模 = S6 外踝臂（浅橙，对齐 SER-II） ----
    cov_y = 11.0
    cov_w, cov_h = 26.0, 7.2
    cov_x_center = stage_xs[1] + stage_w / 2
    cov_x = cov_x_center - cov_w / 2

    # 引线 SER-II 盒底 → 浅橙盒顶
    arr = FancyArrowPatch(
        (cov_x_center, stage_y_bot),
        (cov_x_center, cov_y + cov_h),
        arrowstyle="-|>", mutation_scale=14, linewidth=1.2,
        color="#0a0a0a", shrinkA=0, shrinkB=1, zorder=4,
    )
    ax.add_patch(arr)

    cov_box = FancyBboxPatch(
        (cov_x, cov_y), cov_w, cov_h,
        boxstyle="round,pad=0.3,rounding_size=0.8",
        linewidth=1.0, edgecolor=C_COVER_EDGE, facecolor=C_COVER_FILL,
        zorder=2,
    )
    ax.add_patch(cov_box)
    ax.text(
        cov_x + cov_w / 2, cov_y + cov_h - 2.3,
        L("本文建模覆盖", "Modeled here"),
        ha="center", va="center", fontsize=9.2, fontweight="bold",
        color="#101010", zorder=3,
    )
    ax.text(
        cov_x + cov_w / 2, cov_y + cov_h / 2 - 1.6,
        L("= 外踝臂：外旋 → 腓骨远端越阈\n（S6，§3.4 + 图 7）",
          "= lateral arm: ER → distal fibula\nover-threshold (S6, §3.4 + Fig.7)"),
        ha="center", va="center", fontsize=8.0, color="#222222",
        zorder=3,
    )

    # ---- 覆盖度注解 2：文献支撑 = 后踝臂 + 内踝臂（浅灰，对齐 SER-III+IV） ----
    lit_y = 11.0
    lit_w, lit_h = 56.0, 7.2
    lit_x_center = (stage_xs[2] + stage_w / 2 + stage_xs[3] + stage_w / 2) / 2
    lit_x = lit_x_center - lit_w / 2

    # 引线 SER-III 盒底 + SER-IV 盒底 → 浅灰盒顶
    arr = FancyArrowPatch(
        (stage_xs[2] + stage_w * 0.72, stage_y_bot),
        (lit_x + lit_w * 0.20, lit_y + lit_h),
        arrowstyle="-|>", mutation_scale=14, linewidth=1.2,
        color="#0a0a0a", shrinkA=0, shrinkB=1, zorder=4,
    )
    ax.add_patch(arr)
    arr = FancyArrowPatch(
        (stage_xs[3] + stage_w * 0.28, stage_y_bot),
        (lit_x + lit_w * 0.80, lit_y + lit_h),
        arrowstyle="-|>", mutation_scale=14, linewidth=1.2,
        color="#0a0a0a", shrinkA=0, shrinkB=1, zorder=4,
    )
    ax.add_patch(arr)

    lit_box = FancyBboxPatch(
        (lit_x, lit_y), lit_w, lit_h,
        boxstyle="round,pad=0.3,rounding_size=0.8",
        linewidth=0.9, edgecolor=C_LIT_EDGE, facecolor=C_LIT_FILL,
        zorder=2,
    )
    ax.add_patch(lit_box)
    ax.text(
        lit_x + lit_w / 2, lit_y + lit_h - 2.3,
        L("后踝臂 / 内踝臂：文献支撑，未建模",
          "Posterior / medial arms: literature only, not modeled"),
        ha="center", va="center", fontsize=9.2, fontweight="bold",
        color="#101010", zorder=3,
    )
    ax.text(
        lit_x + lit_w / 2, lit_y + lit_h / 2 - 1.7,
        L("[26] Haraguchi & Armiger 2020 + [27] Zhang 等 2022（详见各阶段下方）",
          "[26] Haraguchi & Armiger 2020 + [27] Zhang et al. 2022 (see below)"),
        ha="center", va="center", fontsize=7.6, color="#444444",
        zorder=3,
    )

    # ---- 底部图题 ----
    title_y = 4.5
    ax.text(
        60.0, title_y,
        L("图 8  三踝骨折链式机制四阶段示意（Lauge-Hansen 旋后-外旋 SER；外踝臂由 §3.4 / 图 7 覆盖）",
          "Fig.8  Trimalleolar chain mechanism (Lauge-Hansen SER; "
          "lateral arm covered by §3.4 / Fig.7)"),
        ha="center", va="center", fontsize=10.8, fontweight="bold",
        color="#101010", zorder=3,
    )

    # ---- 底部脚注 ----
    note_w, note_h = 116.0, 3.5
    note_x, note_y = 2.0, 0.3
    note_box = FancyBboxPatch(
        (note_x, note_y), note_w, note_h,
        boxstyle="round,pad=0.3,rounding_size=0.6",
        linewidth=0.7, edgecolor=C_NOTE_EDGE, facecolor=C_NOTE_FILL,
        zorder=2,
    )
    ax.add_patch(note_box)
    ax.text(
        note_x + note_w / 2, note_y + note_h / 2,
        L("示意图（非数据图）—— 用于说明三踝骨折的旋转主导机制；定量结果见图 4（旋后 / R1，冠状面）与图 7（外旋 / S6，横截面）。",
          "Schematic only (not a data figure) — illustrates the rotation-driven mechanism "
          "of trimalleolar fractures; quantitative results: Fig.4 (R1 / supination, "
          "coronal) and Fig.7 (S6 / external rotation, transverse)."),
        ha="center", va="center", fontsize=8.0, color="#444444", style="italic",
        zorder=3,
    )

    fig.tight_layout(pad=0.4)
    out = FIGDIR / "fig8_trimalleolar_ser.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[ok] {out.name}  ({out.stat().st_size} bytes)")


# =========================================================================== #
# fig9 — 距骨骨折机制（矢状面示意图：胫–距–跟载荷柱 + 强制背屈 + 前唇楔入）
# =========================================================================== #
def fig9_talus_mechanism() -> None:
    r"""论文图 9：抱石落地（强制背屈 + 轴向撞击）致距骨颈骨折的力学机制示意图。

    与图 1 同款视觉语言（圆角浅填色盒 + FancyArrow + 顶部横幅 + 底部脚注）；
    左侧为矢状面踝足骨块几何（胫骨 / 距骨 / 跟骨 + 前足背屈 + 胫骨穹前唇楔入距骨颈），
    右侧为三条结论盒（悬臂 / 阈值 / 软垫），左下为文献来源框，最底部为示意图脚注。
    """
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(11.6, 6.8))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 70)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- 配色（与图 1 同款三段：阶段蓝 / 判据橙 / 灰横幅） ----
    C_TIB_FILL = "#eaf2fb"        # 胫骨 浅蓝
    C_TAL_FILL = "#fdebd3"        # 距骨 浅橙（焦点）
    C_CAL_FILL = "#f4f4f4"        # 跟骨 浅灰
    C_EDGE = "#0a0a0a"
    C_WEDGE = "#d62728"           # 楔入 / 骨折 红色高亮
    C_LOAD = "#0a0a0a"
    C_GROUND = "#444444"
    C_BANNER_FILL = "#f4f4f4"
    C_BANNER_EDGE = "#7a7a7a"
    C_NOTE_FILL = "#fbfbfb"
    C_NOTE_EDGE = "#bdbdbd"
    C_ANN_EDGE = "#7a7a7a"

    # ============ 顶部横幅：图题 ============
    title_box = FancyBboxPatch(
        (2, 64.0), 96, 5.2,
        boxstyle="round,pad=0.3,rounding_size=0.8",
        linewidth=0.8, edgecolor=C_BANNER_EDGE, facecolor=C_BANNER_FILL, zorder=2,
    )
    ax.add_patch(title_box)
    ax.text(
        50, 66.6,
        L(
            "图 9  距骨骨折机制（矢状面）：胫–距–跟载荷柱 + 强制背屈 + 胫骨穹前唇楔入距骨颈",
            "Fig.9  Talus fracture mechanism (sagittal): tibia-talus-calcaneus load column + forced dorsiflexion + anterior tibial lip wedges talar neck",
        ),
        ha="center", va="center", fontsize=11, fontweight="bold", color="#101010", zorder=3,
    )

    # ============ 左侧：矢状面踝足骨块示意 ============
    # 坐标：
    #   地面 y=20
    #   跟骨 y=20..30（梯形：后端宽）
    #   距骨 y=30..38（梯形带颈缩：体大、颈缩、头圆）
    #   胫骨 y=38..56（长方形）
    #   前唇 y=38..41.5（胫骨底前缘凸出）
    #   载荷箭头 y=57..63
    #   前足自距骨头 (33, 37) 向上倾斜至 (48, 47)，背屈

    # 1. 地面（粗实线 + 阴影短线）
    ax.add_line(Line2D([4, 50], [20, 20], color=C_GROUND, lw=2.2, zorder=2))
    for x in range(6, 50, 3):
        ax.add_line(Line2D([x, x - 1.6], [20, 19.0], color=C_GROUND, lw=0.5, alpha=0.5, zorder=1))

    # 2. 跟骨（梯形）
    cal_pts = np.array([[10, 20], [32, 20], [30, 30], [12, 30]])
    ax.add_patch(Polygon(cal_pts, closed=True, facecolor=C_CAL_FILL, edgecolor=C_EDGE, lw=1.2, zorder=3))
    ax.text(11.5, 25, L("跟骨", "Calcaneus"),
            ha="left", va="center", fontsize=10, fontweight="bold", zorder=4)

    # 3. 距骨（梯形 + 颈缩）
    tal_pts = np.array([
        [12, 30], [30, 30],   # 底（与跟骨接合面）
        [33, 36],              # 前下（头）
        [29, 38],              # 前上（头顶/颈前）
        [16, 38],              # 后上（体顶）
        [14, 34],              # 后中（体侧）
    ])
    ax.add_patch(Polygon(tal_pts, closed=True, facecolor=C_TAL_FILL, edgecolor=C_EDGE, lw=1.4, zorder=4))
    # 距骨子区标签（体 / 颈 / 头）
    ax.text(15.5, 33.8, L("体", "Body"), ha="left", va="center", fontsize=7.2, color="#7f3d00", zorder=5)
    ax.text(22.5, 37.3, L("颈", "Neck"), ha="center", va="center", fontsize=7.5, fontweight="bold", color="#a04000", zorder=5)
    ax.text(30.0, 36.0, L("头", "Head"), ha="left", va="bottom", fontsize=7.2, color="#7f3d00", zorder=5)

    # 4. 胫骨（长方形）
    ax.add_patch(mpatches.Rectangle((18, 38), 10, 18, facecolor=C_TIB_FILL, edgecolor=C_EDGE, lw=1.2, zorder=3))
    ax.text(23, 47, L("胫骨", "Tibia"),
            ha="center", va="center", fontsize=11, fontweight="bold", zorder=4)

    # 5. 胫骨穹前唇（凸出）
    lip_pts = np.array([[28, 38], [33, 38], [33, 41.5], [28, 41.5]])
    ax.add_patch(Polygon(lip_pts, closed=True, facecolor=C_TIB_FILL, edgecolor=C_EDGE, lw=1.2, zorder=4))
    ax.text(34, 39.5, L("胫骨穹\n前唇", "Anterior tibial lip"),
            ha="left", va="center", fontsize=8, color="#333333", zorder=5)

    # 6. 楔入红色阴影 + 红色箭头
    wedge_pts = np.array([
        [28, 41.5],   # 前唇下左
        [33, 41.5],   # 前唇下右
        [27, 38],     # 距骨颈前上
        [22, 36],     # 距骨颈中部
    ])
    ax.add_patch(Polygon(wedge_pts, closed=True, facecolor=C_WEDGE, alpha=0.30,
                         edgecolor=C_WEDGE, lw=1.0, zorder=5))
    ax.add_patch(FancyArrowPatch(
        (31, 42.5), (23.5, 36.5),
        arrowstyle="-|>", mutation_scale=16, linewidth=1.6,
        color=C_WEDGE, zorder=6, shrinkA=0, shrinkB=1,
    ))
    ax.text(35, 43.8,
            L("前唇楔入距骨颈", "Plafond lip wedges talar neck"),
            ha="left", va="center", fontsize=8, color=C_WEDGE, fontweight="bold", zorder=6)

    # 7. 骨折符号（距骨颈处，红 ✕）
    fx_x, fx_y = 22.5, 35.5
    ax.add_line(Line2D([fx_x - 1.5, fx_x + 1.5], [fx_y - 1.5, fx_y + 1.5], color=C_WEDGE, lw=2.6, zorder=7))
    ax.add_line(Line2D([fx_x - 1.5, fx_x + 1.5], [fx_y + 1.5, fx_y - 1.5], color=C_WEDGE, lw=2.6, zorder=7))
    ax.text(2.5, 35.5,
            L("★ 距骨颈骨折\n（飞行员距骨）",
              "★ Talar neck fracture\n(aviator's astragalus)"),
            ha="left", va="center", fontsize=9, color=C_WEDGE, fontweight="bold", zorder=7)
    ax.add_line(Line2D([10.5, fx_x - 1.5], [35.5, fx_y], color=C_WEDGE, lw=0.8, ls="--", zorder=6))

    # 8. 前足（跖骨）：自距骨头向上倾斜（背屈）
    meta_x0, meta_y0 = 33, 37
    meta_x1, meta_y1 = 48, 47
    ax.add_line(Line2D([meta_x0, meta_x1], [meta_y0, meta_y1], color=C_EDGE, lw=2.2, zorder=3))
    for i in range(4):
        t = i / 3
        cx = meta_x0 + (meta_x1 - meta_x0) * t
        cy = meta_y0 + (meta_y1 - meta_y0) * t
        ax.add_patch(mpatches.Circle((cx, cy), 1.0, facecolor="#ffffff",
                                     edgecolor=C_EDGE, lw=0.9, zorder=4))
    # 脚趾
    ax.add_patch(mpatches.Circle((meta_x1 + 1.6, meta_y1 + 0.6), 1.1,
                                 facecolor="#ffffff", edgecolor=C_EDGE, lw=0.9, zorder=4))
    ax.text(50, 48, L("前足\n（背屈）", "Forefoot\n(dorsiflex)"),
            ha="left", va="bottom", fontsize=8, color="#333333", zorder=5)

    # 9. 背屈角弧（踝关节处：胫骨轴向上 vs 足轴指向足尖）
    arc_c = np.array([22, 38])  # 踝关节中心
    foot_dir = np.array([meta_x1 - arc_c[0], meta_y1 - arc_c[1]])
    foot_dir_n = foot_dir / np.linalg.norm(foot_dir)
    foot_angle = np.degrees(np.arctan2(foot_dir_n[1], foot_dir_n[0]))   # 足轴（向上偏右）
    tib_angle = 90.0                                                  # 胫骨轴（向上）
    theta_lo, theta_hi = sorted([foot_angle, tib_angle])
    ax.add_patch(mpatches.Arc(arc_c, 7.0, 7.0, angle=0,
                              theta1=theta_lo, theta2=theta_hi,
                              color=C_EDGE, lw=1.0, zorder=6))
    ax.text(17.0, 41.5, L("背屈角 β", "Dorsiflex β"),
            ha="center", va="center", fontsize=8, fontweight="bold", color="#7f2704", zorder=6)

    # 10. 载荷柱箭头（自上向下，进入胫骨顶部）
    # 顶部横杠
    ax.add_line(Line2D([18, 28], [63, 63], color=C_LOAD, lw=2.0, zorder=5))
    for cx in [20, 23, 26]:
        ax.add_patch(FancyArrowPatch(
            (cx, 62.8), (cx, 57.2),
            arrowstyle="-|>", mutation_scale=14, linewidth=1.3,
            color=C_LOAD, zorder=5,
        ))
    ax.text(35, 60.5,
            L("轴向载荷 F\n≈ 6–10 kN 骨折阈值", "Axial load F\n≈ 6–10 kN fracture threshold"),
            ha="left", va="center", fontsize=9, fontweight="bold", zorder=5)

    # 11. 胫–距–跟 载荷柱 左侧标签盒
    ax.text(7, 46,
            L("胫–距–跟\n载荷柱", "Tibia–Talus–\nCalcaneus\nload column"),
            ha="center", va="center", fontsize=8.5, fontweight="bold", color="#444444", zorder=5,
            bbox=dict(boxstyle="round,pad=0.35", fc="#fbfbfb", ec="#bbbbbb", lw=0.7))

    # ============ 右侧：3 个结论盒 + 侧注 ============
    ann_x, ann_w = 52, 46

    # Box 1 — 悬臂（最关键）
    b1_y, b1_h = 47, 15
    ax.add_patch(FancyBboxPatch(
        (ann_x, b1_y), ann_w, b1_h,
        boxstyle="round,pad=0.4,rounding_size=1.2",
        linewidth=1.0, edgecolor=C_ANN_EDGE, facecolor="#fff5e6", zorder=2,
    ))
    ax.text(ann_x + ann_w / 2, b1_y + b1_h - 2.2,
            L("① 距骨 = 胫–距–跟载荷柱的唯一结构悬臂",
              "Talus = sole structural cantilever of the tibia-talus-calcaneus load column"),
            ha="center", va="center", fontsize=10.5, fontweight="bold", color="#7f3d00", zorder=3)
    ax.text(ann_x + ann_w / 2, b1_y + 4.6,
            L("载荷必经胫骨 → 距骨颈 → 跟骨；\n"
              "距骨颈（松质弱区）是悬臂根部，无对偶支撑。",
              "Load must pass tibia → talar neck → calcaneus;\n"
              "talar neck (trabecular weak zone) = cantilever root, unbraced."),
            ha="center", va="center", fontsize=8.5, color="#333333", zorder=3)

    # Box 2 — 阈值
    b2_y, b2_h = 30.5, 14.5
    ax.add_patch(FancyBboxPatch(
        (ann_x, b2_y), ann_w, b2_h,
        boxstyle="round,pad=0.4,rounding_size=1.2",
        linewidth=1.0, edgecolor=C_ANN_EDGE, facecolor="#ffe9e9", zorder=2,
    ))
    ax.text(ann_x + ann_w / 2, b2_y + b2_h - 2.2,
            L("② 骨折阈值 ≈ 6–10 kN   [14, 25, 31]",
              "Fracture threshold ≈ 6–10 kN   [14, 25, 31]"),
            ha="center", va="center", fontsize=10.5, fontweight="bold", color=C_WEDGE, zorder=3)
    ax.text(ann_x + ann_w / 2, b2_y + 4.6,
            L("抱石落地峰值力常落入或越过此带；\n"
              "距骨颈为剪切屈服点（图 7 / 8 腓骨链互补）。",
              "Bouldering peak force often reaches/exceeds this band;\n"
              "talar neck = shear-yield point (complement to fibula chain in figs 7 / 8)."),
            ha="center", va="center", fontsize=8.5, color="#333333", zorder=3)

    # Box 3 — 软垫
    b3_y, b3_h = 14, 14.5
    ax.add_patch(FancyBboxPatch(
        (ann_x, b3_y), ann_w, b3_h,
        boxstyle="round,pad=0.4,rounding_size=1.2",
        linewidth=1.0, edgecolor=C_ANN_EDGE, facecolor="#fbfbfb", zorder=2,
    ))
    ax.text(ann_x + ann_w / 2, b3_y + b3_h - 2.2,
            L("③ 软垫不能屏蔽距骨",
              "Pads do not shield talus"),
            ha="center", va="center", fontsize=10.5, fontweight="bold", color="#101010", zorder=3)
    ax.text(ann_x + ann_w / 2, b3_y + 4.6,
            L("软垫只削低峰值力，不改变足姿（背屈）与载荷路径 →\n"
              "载荷柱几何不变，距骨颈受力比例不变。",
              "Pads only lower peak force; foot posture (dorsiflex) and load path unchanged\n"
              "→ column geometry unchanged, talar neck load ratio unchanged."),
            ha="center", va="center", fontsize=8.5, color="#333333", zorder=3)

    # 侧注（Box 3 下方）
    side_y, side_h = 7, 6
    ax.add_patch(FancyBboxPatch(
        (ann_x, side_y), ann_w, side_h,
        boxstyle="round,pad=0.3,rounding_size=0.8",
        linewidth=0.8, edgecolor="#bbbbbb", facecolor="#fbfbfb", zorder=2,
    ))
    ax.text(ann_x + ann_w / 2, side_y + side_h - 1.6,
            L("侧注：其他距骨骨折", "Side note: other talar fractures"),
            ha="center", va="center", fontsize=8.5, fontweight="bold", color="#444444", zorder=3)
    ax.text(ann_x + ann_w / 2, side_y + 1.6,
            L("距骨体（背屈 + 旋转）≈ 高能坠落；   距骨后突（跖屈）≈ 撞击 / 蹬地。",
              "Talar body (dorsiflex + rotation) ≈ high-energy fall;   "
              "Posterior process (plantarflex) ≈ kickback."),
            ha="center", va="center", fontsize=7.5, color="#444444", zorder=3, style="italic")

    # ============ 左下：文献来源 ============
    src_box = FancyBboxPatch(
        (2, 7), 47, 8,
        boxstyle="round,pad=0.3,rounding_size=0.7",
        linewidth=0.8, edgecolor=C_NOTE_EDGE, facecolor=C_NOTE_FILL, zorder=2,
    )
    ax.add_patch(src_box)
    ax.text(3, 13,
            L("文献来源：", "Sources:"),
            ha="left", va="center", fontsize=8.5, fontweight="bold", color="#101010", zorder=3)
    ax.text(3, 9.8,
            L(
                "[28, 29]  Hawkins 1970 / Canale & Kelly 1978（距骨颈骨折机制）\n"
                "[30]      Peterson 1976（实验：距骨颈 = 载荷柱悬臂）\n"
                "[31]      Wong 2016（FE：距骨剪切屈服）",
                "[28, 29]  Hawkins 1970 / Canale & Kelly 1978 (talar neck fracture mechanism)\n"
                "[30]      Peterson 1976 (exp: talar neck = load-column cantilever)\n"
                "[31]      Wong 2016 (FE: talar shear yielding)",
            ),
            ha="left", va="center", fontsize=7.6, color="#333333", zorder=3)

    # ============ 底部脚注 ============
    ax.text(
        50, 2.5,
        L(
            "示意图（非数据图），用于说明距骨骨折的背屈–轴向撞击机制；"
            "无实测 / 建模数值，仅几何与物理路径示意。",
            "Schematic (non-data); illustrates the dorsiflexion + axial impaction mechanism of talar fracture. "
            "No measured / modeled values — geometry and physical path only.",
        ),
        ha="center", va="bottom", fontsize=8.5, color="#444444", zorder=3, style="italic",
    )

    fig.tight_layout(pad=0.4)
    out = FIGDIR / "fig9_talus_mechanism.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"[ok] {out.name}  ({out.stat().st_size} bytes)")


# =========================================================================== #
def main() -> None:
    print(f"[font] route = {FONT_ROUTE}")
    fig1_pipeline()
    fig2_alignment()
    fig3_single_foot()
    fig4_supination()
    fig5_pad_tradeoff()
    fig6_ligament()
    fig7_external_rotation()
    fig8_trimalleolar_ser()
    fig9_talus_mechanism()
    print(f"[done] 9 figures -> {FIGDIR}")


if __name__ == "__main__":
    main()
