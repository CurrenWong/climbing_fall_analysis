"""Figure: the "pad paradox" for a supinated (inverted) foot — R4 data.

Plots total ankle risk vs pad equivalent stiffness k for four drop heights,
from results/opensim_fe/pad_supination_tradeoff.json (sweep.cells), and marks
the k where risk = 1 (threshold_by_height) and the rigid limit.

Message: total_risk(k) is monotonically DECREASING -> softer pad = more danger;
there is no interior optimum. Needs k >= ~1e5 N/mm to be safe.

Usage: python scripts/opensim_fe/make_pad_paradox_figure.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "results" / "opensim_fe" / "pad_supination_tradeoff.json"
OUT = ROOT / "results" / "opensim_fe" / "figures" / "fig_pad_paradox.png"


def main() -> None:
    j = json.loads(SRC.read_text(encoding="utf-8"))
    cells = j["sweep"]["cells"]
    th = {r["height_m"]: r["k_threshold_n_per_mm"]
          for r in j["sweep"]["threshold_by_height"]}
    heights = sorted({c["height_m"] for c in cells})
    colors = {2.0: "#1f77b4", 3.0: "#2ca02c", 3.5: "#ff7f0e", 4.5: "#c0392b"}

    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=150)
    rigid_pts = {}
    for h in heights:
        cs = [c for c in cells if c["height_m"] == h]
        pts = sorted([(c["k_n_per_mm"], c["total_risk"]) for c in cs
                      if c["k_n_per_mm"] is not None])
        k = np.array([p[0] for p in pts], float)
        r = np.array([p[1] for p in pts], float)
        ax.plot(k, r, "o-", ms=4.5, lw=1.9, color=colors[h], label=f"h = {h:g} m")
        rig = [c["total_risk"] for c in cs if c["k_n_per_mm"] is None]
        if rig:
            rigid_pts[h] = rig[0]
            ax.plot([3.5e6], rig, "s", ms=8, color=colors[h],
                    mec="k", mew=0.6, zorder=6)
        if th.get(h):
            ax.axvline(th[h], color=colors[h], ls=":", lw=1.1, alpha=0.6)

    # rigid asymptote: one stacked label, no per-point collision
    if rigid_pts:
        txt = "   ".join(f"h={h:g}: {rigid_pts[h]:.2f}" for h in sorted(rigid_pts))
        ax.text(3.2e6, 1.55, "刚性极限\n" + txt, fontsize=8.2, color="#14532d",
                ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", fc="#f2f7f2", ec="#9bbf9b", lw=0.7))

    ax.axhline(1.0, color="#333", lw=1.5, ls="--")
    ax.text(55, 1.12, "risk = 1  （越阈 → 判为骨折）", fontsize=9, color="#333")

    # annotate the two extreme regimes
    ax.annotate("软垫（k 小）\n足可继续内翻 → 风险 26–45×",
                xy=(70, 38), xytext=(90, 56), fontsize=9.5, color="#7a2718",
                arrowprops=dict(arrowstyle="->", color="#7a2718", lw=1.2))
    ax.annotate("刚性极限（k→∞）\n风险 0.66–0.78（= 纯轴向）",
                xy=(1.2e6, 0.72), xytext=(6.0e3, 0.30), fontsize=9.5,
                color="#14532d",
                arrowprops=dict(arrowstyle="->", color="#14532d", lw=1.2))
    ax.text(1.0e5, 6.5, "阈值 k_th ≈ 1.0–2.4×10⁵ N/mm\n（需比这更硬才安全 → 真实垫不可达）",
            fontsize=9, color="#333", rotation=0)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(40, 6e6)
    ax.set_ylim(0.2, 80)
    ax.set_xlabel("垫子等效刚度 k  (N/mm)　——→　越右越硬", fontsize=11)
    ax.set_ylabel("踝部总风险  total_risk = σ / σ$_c$", fontsize=11)
    ax.set_title("「垫子悖论」：旋后（内翻）着地时，垫子越软，踝越危险\n"
                 "（数据：pad_supination_tradeoff.json · R4；■ = 刚性极限，点线 = risk=1 阈值）",
                 fontsize=11.5)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=9, loc="upper right")
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    print("wrote", OUT)
    for h in heights:
        print(f"  h={h:g}m  k_th={th.get(h)}  rigid_risk={rigid_pts.get(h):.3f}")


if __name__ == "__main__":
    main()
