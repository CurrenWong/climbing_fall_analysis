"""S2.2 · 画 σ_vm(h)（刚性地面），标注跟骨阈值与首次骨折高度。

读 ``results/opensim_fe/s2_height_sweep_cl{3,4}.json``，输出 ``s2_sigma_h.png``。::

    python scripts/opensim_fe/s2_plot.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Windows 中文字体，避免标题/图例出现方块
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "opensim_fe"
THRESHOLD_MPA = 150.0


def main() -> int:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    plotted = 0
    for tag, style in (("_cl4", "-"), ("_cl3", "--")):
        p = RES / f"s2_height_sweep{tag}.json"
        if not p.is_file():
            continue
        rows = json.loads(p.read_text(encoding="utf-8"))
        label = tag.lstrip("_")
        h = [r["h"] for r in rows]
        ax.plot(h, [r["max"] for r in rows], style, marker="o", ms=3.5, label=f"max  ({label})")
        ax.plot(h, [r["p99"] for r in rows], style, marker="s", ms=3, alpha=0.6, label=f"p99  ({label})")
        ax.plot(h, [r["p95"] for r in rows], style, marker="^", ms=3, alpha=0.45, label=f"p95  ({label})")
        if all("gauge_max" in r for r in rows):
            ax.plot(h, [r["gauge_max"] for r in rows], style, marker="d", ms=3,
                    alpha=0.55, label=f"gauge(R4)  ({label})")
        plotted += 1
    if not plotted:
        print("没有可画的 sweep JSON；先跑 s2_height_sweep.py")
        return 1

    ax.axhline(THRESHOLD_MPA, color="r", lw=1.2, ls=":", label=f"跟骨阈值 {THRESHOLD_MPA:.0f} MPa")
    ax.set_xlabel("drop height h (m)")
    ax.set_ylabel("von Mises (MPa)")
    ax.set_title("S2 · 跟骨 σ_vm(h) · 刚性地面（面压 p=F/A）")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)

    out = RES / "s2_sigma_h.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
