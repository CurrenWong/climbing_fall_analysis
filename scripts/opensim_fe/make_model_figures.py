"""make_model_figures.py -- render the ankle mechanism figures (fig8/fig9)
FROM THE REAL FE GEOMETRY (THUMS AM50 right lower limb, global mm frame).

Read-only over existing meshes; writes only to results/opensim_fe/figures/.

fig8 : 三踝骨折的旋转机制 -- ankle mortise (tibia + fibula + talus), distal
        fibula / lateral malleolus highlighted; external-rotation arrow.
fig9 : 距骨骨折的背屈-轴向撞击机制 -- load column (tibia + talus + calcaneus);
        axial load arrows, talar-neck weak point highlighted.

Chinese labels are composited with matplotlib (Microsoft YaHei) on top of the
pyvista offscreen render, so we keep full CJK control.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pyvista as pv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "temp" / "opensim_fe"
FIG = ROOT / "results" / "opensim_fe" / "figures"

BONES = {
    "tibia": BASE / "fracture_matrix_remesh/parts_extract/tibia_r_union/tibia_r_union_surface.stl",
    "fibula": BASE / "fracture_matrix_remesh/parts_extract/fibula_r_union/fibula_r_union_surface.stl",
    "talus": BASE / "parts/talus_r/talus_r_surface.stl",
    "calcaneus": BASE / "fracture_matrix_remesh/parts_extract/calcaneus_r_union/calcaneus_r_union_surface.stl",
}
BONE_C = {
    "tibia": (0.905, 0.882, 0.827),
    "fibula": (0.882, 0.855, 0.792),
    "talus": (0.855, 0.831, 0.769),
    "calcaneus": (0.827, 0.804, 0.741),
}
HOT = (0.86, 0.20, 0.16)
INK = "#3a2f22"


def cjk_font() -> str:
    for name in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC"):
        try:
            return font_manager.findfont(font_manager.FontProperties(family=name), fallback_to_default=False)
        except Exception:
            continue
    return font_manager.findfont(font_manager.FontProperties())


def load(name: str) -> pv.PolyData:
    return pv.read(str(BONES[name]))


def clip_distal(mesh: pv.PolyData, zmax: float) -> pv.PolyData:
    """Keep only the distal part (z <= zmax) for a clean ankle close-up."""
    return mesh.clip(normal=(0.0, 0.0, -1.0), origin=(0.0, 0.0, zmax), crinkle=False)


def base_plotter(size, bg="white"):
    pl = pv.Plotter(off_screen=True, window_size=size)
    pl.set_background(bg)
    return pl


def hot_blob(center, radius) -> pv.PolyData:
    return pv.Sphere(radius=radius, center=center, theta_resolution=48, phi_resolution=48)


def fig8() -> None:
    tib = clip_distal(load("tibia"), -60.0)
    fib = clip_distal(load("fibula"), -60.0)
    tal = load("talus")

    fib_pts = np.asarray(fib.points)
    distal = fib_pts[fib_pts[:, 2].argsort()][: max(8, len(fib_pts) // 25)].mean(axis=0)
    tal_c = np.asarray(tal.points).mean(axis=0)

    pl = base_plotter((2200, 1700))
    for k, m in (("tibia", tib), ("fibula", fib), ("talus", tal)):
        pl.add_mesh(m, color=BONE_C[k], smooth_shading=True, specular=0.25, ambient=0.40, diffuse=0.78)
    pl.add_mesh(hot_blob(distal, 22.0), color=HOT, opacity=0.42, smooth_shading=True)
    pl.camera.focal_point = (distal + tal_c) / 2.0
    pl.camera_position = "xy"          # start from a standard view
    pl.view_isometric()
    pl.camera.focal_point = (distal + tal_c) / 2.0
    pl.camera.zoom(1.9)
    img = pl.screenshot(return_img=True)
    pl.close()

    fig, ax = plt.subplots(figsize=(img.shape[1] / 200, img.shape[0] / 200), dpi=200)
    ax.imshow(img)
    ax.axis("off")
    fp = font_manager.FontProperties(fname=cjk_font(), size=17, weight="bold")
    fp2 = font_manager.FontProperties(fname=cjk_font(), size=14)
    ax.text(0.5, 0.965, "图 8  三踝骨折的旋转机制（真实 FE 几何：胫骨 + 腓骨 + 距骨）",
            ha="center", va="top", transform=ax.transAxes, fontproperties=fp, color=INK)
    # external-rotation arrow (curved) about the vertical axis
    ax.annotate("", xy=(0.20, 0.34), xytext=(0.30, 0.60), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#b3261e", lw=3.2,
                                connectionstyle="arc3,rad=0.45"))
    ax.text(0.11, 0.60, "脚绕胫骨长轴\n外旋 θ", transform=ax.transAxes, fontproperties=fp2,
            color="#b3261e", ha="left", va="center")
    ax.text(0.62, 0.30, "外踝（腓骨远端）\n斜形骨折热点", transform=ax.transAxes, fontproperties=fp2,
            color="#b3261e", ha="left", va="center")
    ax.plot([0.55, 0.70], [0.40, 0.46], color="#b3261e", lw=1.6, ls="--", transform=ax.transAxes)
    ax.text(0.72, 0.72, "胫骨", transform=ax.transAxes, fontproperties=fp2, color=INK)
    ax.text(0.30, 0.62, "腓骨", transform=ax.transAxes, fontproperties=fp2, color=INK)
    ax.text(0.60, 0.18, "距骨", transform=ax.transAxes, fontproperties=fp2, color=INK)
    ax.text(0.02, 0.03, "真实骨骼表面网格（THUMS AM50 右下肢，global mm）；红=外踝失效热点（示意，非应力场）。",
            transform=ax.transAxes, fontproperties=font_manager.FontProperties(fname=cjk_font(), size=10),
            color="#6b5d4a", ha="left", va="bottom")
    out = FIG / "fig8_trimalleolar_ser.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out, img.shape)


def fig9() -> None:
    tib = clip_distal(load("tibia"), -40.0)
    tal = load("talus")
    cal = load("calcaneus")

    tal_c = np.asarray(tal.points).mean(axis=0)
    tib_d = np.asarray(tib.points)
    tib_low = tib_d[tib_d[:, 2].argsort()][: max(8, len(tib_d) // 20)].mean(axis=0)
    # talar neck approx = anterior-superior of talus
    tal_pts = np.asarray(tal.points)
    neck = tal_pts[tal_pts[:, 2].argsort()][: max(8, len(tal_pts) // 12)].mean(axis=0)

    pl = base_plotter((1900, 1900))
    for k, m in (("tibia", tib), ("talus", tal), ("calcaneus", cal)):
        pl.add_mesh(m, color=BONE_C[k], smooth_shading=True, specular=0.25, ambient=0.40, diffuse=0.78)
    pl.add_mesh(hot_blob(neck, 20.0), color=HOT, opacity=0.45, smooth_shading=True)
    pl.view_xz()
    pl.camera.focal_point = (tib_low + tal_c + np.asarray(cal.points).mean(axis=0)) / 3.0
    pl.camera.zoom(1.55)
    img = pl.screenshot(return_img=True)
    pl.close()

    fig, ax = plt.subplots(figsize=(img.shape[1] / 200, img.shape[0] / 200), dpi=200)
    ax.imshow(img)
    ax.axis("off")
    fp = font_manager.FontProperties(fname=cjk_font(), size=18, weight="bold")
    fp2 = font_manager.FontProperties(fname=cjk_font(), size=15)
    ax.text(0.5, 0.97, "图 9  距骨骨折的背屈–轴向撞击机制（真实 FE 几何：胫骨–距骨–跟骨）",
            ha="center", va="top", transform=ax.transAxes, fontproperties=fp, color=INK)
    # axial load arrows
    for x in (0.30, 0.42, 0.54):
        ax.annotate("", xy=(x, 0.60), xytext=(x, 0.80), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color="#b3261e", lw=3.0))
    ax.text(0.42, 0.84, "轴向载荷 F ≈ 6–10 kN", transform=ax.transAxes, fontproperties=fp2,
            color="#b3261e", ha="center", va="bottom")
    ax.annotate("", xy=(0.62, 0.50), xytext=(0.78, 0.62), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#1f5f8b", lw=2.6,
                                connectionstyle="arc3,rad=0.35"))
    ax.text(0.80, 0.56, "背屈", transform=ax.transAxes, fontproperties=fp2, color="#1f5f8b")
    ax.text(0.30, 0.30, "距骨颈（骨小梁突变\n= 载荷柱弱点）", transform=ax.transAxes, fontproperties=fp2,
            color="#b3261e", ha="left", va="center")
    ax.text(0.60, 0.70, "胫骨", transform=ax.transAxes, fontproperties=fp2, color=INK)
    ax.text(0.55, 0.34, "距骨", transform=ax.transAxes, fontproperties=fp2, color=INK)
    ax.text(0.55, 0.14, "跟骨", transform=ax.transAxes, fontproperties=fp2, color=INK)
    ax.text(0.02, 0.02, "真实骨骼表面网格（THUMS AM50 右下肢，global mm）；红=距骨颈失效热点（示意，非应力场）。",
            transform=ax.transAxes, fontproperties=font_manager.FontProperties(fname=cjk_font(), size=10),
            color="#6b5d4a", ha="left", va="bottom")
    out = FIG / "fig9_talus_mechanism.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out, img.shape)


if __name__ == "__main__":
    pv.OFF_SCREEN = True
    fig8()
    fig9()
