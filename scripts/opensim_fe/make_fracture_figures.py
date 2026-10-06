"""Fracture-mechanism schematics (v9 smoothing), 4 angles, two framings.

fig8 = trimalleolar fracture: red fracture lines on lateral/medial/posterior
       malleolus + external-rotation arc + axial-load arrow.
fig9 = talar fracture: a THIN red fracture line on the talus + axial-impaction
       arrow; dorsiflexion is shown as a 2D curved arrow in blank space
       (a sagittal 3D arc inevitably crosses the tibia).

3D fracture lines are cut from the real mesh (plane slice) so they sit on the
bone and occlude correctly. 2D glyphs are auto-placed in background pixels.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pyvista as pv
import vtk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "temp" / "opensim_fe"
OUT = BASE / "_views_v9_fx22"
BONES = {
    "tibia": BASE / "parts_gmsh/tibia_r/tibia_r_surface_used.stl",
    "fibula": BASE / "fracture_matrix_remesh/parts_extract/fibula_r_union/fibula_r_union_surface.stl",
    "talus": BASE / "parts/talus_r/talus_r_surface.stl",
    "calcaneus": BASE / "fracture_matrix_remesh/parts_extract/calcaneus_r_union/calcaneus_r_union_surface.stl",
    "foot": BASE / "parts/foot_r/foot_r_surface.stl",
}
COLORS = {"tibia": (0.91, 0.89, 0.83), "fibula": (0.86, 0.66, 0.57),
          "talus": (0.80, 0.86, 0.88), "calcaneus": (0.84, 0.81, 0.72),
          "foot": (0.89, 0.85, 0.76)}
ORDER = ["tibia", "fibula", "calcaneus", "foot", "talus"]
VIEWS = [("A 外侧", 0, "lateral"), ("D 前方", 90, "anterior"),
         ("G 内侧", 180, "medial"), ("J 后方", 270, "posterior")]
WIN = (1250, 1500)
FRAMINGS = [("full", 1420.0, 0.0, 0.0), ("near", 780.0, -8.0, 62.0)]
# (name, DIST, focal shift *axis, focal shift *v2) -- near shifts anteriorly so
# the WHOLE foot fits; 4-tuple.

FX = (0.78, 0.08, 0.08)
C_ROT = (0.90, 0.42, 0.05)
C_AX = (0.10, 0.33, 0.70)
C_DORSI = "#d2691e"
C_DORSI3 = (0.82, 0.41, 0.12)

TIBIA_SUBFILTER = "loop"
TIBIA_PASSES = [(200, 0.02), (150, 0.03), (100, 0.05)]
FX_R = 2.7
TAL_R = 1.4
DORSI_DEG = 24.0        # fig9: foot posed in dorsiflexion
EXROT_DEG = 20.0        # fig8: foot posed in external rotation
INV_DEG = 28.0          # fig8: foot also SUPINATED (inverted sole) -- SER


def unit(x):
    return x / np.linalg.norm(x)


def prep(path: Path, passes=None, subdiv=0) -> pv.PolyData:
    m = pv.read(str(path))
    for fn in (lambda x: x.clean(), lambda x: x.triangulate()):
        try:
            m = fn(m)
        except Exception:
            pass
    if subdiv > 0:
        try:
            m = m.subdivide(subdiv, TIBIA_SUBFILTER)
        except Exception as e:
            print("  subdivide failed", e)
    if passes:
        for n_it, pb in passes:
            try:
                m = m.smooth_taubin(n_iter=n_it, pass_band=pb,
                                    boundary_smoothing=True, feature_smoothing=True)
            except Exception as e:
                print("  taubin failed", e)
    return m


def cjk(sz):
    for nm in ("Microsoft YaHei", "SimHei"):
        try:
            return font_manager.FontProperties(
                fname=font_manager.findfont(font_manager.FontProperties(family=nm),
                                            fallback_to_default=False), size=sz)
        except Exception:
            continue
    return font_manager.FontProperties(size=sz)


def ordered_loop(mesh, origin, normal):
    try:
        s = mesh.slice(normal=unit(normal), origin=origin)
    except Exception:
        return np.empty((0, 3))
    P = np.asarray(s.points, float)
    if len(P) < 3:
        return P
    n = unit(normal)
    ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = unit(np.cross(n, ref)); e2 = unit(np.cross(n, e1))
    c = P.mean(0)
    return P[np.argsort(np.arctan2((P - c) @ e2, (P - c) @ e1))]


def offset_out(P, d=0.8):
    c = P.mean(0)
    return P + d * np.array([unit(p - c) for p in P])


def arc_around(P, key, frac=0.30):
    if len(P) < 5:
        return P
    i = int(np.argmax(key(P))); N = len(P); half = max(2, int(N * frac / 2))
    return P[[(i + k) % N for k in range(-half, half + 1)]]


def tube(pts, closed, r=FX_R):
    return pv.lines_from_points(pts, close=closed).tube(radius=r, n_sides=14)


def rot_about(pts, center, axis, deg):
    """Rodrigues rotation of point array about an axis through `center`."""
    a = np.radians(deg); k = unit(axis)
    v = np.asarray(pts, float) - center
    kk = np.broadcast_to(k, v.shape)
    out = (v * np.cos(a) + np.cross(kk, v) * np.sin(a)
           + (v @ k)[:, None] * k[None, :] * (1 - np.cos(a)))
    return center + out


def rotate_mesh(m, center, axis, deg):
    m2 = m.copy()
    m2.points = rot_about(np.asarray(m2.points, float), center, axis, deg)
    return m2


def straight_arrow(start, direction, length):
    return pv.Arrow(start=start, direction=unit(direction), scale=length,
                    tip_length=0.32, tip_radius=0.085,
                    shaft_radius=0.030, tip_resolution=24, shaft_resolution=16)


def curved_arrow(center, normal, r, a0, a1, tube_r=2.9, head=(20.0, 6.5)):
    n = unit(normal)
    ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = unit(np.cross(n, ref)); e2 = unit(np.cross(n, e1))
    t = np.linspace(a0, a1, 72)
    pts = center + r * (np.cos(t)[:, None] * e1 + np.sin(t)[:, None] * e2)
    body = pv.lines_from_points(pts).tube(radius=tube_r, n_sides=14)
    d = unit(pts[-1] - pts[-2])
    cone = pv.Cone(center=pts[-1] + d * head[0] * 0.45, direction=d,
                   height=head[0], radius=head[1], resolution=28)
    return [body, cone], pts


def w2i(renderer, pt, w, h):
    c = vtk.vtkCoordinate()
    c.SetCoordinateSystemToWorld()
    c.SetValue(float(pt[0]), float(pt[1]), float(pt[2]))
    x, y = c.GetComputedDisplayValue(renderer)
    return float(x), float(h - y)


def _blank(img, x0, y0, x1, y1):
    H, W = img.shape[:2]
    x0 = int(max(0, x0)); y0 = int(max(0, y0))
    x1 = int(min(W, x1)); y1 = int(min(H, y1))
    if x1 <= x0 or y1 <= y0:
        return False
    return int(img[y0:y1, x0:x1, :3].min()) > 242


def _ov(a, b):
    return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])


def draw_labels(ax, img, items, w, h, fp_size=20, prefer_right=()):
    """Label each item in the nearest blank patch (background pixel) around its
    anchor, never overlapping another label. Keeps leaders short and off-bone.
    `prefer_right` forces those labels to the +x side of their anchor."""
    fp = cjk(fp_size)
    th = int(fp_size * 2.1)
    placed = []
    for text, (xi, yi) in sorted(items, key=lambda it: it[1][1]):
        tw = int(th * len(text)) + 14
        best = None
        for r in range(int(tw * 0.55), int(max(w, h) * 0.72), 12):
            for adeg in range(0, 360, 10):
                a = np.radians(adeg)
                cx = xi + r * np.cos(a); cy = yi + r * np.sin(a)
                if text in prefer_right and cx < xi + 30:
                    continue
                ha = "left" if cx >= xi else "right"
                x0 = cx if ha == "left" else cx - tw
                rect = (x0 - 4, cy - th / 2 - 4, x0 + tw + 4, cy + th / 2 + 4)
                if not (4 < rect[0] and rect[2] < w - 4 and 4 < rect[1] and rect[3] < h - 4):
                    continue
                if not _blank(img, *rect):
                    continue
                if any(_ov(rect, p) for p in placed):
                    continue
                best = (ha, cx, cy, rect)
                break
            if best:
                break
        if best is None:
            ha = "right" if xi >= w / 2 else "left"
            x0 = (w - 6 - tw) if ha == "right" else 6
            cand_ys = (list(range(int(yi), int(h * 0.05), -14))
                       + list(range(int(yi), int(h * 0.95), 14))) or [int(yi)]
            yy = int(yi)
            for ytest in cand_ys:
                rect = (x0 - 4, ytest - th / 2 - 4, x0 + tw + 4, ytest + th / 2 + 4)
                if _blank(img, *rect) and not any(_ov(rect, p) for p in placed):
                    yy = ytest
                    break
            cx = (x0 + tw) if ha == "right" else x0
            best = (ha, cx, yy, (x0 - 4, yy - th / 2 - 4, x0 + tw + 4, yy + th / 2 + 4))
        ha, cx, cy, rect = best
        placed.append(rect)
        ax.annotate(text, xy=(xi, yi), xytext=(cx, cy), textcoords="data",
                    ha=ha, va="center", fontproperties=fp, color="#111111",
                    arrowprops=dict(arrowstyle="-", color="#888888",
                                    lw=1.1, shrinkA=3, shrinkB=3))
        ax.plot([xi], [yi], marker="o", ms=4.4, color="#111111",
                mec="white", mew=0.8, zorder=6)


def find_blank_spot(img, cx, cy, rmin=110, rmax=430, step=16, patch=14):
    h, w = img.shape[:2]
    for r in range(rmin, rmax, step):
        for adeg in range(0, 360, 10):
            a = np.radians(adeg)
            x = int(cx + r * np.cos(a)); y = int(cy + r * np.sin(a))
            if not (patch + 2 < x < w - patch - 2 and patch + 2 < y < h - patch - 2):
                continue
            p = img[y - patch:y + patch, x - patch:x + patch, :3]
            if p.size and int(p.min()) > 242:
                return x, y
    return int(min(max(cx, 70), w - 70)), int(min(max(cy, 70), h - 70))


def draw_arc_glyph(ax, img, target, text, fp, color=C_DORSI):
    """Curved arrow drawn in a blank area OUTSIDE the bone, dashed leader to the
    anatomical point. Used for dorsiflexion (a sagittal 3D arc would cross the
    tibia)."""
    tx, ty = target
    sx, sy = find_blank_spot(img, tx, ty, rmin=165, rmax=470, patch=36)
    u = np.array([sx - tx, sy - ty], float)
    n = float(np.linalg.norm(u))
    u = u / n if n > 1e-6 else np.array([1.0, 0.0])
    perp = np.array([-u[1], u[0]])
    c = np.array([sx, sy], float)
    arr = FancyArrowPatch(c - perp * 66.0, c + perp * 66.0,
                          connectionstyle="arc3,rad=0.5", arrowstyle="-|>",
                          mutation_scale=32, lw=4.0, color=color,
                          shrinkA=0, shrinkB=0, zorder=8)
    ax.add_patch(arr)
    ax.plot([sx, tx], [sy, ty], color=color, lw=1.3,
            linestyle=(0, (5, 3)), zorder=7)
    ax.text(sx + u[0] * 104, sy + u[1] * 104, text, fontproperties=fp,
            color=color, ha="center", va="center", zorder=8)


def main() -> None:
    pv.OFF_SCREEN = True
    M = {k: prep(p, passes=(TIBIA_PASSES if k == "tibia" else None),
                 subdiv=(1 if k == "tibia" else 0)) for k, p in BONES.items()}

    tib_c, fib_c = M["tibia"].points.mean(0), M["fibula"].points.mean(0)
    tal_c = M["talus"].points.mean(0)
    tp = M["tibia"].points - tib_c
    w, v = np.linalg.eigh(tp.T @ tp)
    axis = unit(v[:, np.argmax(w)])
    if np.dot(axis, tib_c - tal_c) < 0:
        axis = -axis
    ml = unit((fib_c - tib_c) - np.dot(fib_c - tib_c, axis) * axis)
    v2 = unit(np.cross(axis, ml))
    allp = np.vstack([m.points for m in M.values()])
    scene_c = (allp.min(0) + allp.max(0)) / 2

    def pa(P):   return (P - scene_c) @ axis
    def plat(P): return (P - scene_c) @ ml
    def pant(P): return (P - scene_c) @ v2

    T, F, TL = M["tibia"], M["fibula"], M["talus"]
    lo = pa(T.points); dmin = float(lo.min())
    tibial_tip = T.points[int(lo.argmin())]
    fib_tip = F.points[int(pa(F.points).argmin())]
    talus_pt = TL.points.mean(0)
    calc_pt = M["calcaneus"].points.mean(0)
    ankle_c = (tibial_tip + fib_tip) / 2

    # ---- fracture lines -------------------------------------------------
    fx, fx_tal = [], []
    fib_loop = ordered_loop(F, fib_tip + 26 * axis, unit(axis + 0.32 * ml))
    if len(fib_loop) > 4:
        fib_loop = offset_out(fib_loop, 1.3); fx.append(tube(fib_loop, True))
        fib_line_c = fib_loop.mean(0)
    else:
        fib_line_c = fib_tip + 26 * axis
    med_arc = arc_around(ordered_loop(T, tibial_tip + 12 * axis, axis),
                         lambda P: -plat(P), 0.36)
    if len(med_arc) > 3:
        med_arc = offset_out(med_arc, 1.3); fx.append(tube(med_arc, False))
        med_line_c = med_arc.mean(0)
    else:
        med_line_c = tibial_tip + 12 * axis
    post_arc = arc_around(ordered_loop(T, tibial_tip + 18 * axis, axis),
                          lambda P: -pant(P), 0.32)
    if len(post_arc) > 3:
        post_arc = offset_out(post_arc, 1.3); fx.append(tube(post_arc, False))
        post_line_c = post_arc.mean(0)
    else:
        post_line_c = tibial_tip + 18 * axis
    # talus: THIN red line hugging the bone. Fracture plane ~ CORONAL (normal
    # ~ v2) so it reads as a vertical fracture line from the lateral/medial view
    tal_loop = ordered_loop(TL, talus_pt, unit(v2 + 0.12 * axis))
    if len(tal_loop) > 4:
        tal_loop = offset_out(tal_loop, 1.5)
        fx_tal.append(tube(tal_loop, True, r=TAL_R))
        tal_line_c = tal_loop.mean(0)
    else:
        tal_line_c = talus_pt

    # ---- pose the foot in DORSIFLEXION (fig9 only) ----------------------
    # rotate 距骨/跟骨/足 about the dorsiflexion axis (+ml through the talus);
    # the tibia & fibula stay put. fig8 (trimalleolar / external rotation) keeps
    # the neutral pose.
    talus_pt_n, calc_pt_n, tal_line_c_n = talus_pt, calc_pt, tal_line_c
    pivot = talus_pt.copy()
    M_rot = dict(M)
    for k in ("talus", "calcaneus", "foot"):
        M_rot[k] = rotate_mesh(M[k], pivot, ml, DORSI_DEG)
    fx_tal_rot = [rotate_mesh(t, pivot, ml, DORSI_DEG) for t in fx_tal]
    talus_pt_d = M_rot["talus"].points.mean(0)
    calc_pt_d = M_rot["calcaneus"].points.mean(0)
    tal_line_c_d = rot_about(np.asarray([tal_line_c]), pivot, ml, DORSI_DEG)[0]

    # ---- pose the foot in EXTERNAL ROTATION + SUPINATION (fig8) ---------
    # SER mechanism: foot is SUPINATED (inverted sole, 脚心内翻) + externally
    # rotated. Inversion = rotation about the foot's long axis (v2); external
    # rotation = about the shank long axis. Tibia & fibula stay put.
    shank_pt = tib_c + float(np.dot(ankle_c - tib_c, axis)) * axis
    M_ext = dict(M)
    for k in ("talus", "calcaneus", "foot"):
        mm = M[k].copy()
        p = np.asarray(mm.points, float)
        p = rot_about(p, talus_pt_n, v2, INV_DEG)          # supination (inversion)
        p = rot_about(p, shank_pt, axis, -EXROT_DEG)  # external rotation
        mm.points = p
        M_ext[k] = mm
    talus_pt_e = M_ext["talus"].points.mean(0)

    # ---- mechanism arrows (fig8): 内翻 (inversion) + 外旋 (ext. rotation) --
    # 内翻: arc about the FOOT's long axis (v2) at the hindfoot; decreasing
    # angle = sole rolls medially = inversion.
    inv_c = calc_pt_n
    its = np.radians(np.linspace(345, 195, 72))
    ipts = inv_c + 52.0 * (np.cos(its)[:, None] * ml[None, :]
                           + np.sin(its)[:, None] * axis[None, :])
    i_dir = unit(ipts[-1] - ipts[-2])
    inv_arrow = [pv.lines_from_points(ipts).tube(radius=3.6, n_sides=14),
                 pv.Cone(center=ipts[-1] + i_dir * 10.0, direction=i_dir,
                         height=24.0, radius=8.0, resolution=28)]
    inv_pt = ipts[len(ipts) // 2]

    # 外旋: transverse arc in the empty space BELOW the foot; increasing angle
    # = rotation about -axis = external rotation (toe-out), head lateral.
    ex_c = ankle_c + 30 * v2 - 140 * axis
    ets = np.radians(np.linspace(-65, 65, 72))
    epts = ex_c + 62.0 * (np.cos(ets)[:, None] * v2[None, :]
                          + np.sin(ets)[:, None] * ml[None, :])
    e_dir = unit(epts[-1] - epts[-2])
    rot_arrow = [pv.lines_from_points(epts).tube(radius=3.6, n_sides=14),
                 pv.Cone(center=epts[-1] + e_dir * 10.0, direction=e_dir,
                         height=24.0, radius=8.0, resolution=28)]
    rot_pt = epts[len(epts) // 2]

    ax_arrow = straight_arrow(ankle_c + 108 * axis - 40 * v2, -axis, 74.0)
    ax_mid = ankle_c + 108 * axis - 40 * v2 - axis * 37.0

    # dorsiflexion: a REAL 3D curved arrow whose axis is +ml (the true
    # dorsiflexion axis for the right foot). Arc lives in the (v2, axis)
    # sagittal plane; increasing angle rotates v2->axis = toes up = correct.
    dcen = ankle_c + 56 * axis - 78 * ml + 46 * v2
    dtt = np.radians(np.linspace(-75, 45, 72))
    dpts = dcen + 40.0 * (np.cos(dtt)[:, None] * v2[None, :]
                          + np.sin(dtt)[:, None] * axis[None, :])
    d_dir = unit(dpts[-1] - dpts[-2])
    dorsi_arrow = [pv.lines_from_points(dpts).tube(radius=3.0, n_sides=14),
                   pv.Cone(center=dpts[-1] + d_dir * 8.0, direction=d_dir,
                           height=18.0, radius=6.0, resolution=28)]
    dorsi_pt = dpts[len(dpts) // 2]

    fig8 = [("外踝骨折", fib_line_c), ("内踝骨折", med_line_c),
            ("后踝骨折", post_line_c), ("内翻", inv_pt), ("外旋", rot_pt),
            ("轴向载荷", ax_mid), ("距骨", talus_pt_e)]
    fig9 = [("距骨骨折", tal_line_c_d), ("胫骨远端", tibial_tip + 26 * axis),
            ("跟骨", calc_pt_d), ("轴向撞击", ax_mid), ("背屈", dorsi_pt)]

    SCENES = {
        "fig8": ([(m, FX) for m in fx] + [(m, C_ROT) for m in inv_arrow] + [(m, C_ROT) for m in rot_arrow]
                 + [(ax_arrow, C_AX)]),
        "fig9": ([(m, FX) for m in fx_tal_rot] + [(ax_arrow, C_AX)]
                 + [(m, C_DORSI3) for m in dorsi_arrow]),
    }
    keys_n = {"外踝骨折": fib_line_c, "内踝骨折": med_line_c,
              "后踝骨折": post_line_c, "距骨": talus_pt_n, "外旋": rot_pt,
              "内翻": inv_pt, "轴向载荷": ax_mid, "距骨骨折": tal_line_c_n,
              "胫骨远端": tibial_tip + 26 * axis, "跟骨": calc_pt_n,
              "轴向撞击": ax_mid, "背屈": dorsi_pt}
    keys_d = dict(keys_n)
    keys_d.update({"距骨": talus_pt_d, "跟骨": calc_pt_d,
                   "距骨骨折": tal_line_c_d})
    keys_e = dict(keys_n)
    keys_e["距骨"] = talus_pt_e

    W, H = WIN
    for fz_name, dist, fshift, vshift in FRAMINGS:
        (OUT / fz_name).mkdir(parents=True, exist_ok=True)
        focal = ankle_c + fshift * axis + vshift * v2
        for tag in ("fig8", "fig9"):
            meshes = M_ext if tag == "fig8" else M_rot
            keys = keys_e if tag == "fig8" else keys_d
            labels = fig8 if tag == "fig8" else fig9
            pr = () if tag == "fig8" else ("轴向撞击",)
            imgs = {}
            for lab, az, name in VIEWS:
                a = np.radians(az)
                d = np.cos(a) * ml + np.sin(a) * v2
                pl = pv.Plotter(off_screen=True, window_size=WIN)
                pl.set_background("white")
                for k in ORDER:
                    pl.add_mesh(meshes[k], color=COLORS[k], smooth_shading=True,
                                specular=0.32, ambient=0.42, diffuse=0.8)
                for meshp, col in SCENES[tag]:
                    pl.add_mesh(meshp, color=col, smooth_shading=True)
                pl.camera.up = axis
                pl.camera.focal_point = focal
                pl.camera.position = focal + unit(d) * dist
                pl.reset_camera_clipping_range()
                img = pl.screenshot(return_img=True)
                ren = pl.renderer
                proj = {k: w2i(ren, p, W, H) for k, p in keys.items()}
                pl.close()
                imgs[name] = (img, proj)
                figp = plt.figure(figsize=(W / 100, H / 100), dpi=100)
                axp = figp.add_axes([0, 0, 1, 1]); axp.axis("off")
                axp.imshow(img)
                draw_labels(axp, img, [(t, proj[t]) for t, _ in labels], W, H, fp_size=17, prefer_right=pr)
                figp.savefig(OUT / fz_name / f"{tag}_{name}.png", dpi=100)
                plt.close(figp)
            fig, axes = plt.subplots(2, 2, figsize=(2 * W / 100, 2 * H / 100), dpi=100)
            for ax, (lab, az, name) in zip(axes.ravel(), VIEWS):
                img, proj = imgs[name]
                ax.imshow(img); ax.set_title(lab, fontproperties=cjk(16)); ax.axis("off")
                draw_labels(ax, img, [(t, proj[t]) for t, _ in labels], W, H, fp_size=17, prefer_right=pr)
            fig.tight_layout()
            fig.savefig(BASE / f"_{tag}_fx_{fz_name}_v22_sheet.png", dpi=100,
                        bbox_inches="tight")
            plt.close(fig)
            print("wrote", BASE / f"_{tag}_fx_{fz_name}_v22_sheet.png")


if __name__ == "__main__":
    main()
