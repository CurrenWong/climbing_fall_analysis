"""Render the FALL simulation: full-body OpenSim dead-drop forward dynamics.

Timeline: free fall (kinematic, v = g t) -> impact (``run_dead_drop`` forward
dynamics, real GRF from ``pad.py``). Bodies are the Rajagopal2015 .vtp geometry,
posed per frame with ``getTransformInGround()``. Frames -> ffmpeg -> MP4 (+ GIF).

Two fixes over the first version
--------------------------------
1. **Foam compresses instead of bones clipping.** The FD has *no contact
   constraint* (only a prescribed GRF), so the body genuinely travels ~16 cm
   downward during the impact. The pad top surface is now drawn at the frame's
   **lowest bone vertex** (``min(0, min_y)``) so the foam visibly compresses
   under the foot and the bones never pass through the slab.
2. **Joint-mode is labelled on screen** so the viewer does not expect flexion
   that the model is not solving for.

Modes
-----
* default    -> S1 rigid leg: every ``LEG_JOINTS`` coordinate locked
                ([Y25] "直立、绷直腿" 刚性柱近似)。
* ``--passive`` -> S3 path (b): ``LEG_FLEX_JOINTS`` unlocked + passive
                stiffness/damping (``PassiveStiffness.nominal_landing()``);
                the joints bend only ~2 deg in a 20 ms impact
                (``NONVERTICAL_S3_POSTURE_REPORT.md``), shown live as a readout.

Usage
-----
  python scripts/opensim_fe/make_fall_animation.py --height 2.0 --on-pad
  python scripts/opensim_fe/make_fall_animation.py --height 2.0 --on-pad --passive
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import opensim as os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from climbing.coupling.opensim_grf import ground_reaction          # noqa: E402
from climbing.coupling.opensim_fall import (                       # noqa: E402
    run_dead_drop,
    PassiveStiffness,
)

GEO = ROOT / "model" / "opensim" / "FullBodyModel-4.0" / "Geometry"
OUTD = ROOT / "temp" / "opensim_fe" / "_fall_frames"
FFMPEG = shutil.which("ffmpeg") or r"D:\Program\ffmpeg-8.0-full_build\bin\ffmpeg.exe"
BONE = (0.90, 0.87, 0.80)
PAD = (0.30, 0.42, 0.62)
GROUND = (0.85, 0.84, 0.80)
PAD_THICK = 0.20            # m — src/climbing/pad.py CrashPad.thickness_m
GROUND_Y = -PAD_THICK       # pad base rests on the floor at this level
WIN = (1000, 1200)

# pad footprint (matches the reference box) + dent-grid resolution
PAD_X = (-0.62, 0.62)
PAD_Z = (-0.45, 0.55)
PAD_NX, PAD_NZ = 56, 44
FOOT_R = 0.065              # m — horizontal influence radius of a foot vertex


def _is_foot(name: str) -> bool:
    return any(s in name for s in ("calcn", "talus", "toes", "foot", "calcaneus"))


def _pad_top_field(foot_pts: list) -> "np.ndarray":
    """(PAD_NX, PAD_NZ) top-surface height: ``0`` where the pad is uncompressed,
    else the lowest foot surface over that cell, clamped to the floor.

    A foot vertex only pulls the pad down within ``FOOT_R`` horizontally, so the
    slab forms a *local* dent under the sole instead of dropping uniformly.
    """
    top = np.zeros((PAD_NX, PAD_NZ))
    if not foot_pts:
        return top
    fp = np.vstack(foot_pts)
    xs = np.linspace(PAD_X[0], PAD_X[1], PAD_NX)
    zs = np.linspace(PAD_Z[0], PAD_Z[1], PAD_NZ)
    Xg, Zg = np.meshgrid(xs, zs, indexing="ij")
    gx, gz = Xg.ravel(), Zg.ravel()
    # squared horizontal distance from every grid cell to every foot vertex
    d2 = (gx[:, None] - fp[None, :, 0]) ** 2 + (gz[:, None] - fp[None, :, 2]) ** 2
    near = d2 <= FOOT_R ** 2
    ys = fp[:, 1]
    out = top.ravel()
    for gi in range(gx.size):
        m = near[gi]
        if m.any():
            y = float(ys[m].min())
            if y < 0.0:
                out[gi] = max(GROUND_Y, y)
    return top


def _pad_solid(top: "np.ndarray") -> "pv.PolyData":
    """Closed solid between the dented top surface and the floor at ``GROUND_Y``."""
    xs = np.linspace(PAD_X[0], PAD_X[1], PAD_NX)
    zs = np.linspace(PAD_Z[0], PAD_Z[1], PAD_NZ)
    Xg, Zg = np.meshgrid(xs, zs, indexing="ij")
    # StructuredGrid wants x varying fastest -> transpose (NX,NZ) -> (NZ,NX)
    xT, zT, tT = Xg.T.ravel(), Zg.T.ravel(), top.T.ravel()
    layer_top = np.column_stack([xT, tT, zT])
    layer_bot = np.column_stack([xT, np.full_like(tT, GROUND_Y), zT])
    g = pv.StructuredGrid()
    g.points = np.vstack([layer_top, layer_bot])
    g.dimensions = (PAD_NX, PAD_NZ, 2)
    return g.extract_surface()



def mat4(T):
    R = np.array([[T.R().get(r, c) for c in range(3)] for r in range(3)])
    p = np.array([T.p().get(0), T.p().get(1), T.p().get(2)])
    return R, p


def body_meshes(model):
    """[(mesh_polydata, body_name)] for every attached mesh, cached once."""
    cache, out = {}, []
    bs = model.getBodySet()
    for i in range(bs.getSize()):
        b = bs.get(i)
        ag = b.getPropertyByName("attached_geometry")
        for k in range(ag.size()):
            g = os.Mesh.safeDownCast(ag.getValueAsObject(k))
            if g is None:
                continue
            fn = g.get_mesh_file()
            if fn not in cache:
                p = GEO / fn
                cache[fn] = pv.read(str(p)).triangulate() if p.exists() else None
            m = cache[fn]
            if m is None:
                continue
            sc = g.get_scale_factors()
            s = np.array([sc.get(0), sc.get(1), sc.get(2)])
            mm = m.copy()
            if not np.allclose(s, 1.0):
                mm.points = np.asarray(mm.points) * s
            out.append((mm, b.getName()))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--height", type=float, default=2.0)
    ap.add_argument("--on-pad", action="store_true")
    ap.add_argument("--mass", type=float, default=77.0)
    ap.add_argument("--n-free", type=int, default=24)
    ap.add_argument("--n-impact", type=int, default=60)
    ap.add_argument("--fps", type=int, default=25)
    ap.add_argument("--passive", action="store_true",
                    help="S3 path (b): unlock joints + passive stiffness (bends ~2 deg)")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    if args.quick:
        args.n_free, args.n_impact = 3, 4

    grf = ground_reaction(height_m=args.height, mass_kg=args.mass,
                          on_pad=args.on_pad)
    v0 = float(grf.impact_speed_ms)
    kw = dict(scale_mass_kg=args.mass)
    if args.passive:
        kw["passive"] = PassiveStiffness.nominal_landing()
    fall = run_dead_drop(grf, **kw)
    model = fall.model
    table = fall.states
    tt = np.asarray(table.getIndependentColumn(), dtype=float)
    print(f"impact v0={v0:.2f} m/s  peak GRF={fall.peak_applied_n/1000:.1f} kN  "
          f"FD t=[{tt[0]*1e3:.1f},{tt[-1]*1e3:.1f}] ms  rows={tt.size}  "
          f"mode={'PASSIVE(unlock)' if args.passive else 'LOCKED(S1 rigid leg)'}")

    state = model.initSystem()
    labels = [table.getColumnLabel(j) for j in range(table.getNumColumns())]
    cols = {lb: np.asarray(table.getDependentColumn(lb).to_numpy(), dtype=float)
            for lb in labels}

    def set_row(tq):
        i = int(np.argmin(np.abs(tt - tq)))
        for lb in labels:
            model.setStateVariableValue(state, lb, float(cols[lb][i]))
        model.realizePosition(state)

    meshes = body_meshes(model)
    print("meshes:", len(meshes))
    cs = model.getCoordinateSet()

    def flex_deg(name):
        return float(np.degrees(cs.get(name).getValue(state)))

    g = 9.81
    T_fall = v0 / g
    OUTD.mkdir(parents=True, exist_ok=True)
    for f in OUTD.glob("f*.png"):
        f.unlink()

    frames = []
    for k in range(args.n_free):
        t = T_fall * k / max(1, args.n_free - 1)
        frames.append((t, args.height - 0.5 * g * t * t, "pose0", None))
    for k in range(args.n_impact):
        tq = tt[0] + (tt[-1] - tt[0]) * k / max(1, args.n_impact - 1)
        frames.append((T_fall + tq, 0.0, "fd", tq))

    joint_tag = ("joints: UNLOCKED + passive stiffness (S3b)"
                 if args.passive else "joints: LOCKED (S1 rigid leg, [Y25])")

    # fixed floor the pad rests on (visual anchor: it must NOT move)
    ground = pv.Box(bounds=(-1.5, 1.5, GROUND_Y - 0.06, GROUND_Y, -1.5, 1.5))
    # uncompressed-pad reference wireframe: makes the compression legible
    pad_ref = pv.Box(bounds=(-0.62, 0.62, GROUND_Y, 0.0, -0.45, 0.55))

    pv.OFF_SCREEN = True
    for idx, (tdisp, dy, src, tq) in enumerate(frames):
        if src == "pose0":
            set_row(float(tt[0]))
        else:
            set_row(float(tq))

        # --- pose every mesh once; collect foot verts for the pad dent ------
        posed, foot_pts = [], []
        for mm, bn in meshes:
            R, p = mat4(model.getBodySet().get(bn).getTransformInGround(state))
            pts = np.asarray(mm.points) @ R.T + p + np.array([0.0, dy, 0.0])
            posed.append(pts)
            if _is_foot(bn):
                foot_pts.append(pts)

        # pad top deforms under the foot (natural local dent); base stays fixed
        top = _pad_top_field(foot_pts)
        max_dent = float(-top.min())
        pad = _pad_solid(top)

        pl = pv.Plotter(off_screen=True, window_size=WIN)
        pl.set_background("white")
        pl.add_mesh(ground, color=GROUND, opacity=1.0)
        pl.add_mesh(pad_ref, style="wireframe", color="#7b8bab",
                    line_width=1, opacity=0.5)
        pl.add_mesh(pad, color=PAD, opacity=0.85)
        for (mm, bn), pts in zip(meshes, posed):
            m2 = mm.copy()
            m2.points = pts
            pl.add_mesh(m2, color=BONE, smooth_shading=True,
                        specular=0.25, ambient=0.5, diffuse=0.75)

        if src == "fd":
            fN = float(np.interp(tq, grf.t_s, grf.f_left_n + grf.f_right_n))
            txt = f"impact  t = {tq*1e3:5.1f} ms   F_GRF = {fN/1e3:5.2f} kN"
        else:
            txt = f"free fall   v = {g*tdisp:4.2f} m/s   drop {args.height:.1f} m"
        pl.add_text(f"{txt}\n{joint_tag}\n"
                    f"knee_r = {flex_deg('knee_angle_r'):+5.2f} deg    "
                    f"ankle_r = {flex_deg('ankle_angle_r'):+5.2f} deg\n"
                    f"pad: {max_dent*100:4.1f} / {PAD_THICK*100:.0f} cm compressed",
                    position="upper_left", font_size=12, color="#222222")
        pl.camera.up = (0, 1, 0)
        pl.camera.focal_point = (0.0, 1.70, 0.0)
        pl.camera.position = (6.0, 2.9, 6.0)
        pl.reset_camera_clipping_range()
        img = pl.screenshot(return_img=True)
        pl.close()

        # composite: 3D frame + GRF(t) panel
        fig = plt.figure(figsize=(15.0, 9.0), dpi=100)
        ax0 = fig.add_axes([0.005, 0.005, 0.615, 0.99]); ax0.axis("off")
        ax0.imshow(img)
        ax1 = fig.add_axes([0.665, 0.16, 0.30, 0.66])
        ax1.plot(grf.t_s * 1e3, (grf.f_left_n + grf.f_right_n) / 1e3,
                 color="#c0392b", lw=2.0)
        if src == "fd":
            ax1.axvline(tq * 1e3, color="#111", lw=1.2)
            ax1.plot([tq * 1e3], [fN / 1e3], "o", color="#c0392b", ms=6)
        ax1.set_xlim(-1, grf.t_s[-1] * 1e3 + 1)
        ax1.set_ylim(0, max(40.0, (grf.f_left_n + grf.f_right_n).max() / 1e3 * 1.12))
        ax1.set_xlabel("impact time (ms)", fontsize=11)
        ax1.set_ylabel("F$_{GRF}$ (kN)", fontsize=11)
        ax1.grid(alpha=0.3)
        ax1.set_title("openSIM forward dynamics · dead-drop", fontsize=11)
        ax1.text(0.03, 0.97, f"h = {args.height:.1f} m   v0 = {v0:.2f} m/s\n"
                             f"peak = {fall.peak_applied_n/1e3:.1f} kN   "
                             f"mass = {args.mass:.0f} kg"
                             + ("   pad" if args.on_pad else "   rigid")
                             + ("\nunlocked + passive (S3b)" if args.passive
                                else "\nrigid leg (S1)"),
                 transform=ax1.transAxes, fontsize=10, va="top")
        fig.savefig(str(OUTD / f"f{idx:04d}.png"))
        plt.close(fig)
        print(f"  frame {idx+1}/{len(frames)}  {txt}  pad_dent={max_dent*100:.1f}cm")

    if args.quick:
        print("quick mode: frames written to", OUTD)
        return
    name = "fall_simulation_passive" if args.passive else "fall_simulation"
    mp4 = ROOT / "results" / "opensim_fe" / "figures" / f"{name}.mp4"
    subprocess.run([FFMPEG, "-y", "-framerate", str(args.fps),
                    "-i", str(OUTD / "f%04d.png"), "-pix_fmt", "yuv420p",
                    "-vf", "scale=1500:900", str(mp4)], check=True)
    gif = mp4.with_suffix(".gif")
    subprocess.run([FFMPEG, "-y", "-i", str(mp4),
                    "-vf", "fps=15,scale=760:-1:flags=lanczos", str(gif)],
                   check=True)
    print("wrote", mp4, "and", gif)


if __name__ == "__main__":
    main()
