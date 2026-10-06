"""Scenario demo animations for the "dangerous posture" section of the paper.

A) ``--scenario ser``   → ``fall_scenario_ser.mp4``
   三踝 / Lauge-Hansen SER posture: single-foot + supination β=30° +
   external rotation θ=25°, on the crash pad, h = 2 m.
   NOTE: the OpenSim model CANNOT represent inversion or axial rotation
   (``subtalar_angle`` / ``hip_rotation`` are locked — the project's documented
   limit), so the foot posture itself is illustrated in ``fig8``; this
   animation carries the *load history* + a live ``M_inv`` / risk readout
   computed with the real ``ankle_supination`` module.

B) ``--scenario talus`` → ``fall_scenario_talus.mp4``
   距骨 / Hawkins posture: single-foot + ankle dorsiflexed +25°, on the pad,
   h = 2 m. Dorsiflexion IS a real (unlock-then-lock) sagittal DOF, so the
   rendered pose really is dorsiflexed.
   NOTE: the project has NO quantitative talus criterion → mechanism-only.

Performance
-----------
Both the main view and the ankle close-up REUSE one ``pv.Plotter`` for the whole
sequence; each frame only rewrites the point arrays (bones + pad top) instead of
rebuilding 81 actors. That is ~5x faster than the first version.

Usage
-----
  python scripts/opensim_fe/make_scenario_animations.py --scenario both
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

import make_fall_animation as mfa                                    # noqa: E402
from climbing.coupling.opensim_grf import ground_reaction            # noqa: E402
from climbing.coupling.opensim_fall import run_dead_drop, LegPose    # noqa: E402
from climbing.coupling.ankle_supination import supination_load       # noqa: E402

H_M, MASS, FPS = 2.0, 77.0, 25
N_FREE, N_IMPACT = 15, 45
MAIN_WIN, INSET_WIN = (900, 1080), (460, 460)
ATFL_FAIL_NM = 21.0          # ankle_ligament.ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM
LEVER_MM = 30.0
PAD_BOUNDS = (-0.62, 0.62, -0.45, 0.55)      # x0,x1,z0,z1 of the pad slab


def build(scenario: str, two_foot: bool = False):
    """Return (grf, grf_rigid, fall, meta).

    ``two_foot=True`` reproduces the R1/S5 report baseline exactly: the default
    ``split=(0.5, 0.5)`` (each foot carries half) with the model's own mass
    (75.337 kg) — i.e. the same F as ``nonvertical_s5_supination.py``.
    """
    split = (0.5, 0.5) if two_foot else (1.0, 0.0)
    mass_kw = {} if two_foot else dict(scale_mass_kg=MASS)
    load_tag = ("双脚各半 split=(0.5,0.5)，质量 75.337 kg（= S5 报告口径）"
                if two_foot else
                "★ 单脚承载 split=(1,0)，质量 77 kg（比 S5 报告口径严酷 ≈2×）")
    if scenario == "ser":
        kw = dict(height_m=H_M, mass_kg=MASS, posture="one-leg-awkward",
                  split=split)
        grf = ground_reaction(on_pad=True, **kw)
        grf_rigid = ground_reaction(on_pad=False, **kw)
        fall = run_dead_drop(grf, **mass_kw)
        meta = dict(
            title="三踝骨折 · SER 姿势（旋后内翻 + 外旋）",
            sub=("内翻 β=30° · 外旋 θ=25°　"
                 "（内翻/轴向旋转不属模型自由度，足部姿态见 fig8）"),
            scenario="ser", load_tag=load_tag,
        )
    elif scenario == "talus":
        kw = dict(height_m=H_M, mass_kg=MASS, posture="feet-first-stiff",
                  split=split)
        grf = ground_reaction(on_pad=True, **kw)
        grf_rigid = ground_reaction(on_pad=False, **kw)
        fall = run_dead_drop(grf, posture=LegPose.flexed_landing(knee_deg=0.0,
                                                                 ankle_deg=25.0),
                             **mass_kw)
        meta = dict(
            title="距骨骨折 · 背屈楔入姿势（Hawkins）",
            sub="踝背屈 +25°（锁定姿势，画面为真实姿态）· 机制演示，无量化判据",
            scenario="talus", load_tag=load_tag,
        )
    else:
        raise ValueError(scenario)
    return grf, grf_rigid, fall, meta


def per_frame_value(grf, scenario, tq):
    """(F_total_N, F_right_N, readout) at impact time tq.

    The M_inv readout uses the **right-foot** force — exactly the R1/S5
    convention (``nonvertical_s5_supination.py`` takes ``max(grf.f_right_n)``),
    so ``--two-foot`` reproduces the report numbers verbatim.
    """
    i = int(np.argmin(np.abs(grf.t_s - tq)))
    f_tot = float(grf.f_right_n[i] + grf.f_left_n[i])
    f_r = float(grf.f_right_n[i])
    if scenario == "ser":
        sl = supination_load((0.0, f_r, 0.0), supination_deg=30.0,
                             lever_arm_mm=LEVER_MM)
        return f_tot, f_r, dict(m_inv_nm=sl.m_inv_nm, risk=sl.m_inv_nm / ATFL_FAIL_NM)
    return f_tot, f_r, {}


def make_plotter(meshes, win, *, ref_wire):
    """One reusable plotter with the pad/ground + a private copy of each bone."""
    pl = pv.Plotter(off_screen=True, window_size=win)
    pl.set_background("white")
    g = pv.Box(bounds=(-1.5, 1.5, mfa.GROUND_Y - 0.06, mfa.GROUND_Y, -1.5, 1.5))
    pl.add_mesh(g, color=mfa.GROUND)
    if ref_wire:
        r = pv.Box(bounds=(PAD_BOUNDS[0], PAD_BOUNDS[1], mfa.GROUND_Y, 0.0,
                           PAD_BOUNDS[2], PAD_BOUNDS[3]))
        pl.add_mesh(r, style="wireframe", color="#7b8bab", line_width=1, opacity=0.5)
    pad_pd = pv.Box(bounds=(PAD_BOUNDS[0], PAD_BOUNDS[1], mfa.GROUND_Y, 0.0,
                            PAD_BOUNDS[2], PAD_BOUNDS[3]))
    pl.add_mesh(pad_pd, color=mfa.PAD, opacity=0.85)
    bones = []
    for mm, _bn in meshes:
        pd = mm.copy()
        pl.add_mesh(pd, color=mfa.BONE, smooth_shading=True,
                    specular=0.25, ambient=0.5, diffuse=0.75)
        bones.append(pd)
    return pl, pad_pd, bones


def render(scenario: str, quick: bool = False, two_foot: bool = False) -> None:
    grf, grf_rigid, fall, meta = build(scenario, two_foot=two_foot)
    model, table = fall.model, fall.states
    tt = np.asarray(table.getIndependentColumn(), dtype=float)
    v0 = float(grf.impact_speed_ms)
    peak_r = float(np.max(grf_rigid.f_left_n + grf_rigid.f_right_n))
    print(f"[{scenario}] v0={v0:.2f} m/s  peak={fall.peak_applied_n/1e3:.1f} kN  "
          f"rigid_peak={peak_r/1e3:.1f} kN")

    state = model.initSystem()
    labels = [table.getColumnLabel(j) for j in range(table.getNumColumns())]
    cols = {lb: np.asarray(table.getDependentColumn(lb).to_numpy(), dtype=float)
            for lb in labels}

    def set_row(tq):
        i = int(np.argmin(np.abs(tt - tq)))
        for lb in labels:
            model.setStateVariableValue(state, lb, float(cols[lb][i]))
        model.realizePosition(state)

    meshes = mfa.body_meshes(model)
    cs = model.getCoordinateSet()
    n_free, n_impact = (3, 4) if quick else (N_FREE, N_IMPACT)

    outd = ROOT / "temp" / "opensim_fe" / f"_scenario_{scenario}{'_twofoot' if two_foot else ''}"
    outd.mkdir(parents=True, exist_ok=True)
    for f in outd.glob("f*.png"):
        f.unlink()

    g = 9.81
    T_fall = v0 / g
    frames = []
    for k in range(n_free):
        t = T_fall * k / max(1, n_free - 1)
        frames.append((H_M - 0.5 * g * t * t, float(g * t), None))
    for k in range(n_impact):
        tq = tt[0] + (tt[-1] - tt[0]) * k / max(1, n_impact - 1)
        frames.append((0.0, None, float(tq)))

    pv.OFF_SCREEN = True
    pl, pad_pd, bones = make_plotter(meshes, MAIN_WIN, ref_wire=True)
    pl_in, pad_in, bones_in = make_plotter(meshes, INSET_WIN, ref_wire=False)
    for p_ in (pl, pl_in):
        p_.camera.up = (0, 1, 0)

    tag = "_twofoot" if two_foot else ""
    mp4 = ROOT / "results" / "opensim_fe" / "figures" / f"fall_scenario_{scenario}{tag}.mp4"
    for idx, (dy, v_free, tq) in enumerate(frames):
        is_impact = tq is not None
        set_row(float(tt[0]) if not is_impact else float(tq))

        posed, lo_y = [], np.inf
        for mm, bn in meshes:
            R, p = mfa.mat4(model.getBodySet().get(bn).getTransformInGround(state))
            pts = np.asarray(mm.points) @ R.T + p + np.array([0.0, dy, 0.0])
            lo_y = min(lo_y, float(pts[:, 1].min()))
            posed.append(pts)
        pad_top = float(min(0.0, max(mfa.GROUND_Y, lo_y)))
        new_pad = pv.Box(bounds=(PAD_BOUNDS[0], PAD_BOUNDS[1],
                                 mfa.GROUND_Y, pad_top, PAD_BOUNDS[2], PAD_BOUNDS[3]))

        # ---- main view (camera fixed) ----
        for pd, pts in zip(bones, posed):
            pd.points = pts
        pad_pd.points = new_pad.points.copy()
        pl.camera.focal_point = (0.0, 1.70, 0.0)
        pl.camera.position = (6.0, 2.9, 6.0)
        pl.reset_camera_clipping_range()
        pl.render()
        img = pl.screenshot(return_img=True)

        # ---- ankle close-up (camera tracks the right talus) ----
        ank = mfa.mat4(model.getBodySet().get("talus_r")
                       .getTransformInGround(state))[1] + np.array([0.0, dy, 0.0])
        for pd, pts in zip(bones_in, posed):
            pd.points = pts
        pad_in.points = new_pad.points.copy()
        pl_in.camera.focal_point = tuple(ank)
        pl_in.camera.position = tuple(ank + np.array([0.42, 0.30, 0.42]))
        pl_in.reset_camera_clipping_range()
        pl_in.render()
        img2 = pl_in.screenshot(return_img=True)

        # ---- text (matplotlib layer -> full CJK support) ----
        if is_impact:
            fN, fR, ro = per_frame_value(grf, scenario, tq)
            lines = [f"t = {tq*1e3:5.1f} ms   F_GRF(总) = {fN/1e3:5.2f} kN"
                     f"   右足 = {fR/1e3:5.2f} kN",
                     meta["load_tag"]]
            if scenario == "ser":
                lines.append(f"M_inv = {ro['m_inv_nm']:6.1f} N·m"
                             f"   = {ro['risk']:5.1f} × 21 N·m")
                lines.append("★ 已超 21 N·m 阈值 → 韧带 / 腓骨失效判据成立"
                             if ro["risk"] > 1 else "未超 21 N·m 阈值")
            else:
                dflex = float(np.degrees(cs.get("ankle_angle_r").getValue(state)))
                lines.append(f"踝背屈角 = {dflex:+5.1f}° （锁定姿势）")
                lines.append("〔注〕机制演示：项目无距骨量化判据")
        else:
            lines = [f"自由落体   v = {v_free:4.2f} m/s   drop {H_M:.1f} m",
                     "〔非 FD〕纯运动学插值，姿态 = 初始位形"]
        lines.append(f"垫子压缩 {abs(pad_top)*100:4.1f} cm 〔渲染约定：FD 无接触约束〕")
        status = "\n".join(lines)

        fig = plt.figure(figsize=(15.0, 8.4), dpi=100)
        ax0 = fig.add_axes([0.005, 0.01, 0.60, 0.98]); ax0.axis("off")
        ax0.imshow(img)
        fig.text(0.015, 0.982, meta["title"], fontsize=13, weight="bold", va="top",
                 bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#dddddd", alpha=0.93))
        fig.text(0.015, 0.944, status, fontsize=10, va="top", linespacing=1.5,
                 bbox=dict(boxstyle="round,pad=0.4", fc="#fbfbfb", ec="#dddddd", alpha=0.93))
        ax2 = fig.add_axes([0.030, 0.072, 0.215, 0.280])
        ax2.imshow(img2); ax2.axis("off")
        ax2.text(0.5, -0.02, "踝部特写（右脚）", transform=ax2.transAxes,
                 ha="center", va="top", fontsize=9.5, color="#333")
        ax1 = fig.add_axes([0.655, 0.17, 0.31, 0.64])
        ax1.plot(grf_rigid.t_s * 1e3, (grf_rigid.f_left_n + grf_rigid.f_right_n) / 1e3,
                 color="#8a8a8a", lw=1.6, ls="--", label="rigid ground")
        ax1.plot(grf.t_s * 1e3, (grf.f_left_n + grf.f_right_n) / 1e3,
                 color="#c0392b", lw=2.2, label="crash pad")
        if is_impact:
            ax1.axvline(tq * 1e3, color="#111", lw=1.2)
        ax1.set_xlim(-1, grf.t_s[-1] * 1e3 + 1)
        ax1.set_ylim(0, max(peak_r, fall.peak_applied_n) / 1e3 * 1.12)
        ax1.set_xlabel("impact time (ms)", fontsize=11)
        ax1.set_ylabel("F$_{GRF}$ (kN)", fontsize=11)
        ax1.grid(alpha=0.3); ax1.legend(fontsize=9, loc="upper left")
        ax1.set_title(f"openSIM FD · {'two-foot (S5)' if two_foot else 'single-foot'} impact",
                      fontsize=11)
        fig.text(0.655, 0.045, meta["sub"], fontsize=9.5, color="#444", wrap=True)
        fig.savefig(str(outd / f"f{idx:04d}.png"))
        plt.close(fig)
        print(f"  frame {idx+1}/{len(frames)}  {lines[0]}")

    pl.close(); pl_in.close()
    if quick:
        print("quick mode: frames in", outd)
        return
    subprocess.run([mfa.FFMPEG, "-y", "-framerate", str(FPS),
                    "-i", str(outd / "f%04d.png"), "-pix_fmt", "yuv420p",
                    "-vf", "scale=1500:840", str(mp4)], check=True)
    gif = mp4.with_suffix(".gif")
    subprocess.run([mfa.FFMPEG, "-y", "-i", str(mp4),
                    "-vf", "fps=14,scale=760:-1:flags=lanczos", str(gif)], check=True)
    print("wrote", mp4, "and", gif)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", choices=["ser", "talus", "both"], default="both")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--two-foot", action="store_true",
                    help="split=(0.5,0.5) + 75.337 kg = 精确复现 R1/S5 报告口径")
    a = ap.parse_args()
    for s in (["ser", "talus"] if a.scenario == "both" else [a.scenario]):
        render(s, quick=a.quick, two_foot=a.two_foot)


if __name__ == "__main__":
    main()
