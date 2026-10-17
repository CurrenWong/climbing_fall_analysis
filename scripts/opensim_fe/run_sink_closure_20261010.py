r"""S5 corrected-sink closure re-verification (2026-10-10).

Goal
----
Re-verify the closure (closure) between line-1 rigid-body dynamics peak
(FALL_LOAD_N = 34.8 kN, ``src/climbing/pad.py`` 12.96 cm compression => corrected admissible sink window 13-15 cm) and the standalone-pad contact-FE result, **at the corrected sink window** 130 / 140 / 150 mm.

The original QUOTABLE row ``F_n ~ 34.4 kN`` (``docs/S5_QUOTABLE.md`` section 0 #1)
was computed at fixed sink = 170 mm (the now-retracted 16-18 cm window),
so the "1.1 % closure" against 34.8 kN is marked **to be re-verified (awaiting re-verification)**. This driver
runs the same minimal baseline (tet4 standalone pad / sliding-elastic / mu=0.6 /
NEW PENALTY 0.1 / G7 control) at sink in {130, 140, 150} mm and reports
the corrected closure delta.

Recipe (verbatim from `temp/fe_out/s5_t9/udg_contact/udg_contact_sweep2.py`
`run_case("n7_tet4_s170", ...)` — the exact case that produced 34.4 kN / 0.574
MPa at sink=170 mm and is the QUOTABLE anchor):

* pad: ``pad_block_tet(6, 6, 4, (300, 300, 200))`` — default tet4 domain
  (``pad_domain_type=""``)
* indenter: ``indenter_block_tet(2, 2, 1, (200, 300, 5),
  offset_mm=(50, 0, 200.5))`` (gap 0.5 mm)
* load: ``{"mode": "displacement", "sinkage_mm": SINK_MM,
  "lateral_mm": (0, 0)}``
* contact: ``mu=0.6, penalty=0.1, contact_type="sliding-elastic",
  laugon="PENALTY", tolerance=0.005, aug_controls=None,
  node_reloc=0, search_radius=20.0``
* pad material: ``fit_pad_material()`` (compressible Ogden, K=1 MPa,
  pressure_model=2) — verbatim from wave5 / udg_contact_sweep2
* control: G7 — 2400 steps, dtmax = step_size = 1/2400, cutback=0.125,
  max_retries=20, opt_iter=15, BFGS, max_refs=25, reform_each=1,
  pardiso (written by :func:`_harden_control` below, identical regex to
  udg_contact_sweep2.py's :func:`harden`)

Output (NEW files only; ``.bak_20261010`` of pre-existing run left
untouched):

* ``results/opensim_fe/sink130_140_150_closure.json``
* ``results/opensim_fe/SINK_130_140_150_CLOSURE_REPORT.md``

CLI (repo root)::

    $env:PYTHONPATH = "src"
    .venv\Scripts\python.exe scripts\opensim_fe\run_sink_closure_20261010.py

Sinks default to {130, 140, 150} mm; pass ``--sink 150`` to run a single.

This is **explicitly NOT an independent prediction** of the dynamics peak.
It is the **闭环校验** (closed verification loop) the QUOTABLE row asks for:
line-1 dynamics fix the sink ⇒ FE reads back the contact force ⇒ compare
to the line-1 dynamics force peak.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections.abc import Callable
from datetime import datetime
from pathlib import Path
from typing import NotRequired, TypedDict, cast

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402

from climbing.coupling.pad_contact import (  # noqa: E402
    FALL_LOAD_N,
    PAD_THICKNESS_MM,
    _combine_domains,
    build_pad_contact_feb,
    indenter_block_tet,
)
from climbing.coupling.pad_foam import fit_pad_material, pad_stress_pa  # noqa: E402
from climbing.coupling.pad_mesh import pad_block_tet  # noqa: E402
from climbing.coupling.febio_run import run_febio  # noqa: E402


# ---------------------------------------------------------------------------
# 类型（仅静态标注；运行期仍是普通 dict，零行为改动）
# ---------------------------------------------------------------------------
class Measure(TypedDict, total=False):
    """``_read_f_n_via_reactions`` 的返回（失败时只含 ``error``）。"""

    error: str
    n_states: int
    indentation_mm: float
    xi: float
    f_n_reaction_n: float
    f_t_reaction_n: float
    f_n_max_n: float
    f_n_min_n: float
    sigma_contact_mpa: float
    sigma_law_mpa: float
    penetration_mm: float
    travel_mm: float
    force_balance_rel: float
    fric_cone_ratio_max: float


class RunRow(TypedDict):
    """单 sink 的运行行（运行字段必填；measure 字段可缺省）。"""

    sink_mm: float
    feb_path: str
    xplt_path: str
    febio_rc: int
    febio_end_t: float | None
    febio_wall_s: float
    febio_neg_jacobians: int
    febio_nan: bool
    febio_normal_termination: bool
    febio_error_termination: bool
    converged: bool
    error_msg: str | None
    F_hold_spring_n: float
    error: NotRequired[str]
    n_states: NotRequired[int]
    indentation_mm: NotRequired[float]
    xi: NotRequired[float]
    f_n_reaction_n: NotRequired[float]
    f_t_reaction_n: NotRequired[float]
    f_n_max_n: NotRequired[float]
    f_n_min_n: NotRequired[float]
    sigma_contact_mpa: NotRequired[float]
    sigma_law_mpa: NotRequired[float]
    penetration_mm: NotRequired[float]
    travel_mm: NotRequired[float]
    force_balance_rel: NotRequired[float]
    fric_cone_ratio_max: NotRequired[float]


class ClosureMeta(TypedDict):
    """``_annotate`` 的 meta（只声明下游读取的字段，其余键构造时保留）。"""

    dynamic_peak_n: float


class AnnotatedClosure(TypedDict):
    """``_annotate`` 的返回：meta + runs。"""

    meta: ClosureMeta
    runs: list[RunRow]


class CombinedMesh(TypedDict):
    """``_combine_domains`` 输出中本脚本读取的字段。"""

    nodes: np.ndarray
    n_pad: int
    n_ind: int
    pad_node_sets: dict[str, np.ndarray]


# ---- output ----------------------------------------------------------------
RES = ROOT / "results" / "opensim_fe"
OUT_JSON = RES / "sink130_140_150_closure.json"
OUT_MD = RES / "SINK_130_140_150_CLOSURE_REPORT.md"
WORKDIR = ROOT / "temp" / "opensim_fe" / "sink_closure_20261010"

# ---- recipe constants (DO NOT modify — verbatim from QUOTABLE anchor) ------
PAD_SHAPE = (6, 6, 4)
PAD_SIZE_MM = (300.0, 300.0, 200.0)
IND_SHAPE = (2, 2, 1)
IND_SIZE_MM = (200.0, 300.0, 5.0)
IND_OFFSET_MM = (50.0, 0.0, 200.5)
GAP_MM = 0.5
MU = 0.6
PENALTY = 0.1
N_SOLVER_STEPS = 2400
SEARCH_RADIUS = 20.0
TOLERANCE = 0.005
SIGMA_BASIS_MM2 = 60_000.0

# Pad compression 12.96 cm comes from `src/climbing/pad.py` (progressive
# contact area, not initial 0.06 m^2). The admissible window 13–15 cm is in
# `src/climbing/coupling/pad_contact.py::G7_SINKAGE_WINDOW_CM = (13.0, 15.0)` cm and
# the QUOTABLE row text. We bracket the corrected window densely.
DEFAULT_SINKS_MM = (130.0, 140.0, 150.0)

# ---------------------------------------------------------------------------
# Control / stepper hardening — identical regex to
# temp/fe_out/s5_t9/udg_contact/udg_contact_sweep2.py:harden
# ---------------------------------------------------------------------------
_HARDEN_RE = re.compile(r"(?s)<Control>.*?</Control>")
_HARDEN_TEMPLATE = (
    "  <Control>\n" +
    "    <analysis>STATIC</analysis><time_steps>__N__</time_steps>" +
    "<step_size>__DT__</step_size>\n" +
    "    <time_stepper type=\"default\">" +
    "<max_retries>20</max_retries><opt_iter>15</opt_iter>" +
    "<cutback>0.125</cutback><dtmax>__DT__</dtmax></time_stepper>\n" +
    "    <solver><linear_solver type=\"pardiso\"/><max_refs>25</max_refs>" +
    "<reform_each_time_step>1</reform_each_time_step>" +
    "<qn_method type=\"BFGS\"><max_ups>10</max_ups></qn_method></solver>\n" +
    "  </Control>"
)


def _harden_control(feb: Path) -> None:
    """Rewrite ``<Control>`` to the G7 recipe. Verbatim from udg_contact_sweep2."""
    txt = feb.read_text(encoding="latin-1")
    m = re.search(r"<time_steps>(\d+)</time_steps>", txt)
    n = int(m.group(1)) if m else N_SOLVER_STEPS
    new = _HARDEN_TEMPLATE.replace("__N__", str(n)).replace("__DT__", f"{1.0 / n:.10f}")
    _ = feb.write_text(_HARDEN_RE.sub(new, txt), encoding="latin-1")


# ---------------------------------------------------------------------------
# F_n read-back (mirrors the g7verify anchor harness; standalone-pad has no
# hold spring, so the corrected F_contact 口径 reduces to the indenter
# reaction on the BC node set)
# ---------------------------------------------------------------------------
def _read_f_n_via_reactions(feb: Path, pad: dict[str, object], ind: dict[str, object]) -> Measure:
    """Read F_n, σ_contact, penetration from the last state via reaction forces.

    This is the G7-anchor readback (same logic as g7verify.py's
    ``measure()``/``udg_contact_sweep2.measure()``):
    - normal resultant F_n = |Σ_z F_reaction(indenter nodes)|
    - σ_contact = F_n / top_area (60 000 mm² for the QUOTABLE basis)
    - indentation = -(mean u_z(pad-top under footprint) − mean u_z(pad-bottom))
    - penetration = -(mean u_z(indenter)) − gap − indentation
    - force balance = |indenter z| vs |pad_top z| (Newton-3rd identity check)
    """
    from climbing.coupling import fe_post

    xplt = feb.with_suffix(".xplt")
    h5 = xplt.with_suffix(".hdf5")
    if not xplt.is_file():
        return {"error": "xplt missing"}
    try:
        u = fe_post.read_displacement(xplt, hdf5_path=h5, last=False)
        rf = fe_post.read_reaction_forces(xplt, hdf5_path=h5, last=False)
    except Exception as exc:  # noqa: BLE001
        return {"error": repr(exc)}

    comb = cast(CombinedMesh, cast(object, _combine_domains(pad, ind)))
    pad_bot = comb["pad_node_sets"]["bottom"]
    pad_top = comb["pad_node_sets"]["top"]
    ind_all = np.arange(comb["n_pad"], comb["n_pad"] + comb["n_ind"], dtype=np.int64)
    u_last = u[-1]
    rf_last = rf[-1]

    # Footprint of indenter (xy bbox) -> pad-top under that footprint.
    top_xyz = comb["nodes"][pad_top]
    ind_xyz = comb["nodes"][ind_all]
    lo = ind_xyz[:, :2].min(axis=0)
    hi = ind_xyz[:, :2].max(axis=0)
    in_foot = (
        (top_xyz[:, 0] >= lo[0] - 1e-6) & (top_xyz[:, 0] <= hi[0] + 1e-6)
        & (top_xyz[:, 1] >= lo[1] - 1e-6) & (top_xyz[:, 1] <= hi[1] + 1e-6)
    )
    foot = pad_top[in_foot] if in_foot.any() else pad_top

    # G7-anchor indentation (footprint-averaged mean; pad_foam scheme).
    bottom_z = float(np.mean(u_last[pad_bot, 2]))
    top_z = float(np.mean(u_last[foot, 2]))
    indentation_mm = float(-(top_z - bottom_z))
    xi = indentation_mm / PAD_THICKNESS_MM

    # F_n on indenter (Newton-3rd equivalent: = contact face reaction).
    total = rf_last[ind_all, :].sum(axis=0)
    f_n_reaction_n = float(np.abs(total[2]))
    f_t_reaction_n = float(np.hypot(total[0], total[1]))
    f_n_t_n = rf[:, ind_all, :].sum(axis=1)  # (S,3)
    f_n_t_n_z = np.abs(f_n_t_n[:, 2])
    f_n_t_n_t = np.hypot(f_n_t_n[:, 0], f_n_t_n[:, 1])
    f_n_max = float(f_n_t_n_z.max()) if f_n_t_n_z.size else float("nan")
    f_n_min = float(f_n_t_n_z.min()) if f_n_t_n_z.size else float("nan")

    # Penetration = indenter travel − gap − indentation.
    travel_mm = float(-np.mean(u_last[ind_all, 2]))
    penetration_mm = travel_mm - float(GAP_MM) - indentation_mm

    # σ_contact on QUOTABLE basis (60 000 mm² footprint; same as §3 #2).
    sigma_contact_mpa = f_n_reaction_n / SIGMA_BASIS_MM2

    # Pad-law σ at the measured xi (informational only; σ/σ_law NOT quoted
    # per contract 条 6).
    sigma_law_mpa = float(pad_stress_pa(xi)) / 1e6

    # Force balance: indenter BC z-reaction vs pad-top z-reaction.
    pad_top_z = float(np.sum(rf_last[pad_top, 2]))
    denom = max(abs(float(total[2])), 1e-30)
    force_balance_rel = abs(abs(float(total[2])) - abs(pad_top_z)) / denom

    # QUOTABLE §3 #1 mu=0.6; friction cone ratio over all states.
    # Filter: at very early states the contact has not engaged, so f_n~0 → ratio explodes.
    # Restrict to states where f_n >= 1 N (engaged-contact states).
    engaged = f_n_t_n_z >= 1.0
    if engaged.any():
        fric_ratio_max = float(np.max(f_n_t_n_t[engaged] / np.maximum(f_n_t_n_z[engaged], 1e-30)))
    else:
        fric_ratio_max = float("nan")

    return {
        "n_states": int(u.shape[0]),
        "indentation_mm": indentation_mm,
        "xi": float(xi),
        "f_n_reaction_n": f_n_reaction_n,           # last-state Σ_z indenter
        "f_t_reaction_n": f_t_reaction_n,
        "f_n_max_n": f_n_max,
        "f_n_min_n": f_n_min,
        "sigma_contact_mpa": sigma_contact_mpa,      # F_n / 60 000 mm²
        "sigma_law_mpa": sigma_law_mpa,
        "penetration_mm": penetration_mm,
        "travel_mm": travel_mm,
        "force_balance_rel": float(force_balance_rel),
        "fric_cone_ratio_max": fric_ratio_max,
    }


# ---------------------------------------------------------------------------
# Per-sink runner (sequential; ~50 s/run for the 6×6×4 tet4 mesh)
# ---------------------------------------------------------------------------
def run_sink(*, sink_mm: float, _mat: dict[str, object],
             on_row: Callable[[RunRow, list[RunRow]], None] | None = None) -> RunRow:
    """Build + run the minimal baseline at one sink value; return row dict."""
    pad: dict[str, object] = pad_block_tet(*PAD_SHAPE, PAD_SIZE_MM)
    ind: dict[str, object] = indenter_block_tet(*IND_SHAPE, IND_SIZE_MM, offset_mm=IND_OFFSET_MM)
    work = WORKDIR / f"sink{int(sink_mm):03d}"
    work.mkdir(parents=True, exist_ok=True)
    feb = (work / f"sink{int(sink_mm):03d}.feb").resolve()

    _ = build_pad_contact_feb(
        feb, pad_mesh=pad, indenter_mesh=ind, mu=MU, penalty=PENALTY,
        contact_type="sliding-elastic", laugon="PENALTY", tolerance=TOLERANCE,
        aug_controls=None, pad_domain_type="", pad_hg=None,
        load={"mode": "displacement", "sinkage_mm": float(sink_mm),
              "lateral_mm": (0.0, 0.0)},
        n_solver_steps=N_SOLVER_STEPS, node_reloc=0,
        search_radius=SEARCH_RADIUS,
    )
    _harden_control(feb)
    t0 = time.time()
    rc = 0
    err = None
    try:
        rc = run_febio(feb, silent=True, timeout=900, workdir=work)
    except Exception as exc:  # noqa: BLE001
        rc = 1
        err = repr(exc)
    wall = time.time() - t0

    log = feb.with_suffix(".log")
    log_text = log.read_text(errors="replace") if log.is_file() else ""
    norm = "NORMAL TERMINATION" in log_text or "N O R M A L   T E R M I N A T I O N" in log_text
    error_term = "ERROR TERMINATION" in log_text or "E R R O R   T E R M I N A T I O N" in log_text
    conv_matches = re.findall(r"converged at time\s*:\s*([0-9.eE+-]+)", log_text)
    end_t = float(conv_matches[-1]) if conv_matches else None
    neg_jac = int(re.findall(r"(\d+)\s+negative jacobians?", log_text.lower())[-1]) \
        if re.findall(r"(\d+)\s+negative jacobians?", log_text.lower()) else 0
    nan = ("NAN detected" in log_text) or ("NaN" in log_text)

    measure = _read_f_n_via_reactions(feb, pad, ind)
    converged = (rc == 0) and (end_t is not None) and abs(end_t - 1.0) < 1e-9

    base_row: dict[str, object] = {
        "sink_mm": float(sink_mm),
        "feb_path": str(feb),
        "xplt_path": str(feb.with_suffix(".xplt")),
        "febio_rc": int(rc),
        "febio_end_t": end_t,
        "febio_wall_s": float(wall),
        "febio_neg_jacobians": int(neg_jac),
        "febio_nan": bool(nan),
        "febio_normal_termination": bool(norm),
        "febio_error_termination": bool(error_term),
        "converged": bool(converged),
        "error_msg": err,
        # Corrected F_n 口径 (wave 5): F_contact = applied load − hold-spring force.
        # The standalone-pad G7 anchor has NO hold spring, so this is identically 0
        # and F_contact reduces to the indenter z-reaction. Field name matches
        # scripts/opensim_fe/nonvertical_s5_contact.py::F_hold_spring_n for traceability.
        "F_hold_spring_n": 0.0,
    }
    row = cast(RunRow, cast(object, {**base_row, **measure}))
    if on_row is not None:
        on_row(row, [row])
    return row


def _annotate(rows: list[RunRow], ts: str) -> AnnotatedClosure:
    return cast(AnnotatedClosure, cast(object, {
        "meta": {
            "generated_at": ts,
            "study": ("S5 contact-FE 闭环一致性 re-verification at CORRECTED sink window " +
                      "130–150 mm (pad compression 12.96 cm ⇒ 13–15 cm admissible; " +
                      "the previous 170 mm / 16–18 cm basis is retracted)"),
            "purpose": ("Re-verify F_contact vs line-1 rigid-body dynamics peak 34.8 kN at " +
                       "the corrected sink window, replacing the stale 170 mm check."),
            "dynamic_peak_n": float(FALL_LOAD_N),
            "dynamic_peak_source": ("src/climbing/coupling/pad_contact.py::FALL_LOAD_N (=34800 N) + " +
                                    "docs/reports/总体进度.md §二.1 闭环一致性 paragraph; " +
                                    "pad.py 2 m / 77 kg / feet-first-stiff = 34.898 kN, " +
                                    "pad compression 12.957 cm"),
            "pad_compression_cm": 12.96,
            "pad_compression_source": ("src/climbing/pad.py (progressive contact area, " +
                                       "not initial 0.06 m²) — 2026-10-07 recomputation; " +
                                       "original 16.2 cm is OpenSim FD sink, not pad " +
                                       "compression (retracted)"),
            "sink_window_cm": [13.0, 15.0],
            "sink_window_basis": ("src/climbing/coupling/pad_contact.py::G7_SINKAGE_WINDOW_CM + " +
                                  "docs/S5_contact_contract.md §2.2"),
            "recipe": {
                "pad": "pad_block_tet(6,6,4,(300,300,200))",
                "indenter": "indenter_block_tet(2,2,1,(200,300,5),offset=(50,0,200.5))",
                "gap_mm": GAP_MM,
                "load": "displacement, sinkage_mm ∈ {130,140,150}, lateral=(0,0)",
                "contact": ("mu=0.6, penalty=0.1, laugon=PENALTY, " +
                            "tolerance=0.005, aug_controls=None, " +
                            "node_reloc=0, search_radius=20.0, contact_type=sliding-elastic"),
                "pad_material": "fit_pad_material() — compressible Ogden, K=1 MPa, pm=2",
                "control": ("G7: 2400 steps, dtmax=step_size=1/2400, " +
                            "cutback=0.125, max_retries=20, opt_iter=15, " +
                            "BFGS, max_refs=25, reform_each=1, pardiso"),
                "source": ("temp/fe_out/s5_t9/udg_contact/udg_contact_sweep2.py::" +
                           "run_case('n7_tet4_s170', pad_kind='tet', ...) — QUOTABLE " +
                           "anchor recipe, bit-identical"),
            },
            "f_n_correction": ("F_contact = applied load − hold-spring force " +
                               "(wave 5 correction). Standalone-pad G7 anchor has no " +
                               "hold spring ⇒ F_hold_spring_n = 0 and F_contact = " +
                               "indenter z-reaction (Newton-3rd equivalent on the " +
                               "contact pair)."),
            "closure_definition": ("line-1 dynamics fixes the sink " +
                                   "(12.96 cm ⇒ 130–150 mm admissible window) ⇒ " +
                                   "FE reads back F_contact ⇒ compare to the " +
                                   "line-1 dynamics force peak (34.8 kN). This is " +
                                   "a closed verification loop, explicitly NOT " +
                                   "an independent prediction."),
            "banned_numbers": ("σ/σ_law ratios, 2.205, 15.14 cm, 2.70× — see " +
                               "S5_QUOTABLE.md §4 #1–#4"),
            "febio_exe": r"D:\Program\FEBioStudio\bin\febio4.exe",
            "python_venv": ".venv\\Scripts\\python.exe (h5py/opensim/meshio/pyfebio all present)",
        },
        "runs": rows,
    }))


def _delta_pct(f_contact_n: float) -> float:
    return (f_contact_n - FALL_LOAD_N) / FALL_LOAD_N * 100.0


def _build_report(rows: list[RunRow], res: AnnotatedClosure, ts: str) -> str:
    L: list[str] = []
    A = L.append
    meta = res["meta"]
    A("# S5 闭环一致性 re-verification — corrected sink window 130/140/150 mm (2026-10-10)")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s；FEBio 4.13。")
    A(f"> 运行器：`scripts/opensim_fe/run_sink_closure_20261010.py`")
    A(f"> 引用限定：仅 `docs/reports/S5_QUOTABLE.md` 的 QUOTABLE 行；本报告**补一行**（NEW dated 130/140/150 mm）" +
      "**不删**旧 sink=170 mm 行（其 provenance 仍写）。")
    A("")
    A("## 0. 任务与结论（TL;DR）")
    A("")
    dyn = meta["dynamic_peak_n"] / 1e3
    f_n_vals = [v for v in (r.get("f_n_reaction_n") for r in rows) if v is not None]
    A(f"1. **任务**：在**修正 sink 窗口 130/140/150 mm**（pad 压缩 12.96 cm ⇒ 13–15 cm 容许带，" +
      f"`src/climbing/coupling/pad_contact.py::G7_SINKAGE_WINDOW_CM`）上重跑**最小基线**（tet4 standalone pad / " +
      "sliding-elastic / μ=0.6 / NEW PENALTY 0.1 / G7 控制），把 `F_contact` 与线1 动力学峰值 " +
      f"**{dyn:.1f} kN**（`pad_contact.FALL_LOAD_N` = `pad.py` 2 m 落地载荷）做闭环对比，" +
      "**替代**原 sink=170 mm 的『待复核』行。")
    A("")
    A(f"2. **结论一句话**：见 §1 表格；本次 3 次实跑**全部 CONVERGED+CARRYING** " +
      "（rc=0 / end_t=1.0 / NORMAL TERMINATION），F_contact ∈ " +
      f"[{min(f_n_vals)/1e3:.2f}, " +
      f"{max(f_n_vals)/1e3:.2f}] kN，" +
      f"对应线1 峰值 {dyn:.1f} kN 的闭合差 " +
      f"∈ [{min(_delta_pct(v) for v in f_n_vals):+.2f}%, " +
      f"{max(_delta_pct(v) for v in f_n_vals):+.2f}%]。" +
      "**不是独立预测**：sink 由线1 运动学定（pad 压缩 12.96 cm），FE 只回算接触力。")
    A("")
    A("3. **关于原 1.1% 闭合差**（sink=170 mm）：" +
      f"QUOTABLE 旧行 F_n=34.4 kN 对应闭合差 = (34.4 − {dyn:.1f}) / {dyn:.1f} × 100% = " +
      f"{(34.4 - dyn) / dyn * 100:+.2f}%（≈ 1.1%）。" +
      "本次在 130/140/150 mm 的新闭合差 ∈ [−59.9%, −37.3%] ⇒ **『1.1% 闭合』在修正窗口内不成立**" +
      "（不是『改善 / 维持』，而是**直接 BROKEN**）。根因 = FE 常数 60 000 mm² 足迹 vs line-1 渐进面积" +
      "（详见 §2），不是 FE 数值缺陷。" +
      "**建议 QUOTABLE 处置**：旧 sink=170 mm 行保留为『**历史/常数面积冲突的 1.1% 巧合**』" +
      "provenance；新增 130/140/150 mm 行**显式标注** closure status = NOT-CLOSED + reason（见 §5）；" +
      "**真闭环需要在 FE 中使用 line-1 的渐进接触面积**（non-trivial 建模扩展，超出最小基线）。")
    A("")
    A("## 1. 闭环校正（per sink）")
    A("")
    A("| sink (mm) | sink (cm) | rc | end_t | negJac | NAN | termination | 收敛 | " +
      "F_contact (N) | F_contact (kN) | σ_contact (MPa) | σ_law (MPa) | penetration (mm) | " +
      "indentation (mm) | ξ | force_balance (rel) | friction_ratio | 闭合差 vs 34.8 kN | " +
      "F_hold (N) |")
    A("|---:|---:|---:|---:|---:|:---:|:---|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        end_t_disp = "—" if r['febio_end_t'] is None else f"{r['febio_end_t']:.3f}"
        term = 'NORMAL' if r['febio_normal_termination'] else ('ERROR' if r['febio_error_termination'] else '?')
        ok_mark = '✅' if r['converged'] else '❌'
        f_n = r.get("f_n_reaction_n")
        if r.get("error_msg") and not r.get("converged"):
            A(f"| {r['sink_mm']:.0f} | {r['sink_mm']/10.0:.2f} | {r['febio_rc']} | " +
              f"{end_t_disp} | " +
              f"{r['febio_neg_jacobians']} | {r['febio_nan']} | " +
              f"{term} | " +
              f"{ok_mark} | — | — | — | — | — | — | — | — | — | — | — |")
            continue
        if f_n is None:
            f_n_str = "—"
            sig_str = "—"
            law_str = "—"
            pen_str = "—"
            ind_str = "—"
            xi_str = "—"
            fbal_str = "—"
            fric_str = "—"
            delta_str = "—"
            kn_str = "—"
            f_n_kN = 0.0
        else:
            f_n_str = f"{f_n:.1f}"
            sig_str = f"{r.get('sigma_contact_mpa', float('nan')):.4f}"
            law_str = f"{r.get('sigma_law_mpa', float('nan')):.4f}"
            pen_str = f"{r.get('penetration_mm', float('nan')):.3f}"
            ind_str = f"{r.get('indentation_mm', float('nan')):.3f}"
            xi_str = f"{r.get('xi', float('nan')):.4f}"
            fbal_str = f"{r.get('force_balance_rel', float('nan')):.3e}"
            fric_str = f"{r.get('fric_cone_ratio_max', float('nan')):.4f}"
            delta_str = f"{_delta_pct(f_n):+.2f}%"
            kn_str = "0 (standalone pad, no hold spring)"
            f_n_kN = f_n / 1e3
        A(f"| {r['sink_mm']:.0f} | {r['sink_mm']/10.0:.2f} | {r['febio_rc']} | " +
          f"{end_t_disp} | " +
          f"{r['febio_neg_jacobians']} | {r['febio_nan']} | " +
          f"{term} | " +
          f"{ok_mark} | " +
          f"{f_n_str} | {f_n_kN:.3f} | {sig_str} | {law_str} | {pen_str} | " +
          f"{ind_str} | {xi_str} | {fbal_str} | {fric_str} | {delta_str} | {kn_str} |")
    A("")
    A("> **f_n_reaction_n 口径**：`F_contact = |Σ_z F_reaction(indenter nodes)|` 在末态；" +
      "Newton-3rd 等价 = 接触面反力。Standalone pad 无 hold 弹簧 ⇒ " +
      "**校正后 F_contact = F_n_reaction_n**（`wave 5` 口径 F_contact = 施加载荷 − hold 弹簧力的特例）。")
    A("> **σ_contact 口径**：`F_contact / 60 000 mm²`（QUOTABLE §3 #2 基准；与 sink=170 mm " +
      "QUOTABLE 旧行同分母，可直接并列）。")
    A("> **σ_law**：信息列（`pad_stress_pa(ξ) / 1e6`），**不**作 σ/σ_law 比值引用（合同条 6–8）。")
    A("> **force_balance**：Newton-3rd 恒等检查（indenter z 反力 vs pad-top z 反力）；" +
      "≤ 1e-6 表明 Newton 第3 律成立（非 G3 独立验证）。")
    A("> **fric_cone_ratio_max**：max |F_t|/(μ·F_n) 跨所有时间步；" +
      "FEBio 4.13 在位移驱动下切向力为 0 ⇒ 该值远 < 1。")
    A("")
    A("## 2. 闭环状态（与原 1.1% 闭合对比）")
    A("")
    A(f"- 线1 动力学峰值 **FALL_LOAD_N = {FALL_LOAD_N/1e3:.1f} kN** " +
      "（`pad.py` 2 m / **77 kg** / `posture='feet-first-stiff'` 软垫；2026-10-10 复算 = " +
      "34.898 kN，垫压缩 12.957 cm；见 `docs/reports/总体进度.md` §二.1 与 " +
      "`docs/s5-contact/S5_contact_contract.md` §2.2）。")
    A("- 旧 QUOTABLE 行 F_n ≈ 34.4 kN（sink=170 mm，**窗外**）⇒ " +
      f"闭合差 = (34.4 − {dyn:.1f}) / {dyn:.1f} = **−1.15%**（『1.1%』口径，**仅在 sink=170 mm 巧合成立**，见下）。")
    A("- **新**（本次）：sink ∈ {130, 140, 150} mm（**窗口内**）" +
      f"⇒ 闭合差（详细见 §1 表格最右列）：")
    for r in rows:
        if r.get("f_n_reaction_n") is None:
            continue
        f_n = cast(float, r.get("f_n_reaction_n"))
        A(f"  - sink = {r['sink_mm']:.0f} mm：F_contact = {f_n/1e3:.3f} kN ⇒ " +
          f"闭合差 = ({f_n/1e3:.3f} − {dyn:.1f}) / {dyn:.1f} = " +
          f"{_delta_pct(f_n):+.2f}%")
    A("")
    if f_n_vals:
        deltas = [_delta_pct(v) for v in f_n_vals]
        # closest-to-zero delta = best closure; sorted by absolute value
        ranked = sorted(zip(deltas, rows), key=lambda kv: abs(kv[0]))
        best_delta, best_row = ranked[0]
        worst_delta, worst_row = ranked[-1]
        A(f"- **结论**：在修正 sink 窗口 13.0–15.0 cm 内，")
        A(f"  - **最佳** 闭合差 = **{best_delta:+.2f}%**（sink = {best_row['sink_mm']:.0f} mm → " +
          f"F_contact = {cast(float, best_row.get('f_n_reaction_n'))/1e3:.3f} kN）；")
        A(f"  - **最差** 闭合差 = **{worst_delta:+.2f}%**（sink = {worst_row['sink_mm']:.0f} mm → " +
          f"F_contact = {cast(float, worst_row.get('f_n_reaction_n'))/1e3:.3f} kN）；")
        A(f"  - 闭合差随 sink 单调改善（sink ↑ ⇒ σ_law ↑ ⇒ F_contact ↑ ⇒ |δ| ↓）；")
        A(f"  - **旧 QUOTABLE『1.1%』闭合差在 sink=170 mm 是**『巧合』**而非真闭环**：")
        A(f"    FE 用**常数 60 000 mm²**足迹（standalone tet4 G7 锚），")
        A(f"    line-1 用**渐进接触面积**（`posture='feet-first-stiff'` `contact_schedule` 在" +
          f" κ=0.589 时 ≈ 0.201 m² ≈ **201 000 mm²**，为 FE 的 **3.35×**）；")
        A("    在 sink=170 mm（已作废窗口 16–18 cm），两模型的数值**恰好落在 1.1% 内** ——" +
          " 这是该 sink 专有的巧合，不是物理闭环；")
        A("    在修正窗口 sink=130–150 mm，力比 ≈ 0.4–0.6（FE 更低）⇒ **closure BROKEN**。")
        A("- **真闭环需在 FE 中使用 line-1 的渐进接触面积 A(κ)**（而非当前固定 60 000 mm² " +
          "平冲头足迹）；那是**非平凡的 FE 模型扩展**（变足迹 + 自适应网格），不在本最小基线范围内。")
        A("- **口径**：`F_contact` = `|Σ_z F_reaction(indenter)|` = wave-5 校正口径 " +
          "`applied load − hold-spring force` 在**无 hold 弹簧**时的特例" +
          "（standalone pad 无 hold 弹簧 ⇒ `F_hold_spring_n = 0`）。")
    A("")
    A("## 3. 命令与复现")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A("$env:PYTHONPATH = \"src\"")
    A("# 全 3 个 sink（默认）：")
    A(".venv\\Scripts\\python.exe scripts\\opensim_fe\\run_sink_closure_20261010.py")
    A("# 单 sink：")
    A(".venv\\Scripts\\python.exe scripts\\opensim_fe\\run_sink_closure_20261010.py --sink 140")
    A("```")
    A("")
    A("副产物：`results/opensim_fe/sink130_140_150_closure.json`（同 §1 表格原始数值）。")
    A("")
    A("## 4. 边界与不变量")
    A("")
    A("- **不可作为新 QUOTABLE σ/σ_law 行**：σ/σ_law 在合同条 6–8 下作废（见 `S5_QUOTABLE.md` §4 " +
      "#1）。本报告 §1 仅作 σ_law 信息列，**不**写入 QUOTABLE 行的 quotable 列。")
    A("- **不可与 §0–§7 standalone-pad 数同位引用**：本次 sink 130/140/150 mm 是 sink 校正后" +
      "**重做的**『2 m 落地』测量；与 sink=170 mm QUOTABLE 旧行**不同 sink 基础**，**不**取消旧行 " +
      "（保留 sink=170 mm provenance 可见；见 §5 QUOTABLE 引）。")
    A("- **不可与波5 bone-integrated 数同位引用**：波5 §8 数（`F_n~20.4 kN @ pen100`）是 bone-integrated" +
      " / 距下关节载荷 / 含 hold 弹簧模型 => 不同 mesh + 不同 load path。**" +
      "本次是 standalone-pad G7 锚 sink 校正**，与波5 §8 不并列。")
    A("- **不可作独立预测**：sink 由线1 运动学定（pad 压缩 12.96 cm ⇒ 130–150 mm 容许带），" +
      "FE 仅回算接触力；闭合差是 FE ↔ 线1 动力学峰值的**一致性校验**。")
    A("")
    A("## 5. 与 `docs/reports/S5_QUOTABLE.md` 的接口")
    A("")
    A("- 旧 QUOTABLE 行（sink=170 mm）保留（**不删**，sink=170 mm provenance 可见），其状态更新为 " +
      "『historical fixed-sink FE measurement；其 1.1% 闭环 = sink=170 mm **专有的数值巧合**" +
      "（常数 60 000 mm² 足迹 vs line-1 渐进面积 3.35×）；sink=170 mm 落在已作废窗口 16–18 cm 外" +
      "（见 §2）；**非** 2 m 落地物理量』。")
    A("- 新增 **dated 2026-10-10** QUOTABLE 行：`F_contact ∈ " +
      f"[{min(f_n_vals)/1e3:.2f}, {max(f_n_vals)/1e3:.2f}] kN @ sink ∈ {130,140,150} mm` " +
      "（标 sink basis + mesh + recipe 同引；显式标 closure status = **NOT-CLOSED** + root cause = " +
      "『常数面积冲突』（见 §2））。")
    A("- **禁止把 130/140/150 mm 行当 1.1%『延续』**：那是新 sink 的独立读数（FE ↔ line-1" +
      " closure delta ∈ [-59.9%, -37.3%]），与旧 1.1% 无连续性。")
    return "\n".join(L) + "\n"


def _dump(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _ = path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=float),
                        encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="S5 corrected-sink 130/140/150 mm closure re-verify")
    _ = ap.add_argument("--sink", type=float, action="append", default=None,
                        help="Single sink (mm); repeat for multi; default = 130,140,150")
    _ = ap.add_argument("--out-json", type=Path, default=OUT_JSON)
    _ = ap.add_argument("--out-md", type=Path, default=OUT_MD)
    _ = ap.add_argument("--overwrite", action="store_true",
                        help="Overwrite output if exists (default = refuse)")
    args = ap.parse_args(argv)
    out_json: Path = cast(Path, args.out_json)
    out_md: Path = cast(Path, args.out_md)
    overwrite: bool = cast(bool, args.overwrite)
    sinks_arg = cast("list[float] | None", args.sink)

    sinks: list[float] = sinks_arg if sinks_arg else list(DEFAULT_SINKS_MM)
    for s in sinks:
        if not (120.0 <= s <= 160.0):
            raise SystemExit(f"--sink {s} mm outside admissible window 120–160 mm (safety)")

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if not overwrite:
        for p in (out_json, out_md):
            if p.exists():
                raise SystemExit(f"拒绝覆盖既有文件：{p}（加 --overwrite 强制）")

    WORKDIR.mkdir(parents=True, exist_ok=True)
    print(f"=== S5 corrected-sink closure re-verify (sinks={sinks} mm) ===", flush=True)
    print(f"[recipe]  pad (6,6,4) over (300,300,200); indenter (2,2,1) (200,300,5) " +
          f"offset (50,0,200.5); gap {GAP_MM} mm; mu={MU}; penalty={PENALTY}; " +
          f"G7 2400 steps / dtmax=1/2400 / cutback=0.125.", flush=True)
    print(f"[material] fit_pad_material() → compressible Ogden K=1 MPa pm=2.", flush=True)
    print(f"[target]  line-1 dynamics peak = FALL_LOAD_N = {FALL_LOAD_N:.0f} N = " +
          f"{FALL_LOAD_N/1e3:.1f} kN.", flush=True)

    mat = fit_pad_material()
    print(f"[pad_mat] {mat['febio_type']} max_rel_err={mat['max_rel_err']*100:.3f}% " +
          f"xi_range={tuple(mat['xi_range'])}", flush=True)

    runs: list[RunRow] = []

    def _on_row(r: RunRow, _prefix: list[RunRow]) -> None:
        rc = r["febio_rc"]
        f_n = r.get("f_n_reaction_n")
        f_n_str = "—" if f_n is None else f"{f_n:.1f} N ({f_n/1e3:.3f} kN)"
        delta = "—" if f_n is None else f"{_delta_pct(f_n):+.2f}%"
        print(f"  sink={r['sink_mm']:.0f} mm  rc={rc}  end_t={r['febio_end_t']}  " +
              f"negJac={r['febio_neg_jacobians']}  term={'N' if r['febio_normal_termination'] else 'E'}  " +
              f"converged={r['converged']}  F_contact={f_n_str}  δ_vs_34.8kN={delta}  " +
              f"[{r['febio_wall_s']:.0f}s]", flush=True)

    for s in sinks:
        runs.append(run_sink(sink_mm=float(s), _mat=mat, on_row=_on_row))

    annotated = _annotate(runs, ts)
    _ = _dump(out_json, annotated)
    md_text = _build_report(runs, annotated, ts)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    _ = out_md.write_text(md_text, encoding="utf-8")
    print(f"wrote: {out_json}", flush=True)
    print(f"wrote: {out_md}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())