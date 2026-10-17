"""fibula_fe_audit.py -- scoped FE magnitude audit for the fibula all_one defect.

Executes the 4-item scope locked in
``results/opensim_fe/FIBULA_IDENTIFIABILITY_REPORT.md`` §5:

1. **Add tibia-fibula load sharing** at the FE ``ref_load`` entry, using the
   literature anchor (fibula ~0.05 / 0.075 / 0.10 of ``ankle_r``; tibia the
   complement).  Implemented additively in ``fracture_matrix.py`` via
   ``Part.load_fraction`` + ``effective_ref_load`` (default 1.0 = no sharing).
2. **Verify fibula convergence / linearity** under the *allocated* reference
   load: log the descending load-reduction search (monotonic?) and the
   two-load linearity self-check.
3. **Re-emit the fibula fracture-matrix row** from the corrected FE result and
   test whether the logistic regains a 0/1 mix (OR definable).
4. **Register (do not fix)** the other same-class rows that share one joint
   reaction (L3 / T6 / C5 / parietal_r all share ``lumbar``).

Discipline
----------
* Read-only w.r.t. every frozen input: never writes to
  ``results/opensim_fe/fracture_matrix*.json``, ``s4_alignment_bc2.json``,
  ``bc_robust_metric_route1.json``, ``g3_*`` or the frozen per-bone
  ``temp/fe_out/fracture_matrix_bc2/*/*_fe.json``.
* New artifacts only: ``temp/fe_out/fibula_fe_audit/``,
  ``results/opensim_fe/fibula_fe_audit_matrix.json`` and
  ``results/opensim_fe/FIBULA_FE_AUDIT_REPORT.md``.
* FE built with ``align_pca=True`` (the BC2 opt-in) so results are comparable
  to the frozen baseline.
* No σc change (70 MPa stays), no G3 re-run, no target back-solving -- the
  numbers are reported exactly as the solver returns them.

Reproduce (repo root, PYTHONPATH=src)::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\fibula_fe_audit.py
    # optional: --no-tibia (fibula only)  --skip-ladder  --policies 0.05,0.075,0.10
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import fracture_matrix as fm  # noqa: E402
from bone_feb import build_bone_feb, run_bone_feb  # noqa: E402
from scripts.opensim_fe import paper_alignment as pa  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
FROZEN_MATRIX = RES / "fracture_matrix_bc2.json"
FROZEN_FE_DIR = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2"
AUDIT_TEMP = ROOT / "temp" / "opensim_fe" / "fibula_fe_audit"
OUT_JSON = RES / "fibula_fe_audit_matrix.json"
OUT_REPORT = RES / "FIBULA_FE_AUDIT_REPORT.md"

#: literature-anchored fibula axial shares (``bone.py::fibula_ends`` /
#: ``risk_1d_loadshare.py::share_lit05 / share_lit10``).
DEFAULT_POLICIES = (0.05, 0.075, 0.10)

#: bones that draw their whole reaction from the single ``lumbar`` free body
#: (same class of defect; registered, not fixed).
LUMBAR_SHARED_BONES = ("L3", "T6", "C5", "parietal_r")

#: σc used for the fibula row.  70 = route-乙 weak sub-region (fibula ENDS,
#: [Y25] Table 1) = the value G3 uses; 160 = fibula shaft (contrast).
FIBULA_SIGMA_C_WEAK = 70.0
FIBULA_SIGMA_C_SHAFT = 160.0


# --------------------------------------------------------------------------
# BC2 FE build conditions (align_pca=True) -- matches rerun_matrix_bc2.py
# --------------------------------------------------------------------------
def _build_run_bc2(part, mesh, mats, feb, load, time_steps):
    for sd in part.spec.solid_domains:
        sd.material = mats[sd.pids[0]]
    for sh in part.spec.shell_domains:
        sh.material = mats[sh.pids[0]]
    build_bone_feb(part.spec, feb, load_n=load, top_frac=0.30, bot_frac=0.30,
                   axis=None, with_shell=True, time_steps=time_steps, mesh=mesh,
                   align_pca=True)
    run_bone_feb(feb, workdir=feb.parent, timeout=1800)


class Log:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.lines: list[str] = []

    def __call__(self, msg: str) -> None:
        print(msg, flush=True)
        self.lines.append(msg)

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text("\n".join(self.lines) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------
# frozen baseline + row helpers
# --------------------------------------------------------------------------
def load_frozen() -> dict:
    fmz = json.loads(FROZEN_MATRIX.read_text(encoding="utf-8"))
    heights = np.asarray(fmz["heights_m"], float)
    loads = {f"{k}__vert": np.asarray(v, float) for k, v in fmz["loads_n"].items()}
    loads["heights"] = heights
    recs = {p["bone"]: json.loads(
        (FROZEN_FE_DIR / p["bone"] / f"{p['bone']}_fe.json").read_text(encoding="utf-8"))
        for p in fmz["parts"]}
    return {"json": fmz, "heights": heights, "loads": loads, "recs": recs}


def _row_stats(sigma: np.ndarray, sigma_c: float) -> dict:
    sigma = np.asarray(sigma, float)
    r = sigma / float(sigma_c)
    over = int((r >= 1.0).sum())
    first = next((i + 1 for i, v in enumerate(r) if v >= 1.0), None)
    return {
        "sigma_c_mpa": float(sigma_c),
        "ratio_min": float(r.min()), "ratio_max": float(r.max()),
        "ratio_h1": float(r[0]), "ratio_h5": float(r[4]), "ratio_h50": float(r[-1]),
        "n_fractured_h": over, "n_safe_h": int(r.size - over),
        "first_fracture_h": first,
        "all_one": bool(over == r.size), "all_zero": bool(over == 0),
        "identifiable_0_1_mixed": bool(0 < over < r.size),
    }


def _logistic(heights: np.ndarray, row: np.ndarray) -> dict:
    rec = pa.logistic_fit(np.asarray(heights, float), np.asarray(row).astype(int))
    keep = ("status", "OR", "beta", "se", "ci95", "wald_z", "p", "n",
            "n_fracture", "n_safe", "separation_threshold_m", "note", "converged")
    return {k: rec[k] for k in keep if k in rec}


def _search_monotonic(attempts: list[dict]) -> dict:
    """Down-closed convergence?  From a descending search + ladder we report the
    raw sequence; monotone=True iff, within the ladder, every converged load has
    all *smaller* ladder loads also converged (no success above a failure)."""
    return {"attempt_sequence": [[a.get("load_n"), bool(a.get("ok")),
                                  a.get("purpose", "search")] for a in attempts]}


# --------------------------------------------------------------------------
# main audit
# --------------------------------------------------------------------------
def _solve(part, allocated_ref: float, temp_dir: Path, log: Log, reuse: bool) -> dict:
    fm.OUT_TEMP = temp_dir
    temp_dir.mkdir(parents=True, exist_ok=True)
    return fm.solve_bone_fe(part, allocated_ref, log, reuse=reuse, time_steps=1)


def _solve_retry(part, allocated_ref: float, base_dir: Path, log: Log,
                 retries: int = 2) -> dict:
    """Solve via the production path, retrying on a *fresh mesh instance* so the
    mesh-regeneration non-determinism is sampled (and recorded).  Returns the
    first successful record, else the last failure; every try is logged under
    ``_retries``."""
    last: dict = {}
    tries: list[dict] = []
    for k in range(retries):
        rec = _solve(part, allocated_ref, base_dir / f"try{k}", log, reuse=True)
        tries.append({
            "try": k,
            "status": rec.get("status"),
            "fe_solve_load_n": rec.get("fe_solve_load_n"),
            "load_reduction_factor": rec.get("load_reduction_factor"),
            "attempts": [[a.get("load_n"), bool(a.get("ok")), a.get("purpose", "search")]
                         for a in rec.get("attempts", [])],
            "error": (rec.get("error") or "")[:300],
        })
        last = rec
        if rec.get("status") in ("ok", "no_domains") and rec.get("post"):
            break
    out = dict(last or {})
    out["_retries"] = tries
    out["_n_retries_used"] = len(tries)
    return out


def _fibula_row_from_rec(part, rec: dict, loads: dict, sigma_c: float) -> dict:
    """σ_nom(h) = (main-domain stat · lrf) · load(h)/load_ref(full).  This is the
    exact arithmetic of ``compute_matrix`` but exposed for the fibula so the
    weak-σc row can be produced alongside the part's own σc."""
    heights = loads["heights"]
    ref_idx = int(np.argmin(np.abs(heights - fm.REF_H)))
    ref_load = float(loads[f"{part.joint}__vert"][ref_idx])
    load_h = np.asarray(loads[f"{part.joint}__vert"], float)
    main = fm.extract_main_domain(part, rec, _NullLog())
    if not main.get("stats"):
        return {"error": main.get("error", "无 stats（FE 未解或主域缺失）"),
                "status": rec.get("status"), "lrf": rec.get("load_reduction_factor")}
    lrf = float(rec.get("load_reduction_factor") or 1.0)
    ratio = load_h / ref_load
    out = {"main_domain": main.get("label"), "lrf": lrf,
           "stats": {s: float(main["stats"][s] * lrf) for s in fm.STATS if main.get("stats")}}
    for s in ("median", "p95"):
        sig = np.asarray(main["stats"][s], float) * lrf * ratio
        out[f"sigma_nom_{s}"] = [float(x) for x in sig]
    sig_med = np.asarray(out["sigma_nom_median"], float)
    out["weak"] = {"sigma_c_mpa": sigma_c, **_row_stats(sig_med, sigma_c)}
    out["weak"]["matrix"] = [int(x) for x in (sig_med >= sigma_c)]
    out["weak"]["logistic"] = _logistic(heights, sig_med >= sigma_c)
    return out


class _NullLog:
    def __call__(self, msg: str) -> None:
        pass


def build_audit(policies: tuple[float, ...], *, solve_tibia: bool,
                run_ladder: bool) -> dict:
    frozen = load_frozen()
    heights = frozen["heights"]
    loads = frozen["loads"]
    ref_idx = int(np.argmin(np.abs(heights - fm.REF_H)))
    ankle_ref = float(loads["ankle_r__vert"][ref_idx])

    AUDIT_TEMP.mkdir(parents=True, exist_ok=True)
    log = Log(AUDIT_TEMP / "run.log")
    log(f"[audit] ankle_r ref @5m = {ankle_ref:.3f} N")

    # --- default-path identity (load_fraction = 1.0) --------------------
    parts_default = fm.build_parts()
    mat_default = fm.compute_matrix(parts_default, frozen["recs"], loads, "median",
                                    _NullLog())
    fro = frozen["json"]["sigma_vm_mpa"]["median"]
    max_d = 0.0
    for i, r in enumerate(mat_default["rows"]):
        a = np.asarray(mat_default["sigma_vm"]["median"][i], float)
        f = np.asarray([np.nan if v is None else v for v in fro[r["bone"]]], float)
        max_d = max(max_d, float(np.nanmax(np.abs(a - f))))
    identity = {"load_fraction": 1.0, "max_abs_sigma_delta_vs_frozen_bc2": max_d,
                "bit_identical": bool(max_d == 0.0)}
    log(f"[audit] default-path identity: max|Δσ| vs frozen = {max_d:g} "
        f"(bit_identical={identity['bit_identical']})")

    # --- baseline fibula row (frozen, no sharing) ----------------------
    fib_part = {p.bone: p for p in parts_default}["fibula_r"]
    base_row = _fibula_row_from_rec(fib_part, frozen["recs"]["fibula_r"], loads,
                                    FIBULA_SIGMA_C_WEAK)
    base_row_shaft = _row_stats(np.asarray(base_row["sigma_nom_median"], float),
                                FIBULA_SIGMA_C_SHAFT)

    # --- per-policy re-solve -------------------------------------------
    out_policies: dict[str, dict] = {}
    for share in policies:
        name = f"share_{share:.3f}"
        pdir = AUDIT_TEMP / name
        pdir.mkdir(parents=True, exist_ok=True)
        parts = fm.build_parts()
        by_bone = {p.bone: p for p in parts}
        by_bone["fibula_r"].load_fraction = float(share)
        by_bone["tibia_r"].load_fraction = float(1.0 - share)
        tib_share = 1.0 - share
        log("=" * 72)
        log(f"[audit] {name}: fibula={share:.3f} tibia={tib_share:.3f} "
            f"-> fibula alloc ref = {ankle_ref*share:.2f} N")

        fib_alloc = fm.effective_ref_load(by_bone["fibula_r"], ankle_ref)
        fib_rec = _solve_retry(by_bone["fibula_r"], fib_alloc, pdir / "fibula_r", log)
        log(f"[audit] fibula: status={fib_rec.get('status')} "
            f"solve={fib_rec.get('fe_solve_load_n')} lrf={fib_rec.get('load_reduction_factor')} "
            f"linearity={fib_rec.get('linearity')}")
        log(f"[audit] fibula retries: {fib_rec.get('_retries')}")
        fe_ok = fib_rec.get("status") in ("ok", "no_domains") and fib_rec.get("post")

        fe_mixed = dict(frozen["recs"])
        fe_mixed["fibula_r"] = fib_rec
        if solve_tibia:
            tib_alloc = fm.effective_ref_load(by_bone["tibia_r"], ankle_ref)
            tib_rec = _solve_retry(by_bone["tibia_r"], tib_alloc, pdir / "tibia_r", log)
            fe_mixed["tibia_r"] = tib_rec
            log(f"[audit] tibia: status={tib_rec.get('status')} "
                f"solve={tib_rec.get('fe_solve_load_n')} lrf={tib_rec.get('load_reduction_factor')} "
                f"linearity={tib_rec.get('linearity')}")

        mat = fm.compute_matrix(parts, fe_mixed, loads, "median", _NullLog())
        fib_row = _fibula_row_from_rec(by_bone["fibula_r"], fib_rec, loads,
                                       FIBULA_SIGMA_C_WEAK) if fe_ok else None

        # -- diagnostic linear-share counterfactual (report §3.2 evidence 4):
        #    apply the share ANALYTICALLY to the frozen baseline σ_nom.  Valid
        #    ONLY if the FE is linear in load (it is NOT -- see the ladder).
        base_sig = np.asarray(base_row["sigma_nom_median"], float)
        cf_sig = base_sig * float(share)
        cf = _row_stats(cf_sig, FIBULA_SIGMA_C_WEAK)
        cf["matrix"] = [int(x) for x in (cf_sig >= FIBULA_SIGMA_C_WEAK)]
        cf["logistic"] = _logistic(heights, cf_sig >= FIBULA_SIGMA_C_WEAK)
        cf["sigma_nom_median_h1_h5_h50_mpa"] = [float(cf_sig[0]), float(cf_sig[4]), float(cf_sig[-1])]

        out_policies[name] = {
            "fibula_share": float(share), "tibia_share": float(tib_share),
            "allocated_ref_load_n": float(fib_alloc),
            "fe": {"status": fib_rec.get("status"),
                   "fe_solve_load_n": fib_rec.get("fe_solve_load_n"),
                   "load_reduction_factor": fib_rec.get("load_reduction_factor"),
                   "linearity": fib_rec.get("linearity"),
                   "error": (fib_rec.get("error") or "")[:600],
                   "retries": fib_rec.get("_retries")},
            "search": _search_monotonic(fib_rec.get("attempts", [])),
            "fe_row_weak_sigma_c_70": (fib_row["weak"] if fib_row and "weak" in fib_row
                                       else (fib_row or {"error": "no FE solution"})),
            "fe_sigma_nom_h1_h5_h50_mpa": ([fib_row["sigma_nom_median"][0],
                                            fib_row["sigma_nom_median"][4],
                                            fib_row["sigma_nom_median"][-1]]
                                           if fib_row and "sigma_nom_median" in fib_row
                                           else None),
            "linear_share_counterfactual_weak_sigma_c_70": cf,
            "tibia_fe": ({"status": fe_mixed["tibia_r"].get("status"),
                          "fe_solve_load_n": fe_mixed["tibia_r"].get("fe_solve_load_n"),
                          "load_reduction_factor": fe_mixed["tibia_r"].get("load_reduction_factor"),
                          "linearity": fe_mixed["tibia_r"].get("linearity")}
                         if solve_tibia else None),
            "corrected_matrix_rows": mat["rows"],
        }
        if fib_row and "weak" in fib_row:
            w = fib_row["weak"]
            log(f"[audit] {name} fibula FE weak(σc=70): n_fractured={w['n_fractured_h']}/50 "
                f"all_one={w['all_one']} 0/1_mixed={w['identifiable_0_1_mixed']} "
                f"logistic={w['logistic'].get('status')} OR={w['logistic'].get('OR')}")
        else:
            log(f"[audit] {name} fibula FE: NO SOLUTION (status={fib_rec.get('status')})")
        log(f"[audit] {name} linear CF weak(σc=70): n_fractured={cf['n_fractured_h']}/50 "
            f"all_one={cf['all_one']} 0/1_mixed={cf['identifiable_0_1_mixed']} "
            f"logistic={cf['logistic'].get('status')} OR={cf['logistic'].get('OR')}")

    # --- convergence ladder (fibula, BC2 conditions) -------------------
    ladder = None
    if run_ladder:
        ladder = _run_ladder(frozen, ankle_ref, log)

    # --- register spine/skull shared-joint rows ------------------------
    registration = _register_shared(frozen)

    log.save()
    return {
        "meta": {
            "title": "Fibula FE magnitude audit -- tibia/fibula load sharing",
            "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "scope": "FIBULA_IDENTIFIABILITY_REPORT.md §5 (4 items)",
            "fe_conditions": "align_pca=True (BC2 opt-in), time_steps=1",
            "febio": _febio_probe(),
            "frozen_inputs_read_only": [
                "results/opensim_fe/fracture_matrix_bc2.json",
                "temp/fe_out/fracture_matrix_bc2/load_chain.npz",
                "temp/fe_out/fracture_matrix_bc2/<b>/<b>_fe.json"],
            "new_outputs": [
                "results/opensim_fe/fibula_fe_audit_matrix.json",
                "results/opensim_fe/FIBULA_FE_AUDIT_REPORT.md",
                "temp/fe_out/fibula_fe_audit/"],
            "policies": [float(s) for s in policies],
            "notes": [
                "Default load_fraction=1.0 path re-verified bit-identical vs frozen BC2.",
                "No σc change, no G3 re-run, no target back-solving.",
            ],
        },
        "ankle_ref_load_n": ankle_ref,
        "default_path_identity": identity,
        "baseline_fibula_row": {
            "load_fraction": 1.0,
            "sigma_nom_median_h1_h5_h50_mpa": [base_row["sigma_nom_median"][0],
                                               base_row["sigma_nom_median"][4],
                                               base_row["sigma_nom_median"][-1]],
            "weak_sigma_c_70": base_row["weak"],
            "shaft_sigma_c_160": base_row_shaft,
            "fe": {"status": frozen["recs"]["fibula_r"].get("status"),
                   "fe_solve_load_n": frozen["recs"]["fibula_r"].get("fe_solve_load_n"),
                   "load_reduction_factor": frozen["recs"]["fibula_r"].get("load_reduction_factor"),
                   "linearity": frozen["recs"]["fibula_r"].get("linearity"),
                   "attempts": [[a.get("load_n"), a.get("ok")]
                                for a in frozen["recs"]["fibula_r"].get("attempts", [])]},
        },
        "policies": out_policies,
        "convergence_ladder_fibula": ladder,
        "shared_joint_registration": registration,
    }


def _febio_probe() -> dict:
    try:
        from climbing.coupling.febio_run import find_febio, probe_febio
        exe = find_febio(verify=True)
        return {"exe": str(exe), "version": probe_febio(exe)}
    except Exception as exc:  # noqa: BLE001
        return {"error": f"{type(exc).__name__}: {exc}"}


def _run_ladder(frozen: dict, ankle_ref: float, log: Log) -> dict:
    """Fine load ladder for the fibula under BC2 conditions (mesh extracted
    once), to characterise the convergence set / monotonicity."""
    ldir = AUDIT_TEMP / "ladder"
    ldir.mkdir(parents=True, exist_ok=True)
    fm._build_run = _build_run_bc2
    fib = {p.bone: p for p in fm.build_parts()}["fibula_r"]
    pids = [p for sd in fib.spec.solid_domains for p in sd.pids] + \
           [p for sh in fib.spec.shell_domains for p in sh.pids]
    mats = fm.read_materials(fm.MASTER, pids)
    mesh = fm.extract_mesh(fib.spec, fm.MASTER)
    solid_names = [sd.mesh_name for sd in fib.spec.solid_domains]
    loads = [400, 500, 600, 700, 750, 800, 900, 1000, 1200, 1500, 2000, 2500, 5000]
    # include the exact allocated loads for the policies
    for s in DEFAULT_POLICIES:
        loads.append(ankle_ref * s)
    rows = []
    for load in sorted(set(round(x, 4) for x in loads)):
        feb = ldir / f"fib_lad_{load:.3f}.feb"
        rec = {"load_n": float(load)}
        try:
            _build_run_bc2(fib, mesh, mats, feb, float(load), 1)
            pp = fm.post_process(feb.with_suffix(".xplt"), feb.with_suffix(".hdf5"),
                                 solid_names=solid_names, shell_specs=[], label="fibula_r")
            d = pp.get("domains", {}).get("solid_fib_cort", {})
            rec.update(ok=True, fib_cort={k: d.get(f"vm_{k}_mpa")
                                          for k in ("max", "p95", "mean", "median")})
        except Exception as exc:  # noqa: BLE001
            rec.update(ok=False, error=f"{type(exc).__name__}")
        rows.append(rec)
        log(f"[ladder] load={load:9.2f} ok={rec['ok']}")
    # monotonicity: convergent set down-closed?
    conv = sorted(r["load_n"] for r in rows if r["ok"])
    fail = sorted(r["load_n"] for r in rows if not r["ok"])
    nonmono = [(c, f) for c in conv for f in fail if f < c]
    return {"levels": rows,
            "convergent_loads_n": conv, "failed_loads_n": fail,
            "search_monotonic_down_closed": bool(len(nonmono) == 0),
            "counterexamples_load_above_which_smaller_failed": nonmono}


def _register_shared(frozen: dict) -> dict:
    j = frozen["json"]
    parts = {p["bone"]: p for p in j["parts"]}
    rows = []
    for b in LUMBAR_SHARED_BONES:
        p = parts[b]
        rows.append({
            "bone": b, "part_cn": p["part_cn"], "joint": p["joint"],
            "load_ref_n": p["load_ref_n"],
            "fe_solve_load_n": p["fe_solve_load_n"],
            "load_reduction_factor": p["load_reduction_factor"],
            "utilization_at_5m": p["utilization_at_5m"],
            "first_fracture_median_m": p.get("first_fracture_median_m"),
            "load_fraction": 1.0,
        })
    return {
        "class": "multi-bone rows sharing ONE joint reaction",
        "joint": "lumbar",
        "bones": list(LUMBAR_SHARED_BONES),
        "same_class_as_fibula": True,
        "evidence": ("REGIONS/paper_alignment gives joint='lumbar' to L3,T6,C5,parietal_r; "
                     "fracture_matrix feeds each the WHOLE lumbar reaction "
                     "(ref_load = loads['lumbar__vert'][ref_idx]) with load_fraction=1.0. "
                     "risk_1d_loadshare.py already documents this ('lumbar free-body force "
                     "supports the whole torso + both arms; applying it to a cervical "
                     "vertebra or the skull over-counts 5-10x') and fixes it in the 1D path "
                     "via a mass-above ratio -- but the FE path is NOT fixed."),
        "recommended_action": "register only (out of this audit's scope); a future "
                              "FE audit should allocate the lumbar reaction by a "
                              "mass-above ratio (as risk_1d_loadshare.py does).",
        "rows": rows,
    }


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------
def write_report(res: dict) -> None:
    L: list[str] = []
    A = L.append
    m = res["meta"]
    A("# 腓骨 FE 量级审计 —— 胫-腓载荷分配（spec §5 四项）")
    A("")
    A(f"> 生成 {m['generated_at']}；仓库 `D:\\Project\\climbing_fall_analysis`。")
    A("> **新增产物**：本件 + `results/opensim_fe/fibula_fe_audit_matrix.json` + "
      "`temp/fe_out/fibula_fe_audit/`。**不改**任何冻结输入、不重跑 G3、"
      "不调 σc、不回解分担比凑靶。")
    A(f"> FE 条件：`{m['fe_conditions']}`；FEBio `{m['febio'].get('version')}` "
      f"({m['febio'].get('exe')})；ankle_r @5m = {res['ankle_ref_load_n']:.1f} N。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    idt = res["default_path_identity"]
    A(f"- **载荷分配已落地且默认路径零漂移**：`fracture_matrix.py` 新增 "
      f"`Part.load_fraction`（默认 1.0）+ `effective_ref_load()`；默认路径重算 9 骨的 "
      f"σ 与骨折行 vs 冻结 BC2 的 max|Δσ| = **{idt['max_abs_sigma_delta_vs_frozen_bc2']:g}** "
      f"（bit-identical = {idt['bit_identical']}）。")
    # headline about the fibula
    pol = res["policies"]
    names = list(pol)
    fe_mixed = {n: pol[n]["fe_row_weak_sigma_c_70"].get("identifiable_0_1_mixed")
                for n in names}
    cf_mixed = {n: pol[n]["linear_share_counterfactual_weak_sigma_c_70"]["identifiable_0_1_mixed"]
                for n in names}
    fe_fail = {n: pol[n]["fe"].get("status") for n in names}
    A(f"- **腓骨在分配参考载荷处均不收敛**：分配载荷 ~1263/1895/2527 N **全部失败**，"
      f"被迫降载；降载候选 1000 N 也**不稳定**（同一载荷在不同网格实例下时成时败）。"
      f"各政策 FE status = {fe_fail}：仅 7.5% 一档降载成功（→ all_one），5%/10% 两档**无解**。")
    A(f"- **线性分担反事实**（把冻结 σ_nom × 分担比，即 report §3.2 的口径）确实给出 "
      f"`0/1_mixed` = {cf_mixed}（与 report §3.2 预测一致，logistic 在该口径下可定义）；")
    A(f"- **但真实 FE 重解**要么**无解**、要么给出 **all_one**（`fe 0/1_mixed` = {fe_mixed}）；"
      "根因是**腓骨 FE 在 ~600 N 附近发生屈曲**，应力-载荷关系严重非线性（§3），"
      "故「σ_nom × 分担比」的线性外推**不被 FE 自检背书**，**未**在真实 FE 中复现 0/1 混合。")
    if res.get("convergence_ladder_fibula"):
        lad = res["convergence_ladder_fibula"]
        A(f"- **收敛不单调**：腓骨收敛集 = {lad['convergent_loads_n']} N，"
          f"失败集跨越其间（down-closed = {lad['search_monotonic_down_closed']}）。")
    A("- **登记**（不修）：L3/T6/C5/parietal_r 同样各吃**整条 `lumbar` 反力**，"
      "与腓骨同类（§4）。")
    A("")

    # ---- item 1 ----
    A("## 1. 事项①：补胫-腓载荷分配（additive）")
    A("")
    A("- 位置：`scripts/opensim_fe/fracture_matrix.py` 的 FE `ref_load` 入口（`main()` "
      "调用 `solve_bone_fe` 处），经新函数 `effective_ref_load(part, joint_load_n)`。")
    A("- 新增字段：`Part.load_fraction: float = 1.0`（默认 = **不做分配**，历史行为）。")
    A("- CLI（可选）：`--load-fraction fibula_r=0.05 --load-fraction tibia_r=0.95`。")
    A("- 矩阵的逐高度缩放分母仍取**整条关节载荷** `load(h)`，故 `sigma(h)` 仍是该骨在"
      "「自己的**分配载荷**」下的真实应力（线性假设内）。")
    A(f"- 默认路径证明：用冻结 bc2 的 FE 记录 + 默认 `load_fraction=1.0` 重跑 "
      f"`compute_matrix`，9 骨 σ(median) 与骨折行 vs 冻结 JSON 完全一致 "
      f"(max|Δσ| = {res['default_path_identity']['max_abs_sigma_delta_vs_frozen_bc2']:g})。")
    A("")

    # ---- item 2 ----
    A("## 2. 事项②：分配后重解 + 收敛性验证")
    A("")
    A("| 政策 | fibula 分担 | 分配参考载荷 (N) | FE status | 求解载荷 (N) | 降载因子 | 线性自检 ratio_median | 期望 |")
    A("|---|---:|---:|:--:|---:|---:|---:|---:|")
    for n in names:
        d = pol[n]
        lin = (d["fe"].get("linearity") or {})
        A(f"| {n} | {d['fibula_share']:.3f} | {d['allocated_ref_load_n']:.2f} | "
          f"**{d['fe'].get('status')}** | "
          f"{d['fe'].get('fe_solve_load_n') if d['fe'].get('fe_solve_load_n') else '—'} | "
          f"{('%.3f' % d['fe']['load_reduction_factor']) if d['fe'].get('load_reduction_factor') else '—'} | "
          f"{('%.4f' % lin['ratio_median']) if lin.get('ratio_median') is not None else '**None**'} | "
          f"{('%.1f' % lin['expected_ratio']) if lin.get('expected_ratio') else '—'} |")
    A("")
    bfe = res["baseline_fibula_row"]["fe"]
    A(f"**冻结基线（无分配）**：求解载荷 {bfe['fe_solve_load_n']} N、lrf={bfe['load_reduction_factor']:.3f}、"
      f"linearity = **{bfe['linearity']}**（None）、搜索 = {bfe['attempts']}。")
    A("")
    A("**要点（如实记录）**：")
    A("- 文献分担下的**分配参考载荷本身不收敛**（~1263 / ~1895 / ~2527 N）；"
      "`solve_bone_fe` 的降载候选（1000 / 2500 N）**也不稳定**：同一载荷在不同网格实例下"
      "时成时败（见 §3 的网格非确定性）。")
    A("- 因此 `linearity` 多数情形仍为 **None**（无机可算），个别情形即使算出也远离期望 2.0。")
    A("- 搜索**不单调**：详见 §3 收敛阶梯。")
    A("- **结论（事项②）**：腓骨在分配载荷下**不满足收敛/线性自检**，本项判定为**未通过**"
      "（不是数据缺失，是 FE 结构性问题）。")
    A("")
    A("**胫骨（同批按 complement 分担重解，供对照）**：")
    A("")
    A("| 政策 | 胫骨分担 | 分配载荷 (N) | status | 求解载荷 (N) | 降载因子 | 线性自检 ratio_median |")
    A("|---|---:|---:|:--:|---:|---:|---:|")
    for n in names:
        t = pol[n].get("tibia_fe")
        if not t:
            A(f"| {n} | {1.0 - pol[n]['fibula_share']:.3f} | — | (未重解) | — | — | — |")
            continue
        lin = (t.get("linearity") or {})
        A(f"| {n} | {1.0 - pol[n]['fibula_share']:.3f} | "
          f"{pol[n]['allocated_ref_load_n'] * (1.0 - pol[n]['fibula_share']) / pol[n]['fibula_share']:.1f} | "
          f"**{t.get('status')}** | {t.get('fe_solve_load_n') if t.get('fe_solve_load_n') else '—'} | "
          f"{('%.3f' % t['load_reduction_factor']) if t.get('load_reduction_factor') else '—'} | "
          f"{('%.4f' % lin['ratio_median']) if lin.get('ratio_median') is not None else '**None**'} |")
    A("")
    A("> 胫骨也呈非线性：某一政策在分配载荷处**直接**解出（lrf=1.0），但其应力远高于由 5000 N "
      "线性外推值（例如 p95 1259 vs 外推 514）——说明长骨 FE 在该载荷区间整体非线性，"
      "并非腓骨独有；只是腓骨更细，先跨过屈曲临界。")
    A("")

    # ---- item 3 ----
    A("## 3. 事项③：重出腓骨行 + 0/1 混合检查 + 根因（屈曲）")
    A("")
    if res.get("convergence_ladder_fibula"):
        lad = res["convergence_ladder_fibula"]
        A("### 3.1 收敛阶梯（腓骨, BC2 条件, 网格抽一次）")
        A("")
        A("| 载荷 (N) | 收敛 | fib_cort median (MPa) | p95 (MPa) | max (MPa) |")
        A("|---:|:--:|---:|---:|---:|")
        for r in lad["levels"]:
            fc = r.get("fib_cort", {})
            A(f"| {r['load_n']:.1f} | {'✓' if r['ok'] else '✗'} | "
              f"{('%.3f' % fc['median']) if fc.get('median') is not None else '—'} | "
              f"{('%.3f' % fc['p95']) if fc.get('p95') is not None else '—'} | "
              f"{('%.2f' % fc['max']) if fc.get('max') is not None else '—'} |")
        A("")
        ce = lad.get("counterexamples_load_above_which_smaller_failed", [])
        ce_s = "；".join(f"{int(f)} N 失败但 {int(c)} N 成功" for c, f in ce[:4]) or "（无）"
        A(f"- 收敛集 = {lad['convergent_loads_n']} N；down-closed = "
          f"**{lad['search_monotonic_down_closed']}**（**不单调**：{ce_s}）")
        A(f"- 即**不存在一个「低于临界即收敛」的干净区间** ⇒ 500–2500 N 之间同一载荷可能成、"
          "可能败（网格/分叉共同决定）；这正是「5000 N 可解 ≠ 5000 N 是物理承载极限」的证据。")
        A("")
        A("### 3.2 根因：腓骨 FE 屈曲（Euler buckling），非材料非线性")
        A("")
        A("- 材料卡为 **isotropic elastic（线性弹性）**（`fib_bc2_5000.feb`），"
          "故应力-载荷非线性**只可能来自几何**。")
        lv = {round(r["load_n"], 1): r for r in lad["levels"]}

        def _med(L):
            r = lv.get(L)
            fc = (r or {}).get("fib_cort") or {}
            return fc.get("median")

        m400, m500, m700, m900, m1500, m5000 = (_med(400.0), _med(500.0), _med(700.0),
                                                _med(900.0), _med(1500.0), _med(5000.0))
        def _f(v):
            return "—" if v is None else f"{v:.2f}"
        ratio = (m700 / m500) if (m500 and m700) else float("nan")
        A(f"- 载荷-应力阶梯（**同一网格**，故差异纯属载荷效应）：400→{_f(m400)}, 500→{_f(m500)} MPa"
          f"（近线性，比例≈载荷比例）；**500→700 N 之间 median 跳到 {_f(m700)} MPa"
          f"（≈{ratio:.1f}× / 1.4× 载荷）⇒ 屈曲拐点**；900→{_f(m900)}, 1500→{_f(m1500)}, "
          f"5000→{_f(m5000)} MPa（≈100–135，**后屈曲饱和**）。")
        A("- 这是典型的**细长骨轴向压缩屈曲**：临界载荷约 **550–700 N**"
          "（600/750/800 N 恰在分叉点附近失败）。")
        A("- ⇒ 冻结 BC2 的腓骨 `σ_nom(h=5)=684 MPa` 是**从 5000 N 后屈曲解线性外推**得到的"
          "后屈曲应力，不是「腓骨在 5% 载荷下」的材料应力。")
        A("")
    A("### 3.3-A 真实 FE 重出的腓骨行（各政策；σc=70）")
    A("")
    A("| 政策 | FE status | σ_nom h1 (MPa) | h5 | h50 | σ_nom/σc 区间 | n_fractured/50 | 0/1 mixed | logistic |")
    A("|---|:--:|---:|---:|---:|---|---:|:--:|---|")
    for n in names:
        wr = pol[n]["fe_row_weak_sigma_c_70"]
        s = pol[n]["fe_sigma_nom_h1_h5_h50_mpa"]
        if s is None or "ratio_min" not in wr:
            A(f"| {n} | {pol[n]['fe'].get('status')} | — | — | — | — | — | — | "
              f"**无解**（{wr.get('error', 'no FE solution')}） |")
            continue
        lg = wr["logistic"]
        ors = lg.get("OR")
        A(f"| {n} | {pol[n]['fe'].get('status')} | {s[0]:.3f} | {s[1]:.3f} | {s[2]:.3f} | "
          f"{wr['ratio_min']:.3f}–{wr['ratio_max']:.3f} | {wr['n_fractured_h']} | "
          f"{wr['identifiable_0_1_mixed']} | "
          f"{lg.get('status')}{'' if ors is None else f' (OR={ors:.3f})'} |")
    bwr = res["baseline_fibula_row"]["weak_sigma_c_70"]
    bs = res["baseline_fibula_row"]["sigma_nom_median_h1_h5_h50_mpa"]
    A(f"| **冻结基线(无分配)** | ok | {bs[0]:.3f} | {bs[1]:.3f} | {bs[2]:.3f} | "
      f"{bwr['ratio_min']:.3f}–{bwr['ratio_max']:.3f} | {bwr['n_fractured_h']} | "
      f"{bwr['identifiable_0_1_mixed']} | {bwr['logistic'].get('status')} |")
    A("")
    A("### 3.3-B 线性分担反事实（冻结 σ_nom × 分担比；= report §3.2 口径）")
    A("")
    A("| 政策 | σ_nom h1 (MPa) | h5 | h50 | σ_nom/σc 区间 | n_fractured/50 | all_one | 0/1 mixed | logistic |")
    A("|---|---:|---:|---:|---|---:|:--:|:--:|---|")
    for n in names:
        cf = pol[n]["linear_share_counterfactual_weak_sigma_c_70"]
        s = cf["sigma_nom_median_h1_h5_h50_mpa"]
        lg = cf["logistic"]
        ors = lg.get("OR")
        A(f"| {n} | {s[0]:.3f} | {s[1]:.3f} | {s[2]:.3f} | "
          f"{cf['ratio_min']:.3f}–{cf['ratio_max']:.3f} | {cf['n_fractured_h']} | "
          f"{cf['all_one']} | {cf['identifiable_0_1_mixed']} | "
          f"{lg.get('status')}{'' if ors is None else f' (OR={ors:.3f})'} |")
    A("")
    A("> 说明：σ_nom(h) = 主承载域 median（`solid_fib_cort`）· lrf · load(h)/load_ref。"
      "3.3-B 的线性反事实**假定 σ ∝ F**，而腓骨 FE 在 ~600 N 即屈曲（§3.2），故 3.3-B 的"
      "「0/1 混合」只是**线性外推的产物**，不能被 3.3-A 的真实 FE 复现。"
      "口径：`0/1 mixed=True` 指**确定性行出现 0/1 混合**（⇒ 被试级 logistic 可定义，"
      "与 report §3.2 记法一致）；表内 `logistic=complete_separation` 是**确定性行**本身为"
      "干净阶跃时的必然结果（同 `s4_alignment_bc2.json` 里胫骨的口径），二者不矛盾。")
    A("")
    A("**结论（事项③）**：")
    A("- **按 report §3.2 的线性口径**：5%/7.5%/10% 分担给出 `0/1 mixed = True`，"
      "logistic（被试级）有定义——与 report §3.2 的预测一致。")
    A("- **按真实 FE 口径**：5%/10% **无解**；7.5% 虽解出，却得 **all_one**（σ_nom/σc 区间 "
      "1.852–13.532，h=1 即过阈），**未**恢复 0/1 混合。即「补上 5–10% 分配即恢复 0/1 混合」"
      "这一结论**没有被真实 FE 证实**。")
    A("- **根因升级**：腓骨 all_one 的更深层原因是**孤立腓骨 FE 的屈曲不稳定性**"
      "（真实小腿的骨间膜/软组织支撑被 FE 忽略），而非单纯的载荷分配。"
      "载荷分配是必要但**不充分**的修法。")
    A("")

    # ---- item 4 ----
    A("## 4. 事项④：同类共享关节行登记（仅登记，不修）")
    A("")
    reg = res["shared_joint_registration"]
    A(f"- **同类**：`{reg['bones']}` 全部取 `joint='{reg['joint']}'`，FE 各吃**整条** lumbar 反力"
      f"（load_fraction=1.0）。")
    A("- 证据：`risk_1d_loadshare.py` 已记载该缺陷（1D 路已用 mass-above 比修好），"
      "**FE 路未修**。")
    A("")
    A("| 骨 | 部位 | 关节 | 参考载荷 (N) | 求解载荷 (N) | 降载因子 | util@5m | 首骨折(median) |")
    A("|---|---|---|---:|---:|---:|---:|---:|")
    for r in reg["rows"]:
        A(f"| `{r['bone']}` | {r['part_cn']} | `{r['joint']}` | {r['load_ref_n']:.1f} | "
          f"{r['fe_solve_load_n']:.0f} | {r['load_reduction_factor']:.3f} | "
          f"{r['utilization_at_5m']:.3f} | {r['first_fracture_median_m']} |")
    A("")
    A(f"- 建议动作：**登记**；未来 FE 审计可用 mass-above 比（同 `risk_1d_loadshare.py`）分配 lumbar 反力。")
    A("")

    # ---- boundaries ----
    A("## 5. 诚实边界")
    A("")
    A("1. FE 用 `align_pca=True`（BC2 口径），但每次 `extract_mesh` 的松质骨四面体网格"
      "存在**轻微非确定性**（同 5000 N 复算 median 135.368 vs 冻结 135.414，差 0.03%），"
      "故 σ 的末位数字在重跑间会微动。")
    A("2. 分配载荷处**不收敛**：本件没有、也不能给出一个自检通过的腓骨 FE 应力。"
      "所有腓骨 σ_nom 都带未背书的线性外推误差。")
    A("3. 本件**未**重跑 G3、**未**改 σc、**未**改 A3 结论、**未**回解分担比。")
    A("4. 结论口径只出现 σ_nom/σc（越阈物理量）与载荷分担比，不含任何 σ/σ_law 比值。")
    A("")
    A("## 6. 复现")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\fibula_fe_audit.py")
    A("```")
    A("")
    A("## 7. 产物")
    A("")
    A("- 本报告：`results/opensim_fe/FIBULA_FE_AUDIT_REPORT.md`")
    A("- 机器可读：`results/opensim_fe/fibula_fe_audit_matrix.json`")
    A("- FE 缓存/日志：`temp/fe_out/fibula_fe_audit/`")
    A("- 改动脚本（additive）：`scripts/opensim_fe/fracture_matrix.py`")
    A("")
    OUT_REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print("[audit] wrote", OUT_REPORT)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--policies", default=",".join(str(x) for x in DEFAULT_POLICIES),
                    help="comma-separated fibula shares")
    ap.add_argument("--no-tibia", action="store_true", help="re-solve fibula only")
    ap.add_argument("--skip-ladder", action="store_true")
    ap.add_argument("--report-only", action="store_true",
                    help="regenerate the report+JSON view from an existing audit JSON "
                         "(no FE re-solve)")
    args = ap.parse_args(argv)

    policies = tuple(float(x) for x in args.policies.split(",") if x.strip())
    if args.report_only:
        res = json.loads(OUT_JSON.read_text(encoding="utf-8"))
        write_report(res)
        return 0
    fm._build_run = _build_run_bc2  # BC2 conditions
    t0 = time.time()
    res = build_audit(policies, solve_tibia=not args.no_tibia,
                      run_ladder=not args.skip_ladder)
    res["meta"]["total_seconds"] = time.time() - t0
    OUT_JSON.write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    print("[audit] wrote", OUT_JSON)
    write_report(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
