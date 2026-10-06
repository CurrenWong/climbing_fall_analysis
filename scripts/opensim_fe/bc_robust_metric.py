"""bc_robust_metric.py -- BC-robust cross-check of the 9-bone relative fracture ranking.

Additive diagnostic only.  Does **not** modify any existing module and does not
re-run FEBio for the general case: it re-analyses artifacts that already exist.

Motivation
----------
``results/opensim_fe/MESH_FIX_REPORT.md`` diagnosed that the 9x50 sigma_vm matrix
is **boundary-condition limited**, not mesh limited: the oblique global load is
applied along the bbox longest axis (tibia/fibula tilted ~42 deg -> bending), and
the end bands are fully clamped.  Because the artifact is bone-shape dependent it
can corrupt the *relative* ranking that S4 (Logistic + hierarchical clustering)
is supposed to consume.  This script tests that relative pattern by two
independent, BC-robust routes.

Route 1 (FE-free, 1D strength cross-section)
    risk(h) = F_joint(h) / (A_section * sigma_c)
    * F_joint(h)   : ``fracture_matrix.json`` -> ``loads_n`` (per-height, measured).
    * sigma_c      : ``climbing.bone.MATERIAL_STRENGTH_MPA`` (compression), mapped
                     per bone (see ``SIGMA_KEY``).
    * A_section    : solid-domain volume / bbox-longest-axis length, read from the
                     per-bone ``temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5``.
                     For the 4 limb bones this is the CORT solid (reproduces the
                     measured V_cort of ``MESH_FIX_REPORT.md`` section 2.3 exactly).
                     For the 5 axial bones (CORT is a shell, no solid) it is the
                     **estimated** load-bearing body cross-section
                     (SPON body / diploe).  Clearly labelled.
    First-fracture height = smallest h with risk(h) >= 1 (else > 50 m).

Route 2 (existing FE outputs, best-effort BC-robust stress)
    Recompute sigma_vm statistics from the stored per-element stress, EXCLUDING the
    clamped end bands (top/bottom 30% along the loading axis) plus a volume-weighted
    average.  Sources: clean remesh hdf5 (4 limb bones + 3 axial diagnostics) and the
    original 9-bone fracture-matrix hdf5 (re-analysis).  The stored median is also
    reported.  No FEBio rerun.

Outputs (never overwrites existing artifacts):
    results/opensim_fe/bc_robust_metric_route1.json
    results/opensim_fe/bc_robust_metric_route2.json
    results/opensim_fe/bc_robust_metric_result.json   (both + cross-check)

Usage (repo root)::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\bc_robust_metric.py
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from climbing.bone import MATERIAL_STRENGTH_MPA  # noqa: E402  (PYTHONPATH=src)

ROOT = _HERE.parents[1]
FM_JSON = ROOT / "results" / "opensim_fe" / "fracture_matrix.json"
FM_DIR = ROOT / "temp" / "opensim_fe" / "fracture_matrix"
RM_DIR = ROOT / "temp" / "opensim_fe" / "fracture_matrix_remesh"

# --------------------------------------------------------------------------
# hex8 -> 6 tets (main diagonal 0-6).  Verified to reproduce the measured CORT
# element volumes of MESH_FIX_REPORT.md section 2.3 to 4 significant figures.
# --------------------------------------------------------------------------
_HEX_TETS = ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
             (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6))

# --------------------------------------------------------------------------
# Per-bone registry
#   sigma_key  -> MATERIAL_STRENGTH_MPA key (compression = 2nd tuple element)
#   joint      -> fracture_matrix.json['loads_n'] key feeding this bone
#   solid_dom  -> hdf5 domain used for A_section (CORT for limbs, body for axial)
#   a_source   -> measured / estimated provenance (documented in the report)
# --------------------------------------------------------------------------
BONES: list[dict] = [
    dict(bone="calcaneus_r", part_cn="足部", joint="subtalar_r", sigma_key="calcaneus",
         solid_dom="calc_cort", a_kind="cort", a_source="measured(THUMS CORT hex)"),
    dict(bone="tibia_r", part_cn="胫骨", joint="ankle_r", sigma_key="tibia_shaft",
         solid_dom="tibia_cort", a_kind="cort", a_source="measured(THUMS CORT hex)"),
    dict(bone="fibula_r", part_cn="腓骨", joint="ankle_r", sigma_key="fibula_shaft",
         solid_dom="fib_cort", a_kind="cort", a_source="measured(THUMS CORT hex)"),
    dict(bone="femur_r", part_cn="股骨", joint="hip_r", sigma_key="femur_shaft",
         solid_dom="fem_cort", a_kind="cort", a_source="measured(THUMS CORT hex)"),
    dict(bone="R_HIPBONE", part_cn="骨盆", joint="hip_r", sigma_key="pelvis",
         solid_dom="hip_spon", a_kind="body", a_source="estimated(SPON body, CORT is shell)"),
    dict(bone="L3", part_cn="腰椎", joint="lumbar", sigma_key="spine",
         solid_dom="l3_spon", a_kind="body", a_source="estimated(SPON body, CORT is shell)"),
    dict(bone="T6", part_cn="胸椎", joint="lumbar", sigma_key="spine",
         solid_dom="t6_spon", a_kind="body", a_source="estimated(SPON body, CORT is shell)"),
    dict(bone="C5", part_cn="颈椎", joint="lumbar", sigma_key="spine",
         solid_dom="c5_spon", a_kind="body", a_source="estimated(SPON body, CORT is shell)"),
    dict(bone="parietal_r", part_cn="颅骨", joint="lumbar", sigma_key="skull",
         solid_dom="par_diploe", a_kind="body", a_source="estimated(diploe solid, CORT is shell)"),
]

# Route-2 clean-remesh registry (converged solves only; fibula did not converge).
RM_OBLIQUE = [
    dict(bone="calcaneus_r", part_cn="足部", joint="subtalar_r", sigma_key="calcaneus",
         hdf5="calcaneus_r/calcaneus_r_25354N.hdf5", dom="calcaneus_r_solid",
         solve_load_n=25354.0, ref_load_n=25354.0),
    dict(bone="tibia_r", part_cn="胫骨", joint="ankle_r", sigma_key="tibia_shaft",
         hdf5="tibia_r/tibia_r_10000N.hdf5", dom="tibia_r_solid",
         solve_load_n=10000.0, ref_load_n=25266.0),
    dict(bone="femur_r", part_cn="股骨", joint="hip_r", sigma_key="femur_shaft",
         hdf5="femur_r/femur_r_15228N.hdf5", dom="femur_r_solid",
         solve_load_n=15228.0, ref_load_n=15228.0),
]
RM_AXIAL = [
    dict(bone="calcaneus_r", part_cn="足部", joint="subtalar_r", sigma_key="calcaneus",
         hdf5="calcaneus_r_axial/calcaneus_r_axial_25354N.hdf5", dom="calcaneus_r_solid",
         solve_load_n=25354.0, ref_load_n=25354.0),
    dict(bone="tibia_r", part_cn="胫骨", joint="ankle_r", sigma_key="tibia_shaft",
         hdf5="tibia_r_axial/tibia_r_axial_25266N.hdf5", dom="tibia_r_solid",
         solve_load_n=25266.0, ref_load_n=25266.0),
    dict(bone="femur_r", part_cn="股骨", joint="hip_r", sigma_key="femur_shaft",
         hdf5="femur_r_axial/femur_r_axial_15228N.hdf5", dom="femur_r_solid",
         solve_load_n=15228.0, ref_load_n=15228.0),
]


# --------------------------------------------------------------------------
# small numeric helpers
# --------------------------------------------------------------------------
def _tet_vol(P: np.ndarray) -> np.ndarray:
    return np.abs(np.einsum("ij,ij->i",
                            np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                            P[:, 3] - P[:, 0])) / 6.0


def _hex_vol(P: np.ndarray) -> np.ndarray:
    v = np.zeros(len(P))
    for a, b, c, d in _HEX_TETS:
        v += np.abs(np.einsum("ij,ij->i",
                              np.cross(P[:, b] - P[:, a], P[:, c] - P[:, a]),
                              P[:, d] - P[:, a])) / 6.0
    return v


def _von_mises_voigt(s: np.ndarray) -> np.ndarray:
    xx, yy, zz, xy, yz, xz = (s[:, i] for i in range(6))
    return np.sqrt(0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2)
                   + 3.0 * (xy ** 2 + yz ** 2 + xz ** 2))


def _read_domain(path: Path, dom: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (coords (N,3), conn (n,k) 0-based node ids, stress (n,6))."""
    with h5py.File(path, "r") as h:
        conn = h["meshes/0/domains"][dom][:]
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(np.float64)
        sid = sorted(h["states"].keys(), key=int)[-1]
        s = np.asarray(h["states"][sid]["element_data"]["stress"][dom][:], dtype=np.float64)
    return coords, conn, s


def _centroids_weights(coords: np.ndarray, conn: np.ndarray, kind: str):
    """Return (centroids (n,3), weights (n,)) where weight is volume (solid) or area (shell)."""
    ncol = conn.shape[1] - 1
    P = coords[conn[:, 1:]]
    cent = P.mean(axis=1)
    if kind == "solid":
        if ncol == 4:
            w = _tet_vol(P)
        elif ncol == 8:
            w = _hex_vol(P)
        else:
            raise ValueError(f"unsupported solid element width {ncol}")
        return cent, w
    # shell: polygon area (fan), tolerant of repeated/degenerate nodes
    w = np.zeros(len(conn))
    for i, pid in enumerate(conn[:, 1:]):
        poly = coords[pid]
        # drop consecutive duplicates / keep unique
        keep = [0]
        for j in range(1, len(poly)):
            if np.linalg.norm(poly[j] - poly[keep[-1]]) > 1e-9:
                keep.append(j)
        if len(keep) < 3:
            continue
        q = poly[keep]
        area = 0.0
        for j in range(1, len(q) - 1):
            area += np.linalg.norm(np.cross(q[j] - q[0], q[j + 1] - q[0])) / 2.0
        w[i] = area
    return cent, w


def _domain_stats(path: Path, dom: str, kind: str, scale: float,
                  sigma_nom: float, excl_frac: float = 0.30) -> dict:
    """BC-robust per-domain statistics, end bands (top/bottom excl_frac) excluded."""
    coords, conn, s = _read_domain(path, dom)
    vm = _von_mises_voigt(s) * scale
    cent, w = _centroids_weights(coords, conn, kind)
    lo, hi = coords.min(0), coords.max(0)
    span = hi - lo
    ax = int(np.argmax(span))
    top = cent[:, ax] >= hi[ax] - excl_frac * span[ax]
    bot = cent[:, ax] <= lo[ax] + excl_frac * span[ax]
    mid = ~(top | bot)
    if not mid.any():
        # tiny / thin domain: the whole domain lies in the end bands; keep all
        # elements (degenerate exclusion) and flag it via mid_n == n_elem.
        mid = np.ones_like(mid, dtype=bool)
    wsum = float(w.sum())
    volavg = float((vm * w).sum() / wsum) if wsum > 0 else float("nan")
    return dict(
        n_elem=int(len(vm)), axis=int(ax), axis_label="XYZ"[ax],
        bbox_span_mm=[float(x) for x in span],
        all_p95=float(np.percentile(vm, 95)),
        all_median=float(np.median(vm)),
        all_mean=float(vm.mean()),
        mid_n=int(mid.sum()),
        mid_p95=float(np.percentile(vm[mid], 95)),
        mid_median=float(np.median(vm[mid])),
        mid_mean=float(vm[mid].mean()),
        volume_avg=volavg,
        weight_sum=wsum,
        sigma_nom_mpa=float(sigma_nom),
        all_p95_over_nom=float(np.percentile(vm, 95) / sigma_nom) if sigma_nom else None,
        mid_p95_over_nom=float(np.percentile(vm[mid], 95) / sigma_nom) if sigma_nom else None,
        volume_avg_over_nom=float(volavg / sigma_nom) if sigma_nom else None,
        all_median_over_nom=float(np.median(vm) / sigma_nom) if sigma_nom else None,
    )


# --------------------------------------------------------------------------
# Route 1
# --------------------------------------------------------------------------
def route1(loads: dict, heights: list[float]) -> dict:
    rows = []
    for b in BONES:
        hdf5 = FM_DIR / b["bone"] / f"{b['bone']}.hdf5"
        coords, conn, _ = _read_domain(hdf5, b["solid_dom"])
        kind = "solid"
        _, w = _centroids_weights(coords, conn, kind)
        V = float(w.sum())
        span = coords.max(0) - coords.min(0)
        L = float(span.max())
        A = V / L
        _, sigma_c = MATERIAL_STRENGTH_MPA[b["sigma_key"]]
        F = loads[b["joint"]]
        risk = [float(f) / (A * sigma_c) for f in F]
        h1 = next((heights[i] for i, r in enumerate(risk) if r >= 1.0), None)
        thr = A * sigma_c
        rows.append(dict(
            bone=b["bone"], part_cn=b["part_cn"], joint=b["joint"],
            sigma_key=b["sigma_key"], sigma_c_mpa=float(sigma_c),
            solid_dom=b["solid_dom"], a_kind=b["a_kind"], a_source=b["a_source"],
            V_solid_mm3=V, L_bbox_mm=L,
            A_section_mm2=A, force_threshold_n=thr,
            F_at_5m_n=float(F[4]), risk_at_5m=risk[4],
            first_fracture_h=None if h1 is None else int(h1),
            risk_by_height=risk,
        ))
    # rank: earliest fracture first, then higher risk@5m (bones that never fail last)
    order = sorted(rows, key=lambda r: (r["first_fracture_h"] if r["first_fracture_h"]
                                        else 10**9, -r["risk_at_5m"]))
    for i, r in enumerate(order, 1):
        r["route1_rank"] = i
    return dict(
        method="risk(h) = F_joint(h) / (A_section * sigma_c); first h with risk>=1",
        sigma_c_source="src/climbing/bone.py MATERIAL_STRENGTH_MPA (compression)",
        a_section_source=("V_solid_domain / bbox-longest-axis from "
                          "temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5"),
        a_caveat=("limb bones = measured THUMS CORT solid; axial bones = estimated "
                  "body cross-section (SPON/diploe) because THUMS CORT is a shell"),
        rows=rows,
        ranking=[{"rank": r["route1_rank"], "bone": r["bone"],
                  "first_fracture_h": r["first_fracture_h"],
                  "risk_at_5m": r["risk_at_5m"]} for r in order],
    )


# --------------------------------------------------------------------------
# Route 2
# --------------------------------------------------------------------------
def _nom_from_a(a_section: float, ref_load_n: float) -> float:
    return ref_load_n / a_section if a_section else float("nan")


def route2(loads: dict, heights: list[float], a_by_bone: dict,
           matrix_parts: dict) -> dict:
    out: dict = {"excl_frac": 0.30, "clean_remesh_oblique": [], "clean_remesh_axial": [],
                 "original_mesh_all9": [],
                 "note": ("clean-remesh stress is the raw solved value; sigma stored "
                          "at the solve load, scaled here by ref/solve to the reference load. "
                          "The original-mesh matrix stored sigma_vm_ref WITHOUT applying "
                          "load_reduction_factor -- reported as stored_p95 for comparison.")}

    for spec in RM_OBLIQUE + RM_AXIAL:
        p = RM_DIR / spec["hdf5"]
        a = a_by_bone[spec["bone"]]
        nom = _nom_from_a(a, spec["ref_load_n"])
        scale = spec["ref_load_n"] / spec["solve_load_n"]
        st = _domain_stats(p, spec["dom"], "solid", scale, nom, 0.30)
        mp = matrix_parts.get(spec["bone"], {})
        rec = dict(bone=spec["bone"], part_cn=spec["part_cn"], joint=spec["joint"],
                   sigma_key=spec["sigma_key"], sigma_c_mpa=float(MATERIAL_STRENGTH_MPA[spec["sigma_key"]][1]),
                   source=spec["hdf5"], solve_load_n=spec["solve_load_n"],
                   ref_load_n=spec["ref_load_n"], scale=scale,
                   A_section_mm2=a,
                   stored_p95=mp.get("sigma_vm_ref_mpa", {}).get("p95"),
                   stored_median=mp.get("sigma_vm_ref_mpa", {}).get("median"),
                   **st)
        (out["clean_remesh_oblique"] if spec in RM_OBLIQUE
         else out["clean_remesh_axial"]).append(rec)

    # original 9-bone re-analysis (per-bone main domain; multi-shell -> max per stat)
    for b in BONES:
        bone = b["bone"]
        fej = json.loads((FM_DIR / bone / f"{bone}_fe.json").read_text(encoding="utf-8"))
        mesh = fej["mesh"]
        doms = fej["post"]["domains"]
        ref = float(fej.get("ref_load_n", 1.0))
        solve = float(fej.get("fe_solve_load_n", ref))
        scale = ref / solve
        a = a_by_bone[bone]
        nom = _nom_from_a(a, ref)
        c = MATERIAL_STRENGTH_MPA[b["sigma_key"]][1]
        # candidate domains = the stat-bearing main domain(s): CORT solid for limbs, all shells for axial
        cand = []
        if mesh.get("by_shell"):
            for name in mesh["by_shell"]:
                key = f"shell_{name}"
                if key in doms:
                    cand.append((f"shell_{name}", doms[key]["stress_key"].strip("'"), "shell"))
        else:
            # no shell -> solid CORT domain
            for name in mesh.get("by_solid", {}):
                if name in ("calc_cort", "tibia_cort", "fib_cort", "fem_cort"):
                    cand.append((f"solid_{name}", name, "solid"))
        recs = []
        for label, h5dom, kind in cand:
            st = _domain_stats(FM_DIR / bone / f"{bone}.hdf5", h5dom, kind, scale, nom, 0.30)
            recs.append(dict(domain=label, hdf5_domain=h5dom, **st))
        # matrix rule: for p95/median take the max across candidate domains
        best_p95 = max(r["all_p95"] for r in recs)
        best_p95_mid = max(r["mid_p95"] for r in recs)
        best_med = max(r["all_median"] for r in recs)
        best_mid_med = max(r["mid_median"] for r in recs)
        best_volavg = max(r["volume_avg"] for r in recs)
        mp = matrix_parts.get(bone, {})
        stored_p95 = mp.get("sigma_vm_ref_mpa", {}).get("p95")
        stored_median = mp.get("sigma_vm_ref_mpa", {}).get("median")
        out["original_mesh_all9"].append(dict(
            bone=bone, part_cn=b["part_cn"], joint=b["joint"], sigma_key=b["sigma_key"],
            sigma_c_mpa=float(c), ref_load_n=ref, solve_load_n=solve, scale=scale,
            A_section_mm2=a, sigma_nom_mpa=nom, n_candidate_domains=len(recs),
            all_p95=best_p95, all_p95_solved=best_p95 / scale,
            mid_p95=best_p95_mid, all_median=best_med, mid_median=best_mid_med,
            volume_avg=best_volavg,
            stored_p95=stored_p95, stored_median=stored_median,
            all_p95_over_nom=best_p95 / nom if nom else None,
            mid_p95_over_nom=best_p95_mid / nom if nom else None,
            volume_avg_over_nom=best_volavg / nom if nom else None,
            all_median_over_nom=best_med / nom if nom else None,
            all_p95_over_sigc=best_p95 / c,
            mid_p95_over_sigc=best_p95_mid / c,
            volume_avg_over_sigc=best_volavg / c,
            all_median_over_sigc=best_med / c,
            stored_p95_over_sigc=(stored_p95 / c) if stored_p95 is not None else None,
            domains=recs,
        ))

    # ranking of the original-9 by the two BC-robust-ish FE metrics
    def _rank(rows, key):
        order = sorted(rows, key=lambda r: -(r[key] or 0.0))
        return [{"rank": i, "bone": r["bone"], "value": r[key]} for i, r in enumerate(order, 1)]
    out["ranking_original_by_volume_avg_over_sigc"] = _rank(
        out["original_mesh_all9"], "volume_avg_over_sigc")
    out["ranking_original_by_median_over_sigc"] = _rank(
        out["original_mesh_all9"], "all_median_over_sigc")
    out["ranking_original_by_scaled_p95_over_sigc"] = _rank(
        out["original_mesh_all9"], "all_p95_over_sigc")
    out["ranking_original_by_stored_p95_over_sigc"] = _rank(
        out["original_mesh_all9"], "stored_p95_over_sigc")

    # fibula clean remesh = FAILED (all loads diverged); record explicitly
    out["fibula_remesh_status"] = (
        "FAILED: all candidate loads (25266/10000/5000/2500/1000 N) diverged "
        "(negative jacobians, slender oblique bar buckling). No BC-robust sigma available."
    )
    return out


# --------------------------------------------------------------------------
# cross-check
# --------------------------------------------------------------------------
def _ranks(vals: list[float | None]) -> list[float]:
    """1-based ranks, None pushed to the end (ties -> average)."""
    idx = sorted(range(len(vals)), key=lambda i: (vals[i] is None, vals[i]))
    r = [0.0] * len(vals)
    pos = 0
    while pos < len(idx):
        j = pos
        if vals[idx[pos]] is None:
            for k in range(pos, len(idx)):
                r[idx[k]] = float(len(idx))
            break
        while j + 1 < len(idx) and vals[idx[j + 1]] == vals[idx[pos]]:
            j += 1
        avg = (pos + j) / 2.0 + 1.0
        for k in range(pos, j + 1):
            r[idx[k]] = avg
        pos = j + 1
    return r


def _spearman(a: list[float], b: list[float]) -> float:
    ra, rb = np.array(_ranks(a)), np.array(_ranks(b))
    ra, rb = ra - ra.mean(), rb - rb.mean()
    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


def _rank_tbl(bones: list[str], vals: list[float], desc: bool) -> list[dict]:
    """Return [{rank, bone, value}] ordered by value (desc=True -> largest first)."""
    order = sorted(range(len(bones)), key=lambda i: (-(vals[i]) if desc else vals[i]))
    return [{"rank": k, "bone": bones[i], "value": vals[i]} for k, i in enumerate(order, 1)]


def cross_check(r1: dict, r2: dict, matrix: dict) -> dict:
    heights = matrix["heights_m"]
    # matrix p95-based first-fracture height per bone (from the JSON itself)
    by_bone = {p["bone"]: p for p in matrix["parts"]}
    # Route 1 vs matrix (p95)
    bones = [r["bone"] for r in r1["rows"]]
    r1_h = [r["first_fracture_h"] for r in r1["rows"]]
    mx_h = [by_bone[b]["first_fracture_p95_m"] for b in bones]
    mx_med_h = [by_bone[b]["first_fracture_median_m"] for b in bones]
    # Route 2 original-mesh volume-avg normalized -> first fracture height
    r2orig = {r["bone"]: r for r in r2["original_mesh_all9"]}
    r2_vol_ratio = [r2orig[b]["volume_avg_over_sigc"] for b in bones]
    r2_med_ratio = [r2orig[b]["all_median_over_sigc"] for b in bones]
    r2_p95_ratio = [r2orig[b]["all_p95_over_sigc"] for b in bones]
    mx_util = [by_bone[b]["utilization_at_5m"] for b in bones]
    # Route 1 normalized risk@5m (proxy; first-fracture-height None -> inf)
    r1_risk = [r["risk_at_5m"] for r in r1["rows"]]

    return dict(
        note=("monotone agreement measured by Spearman rho over the 9 bones "
              "(1 = same order, -1 = reversed; None/never-fail ranked last)"),
        spearman_route1_vs_matrix_p95_risk5m=_spearman(r1_risk, mx_util),
        spearman_route1_vs_route2_volume_avg=_spearman(r1_risk, r2_vol_ratio),
        spearman_route1_vs_route2_median=_spearman(r1_risk, r2_med_ratio),
        spearman_route1_vs_route2_scaled_p95=_spearman(r1_risk, r2_p95_ratio),
        spearman_matrix_p95_vs_route2_scaled_p95=_spearman(mx_util, r2_p95_ratio),
        rank_route1=_rank_tbl(bones, r1_risk, True),
        rank_matrix_p95=_rank_tbl(bones, mx_util, True),
        rank_route2_volume_avg=_rank_tbl(bones, r2_vol_ratio, True),
        rank_route2_median=_rank_tbl(bones, r2_med_ratio, True),
        rank_route2_scaled_p95=_rank_tbl(bones, r2_p95_ratio, True),
        route1_first_fracture_h={b: h for b, h in zip(bones, r1_h)},
        matrix_first_fracture_p95_h={b: h for b, h in zip(bones, mx_h)},
        matrix_first_fracture_median_h={b: h for b, h in zip(bones, mx_med_h)},
        route2_original_volume_avg_over_nom={b: v for b, v in zip(bones, r2_vol_ratio)},
        route2_original_median_over_nom={b: v for b, v in zip(bones, r2_med_ratio)},
        route1_risk_at_5m={b: v for b, v in zip(bones, r1_risk)},
    )


# --------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(ROOT / "results" / "opensim_fe"))
    ap.add_argument("--excl-frac", type=float, default=0.30)
    args = ap.parse_args()

    matrix = json.loads(FM_JSON.read_text(encoding="utf-8"))
    loads = matrix["loads_n"]
    heights = [float(h) for h in matrix["heights_m"]]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    r1 = route1(loads, heights)
    a_by_bone = {r["bone"]: r["A_section_mm2"] for r in r1["rows"]}
    matrix_parts = {p["bone"]: p for p in matrix["parts"]}
    r2 = route2(loads, heights, a_by_bone, matrix_parts)
    xc = cross_check(r1, r2, matrix)

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    (out_dir / "bc_robust_metric_route1.json").write_text(
        json.dumps({"generated_at": ts, **r1}, indent=2, ensure_ascii=False), encoding="utf-8")
    (out_dir / "bc_robust_metric_route2.json").write_text(
        json.dumps({"generated_at": ts, **r2}, indent=2, ensure_ascii=False), encoding="utf-8")
    (out_dir / "bc_robust_metric_result.json").write_text(
        json.dumps({"generated_at": ts, "route1": r1, "route2": r2, "cross_check": xc},
                   indent=2, ensure_ascii=False), encoding="utf-8")

    # ---- console report ----
    print("=" * 100)
    print("Route 1 (FE-free 1D strength cross-section): risk(h)=F_joint(h)/(A*sigma_c)")
    print(f"{'bone':13s}{'dom':12s}{'V(mm3)':>11s}{'L':>7s}{'A(mm2)':>9s}"
          f"{'sc':>5s}{'F@5m(N)':>9s}{'thr(kN)':>8s}{'risk@5m':>9s}{'h1st':>6s}")
    for r in sorted(r1["rows"], key=lambda x: x["route1_rank"]):
        print(f"{r['bone']:13s}{r['solid_dom']:12s}{r['V_solid_mm3']:11.1f}{r['L_bbox_mm']:7.1f}"
              f"{r['A_section_mm2']:9.2f}{r['sigma_c_mpa']:5.0f}{r['F_at_5m_n']:9.0f}"
              f"{r['force_threshold_n']/1e3:8.2f}{r['risk_at_5m']:9.4f}"
              f"{str(r['first_fracture_h']):>6s}")
    print("\nRoute 2 (clean remesh, oblique main BC; mid = end bands 30% excluded):")
    print(f"{'bone':13s}{'stored?':>8s}{'all_p95':>10s}{'mid_p95':>10s}{'volavg':>10s}"
          f"{'all_med':>10s}{'mid_med':>10s}{'nom':>8s}")
    for r in r2["clean_remesh_oblique"]:
        print(f"{r['bone']:13s}{'-':>8s}{r['all_p95']:10.2f}{r['mid_p95']:10.2f}"
              f"{r['volume_avg']:10.2f}{r['all_median']:10.2f}{r['mid_median']:10.2f}"
              f"{r['sigma_nom_mpa']:8.2f}")
    print("  axial diagnostic (PCA-aligned load):")
    for r in r2["clean_remesh_axial"]:
        print(f"{r['bone']+'_axial':18s}{'':>3s}{r['all_p95']:10.2f}{r['mid_p95']:10.2f}"
              f"{r['volume_avg']:10.2f}{r['all_median']:10.2f}{r['mid_median']:10.2f}"
              f"{r['sigma_nom_mpa']:8.2f}")
    print("\nRoute 2b (original 9-bone re-analysis, max-across-candidate-domain;")
    print("         all_p95 = raw solved x (ref/solve); stored_p95 = the matrix JSON value):")
    print(f"{'bone':13s}{'scale':>7s}{'all_p95':>10s}{'mid_p95':>10s}{'volavg':>10s}"
          f"{'all_med':>10s}{'stored_p95':>11s}{'vol/sc':>8s}{'med/sc':>8s}")
    for r in r2["original_mesh_all9"]:
        sp = r["stored_p95"] if r["stored_p95"] is not None else float("nan")
        print(f"{r['bone']:13s}{r['scale']:7.3f}{r['all_p95']:10.2f}{r['mid_p95']:10.2f}"
              f"{r['volume_avg']:10.2f}{r['all_median']:10.2f}{sp:11.2f}"
              f"{r['volume_avg_over_sigc']:8.2f}{r['all_median_over_sigc']:8.2f}")

    print("\nRankings by stress/sigma_c usage (1 = highest risk):")
    print("  matrix p95 (stored):", " > ".join(x["bone"] for x in xc["rank_matrix_p95"])[:200])
    print("  Route 1 (1D nominal):", " > ".join(x["bone"] for x in xc["rank_route1"]))
    print("  Route 2 (vol-avg/sc):", " > ".join(x["bone"] for x in xc["rank_route2_volume_avg"]))
    print("  Route 2 (median/sc) :", " > ".join(x["bone"] for x in xc["rank_route2_median"]))
    print("  Route 2 (scaled p95):", " > ".join(x["bone"] for x in xc["rank_route2_scaled_p95"]))
    print(f"  Spearman Route1 vs matrix-p95 = {xc['spearman_route1_vs_matrix_p95_risk5m']:.3f}"
          f" | Route1 vs Route2-volavg = {xc['spearman_route1_vs_route2_volume_avg']:.3f}"
          f" | Route1 vs Route2-median = {xc['spearman_route1_vs_route2_median']:.3f}"
          f" | Route1 vs Route2-scaledp95 = {xc['spearman_route1_vs_route2_scaled_p95']:.3f}"
          f" | matrix vs Route2-scaledp95 = {xc['spearman_matrix_p95_vs_route2_scaled_p95']:.3f}")
    print("=" * 100)
    print(f"wrote: {out_dir / 'bc_robust_metric_result.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
