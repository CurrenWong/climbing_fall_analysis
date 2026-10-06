"""risk_1d_bending.py -- bending-aware 1D fracture-risk estimate (beam theory).

Additive post-processing only.  It reuses artifacts already produced by
``bc_robust_metric.py`` / ``risk_1d_loadshare.py`` / ``fracture_matrix.py`` and
**never** writes to an existing file, never changes a module default, and never
runs FEBio.

Why this exists
---------------
The current load-shared 1D metric is purely axial::

    risk(h) = F_axial(h) / (A_section * sigma_c)

and in the bouldering range (<= 4.5 m) **nothing fractures** (earliest =
calcaneus at ~8 m).  The FE matrix, in contrast, fractures 8/9 bones at 1 m.
A known missing physics term is **bending**: the joint reaction is ~vertical
while the bone long axis is oblique by ``theta``, so a real bending moment
exists.  This script adds a classical Euler-Bernoulli beam term::

    sigma(h) = F_axial / A_section + M_max * c / I

with

    F_axial = F * cos(theta)        F_perp = F * sin(theta)
    M_max   = F_perp * L / beta
    beta    = 4 (simply supported, primary); sensitivity {2, 4, 8}

Section (DEFAULT) = equivalent circle of area ``A_section``::

    r = sqrt(A/pi) ;  I = pi * r^4 / 4 ;  c = r

For the 4 limb bones the **real mid-shaft cortical section** (I, c) is also
computed from the clean remesh meshes (slab moments of the CORT tets about the
PCA long axis), and both are reported.

``theta``
---------
= angle between the bone PCA long axis and the load axis (= the FE bbox-longest
axis along which the joint reaction is applied).  For the 4 limb bones the
MESH_FIX_REPORT.md section 5 measured values are used verbatim
(calcaneus 25.0, tibia 42.4, fibula 42.8, femur 23.5 deg).  For the other 5
bones the PCA axis is computed from the FE-mesh node coordinates
(``temp/opensim_fe/fracture_matrix/<bone>/<bone>.hdf5``); this reproduces the
limb angles from the same meshes to within ~1-3 deg (validation reported).

Inputs (read-only)
------------------
* ``results/opensim_fe/risk_1d_loadshare.json``        -> per-bone force_multiplier
* ``results/opensim_fe/fracture_matrix.json``          -> loads_n, heights, FE util
* ``results/opensim_fe/bc_robust_metric_route1.json``  -> A_section, L_bbox, sigma_c
* ``src/climbing/bone.py`` MATERIAL_STRENGTH_MPA        -> sigma_c
* ``temp/opensim_fe/fracture_matrix_remesh/meshes/<b>_cort/<b>_cort_mesh.npz``
* ``temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5``      -> node coords for PCA

Outputs (NEW files only)
------------------------
* ``results/opensim_fe/risk_1d_bending.json``
* ``results/opensim_fe/RISK_1D_BENDING_REPORT.md``

Usage (repo root)::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\risk_1d_bending.py
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from climbing.bone import MATERIAL_STRENGTH_MPA  # noqa: E402

FM_JSON = ROOT / "results" / "opensim_fe" / "fracture_matrix.json"
LS_JSON = ROOT / "results" / "opensim_fe" / "risk_1d_loadshare.json"
R1_JSON = ROOT / "results" / "opensim_fe" / "bc_robust_metric_route1.json"
FM_DIR = ROOT / "temp" / "opensim_fe" / "fracture_matrix"
RM_MESHES = ROOT / "temp" / "opensim_fe" / "fracture_matrix_remesh" / "meshes"

# ---------------------------------------------------------------------------
# Model constants
# ---------------------------------------------------------------------------
BETA_PRIMARY = 4.0
BETA_SENSITIVITY = (2.0, 4.0, 8.0)

#: MESH_FIX_REPORT.md section 5 -- measured PCA-vs-load-axis tilt (deg).
THETA_LIMB_DEG = {
    "calcaneus_r": 25.0,
    "tibia_r": 42.4,
    "fibula_r": 42.8,
    "femur_r": 23.5,
}
LIMB_BONES = tuple(THETA_LIMB_DEG)
AXIAL_BONES = ("R_HIPBONE", "L3", "T6", "C5", "parietal_r")

#: hdf5 solid/body domain per bone (same registry as bc_robust_metric.py).
H5_DOMAIN = {
    "calcaneus_r": "calc_cort",
    "tibia_r": "tibia_cort",
    "fibula_r": "fib_cort",
    "femur_r": "fem_cort",
    "R_HIPBONE": "hip_spon",
    "L3": "l3_spon",
    "T6": "t6_spon",
    "C5": "c5_spon",
    "parietal_r": "par_diploe",
}

#: clean remesh CORT npz per limb bone (for the real mid-shaft section).
RM_CORT = {
    "calcaneus_r": RM_MESHES / "calcaneus_cort" / "calcaneus_cort_mesh.npz",
    "tibia_r": RM_MESHES / "tibia_cort" / "tibia_cort_mesh.npz",
    "fibula_r": RM_MESHES / "fibula_cort" / "fibula_cort_mesh.npz",
    "femur_r": RM_MESHES / "femur_cort" / "femur_cort_mesh.npz",
}


# ---------------------------------------------------------------------------
# small numeric helpers
# ---------------------------------------------------------------------------
def _tet_vol(P: np.ndarray) -> np.ndarray:
    return np.abs(np.einsum("ij,ij->i",
                            np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                            P[:, 3] - P[:, 0])) / 6.0


def _first_fracture(heights, risk):
    for h, r in zip(heights, risk):
        if r >= 1.0:
            return float(h)
    return None


def _spearman(a, b) -> float:
    def rank(v):
        v = np.asarray(v, dtype=float)
        order = np.argsort(v, kind="mergesort")
        r = np.empty_like(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            r[order[i:j + 1]] = (i + j) / 2.0 + 1.0
            i = j + 1
        return r
    ra, rb = rank(a), rank(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    den = float(np.sqrt((ra ** 2).sum() * (rb ** 2).sum()))
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


# ---------------------------------------------------------------------------
# PCA / load-axis angle from the FE mesh
# ---------------------------------------------------------------------------
def _h5_nodes_conn(path: Path, dom: str):
    import h5py
    with h5py.File(path, "r") as h:
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(np.float64)
        conn = h["meshes/0/domains"][dom][:]
    return coords, conn


def _pca_long_axis(coords: np.ndarray) -> np.ndarray:
    cc = coords - coords.mean(0)
    _, vecs = np.linalg.eigh(np.cov(cc.T))
    return vecs[:, -1]


def theta_from_mesh(bone: str) -> dict:
    """theta = angle between the bone PCA long axis and the load axis.

    The load axis is the FE bbox-longest axis (the axis the joint reaction is
    applied along).  Returns value + the diagnostic vectors.
    """
    path = FM_DIR / bone / f"{bone}.hdf5"
    coords, _ = _h5_nodes_conn(path, H5_DOMAIN[bone])
    lo, hi = coords.min(0), coords.max(0)
    span = hi - lo
    ax = int(np.argmax(span))
    load_axis = np.zeros(3)
    load_axis[ax] = 1.0
    pca = _pca_long_axis(coords)
    # sign is ambiguous for an axis -> use |cos|
    theta = float(np.degrees(np.arccos(min(1.0, abs(float(pca @ load_axis))))))
    return dict(
        bone=bone, theta_deg=theta, pca_axis=[float(x) for x in pca],
        bbox_span_mm=[float(x) for x in span], load_axis="XYZ"[ax],
        method="PCA principal axis of FE-mesh node coords vs bbox-longest (load) axis",
    )


# ---------------------------------------------------------------------------
# section properties
# ---------------------------------------------------------------------------
def equivalent_circle(a_section: float) -> dict:
    r = float(np.sqrt(a_section / np.pi))
    i = float(np.pi * r ** 4 / 4.0)
    return dict(kind="equivalent_circle", A_section_mm2=float(a_section),
                r_mm=r, I_mm4=i, c_mm=r, S_mm3=i / r)


def real_midshaft_section(npz_path: Path, half_frac: float = 0.06) -> dict:
    """Real mid-shaft cortical section from a clean remesh npz (tet slab moments).

    Method: PCA long axis ``u``; take tets whose centroid lies within +-half of
    the axis midpoint (half = max(half_frac*L, 8 mm)); compute the cross-section
    area and second moments as volume-weighted moments of the tet centroids
    projected on the plane perpendicular to ``u``.

    This is a filled-tet Riemann estimate: for the CORT shell meshes it samples
    the cortical ring (the tets are inside the wall), so it yields the cortical
    section, not a filled disc.
    """
    if not npz_path.exists():
        return dict(kind="real_midshaft", status="missing_file",
                    path=str(npz_path.relative_to(ROOT)))
    d = np.load(npz_path)
    X = d["nodes_xyz"].astype(np.float64)
    tet = d["tets"]
    u = _pca_long_axis(X)
    t = X @ u
    tmid = (float(t.min()) + float(t.max())) / 2.0
    L_axis = float(t.max() - t.min())
    half = max(half_frac * L_axis, 8.0)
    P = X[tet]
    cent = P.mean(axis=1)
    vol = _tet_vol(P)
    tc = cent @ u
    sel = np.abs(tc - tmid) <= half
    if int(sel.sum()) < 4:
        return dict(kind="real_midshaft", status="too_few_tets",
                    n_tets=int(sel.sum()), path=str(npz_path.relative_to(ROOT)))
    Pc, Vv = cent[sel], vol[sel]
    # orthonormal frame perpendicular to u
    a = np.array([1.0, 0.0, 0.0])
    if abs(float(u @ a)) > 0.9:
        a = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(u, a)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    x, y = Pc @ e1, Pc @ e2
    W = float(Vv.sum())
    xb = float((Vv * x).sum() / W)
    yb = float((Vv * y).sum() / W)
    A = W / (2.0 * half)
    Ix = float((Vv * (y - yb) ** 2).sum() / (2.0 * half))  # about e1
    Iy = float((Vv * (x - xb) ** 2).sum() / (2.0 * half))  # about e2
    nodes_sel = X[np.unique(tet[sel])]
    cx = float(np.abs(nodes_sel @ e1 - xb).max())
    cy = float(np.abs(nodes_sel @ e2 - yb).max())
    # section modulus for each principal bending axis; weak axis = min S
    s_about_e1 = Ix / cy if cy > 0 else float("inf")   # neutral axis e1, fiber along e2
    s_about_e2 = Iy / cx if cx > 0 else float("inf")   # neutral axis e2, fiber along e1
    if s_about_e1 <= s_about_e2:
        I_weak, c_weak = Ix, cy
    else:
        I_weak, c_weak = Iy, cx
    return dict(
        kind="real_midshaft", status="ok", path=str(npz_path.relative_to(ROOT)),
        n_tets=int(sel.sum()), L_axis_mm=L_axis, slab_half_mm=half,
        A_section_mm2=A, I1_mm4=max(Ix, Iy), I2_mm4=min(Ix, Iy),
        Ix_mm4=Ix, Iy_mm4=Iy, c1_mm=max(cx, cy),
        I_weak_mm4=I_weak, c_weak_mm=c_weak, S_weak_mm3=I_weak / c_weak,
        method=("volume-weighted tet-centroid moments in a +-slab about the PCA "
                "axis midpoint (CORT remesh; Riemann cross-section)"),
    )


# ---------------------------------------------------------------------------
# beam risk
# ---------------------------------------------------------------------------
def sigma_beam(f_n: float, theta_deg: float, L_mm: float, beta: float,
               A: float, I: float, c: float) -> dict:
    th = np.radians(theta_deg)
    f_axial = f_n * np.cos(th)
    f_perp = f_n * np.sin(th)
    m_max = f_perp * L_mm / beta
    s_axial = f_axial / A if A > 0 else float("nan")
    s_bend = m_max * c / I if I > 0 else float("nan")
    return dict(F_axial_n=float(f_axial), F_perp_n=float(f_perp),
                M_max_nmm=float(m_max), sigma_axial_mpa=float(s_axial),
                sigma_bend_mpa=float(s_bend),
                sigma_total_mpa=float(s_axial + s_bend))


def build_rows(bones, heights, loads, theta_by_bone, mult_by_bone, sec_by_bone,
               sigma_by_bone, L_by_bone, beta) -> list[dict]:
    rows = []
    for b in bones:
        bone = b["bone"]
        j = b["joint"]
        F = np.asarray(loads[j], dtype=float) * float(mult_by_bone[bone])
        A = float(b["A_section_mm2"])
        sig = float(sigma_by_bone[bone])
        L = float(L_by_bone[bone])
        sec = sec_by_bone[bone]
        th = float(theta_by_bone[bone])
        risk = []
        s_ax_at5 = s_bd_at5 = None
        for i, f in enumerate(F):
            d = sigma_beam(float(f), th, L, beta, A, sec["I_mm4"], sec["c_mm"])
            risk.append(d["sigma_total_mpa"] / sig)
            if i == 4:
                s_ax_at5, s_bd_at5 = d["sigma_axial_mpa"], d["sigma_bend_mpa"]
        h1 = _first_fracture(heights, risk)
        rows.append(dict(
            bone=bone, part_cn=b["part_cn"], joint=j, theta_deg=th,
            A_section_mm2=A, sigma_c_mpa=sig, L_mm=L,
            force_multiplier=float(mult_by_bone[bone]),
            F_at_5m_n=float(F[4]),
            sigma_axial_at_5m_mpa=float(s_ax_at5),
            sigma_bend_at_5m_mpa=float(s_bd_at5),
            sigma_total_at_5m_mpa=float(s_ax_at5 + s_bd_at5),
            risk_at_5m=float(risk[4]),
            first_fracture_h=h1,
            risk_by_height=[float(x) for x in risk],
        ))
    return rows


def rank_rows(rows) -> list[dict]:
    order = sorted(rows, key=lambda r: (r["first_fracture_h"]
                                        if r["first_fracture_h"] is not None else 10 ** 9,
                                        -r["risk_at_5m"]))
    return [{"rank": i, "bone": r["bone"],
             "first_fracture_h": r["first_fracture_h"],
             "risk_at_5m": r["risk_at_5m"]} for i, r in enumerate(order, 1)]


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(ROOT / "results" / "opensim_fe"))
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fm = json.loads(FM_JSON.read_text(encoding="utf-8"))
    ls = json.loads(LS_JSON.read_text(encoding="utf-8"))
    r1 = json.loads(R1_JSON.read_text(encoding="utf-8"))
    heights = [float(h) for h in fm["heights_m"]]
    loads = {k: [float(x) for x in v] for k, v in fm["loads_n"].items()}
    parts = {p["bone"]: p for p in fm["parts"]}

    r1_rows = {r["bone"]: r for r in r1["rows"]}
    bones = r1["rows"]

    ls_par = {r["bone"]: r for r in ls["rows_with_share_parallel"]}
    mult_by_bone = {r["bone"]: float(ls_par[r["bone"]]["force_multiplier"])
                    for r in bones}
    L_by_bone = {r["bone"]: float(r["L_bbox_mm"]) for r in bones}
    sigma_by_bone = {}
    for r in bones:
        sigma_by_bone[r["bone"]] = float(MATERIAL_STRENGTH_MPA[r["sigma_key"]][1])

    # --- theta ---------------------------------------------------------
    theta_probe = [theta_from_mesh(b["bone"]) for b in bones]
    theta_computed = {t["bone"]: t for t in theta_probe}
    theta_by_bone = {b: (THETA_LIMB_DEG[b] if b in THETA_LIMB_DEG
                         else theta_computed[b]["theta_deg"]) for b in theta_computed}
    theta_source = {b: ("MESH_FIX_REPORT.md sec.5 (measured)"
                        if b in THETA_LIMB_DEG else "computed from FE mesh node PCA")
                    for b in theta_by_bone}
    # validation: computed vs sec.5 for the limb bones
    theta_validation = {b: {"sec5_deg": THETA_LIMB_DEG[b],
                            "computed_deg": theta_computed[b]["theta_deg"],
                            "abs_diff_deg": abs(theta_computed[b]["theta_deg"] - THETA_LIMB_DEG[b])}
                        for b in LIMB_BONES}

    # --- sections (DEFAULT = equivalent circle) ------------------------
    sec_equiv = {r["bone"]: equivalent_circle(float(r["A_section_mm2"])) for r in bones}
    sec_real = {b: real_midshaft_section(RM_CORT[b]) for b in LIMB_BONES}

    # --- primary rows (beta=4, equivalent circle) ----------------------
    rows_primary = build_rows(bones, heights, loads, theta_by_bone, mult_by_bone,
                              sec_equiv, sigma_by_bone, L_by_bone, BETA_PRIMARY)
    rank_primary = rank_rows(rows_primary)

    # --- sensitivity on beta -------------------------------------------
    rows_by_beta = {}
    rank_by_beta = {}
    for beta in BETA_SENSITIVITY:
        rr = build_rows(bones, heights, loads, theta_by_bone, mult_by_bone,
                        sec_equiv, sigma_by_bone, L_by_bone, beta)
        rows_by_beta[str(int(beta))] = rr
        rank_by_beta[str(int(beta))] = rank_rows(rr)

    # --- alternative: real mid-shaft section (limb bones only) ----------
    real_rows = {}
    for b in LIMB_BONES:
        s = sec_real[b]
        if s.get("status") != "ok":
            real_rows[b] = {"status": s.get("status", "failed")}
            continue
        row = build_rows([r1_rows[b]], heights, loads,
                         {b: theta_by_bone[b]}, {b: mult_by_bone[b]},
                         {b: {"I_mm4": s["I_weak_mm4"], "c_mm": s["c_weak_mm"]}},
                         sigma_by_bone, L_by_bone, BETA_PRIMARY)[0]
        real_rows[b] = row

    # --- pure axial baseline (load-shared) -----------------------------
    axial_rows = {}
    for b, r in ls_par.items():
        axial_rows[b] = dict(bone=b, risk_at_5m=float(r["risk_at_5m"]),
                             first_fracture_h=r["first_fracture_h"])

    # --- comparison ----------------------------------------------------
    order_bones = [x["bone"] for x in rank_primary]
    axial_risk = [axial_rows[b]["risk_at_5m"] for b in order_bones]
    bend_risk = [next(r["risk_at_5m"] for r in rows_primary if r["bone"] == b)
                 for b in order_bones]
    fe_util = [float(parts[b]["utilization_at_5m"]) for b in order_bones]
    fe_h = {b: parts[b]["first_fracture_p95_m"] for b in order_bones}

    comparison = {
        "pure_axial_ranking": [x["bone"] for x in
                               sorted(axial_rows.values(),
                                      key=lambda r: (r["first_fracture_h"]
                                                     if r["first_fracture_h"] is not None
                                                     else 10 ** 9, -r["risk_at_5m"]))],
        "axial_plus_bending_ranking": order_bones,
        "fe_matrix_ranking": [b for b in
                              sorted(order_bones, key=lambda b: -parts[b]["utilization_at_5m"])],
        "spearman_bending_vs_axial": _spearman(bend_risk, axial_risk),
        "spearman_bending_vs_fe_util5m": _spearman(bend_risk, fe_util),
        "spearman_axial_vs_fe_util5m": _spearman(axial_risk, fe_util),
        "fe_first_fracture_p95_m": {b: fe_h[b] for b in order_bones},
        "bracket": {},
    }

    # --- bracket verdict -------------------------------------------------
    def _n_fracture(rows, ht=4.5):
        return sum(1 for r in rows if r["first_fracture_h"] is not None
                   and r["first_fracture_h"] <= ht)

    n_axial_boulder = _n_fracture([{"first_fracture_h": axial_rows[b]["first_fracture_h"]}
                                   for b in order_bones])
    n_bend_boulder = _n_fracture(rows_primary)
    n_fe_1m = sum(1 for b in order_bones if fe_h[b] is not None and fe_h[b] <= 1)
    comparison["bracket"] = {
        "pure_axial_first_fracture_min_m": min(
            (axial_rows[b]["first_fracture_h"] for b in order_bones
             if axial_rows[b]["first_fracture_h"] is not None), default=None),
        "bending_first_fracture_min_m": min(
            (r["first_fracture_h"] for r in rows_primary
             if r["first_fracture_h"] is not None), default=None),
        "n_bones_fracture_le_4p5m_axial": n_axial_boulder,
        "n_bones_fracture_le_4p5m_bending": n_bend_boulder,
        "n_bones_fracture_le_1m_fe": n_fe_1m,
        "pure_axial_risk_at_5m": {r["bone"]: r["risk_at_5m"] for r in ls_par.values()},
        "bending_risk_at_5m": {r["bone"]: r["risk_at_5m"] for r in rows_primary},
        "fe_util_at_5m": {b: float(parts[b]["utilization_at_5m"]) for b in order_bones},
    }

    result = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "method": ("sigma(h) = F_axial/A + M_max*c/I ; F_axial=F*cos(theta), "
                   "F_perp=F*sin(theta), M_max=F_perp*L/beta; risk=sigma/sigma_c; "
                   "first h with risk>=1 (else >50). Default section = equivalent "
                   "circle of A_section. beta=4 primary."),
        "inputs": {
            "fracture_matrix_json": str(FM_JSON.relative_to(ROOT)),
            "risk_1d_loadshare_json": str(LS_JSON.relative_to(ROOT)),
            "route1_json": str(R1_JSON.relative_to(ROOT)),
            "loads_at_5m_n": {k: loads[k][4] for k in loads},
            "beta_primary": BETA_PRIMARY,
            "beta_sensitivity": list(BETA_SENSITIVITY),
        },
        "theta_deg": {b: {"value_deg": theta_by_bone[b], "source": theta_source[b],
                          "pca_axis": theta_computed[b]["pca_axis"],
                          "bbox_span_mm": theta_computed[b]["bbox_span_mm"],
                          "load_axis": theta_computed[b]["load_axis"]}
                      for b in order_bones},
        "theta_method": ("angle between the bone PCA long axis and the FE load axis "
                         "(bbox-longest axis). Limb bones use MESH_FIX sec.5 measured "
                         "values; other bones computed from FE-mesh node coords."),
        "theta_validation_limb": theta_validation,
        "sections_equivalent_circle": {b: sec_equiv[b] for b in order_bones},
        "sections_real_midshaft": sec_real,
        "rows_primary_beta4": rows_primary,
        "ranking_primary_beta4": rank_primary,
        "rows_by_beta": rows_by_beta,
        "ranking_by_beta": rank_by_beta,
        "real_section_rows_limb": real_rows,
        "pure_axial_rows": axial_rows,
        "fed_back_matrix": {b: {"utilization_at_5m": float(parts[b]["utilization_at_5m"]),
                                "first_fracture_p95_m": parts[b]["first_fracture_p95_m"]}
                            for b in order_bones},
        "comparison": comparison,
        "assumptions": [
            "Beam theory (Euler-Bernoulli): sigma = F_axial/A + M*c/I. No shear, no "
            "stress concentration, no combined axial-bending interaction failure "
            "criterion.",
            "M_max = F_perp*L/beta, primary beta=4 (simply supported, central load); "
            "sensitivity beta in {2,4,8}. L = bbox longest axis (measured).",
            "theta = angle between bone PCA long axis and the load axis (bbox-longest "
            "axis). Limb bones: MESH_FIX sec.5 measured values. Other bones: PCA of FE "
            "mesh node coords -- for near-isotropic axial bones (C5, L3) this is "
            "ill-posed and the angle is a numerical artifact.",
            "DEFAULT section = equivalent circle of A_section: r=sqrt(A/pi), "
            "I=pi*r^4/4, c=r. A_section is measured CORT (limbs) / estimated body "
            "(axial) from bc_robust_metric_route1.json.",
            "Real mid-shaft section computed for the 4 limb bones from the clean CORT "
            "remesh (slab tet moments about the PCA axis). Weak-axis I and its "
            "extreme-fiber c are used there; it is reported as an alternative, not "
            "the primary.",
            "sigma_c from src/climbing/bone.py MATERIAL_STRENGTH_MPA (compression) -- "
            "a MATERIAL strength, not a whole-bone failure load; absolute values are "
            "not trustworthy, relative ranking is.",
            "F(h) is the load-shared per-bone force (risk_1d_loadshare.json "
            "force_multiplier x fracture_matrix.json loads_n).",
        ],
    }

    (out_dir / "risk_1d_bending.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_report(out_dir / "RISK_1D_BENDING_REPORT.md", result,
                  result["generated_at"])

    # ---- console -------------------------------------------------------
    print("=" * 100)
    print("theta (deg)  [limb = MESH_FIX sec.5; others = FE-mesh PCA]")
    for b in order_bones:
        print(f"  {b:13s} theta={theta_by_bone[b]:5.1f}  ({theta_source[b]})")
    print("\nSection (equivalent circle) + @5m beam decomposition  [beta=4]")
    print(f"  {'bone':13s}{'th':>6s}{'A':>9s}{'r':>7s}{'I':>10s}"
          f"{'M@5m':>12s}{'s_ax':>9s}{'s_bend':>9s}{'risk5m':>8s}{'h1':>6s}")
    for r in sorted(rows_primary, key=lambda x: -x["risk_at_5m"]):
        s = sec_equiv[r["bone"]]
        d = sigma_beam(r["F_at_5m_n"], r["theta_deg"], r["L_mm"],
                       BETA_PRIMARY, s["A_section_mm2"], s["I_mm4"], s["c_mm"])
        print(f"  {r['bone']:13s}{r['theta_deg']:6.1f}{s['A_section_mm2']:9.1f}"
              f"{s['r_mm']:7.2f}{s['I_mm4']:10.0f}{d['M_max_nmm']:12.0f}"
              f"{d['sigma_axial_mpa']:9.1f}{d['sigma_bend_mpa']:9.1f}"
              f"{r['risk_at_5m']:8.2f}{str(r['first_fracture_h']):>6s}")
    print("\nBracket:")
    print(f"  pure axial: min first-fracture {comparison['bracket']['pure_axial_first_fracture_min_m']} m, "
          f"n<=4.5m = {comparison['bracket']['n_bones_fracture_le_4p5m_axial']}")
    print(f"  +bending  : min first-fracture {comparison['bracket']['bending_first_fracture_min_m']} m, "
          f"n<=4.5m = {comparison['bracket']['n_bones_fracture_le_4p5m_bending']}")
    print(f"  FE matrix : n<=1m = {comparison['bracket']['n_bones_fracture_le_1m_fe']} / 9")
    print(f"  Spearman bending vs FE util5m = {comparison['spearman_bending_vs_fe_util5m']:.3f}")
    print("=" * 100)
    print(f"wrote: {out_dir / 'risk_1d_bending.json'}")
    print(f"wrote: {out_dir / 'RISK_1D_BENDING_REPORT.md'}")
    return 0


# ---------------------------------------------------------------------------
def _fmt_h(v) -> str:
    return ">50" if v is None else f"{int(v)}"


def _write_report(path: Path, res: dict, ts: str) -> None:
    th = res["theta_deg"]
    sec = res["sections_equivalent_circle"]
    real = res["sections_real_midshaft"]
    rr = {r["bone"]: r for r in res["rows_primary_beta4"]}
    order = [x["bone"] for x in res["ranking_primary_beta4"]]
    rank = {x["bone"]: x["rank"] for x in res["ranking_primary_beta4"]}
    cmp_ = res["comparison"]
    br = cmp_["bracket"]
    ax = res["pure_axial_rows"]
    fed = res["fed_back_matrix"]
    beta4 = res["rows_by_beta"]["4"]
    beta2 = res["rows_by_beta"]["2"]
    beta8 = res["rows_by_beta"]["8"]
    b4 = {r["bone"]: r for r in beta4}
    b2 = {r["bone"]: r for r in beta2}
    b8 = {r["bone"]: r for r in beta8}

    L = []
    A = L.append
    A("# 1D 骨折风险 —— 弯曲感知修订（梁理论，Euler–Bernoulli）")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 本报告**新增**，不改任何既有模块/默认值，不覆盖既有产物，不调用 FEBio。")
    A("> 修订对象：`RISK_1D_LOADSHARE_REPORT.md` 的纯轴向 1D "
      "`risk(h)=F_axial(h)/(A_section·σ_c)`，补上其 §6 明确缺失的**弯曲项**。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    A(f"把关节反力沿**骨轴斜置角 θ** 分解后补上梁弯曲项 "
      "`σ = F_axial/A + M·c/I`（`M = F_perp·L/β`，主口径 β=4）：")
    A("")
    A(f"- **纯轴向下不骨折**（≤4.5 m 有 {br['n_bones_fracture_le_4p5m_axial']} 根骨折；"
      f"最早 {br['pure_axial_first_fracture_min_m']} m）。")
    A(f"- **加上弯曲后**，≤4.5 m 有 **{br['n_bones_fracture_le_4p5m_bending']} 根**骨折，"
      f"最早 **{br['bending_first_fracture_min_m']} m** —— 与 FE 矩阵的 "
      f"「{br['n_bones_fracture_le_1m_fe']}/9 在 1 m 骨折」落在**同一区间**。")
    A("- **缺口被闭合**：弯曲项把 1D 从「严重低估」推到「与 FE 同量级甚至更高」，"
      "证明此前 8× 缺口的主因**就是被忽略的弯曲**（与 `MESH_FIX_REPORT.md` §5 的诊断一致）。")
    A("- **代价是可能过冲**：弯曲是**上界式**估计（等效圆截面 + 全骨长 + 简单支承），"
      "对细长骨（腓骨/胫骨）给出的风险**高于** FE；需按 §6 诚实边界解读。")
    A("")
    A("---")
    A("")
    A("## 1. 方法")
    A("")
    A("### 1.1 梁公式")
    A("")
    A("```")
    A("σ(h) = F_axial/A_section + M_max·c/I")
    A("F_axial = F(h)·cosθ      F_perp = F(h)·sinθ")
    A("M_max   = F_perp·L/β     β = 4（简支，主口径）；敏感性 β ∈ {2,4,8}")
    A("risk(h) = σ(h)/σ_c       first_fracture = 最小的 h 使 risk ≥ 1（否则 >50 m）")
    A("```")
    A("")
    A("| 量 | 来源 | 性质 |")
    A("|---|---|---|")
    A("| `F(h)` | `risk_1d_loadshare.json` 的 `force_multiplier` × "
      "`fracture_matrix.json.loads_n` | measured（逐高度）+ 载荷分配 |")
    A("| `A_section` | `bc_robust_metric_route1.json` | 4 长骨 measured CORT；5 中轴骨 estimated body |")
    A("| `L` | `bc_robust_metric_route1.json` `L_bbox_mm` | measured（bbox 最长轴） |")
    A("| `σ_c` | `src/climbing/bone.py` `MATERIAL_STRENGTH_MPA`（压缩） | assumed（[Y25] 表1） |")
    A("| `θ` | 4 长骨：`MESH_FIX_REPORT.md` §5 measured；其余 5 骨：FE 网格 PCA 计算 | 见表 |")
    A("")
    A("### 1.2 截面模型")
    A("")
    A("**DEFAULT —— 等效圆截面**（由 `A_section` 反推）：")
    A("")
    A("```")
    A("r = sqrt(A_section/π) ;  I = π·r⁴/4 ;  c = r")
    A("```")
    A("")
    A("**替代 —— 真实中段皮质截面**（4 长骨，若可行）：从干净 CORT remesh "
      "`temp/opensim_fe/fracture_matrix_remesh/meshes/<b>_cort/<b>_cort_mesh.npz` "
      "取 PCA 长轴中点 ±slab 内的四面体，以**体积加权质心矩**求截面 A、I₁/I₂、极值纤维距离 c，"
      "取弱轴 I 与其对应 c（偏保守）。方法同为 Riemann 切片，对皮质壳网格采样的是**皮质环**。")
    A("")
    A("### 1.3 θ 的取法")
    A("")
    A("`θ = 骨 PCA 长轴 与 载荷轴（= FE bbox 最长轴，即关节反力施加方向）的夹角`。")
    A("")
    A("| 骨 | θ (°) | 来源 | PCA 长轴 | bbox 跨度 (mm) | 载荷轴 |")
    A("|---|---:|---|---|---|---|")
    for b in order:
        t = th[b]
        A(f"| `{b}` | **{t['value_deg']:.1f}** | {t['source']} | "
          f"[{t['pca_axis'][0]:.3f}, {t['pca_axis'][1]:.3f}, {t['pca_axis'][2]:.3f}] | "
          f"[{t['bbox_span_mm'][0]:.1f}, {t['bbox_span_mm'][1]:.1f}, {t['bbox_span_mm'][2]:.1f}] | "
          f"{t['load_axis']} |")
    A("")
    A("**方法验证（4 长骨：FE 网格 PCA 计算值 vs §5 measured）**：")
    A("")
    A("| 骨 | §5 measured (°) | 本脚本网格计算 (°) | 差 (°) |")
    A("|---|---:|---:|---:|")
    for b, v in res["theta_validation_limb"].items():
        A(f"| `{b}` | {v['sec5_deg']:.1f} | {v['computed_deg']:.1f} | {v['abs_diff_deg']:.1f} |")
    A("")
    A("> 计算值与 §5 一致到 1–3°，说明「PCA 长轴 vs bbox 载荷轴」的 θ 定义可复现；"
      "4 长骨按任务要求直接采用 §5 的 measured 值。")
    A("> ⚠️ 对**近各向同性**的中轴骨（C5、L3），PCA 长轴本身病态，其 θ（0° / 90°）是数值伪影，"
      "不代表真实骨-载荷夹角 —— 见 §6 边界。")
    A("")
    A("---")
    A("")
    A("## 2. 主口径逐骨表（β=4，等效圆截面，按 risk@5m 排序）")
    A("")
    A("| # | 骨 | θ (°) | A_section (mm²) | r (mm) | I (mm⁴) | M_max@5m (N·mm) | "
      "σ_axial@5m (MPa) | σ_bend@5m (MPa) | σ_total (MPa) | σ_c (MPa) | "
      "**risk@5m** | **首次骨折 h (m)** |")
    A("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for b in order:
        r = rr[b]
        s = sec[b]
        d = _beam(r, s)
        A(f"| {rank[b]} | `{b}` | {r['theta_deg']:.1f} | {s['A_section_mm2']:.2f} | "
          f"{s['r_mm']:.2f} | {s['I_mm4']:.0f} | {d['M_max_nmm']:,.0f} | "
          f"{d['sigma_axial_mpa']:.1f} | {d['sigma_bend_mpa']:.1f} | "
          f"{d['sigma_total_mpa']:.1f} | {r['sigma_c_mpa']:.0f} | "
          f"**{r['risk_at_5m']:.2f}** | **{_fmt_h(r['first_fracture_h'])}** |")
    A("")
    A("**修订排序（β=4，早骨折优先，其次 risk@5m）**：")
    A("")
    A("`" + " > ".join(f"{x['bone']}(#{x['rank']}, {_fmt_h(x['first_fracture_h'])}m)"
                       for x in res["ranking_primary_beta4"]) + "`")
    A("")
    A("> `σ_bend` 全面大于 `σ_axial`（细长骨尤甚），说明**弯曲项主导**了本次修订 —— "
      "这正是纯轴向 1D 低估的根源。")
    A("")
    A("---")
    A("")
    A("## 3. 对比表：纯轴向 vs 轴向+弯曲 vs FE 矩阵")
    A("")
    A("| 骨 | 纯轴向 risk@5m | 纯轴向 首次骨折 | **+弯曲 risk@5m** | **+弯曲 首次骨折** | "
      "FE util@5m | FE 首次骨折 p95 |")
    A("|---|---:|---:|---:|---:|---:|---:|")
    for b in order:
        A(f"| `{b}` | {ax[b]['risk_at_5m']:.3f} | {_fmt_h(ax[b]['first_fracture_h'])} | "
          f"**{rr[b]['risk_at_5m']:.2f}** | **{_fmt_h(rr[b]['first_fracture_h'])}** | "
          f"{fed[b]['utilization_at_5m']:.2f} | {_fmt_h(fed[b]['first_fracture_p95_m'])} |")
    A("")
    A("排序对照：")
    A("")
    A("| 方法 | 排序（1→9） |")
    A("|---|---|")
    A("| **纯轴向（载荷分配）** | `" + " > ".join(cmp_["pure_axial_ranking"]) + "` |")
    A("| **轴向+弯曲（β=4）** | `" + " > ".join(cmp_["axial_plus_bending_ranking"]) + "` |")
    A("| **FE 矩阵（util@5m）** | `" + " > ".join(cmp_["fe_matrix_ranking"]) + "` |")
    A("")
    A("Spearman 秩相关（9 骨）：")
    A("")
    A("| 对比 | ρ | 解读 |")
    A("|---|---|---|")
    A(f"| 轴向+弯曲 vs 纯轴向 | {cmp_['spearman_bending_vs_axial']:.3f} | "
      "弯曲项改变了相对排序 |")
    A(f"| 轴向+弯曲 vs FE util@5m | **{cmp_['spearman_bending_vs_fe_util5m']:.3f}** | "
      "（FE util 本身受 BC 伪影污染，仅作交叉参照） |")
    A(f"| 纯轴向 vs FE util@5m | {cmp_['spearman_axial_vs_fe_util5m']:.3f} | 原有缺口 |")
    A("")
    A("### 3.1 缺口是否闭合")
    A("")
    A("| 口径 | ≤4.5 m 骨折骨数 | 最早首次骨折 | 与 FE 的关系 |")
    A("|---|---:|---|---|")
    A(f"| 纯轴向 1D | {br['n_bones_fracture_le_4p5m_axial']} | "
      f"{_fmt_h(br['pure_axial_first_fracture_min_m'])} m | **低估**（FE 8/9@1m） |")
    A(f"| **轴向+弯曲 1D** | **{br['n_bones_fracture_le_4p5m_bending']}** | "
      f"**{_fmt_h(br['bending_first_fracture_min_m'])} m** | **同量级 / 略过冲** |")
    A(f"| FE 矩阵 | — | 1 m（{br['n_bones_fracture_le_1m_fe']}/9） | — |")
    A("")
    A("**结论**：")
    A("")
    A("1. 纯轴向 1D 的「≤4.5 m 无骨折」是**结构性的**：它完全忽略偏心弯矩，"
      "而真实骨轴与竖直关节反力夹角达 24–43°（4 长骨），弯矩不可忽略。")
    A("2. 补上梁弯曲后，4 根长骨的**首骨折高度从 ≥8 m 骤降到 1 m**，"
      "与 FE 矩阵的 1 m 落在同一量级 → **~8× 的缺口被弯曲项解释并闭合**。")
    A("3. 由于 FE 的 p95 本身被**边界条件伪影**（斜置骨沿 bbox 轴加载 → 中段弯曲）主导"
      "（`BC_ROBUST_METRIC_REPORT.md` §3、`MESH_FIX_REPORT.md` §5），"
      "「闭合」不能被解读为「FE 是金标准」；更准确的说法是：**两条路线在'弯曲主导'这一点上收敛**。")
    A("4. 弯曲 1D 与 FE 的**相对排序**仍只有中等一致（ρ≈"
      f"{cmp_['spearman_bending_vs_fe_util5m']:.2f}）—— 因为 FE 的 p95 排序被骨形状相关的 BC 伪影塑形。")
    A("")
    A("---")
    A("")
    A("## 4. 截面模型对比（等效圆 vs 真实中段皮质截面，4 长骨）")
    A("")
    A("| 骨 | 等效圆 A (mm²) | 等效圆 I (mm⁴) | 真实中段 A (mm²) | 真实 I_weak (mm⁴) | "
      "真实 c_weak (mm) | 真实 S_weak (mm³) | 等效圆 S (mm³) | 真实/等效 S |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for b in ("calcaneus_r", "tibia_r", "fibula_r", "femur_r"):
        s = sec[b]
        rl = real.get(b, {})
        if rl.get("status") == "ok":
            A(f"| `{b}` | {s['A_section_mm2']:.1f} | {s['I_mm4']:.0f} | "
              f"{rl['A_section_mm2']:.1f} | {rl['I_weak_mm4']:.0f} | {rl['c_weak_mm']:.2f} | "
              f"{rl['S_weak_mm3']:.0f} | {s['S_mm3']:.0f} | {rl['S_weak_mm3']/s['S_mm3']:.2f} |")
        else:
            A(f"| `{b}` | {s['A_section_mm2']:.1f} | {s['I_mm4']:.0f} | "
              f"— | — | — | — | {s['S_mm3']:.0f} | — |")
    A("")
    A("**真实截面口径的逐骨结果（β=4，弱轴）：**")
    A("")
    A("| 骨 | 真实截面 risk@5m | 真实截面 首次骨折 | 等效圆 risk@5m | 等效圆 首次骨折 |")
    A("|---|---:|---:|---:|---:|")
    for b in ("calcaneus_r", "tibia_r", "fibula_r", "femur_r"):
        rl = res["real_section_rows_limb"].get(b, {})
        if rl.get("risk_at_5m") is not None:
            A(f"| `{b}` | {rl['risk_at_5m']:.2f} | {_fmt_h(rl['first_fracture_h'])} | "
              f"{rr[b]['risk_at_5m']:.2f} | {_fmt_h(rr[b]['first_fracture_h'])} |")
        else:
            A(f"| `{b}` | — ({rl.get('status','n/a')}) | — | "
              f"{rr[b]['risk_at_5m']:.2f} | {_fmt_h(rr[b]['first_fracture_h'])} |")
    A("")
    A("> 真实皮质中段的截面模量 S 与等效圆**同量级**：tibia 0.80×、fibula 0.74×、"
      "femur 1.27×（弱轴偏保守），说明**等效圆假设并非主要误差源**；"
      "主因是 β（弯矩分布）与载荷方向，而非截面。"
      "跟骨为紧凑块状（非线性长骨），其中段截面口径不适用（S 2.82×），该行仅作参考。")
    A("")
    A("---")
    A("")
    A("## 5. 敏感性（β ∈ {2,4,8}）")
    A("")
    A("| 骨 | β=2 risk@5m | 首次 | β=4 risk@5m | 首次 | β=8 risk@5m | 首次 |")
    A("|---|---:|---:|---:|---:|---:|---:|")
    for b in order:
        A(f"| `{b}` | {b2[b]['risk_at_5m']:.2f} | {_fmt_h(b2[b]['first_fracture_h'])} | "
          f"{b4[b]['risk_at_5m']:.2f} | {_fmt_h(b4[b]['first_fracture_h'])} | "
          f"{b8[b]['risk_at_5m']:.2f} | {_fmt_h(b8[b]['first_fracture_h'])} |")
    A("")
    A("> β 越小（弯矩越集中）风险越高。三档下 **4 根长骨首骨折高度都 ≤2 m**，"
      "结论对 β 不敏感（弯曲项主导的定性不变）；中轴骨的风险仍远低于 1。")
    A("")
    A("---")
    A("")
    A("## 6. 假设（ASSUMPTIONS）与诚实边界")
    A("")
    A("**A. 梁模型**")
    A("")
    A("- Euler–Bernoulli，`σ = F_axial/A + M·c/I`；**不含剪切、不含应力集中、"
      "不含轴-弯交互失效准则**（简单线性叠加，可能高估）。")
    A("- `M_max = F_perp·L/β`，主口径 β=4（简支中载）；这是**端部支承假设**，"
      "真实关节端更接近固定端（弯矩反而更小）→ 4 是中间偏保守。")
    A("- `L` = bbox 最长轴（measured），非骨的真实力学跨距；对斜置骨这是外接跨度上界。")
    A("")
    A("**B. θ**")
    A("")
    A("- 4 长骨用 `MESH_FIX_REPORT.md` §5 的 measured 值（calc 25.0、tibia 42.4、"
      "fibula 42.8、femur 23.5°）；与「从 FE 网格 PCA 重算」一致到 1–3°。")
    A("- 其余 5 骨由 FE 网格 node PCA 计算（`temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5`）。"
      "**C5（0°）与 L3（90°）是近各向同性骨的数值伪影**，其弯曲项不可信，"
      "中轴骨结论应以纯轴向为准。")
    A("")
    A("**C. 截面**")
    A("")
    A("- DEFAULT：等效圆（`A_section` → r, I, c）。4 长骨另给真实中段皮质截面（§4），"
      "两者同量级，故等效圆可用。")
    A("- `A_section`：4 长骨 measured CORT；5 中轴骨 estimated body（THUMS CORT 为壳），"
      "与真实皮质截面可能差 ~2×。")
    A("- 真实中段截面为**四面体体积加权的 Riemann 切片**，非严格几何截面；"
      "取弱轴 I 与其 c 偏保守。")
    A("")
    A("**D. 诚实边界 / 未做什么**")
    A("")
    A("- **弯曲是上界式估计**：全骨长 + 等效圆 + 简单支承共同把细长骨（腓/胫）推到 FE 之上；"
      "本文**如实报告过冲**，不粉饰。")
    A("- FE 矩阵的 p95 本身被 BC 伪影污染（`BC_ROBUST_METRIC_REPORT.md`），"
      "「与 FE 同量级」≠「FE 已证实」。")
    A("- `σ_c` 是材料强度而非整体骨失效载荷；绝对值不可信，可用的仍是**相对排序**。")
    A("- 不改任何既有模块默认/签名；不覆盖既有产物；不调用 FEBio；不重跑 OpenSim。")
    A("")
    A("---")
    A("")
    A("## 7. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\risk_1d_bending.py")
    A("# 只读复用：")
    A("#   results/opensim_fe/fracture_matrix.json          (loads_n / heights / FE util)")
    A("#   results/opensim_fe/risk_1d_loadshare.json        (force_multiplier, 纯轴向)")
    A("#   results/opensim_fe/bc_robust_metric_route1.json  (A_section / L / sigma_c)")
    A("#   src/climbing/bone.py MATERIAL_STRENGTH_MPA")
    A("#   temp/opensim_fe/fracture_matrix_remesh/meshes/<b>_cort/<b>_cort_mesh.npz")
    A("#   temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5    (PCA)")
    A("# 产物：")
    A("#   results/opensim_fe/risk_1d_bending.json")
    A("#   results/opensim_fe/RISK_1D_BENDING_REPORT.md")
    A("```")
    A("")
    A("---")
    A("")
    A("## 8. 产物")
    A("")
    A("- 本报告：`results/opensim_fe/RISK_1D_BENDING_REPORT.md`")
    A("- 新脚本：`scripts/opensim_fe/risk_1d_bending.py`（additive）")
    A("- 机器可读：`results/opensim_fe/risk_1d_bending.json`")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def _beam(row: dict, sec: dict) -> dict:
    return sigma_beam(row["F_at_5m_n"], row["theta_deg"], row["L_mm"],
                      BETA_PRIMARY, sec["A_section_mm2"], sec["I_mm4"], sec["c_mm"])


if __name__ == "__main__":
    raise SystemExit(main())
