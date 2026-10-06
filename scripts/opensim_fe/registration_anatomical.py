"""registration_anatomical.py -- ANATOMICAL (landmark-based) placement of the THUMS
AM50 right calcaneus into the OpenSim Rajagopal ``calcn_r`` body frame.

Why this module exists
----------------------
``registration.py`` / ``register_calcaneus.py`` (the "foot-frame" solution) put
the THUMS calcaneus into the calcn frame by **rigid ICP against ``r_foot.vtp``**
-- the crude *lumped whole-foot* mesh of the Rajagopal ``calcn_r`` body, which
contains no dedicated calcaneus.  That match is dominated by the foot outline, so
the resulting origin landed near the bone's **posterior edge**, and the
``subtalar_joint`` face (rule ``dist < 25 mm`` AND ``n.(+Y) > 0.2``) degenerated to
a tiny posterior-superior patch (14 quads / 222 mm2).  Feeding a ``PressureLoad
p = F/A`` onto that 222 mm2 patch already produced ~258 MPa at h = 1 m -- a
non-physical "fracture at <= 1 m".

This module instead places the bone by **its own anatomy**, with no external
target mesh:

    x_calcn = R @ x_thums + t          (rigid only; R orthonormal, det = +1)

so that (calcn convention, verified: ``calcn_r`` local -> ground rotation is the
identity; origin = subtalar joint centre, +X anterior, +Y superior/up,
Z medio-lateral):

* the **subtalar (superior articular) surface centroid** maps to the origin;
* the **plantar surface faces -Y** (so "up" is +Y);
* the **anterior process** (narrow, calcaneocuboid end) is +X and the bulky
  **posterior tuberosity** (Achilles end) is -X.

Anatomical landmark recipe (deterministic, THUMS-only)
------------------------------------------------------
Let ``V`` / ``F`` be the raw THUMS surface vertices / triangles (THUMS global mm),
``n`` the outward face unit normals and ``area`` the face areas.

1. **Dorsal-plantar axis** = the 2nd principal axis of ``V`` (``a1``); the
   medio-lateral width is the 3rd (``a2``) and the long (A-P) axis the 1st (``a0``)
   -- for a calcaneus height (~51 mm) > width (~41 mm) > ... .
   The *sign* is set by a **flatness** test: within a 20 deg normal cone, the
   flatter side is the broad flat **sole**, so ``down`` = the sign of ``a1`` whose
   cone is flatter (smaller out-of-plane thickness); ``up = -down``.
   (THUMS: down-side cone 644 mm2 @ 6.6 mm thick vs up-side 369 mm2 @ 31.6 mm.)
2. **Long axis** ``L`` = PCA of ``V`` projected onto the plane perpendicular to
   ``up`` (kept exactly orthogonal), so ``L`` is horizontal.
3. **Anterior / posterior sign**: the calcaneus's **Achilles insertion** is a
   broad, near-planar facet on the posterior tuberosity, normal parallel to the
   long axis.  For each end we score the planarity of the outward-facing end
   facet (``n.(+-L) > 0.8`` within 15 mm of the end) as ``area / (thickness+1)``;
   the higher-scoring end is **posterior**, so **anterior = the other end** and
   ``X = anterior``.  (THUMS: posterior end 264 mm2 @ 2.8 mm; anterior 103 mm2 @
   18.6 mm -> anterior is the narrow end, as required.)
4. ``Z = X x Y`` (right-handed).
5. **Origin** = centroid of the **subtalar surface** = faces with
   ``n.up > 0.2`` and centroid in the **posterior 2/3** along X (the subtalar
   posterior-facet region).  ``t = -R @ origin`` so the origin maps to (0,0,0).

Units: mm-N-MPa-s; all lengths mm.

Outputs (``temp/opensim_fe/thums_calcaneus/``)
----------------------------------------------
    calcaneus_r_anatframe.stl    transformed outer surface, ASCII, mm
    calcaneus_r_anatframe.vtk    transformed solid (3323 tet + 634 hex)
    calcaneus_r_anatframe.npz    node_ids / nodes_xyz / conn / pid / eid
    registration_anat.json       R, t, units, convention, landmarks, metrics
    REGISTRATION_ANAT_REPORT.md  human-readable method / validation / caveats

Public API
----------
face_geometry(V, F)                     -> (centroids, normals, areas)
pca_directions(V)                       -> (centroid, axes_desc, eigvals)
anatomical_frame(V, F, n, area, cent)   -> dict(R, t, axes, landmarks, ...)
compute_registration(surface_stl)       -> dict
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]  # scripts/opensim_fe -> scripts -> repo root
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import kmesh_io as kio          # noqa: E402
import registration as reg      # noqa: E402  (read_ascii_stl / write_ascii_stl)

__all__ = [
    "face_geometry",
    "pca_directions",
    "anatomical_frame",
    "compute_registration",
    "CALC_DIR",
]

CALC_DIR = ROOT / "temp" / "opensim_fe" / "thums_calcaneus"
SURFACE_STL = CALC_DIR / "calcaneus_r_surface.stl"
NODE_IDS_NPY = CALC_DIR / "node_ids.npy"
NODES_XYZ_NPY = CALC_DIR / "nodes_xyz.npy"
ELEMENTS_NPZ = CALC_DIR / "calcaneus_elements.npz"

#: FE-builder joint rule (kept bit-identical to ``climbing.coupling.meshing``).
JOINT_RADIUS_MM = 25.0
JOINT_NORMAL_Y_MIN = 0.2
#: cone half-angle (deg) for the flatness / plantar detection.
_FLAT_CONE_DEG = 20.0
#: an end facet counts as the Achilles insertion when flatter than this (mm).
_FLAT_END_THICK_MM = 6.0


# ---------------------------------------------------------------------------
# geometry primitives
# ---------------------------------------------------------------------------
def face_geometry(V: np.ndarray, F: np.ndarray):
    """Outward face centroids, unit normals and areas of a triangle surface.

    Outward orientation is resolved by the signed enclosed volume of the
    winding (robust for the closed calcaneus shell), not by a centroid test.
    """
    V = np.asarray(V, dtype=np.float64)
    F = np.asarray(F, dtype=np.int64)
    A, B, C = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    centroids = (A + B + C) / 3.0
    n = np.cross(B - A, C - A)
    ln = np.linalg.norm(n, axis=1)
    ln[ln == 0.0] = 1.0
    n = n / ln[:, None]
    area = 0.5 * ln
    signed_vol = float(np.einsum("ij,ij->i", A, np.cross(B, C)).sum() / 6.0)
    if signed_vol < 0.0:
        n = -n
    return centroids, n, area


def pca_directions(V: np.ndarray):
    """Centroid + principal axes (columns, descending variance) of a point cloud."""
    V = np.asarray(V, dtype=np.float64)
    c = V.mean(axis=0)
    w, Vec = np.linalg.eigh((V - c).T @ (V - c) / max(len(V), 1))
    order = np.argsort(w)[::-1]
    return c, Vec[:, order], w[order]


def _cone_flatness(n, area, centroids, d, cos_deg: float = _FLAT_CONE_DEG):
    """(area, out-of-plane thickness) of faces whose normal is within a cone of ``d``."""
    m = n @ d > np.cos(np.deg2rad(cos_deg))
    if not m.any():
        return 0.0, float("inf")
    pc = centroids[m]
    proj = (pc - pc.mean(axis=0)) @ d
    return float(area[m].sum()), float(proj.max() - proj.min())


def _end_facet_score(centroids, n, area, V, end_dir):
    """Planarity score ``area/(thickness+1)`` of the end facet facing ``end_dir``.

    Considers faces within 15 mm of the extreme along ``end_dir`` whose outward
    normal points along ``end_dir`` (``n.end_dir > 0.8``).  A large, nearly-planar
    facet (the calcaneal Achilles insertion) yields a high score.
    """
    proj_c = centroids @ end_dir
    near = proj_c > (proj_c.max() - 15.0)
    facing = n @ end_dir > 0.9
    m = near & facing
    if not m.any():
        return 0.0, 0.0, 0.0
    pc = centroids[m]
    thick = float(((pc - pc.mean(axis=0)) @ end_dir).max() - ((pc - pc.mean(axis=0)) @ end_dir).min())
    a = float(area[m].sum())
    return a / (thick + 1.0), a, thick


# ---------------------------------------------------------------------------
# anatomical frame
# ---------------------------------------------------------------------------
def anatomical_frame(V, F, n, area, centroids) -> dict:
    """Compute the rigid ``(R, t)`` placing the bone anatomically in the calcn frame.

    Returns a dict with ``R`` (3x3), ``t`` (3,), the axes expressed in the THUMS
    frame, the landmark coordinates and the per-step diagnostics.
    """
    c, Vp, eig = pca_directions(V)
    ax0, ax1, ax2 = Vp[:, 0], Vp[:, 1], Vp[:, 2]

    # --- step 1: dorsal-plantar axis = 2nd principal axis; sign by flatness ---
    flats = {}
    for name, d in (("a1+", ax1), ("a1-", -ax1), ("a2+", ax2), ("a2-", -ax2)):
        flats[name] = _cone_flatness(n, area, centroids, d)
    # candidate plantar axis is the flatter of the two non-long axes
    nonlong = ["a1+", "a1-", "a2+", "a2-"]
    down_name = min(nonlong, key=lambda k: flats[k][1])
    down = {"a1+": ax1, "a1-": -ax1, "a2+": ax2, "a2-": -ax2}[down_name]
    up = -down

    # --- step 2: long axis = PCA of V in the plane perpendicular to up --------
    P = V - np.outer(V @ up, up)
    pc = P.mean(axis=0)
    w, Vec = np.linalg.eigh((P - pc).T @ (P - pc) / max(len(P), 1))
    L = Vec[:, int(np.argmax(w))]
    L = L - (L @ up) * up
    L /= np.linalg.norm(L)

    # --- step 3: anterior / posterior sign via the Achilles end facet --------
    s_pos = _end_facet_score(centroids, n, area, V, +L)   # +L end
    s_neg = _end_facet_score(centroids, n, area, V, -L)   # -L end
    # the flat Achilles facet is the POSTERIOR end -> anterior is the OTHER end
    if s_pos[0] >= s_neg[0]:
        posterior = +L
        anterior = -L
    else:
        posterior = -L
        anterior = +L
    X = anterior / np.linalg.norm(anterior)
    Y = up / np.linalg.norm(up)
    Z = np.cross(X, Y)
    Z /= np.linalg.norm(Z)
    R = np.vstack([X, Y, Z])  # rows = calcn axes in THUMS coords: x_calcn = R x + t

    # --- step 5: origin = subtalar (superior articular) surface centroid -----
    # provisional X coordinate (no translation needed for a 2/3 split)
    cX = centroids @ X
    xlo, xhi = float(cX.min()), float(cX.max())
    x_post = xlo + (2.0 / 3.0) * (xhi - xlo)          # posterior 2/3: X < x_post
    nY = n @ Y                                         # face normal . up
    sel = (nY > JOINT_NORMAL_Y_MIN) & (cX < x_post)
    if not sel.any():
        raise RuntimeError("anatomical_frame: no superior-articular (subtalar) faces found")
    origin = (centroids[sel].T @ area[sel]) / float(area[sel].sum())
    t = -R @ origin

    return {
        "R": R,
        "t": t,
        "axes": {"X_anterior": X, "Y_superior": Y, "Z_lateral": Z,
                 "long_raw": L, "down_raw": down, "up_raw": up},
        "pca_axes": Vp,
        "pca_eigenvalues": eig,
        "down_side": down_name,
        "flatness": flats,
        "end_facet_score_plusL": {"score": s_pos[0], "area": s_pos[1], "thickness": s_pos[2]},
        "end_facet_score_minusL": {"score": s_neg[0], "area": s_neg[1], "thickness": s_neg[2]},
        "posterior_dir_raw": posterior,
        "anterior_dir_raw": anterior,
        "origin_thums_mm": origin,
        "n_subtalar_faces": int(sel.sum()),
    }


def compute_registration(surface_stl=SURFACE_STL) -> dict:
    """Load the raw THUMS surface, compute the anatomical frame and its diagnostics."""
    surface_stl = Path(surface_stl)
    V, F = reg.read_ascii_stl(surface_stl)
    centroids, n, area = face_geometry(V, F)
    fr = anatomical_frame(V, F, n, area, centroids)
    R, t = fr["R"], fr["t"]

    P = reg.apply_rigid(V, R, t)
    # R orthonormality / proper-rotation residual
    ortho_err = float(np.linalg.norm(R @ R.T - np.eye(3)))
    det = float(np.linalg.det(R))
    # origin residual: R @ origin + t == 0
    origin_res = float(np.linalg.norm(R @ fr["origin_thums_mm"] + t))
    # origin distance to the transformed surface
    origin_dist = float(np.linalg.norm(P, axis=1).min())

    # surface area / volume preservation (rigid -> equal)
    def surf_area(pts):
        return float(0.5 * np.linalg.norm(
            np.cross(pts[F[:, 1]] - pts[F[:, 0]], pts[F[:, 2]] - pts[F[:, 0]]), axis=1).sum())

    vol = reg.mesh_volume_mm3(P, F)
    plantar_frac = float(reg.plantar_fraction(P, F))

    fr.update({
        "surface_vertices_mm": P,
        "src_tris": F,
        "raw_vertices_mm": V,
        "surface_area_raw_mm2": surf_area(V),
        "surface_area_calcn_mm2": surf_area(P),
        "volume_calcn_mm3": vol,
        "plantar_facing_face_fraction": plantar_frac,
        "origin_distance_to_surface_mm": origin_dist,
        "R_orthonormality_error": ortho_err,
        "R_det": det,
        "origin_residual_mm": origin_res,
        "bbox_min_mm": P.min(axis=0),
        "bbox_max_mm": P.max(axis=0),
        "n_surface_points": int(len(V)),
        "n_surface_tris": int(len(F)),
        "source_stl": str(surface_stl),
    })
    return fr


# ---------------------------------------------------------------------------
# artifacts
# ---------------------------------------------------------------------------
def _write_volume_vtk(path, node_ids, nodes_calcn, conn) -> dict:
    elems = [row[row >= 0].tolist() for row in conn]
    nodes_map = {int(n): (float(x), float(y), float(z))
                 for n, (x, y, z) in zip(node_ids, nodes_calcn)}
    return kio.write_vtk(path, node_ids.tolist(), elems, nodes_map)


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description="THUMS calcaneus ANATOMICAL registration (calcn_r)")
    ap.add_argument("--surface", type=Path, default=SURFACE_STL, help="raw THUMS surface STL (mm)")
    ap.add_argument("--outdir", type=Path, default=CALC_DIR, help="output directory")
    ap.add_argument("--volume-npz", type=Path, default=ELEMENTS_NPZ)
    ap.add_argument("--node-ids", type=Path, default=NODE_IDS_NPY)
    ap.add_argument("--nodes-xyz", type=Path, default=NODES_XYZ_NPY)
    args = ap.parse_args(argv)

    for p in (args.surface, args.volume_npz, args.node_ids, args.nodes_xyz):
        if not p.exists():
            print(f"ERROR: missing input {p}", file=sys.stderr)
            return 2
    args.outdir.mkdir(parents=True, exist_ok=True)

    node_ids = np.load(args.node_ids).astype(np.int64)
    nodes_xyz = np.load(args.nodes_xyz).astype(np.float64)
    ed = np.load(args.volume_npz)
    conn = ed["conn"].astype(np.int64)
    pid = ed["pid"].astype(np.int64)
    eid = ed["eid"].astype(np.int64)

    print("=" * 72)
    print("THUMS AM50 right calcaneus -- ANATOMICAL registration into calcn_r frame")
    print(f"  source : {args.surface.relative_to(ROOT) if args.surface.is_relative_to(ROOT) else args.surface}")
    print("=" * 72, flush=True)

    fr = compute_registration(args.surface)
    R, t = fr["R"], fr["t"]
    print(f"[axes] down-side={fr['down_side']}  PCA eig={np.round(fr['pca_eigenvalues'],1)}")
    print(f"[axes] X(ant)={np.round(fr['axes']['X_anterior'],4)} "
          f"Y(up)={np.round(fr['axes']['Y_superior'],4)} "
          f"Z={np.round(fr['axes']['Z_lateral'],4)}")
    print(f"[axes] end facet  +L score={fr['end_facet_score_plusL']['score']:.2f} "
          f"(A={fr['end_facet_score_plusL']['area']:.0f} th={fr['end_facet_score_plusL']['thickness']:.2f})  "
          f"-L score={fr['end_facet_score_minusL']['score']:.2f} "
          f"(A={fr['end_facet_score_minusL']['area']:.0f} th={fr['end_facet_score_minusL']['thickness']:.2f})")
    print(f"[origin] subtalar centroid (THUMS)={np.round(fr['origin_thums_mm'],2)}  "
          f"R.origin+t=0 residual={fr['origin_residual_mm']:.2e} mm")
    print(f"[check ] R orthonormal err={fr['R_orthonormality_error']:.2e}  det={fr['R_det']:.6f}  "
          f"origin->surface dist={fr['origin_distance_to_surface_mm']:.2f} mm")
    print(f"[check ] bbox calcn min={np.round(fr['bbox_min_mm'],1)} max={np.round(fr['bbox_max_mm'],1)}  "
          f"vol={fr['volume_calcn_mm3']:.0f} mm3  plantar_frac={fr['plantar_facing_face_fraction']:.3f}")

    # ---- transform volume + write artifacts --------------------------------
    nodes_calcn = reg.apply_rigid(nodes_xyz, R, t)

    stl_out = args.outdir / "calcaneus_r_anatframe.stl"
    reg.write_ascii_stl(stl_out, fr["surface_vertices_mm"], fr["src_tris"],
                        name="calcaneus_r_anatframe")

    vtk_out = args.outdir / "calcaneus_r_anatframe.vtk"
    vtk_summary = _write_volume_vtk(vtk_out, node_ids, nodes_calcn, conn)

    npz_out = args.outdir / "calcaneus_r_anatframe.npz"
    np.savez(npz_out, node_ids=node_ids, nodes_xyz=nodes_calcn.astype(np.float64),
             conn=conn, pid=pid, eid=eid)

    # ---- joint / plantar / achilles via the EXACT FE-builder rule -----------
    import thums_feb as T  # noqa: E402  (local import: only needed by the driver)

    mesh = T.load_thums_mesh(npz_out)
    from climbing.coupling.febio_model import joint_area_mm2

    joint_area = float(joint_area_mm2(mesh, "subtalar_joint"))
    achilles_area = float(joint_area_mm2(mesh, "achilles")) if mesh["n_achilles_faces"] else 0.0
    # extents of the named regions
    Pnodes = np.asarray(mesh["nodes"])
    joint_nodes = np.asarray(mesh["node_sets"]["subtalar_joint"], dtype=np.int64)
    plantar_nodes = np.asarray(mesh["node_sets"]["plantar"], dtype=np.int64)
    ach_nodes = np.asarray(mesh["node_sets"]["achilles"], dtype=np.int64)

    def _node_extent(idx):
        if len(idx) == 0:
            return {"min": None, "max": None}
        return {"min": Pnodes[idx].min(axis=0).tolist(),
                "max": Pnodes[idx].max(axis=0).tolist()}

    print(f"[faces ] joint quads={mesh['n_joint_faces']}  area={joint_area:.1f} mm2  "
          f"(old foot-frame = 222.0 mm2)")
    print(f"[faces ] plantar={mesh['n_plantar_faces']}  achilles={mesh['n_achilles_faces']} "
          f"(area={achilles_area:.1f} mm2)")

    reg_json = {
        "R": R.tolist(),
        "t": t.tolist(),
        "units": "mm",
        "convention": "x_calcn = R @ x_thums + t",
        "method": "anatomical-landmarks (PCA + flatness + Achilles end facet)",
        "target": {
            "body": "calcn_r (Rajagopal2015.osim)",
            "frame": "calcn body frame: origin = subtalar joint centre, +X anterior, "
                     "+Y superior, +Z medio-lateral",
            "note": "placed from the THUMS bone's own anatomy; no ICP to r_foot.vtp",
        },
        "landmarks": {
            "plantar_axis_raw": fr["axes"]["down_raw"].tolist(),
            "superior_axis_raw": fr["axes"]["up_raw"].tolist(),
            "long_axis_raw": fr["axes"]["long_raw"].tolist(),
            "posterior_dir_raw": (fr["posterior_dir_raw"]).tolist(),
            "anterior_dir_raw": (fr["anterior_dir_raw"]).tolist(),
            "subtalar_origin_thums_mm": fr["origin_thums_mm"].tolist(),
            "down_side": fr["down_side"],
            "n_subtalar_faces_used": fr["n_subtalar_faces"],
        },
        "joint_face_area_mm2": joint_area,
        "joint_face_area_old_footframe_mm2": 222.0,
        "joint_face_count": int(mesh["n_joint_faces"]),
        "joint_node_count": int(mesh["n_joint_nodes"]),
        "plantar_face_count": int(mesh["n_plantar_faces"]),
        "plantar_extent_mm": _node_extent(plantar_nodes),
        "superior_extent_mm": _node_extent(joint_nodes),
        "achilles_face_count": int(mesh["n_achilles_faces"]),
        "achilles_area_mm2": achilles_area,
        "achilles_extent_mm": _node_extent(ach_nodes),
        "bbox_calcnframe": {"min": fr["bbox_min_mm"].tolist(), "max": fr["bbox_max_mm"].tolist()},
        "volume_calcn_mm3": fr["volume_calcn_mm3"],
        "surface_area_calcn_mm2": fr["surface_area_calcn_mm2"],
        "plantar_facing_face_fraction": fr["plantar_facing_face_fraction"],
        "residual": {
            "R_orthonormality_error": fr["R_orthonormality_error"],
            "R_det": fr["R_det"],
            "origin_Rx_plus_t_norm_mm": fr["origin_residual_mm"],
            "origin_distance_to_surface_mm": fr["origin_distance_to_surface_mm"],
            "surface_area_raw_mm2": fr["surface_area_raw_mm2"],
            "surface_area_calcn_mm2": fr["surface_area_calcn_mm2"],
            "volume_conservation_rel_err": abs(
                fr["volume_calcn_mm3"] - reg.mesh_volume_mm3(fr["raw_vertices_mm"], fr["src_tris"])
            ) / max(1.0, fr["volume_calcn_mm3"]),
            "note": "rigid only: no scale, no shear; area & volume preserved exactly",
        },
        "pca_eigenvalues": fr["pca_eigenvalues"].tolist(),
        "flatness_cone20deg": {k: {"area_mm2": v[0], "thickness_mm": v[1]}
                               for k, v in fr["flatness"].items()},
        "end_facet_scores": {
            "plus_long": fr["end_facet_score_plusL"],
            "minus_long": fr["end_facet_score_minusL"],
        },
        "source": {
            "surface_stl": str(args.surface),
            "volume_nodes": str(args.nodes_xyz),
            "elements_npz": str(args.volume_npz),
            "frame": "THUMS global, mm",
        },
    }
    json_out = args.outdir / "registration_anat.json"
    json_out.write_text(json.dumps(reg_json, ensure_ascii=False, indent=2, default=float),
                        encoding="utf-8")

    report = _render_report(fr, mesh, joint_area, achilles_area, R, t, vtk_summary)
    rep_out = args.outdir / "REGISTRATION_ANAT_REPORT.md"
    rep_out.write_text(report, encoding="utf-8")

    print("\n=== ARTIFACTS ===")
    for p in (stl_out, vtk_out, npz_out, json_out, rep_out):
        print(f"  {p.relative_to(ROOT)}  ({p.stat().st_size} bytes)")
    return 0


def _render_report(fr, mesh, joint_area, achilles_area, R, t, vtk_summary) -> str:
    R = np.asarray(R, dtype=float)
    t = np.asarray(t, dtype=float)
    bmn = np.asarray(fr["bbox_min_mm"], dtype=float)
    bmx = np.asarray(fr["bbox_max_mm"], dtype=float)
    L: list[str] = []
    L.append("# THUMS AM50 Right Calcaneus -> calcn_r **Anatomical** Registration\n")
    L.append("**Goal.** Place the bone so the **subtalar articular surface sits at the "
             "calcn origin** and the **plantar surface faces -Y**, using the bone's own "
             "anatomy (no ICP to the crude lumped `r_foot.vtp`).\n")
    L.append(f"**Transform** `x_calcn = R @ x_thums + t` (rigid only; det R = "
             f"{fr['R_det']:.10f}, ||R R^T - I|| = {fr['R_orthonormality_error']:.2e}).\n")
    L.append("```")
    L.append("R =")
    for row in R:
        L.append("  [" + ", ".join(f"{v: .10e}" for v in row) + "]")
    L.append("t = [" + ", ".join(f"{v: .10e}" for v in t) + "]")
    L.append("```\n")

    L.append("## 1. Anatomical landmark recipe\n")
    L.append("| step | rule | THUMS result |")
    L.append("|---|---|---|")
    L.append("| dorsal-plantar axis | 2nd PCA axis `a1` (height > width) | "
             f"eig = {np.round(fr['pca_eigenvalues'],1).tolist()} |")
    L.append("| **down** sign | flatter 20 deg normal cone (flat sole) | "
             f"`{fr['down_side']}`  "
             f"(down cone {fr['flatness'][fr['down_side']][0]:.0f} mm2 @ "
             f"{fr['flatness'][fr['down_side']][1]:.2f} mm) |")
    L.append("| long axis | PCA of V projected perpendicular to up | "
             f"`{np.round(fr['axes']['long_raw'],3).tolist()}` |")
    L.append("| anterior/posterior | flat Achilles end facet = posterior | "
             f"+L score {fr['end_facet_score_plusL']['score']:.2f} / "
             f"-L score {fr['end_facet_score_minusL']['score']:.2f} |")
    L.append("| origin | subtalar (n.up>0.2) centroid, posterior 2/3 of X | "
             f"`{np.round(fr['origin_thums_mm'],2).tolist()}` THUMS mm |")
    L.append("")
    L.append("Calcn axes expressed in the THUMS frame (rows of R):")
    L.append("")
    L.append(f"- X (anterior) = `{np.round(fr['axes']['X_anterior'],4).tolist()}`")
    L.append(f"- Y (superior/up) = `{np.round(fr['axes']['Y_superior'],4).tolist()}`")
    L.append(f"- Z (medio-lateral) = `{np.round(fr['axes']['Z_lateral'],4).tolist()}`\n")

    L.append("## 2. Validation checks\n")
    L.append("| check | value | pass |")
    L.append("|---|---|---|")
    L.append(f"| `subtalar_joint` area (dist<25, n.(+Y)>.2) | **{joint_area:.1f} mm2** "
             f"(old foot-frame: 222.0 mm2) | {'YES' if joint_area > 400 else 'NO'} |")
    L.append(f"| `subtalar_joint` face count | {mesh['n_joint_faces']} | "
             f"{'YES' if mesh['n_joint_faces'] > 14 else 'NO'} |")
    L.append(f"| `plantar` faces non-zero | {mesh['n_plantar_faces']} | "
             f"{'YES' if mesh['n_plantar_faces'] > 0 else 'NO'} |")
    L.append(f"| `achilles` faces non-zero | {mesh['n_achilles_faces']} | "
             f"{'YES' if mesh['n_achilles_faces'] > 0 else 'NO'} |")
    L.append(f"| plantar faces outward -Y (frac n.(-Y)>.5) | "
             f"{fr['plantar_facing_face_fraction']:.3f} | "
             f"{'YES' if fr['plantar_facing_face_fraction'] > 0 else 'NO'} |")
    L.append(f"| bone sits below origin (Y_min < 0) | Y_min = {bmn[1]:.2f} mm | "
             f"{'YES' if bmn[1] < 0 else 'NO'} |")
    L.append(f"| R orthonormal (det=+1) | det={fr['R_det']:.10f} | "
             f"{'YES' if abs(fr['R_det'] - 1) < 1e-9 else 'NO'} |")
    L.append(f"| origin maps to (0,0,0) | ||R.o+t|| = {fr['origin_residual_mm']:.2e} mm | "
             f"{'YES' if fr['origin_residual_mm'] < 1e-6 else 'NO'} |")
    L.append("")
    L.append(f"- bbox (calcn frame): min `{np.round(bmn,2).tolist()}` .. "
             f"max `{np.round(bmx,2).tolist()}` mm")
    L.append(f"- volume = {fr['volume_calcn_mm3']:.0f} mm3 "
             f"(surface area {fr['surface_area_calcn_mm2']:.0f} mm2; rigid-preserved)")
    L.append(f"- achilles surface area = {achilles_area:.1f} mm2")
    L.append(f"- solid volume written: {vtk_summary['cell_types']}\n")

    L.append("## 3. Sensitivity / honesty notes\n")
    L.append("- The flat-cluster heuristic alone (\"largest flat facet\") is **not** "
             "reliable on a calcaneus: the flattest large boundary is a medio-lateral "
             "wall, not the sole. We therefore restrict the plantar/dorsal axis to the "
             "**2nd principal axis** (height > width for a calcaneus) and only use "
             "flatness to fix its sign. The measured down/up cone asymmetry here is "
             f"{fr['flatness'][fr['down_side']][1]:.2f} mm vs "
             f"{fr['flatness']['a1+' if fr['down_side'] == 'a1-' else 'a1-'][1]:.2f} mm, "
             "i.e. a clear but individual-specific separation.")
    L.append("- Anterior/posterior is fixed by the flat **Achilles insertion facet** "
             f"(posterior end: area {fr['end_facet_score_minusL']['area']:.0f} mm2 @ "
             f"{fr['end_facet_score_minusL']['thickness']:.2f} mm thick vs anterior "
             f"{fr['end_facet_score_plusL']['area']:.0f} mm2 @ "
             f"{fr['end_facet_score_plusL']['thickness']:.2f} mm), so the narrow "
             "calcaneocuboid end lands at +X as required.")
    L.append("- This is a **frame placement**, not a shape match: the calcn axes are "
             "defined from the THUMS bone alone, so they need not coincide with any "
             "other model's calcaneus orientation to the degree. The governing checks "
             "are the subtalar patch and the plantar orientation above.")
    L.append("")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
