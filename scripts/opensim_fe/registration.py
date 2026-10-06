"""registration.py -- Rigid registration of the THUMS AM50 right calcaneus into the
OpenSim Rajagopal ``calcn_r`` body frame (the "calcn frame").

Option C of the OpenSim->FEBio calcaneus study: transplant the real AM50 THUMS
right-calcaneus geometry into the OpenSim/calcn load space.  The OpenSim loads
(subtalar joint surface load + Achilles traction + plantar fixation) live in the
``calcn_r`` body frame, so the bone MUST be expressed in that frame for the load
space to be consistent.

Frame
-----
Target = ``r_foot.vtp`` (the ``calcn_r`` body geometry in the Rajagopal model).
Its coordinates are METERS in the calcn body frame; the origin is the subtalar
joint centre.  Multiplying by 1000 gives millimetres.  The calcn axes at default
stance are approximately::

    +X  anterior
    +Y  superior (up)
    +Z  medio-lateral

so the posterior tuberosity points toward -X and the plantar surface toward -Y.

Units: THUMS source and the target are both handled in **millimetres** here.

Rigid only.  The transform is (R 3x3, t 3) with ``x_calcn = R @ x_thums + t``;
no scaling and no shear are ever applied.

Pipeline (see :func:`register_to_calcn_frame`)
-----------------------------------------------
1. target -> mm via :func:`load_calcn_target`;
2. initial alignment via centroid + principal axes with sign resolution
   (:func:`pca_align`, plus a small deterministic search over the proper
   axis/sign variants);
3. rigid refinement against the target surface via vtk
   ``vtkIterativeClosestPointTransform`` (:func:`icp_refine`);
4. an independent landmark-based Kabsch fit (:func:`landmark_fit`) is computed as
   well; the candidate with the better physical validation is selected;
5. the caller applies the transform to both the surface and the solid volume.

All routines are deterministic: no random numbers anywhere.

Public API
----------
load_calcn_target(vtp_path)                     -> (points_mm (N,3), faces)
pca_align(src_pts, dst_pts)                     -> (R, t)
icp_refine(src_pts, dst_pts, R0, t0, ...)       -> (R, t, rms_mm, n_iter)
apply_rigid(pts, R, t)                          -> pts
register_to_calcn_frame(src_stl, target_vtp)    -> dict
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np

# --------------------------------------------------------------------------- #
#  module constants
# --------------------------------------------------------------------------- #

#: calcn-frame anatomical axes (default stance)
AXIS_ANTERIOR = np.array([1.0, 0.0, 0.0])
AXIS_SUPERIOR = np.array([0.0, 1.0, 0.0])
AXIS_LATERAL = np.array([0.0, 0.0, 1.0])

#: OpenSim Achilles-insertion landmark (mm).  Used ONLY as a soft validation
#: check -- see the honesty note in the report.  It is expressed in the frame
#: the earlier computation used (see REGISTRATION_REPORT.md).
ACHILLES_LANDMARK_MM = np.array([-110.6, 46.1, 78.6])

#: Plausible calcaneus volume range (mm^3).
VOLUME_RANGE_MM3 = (4.0e4, 1.3e5)

#: A fit is called "poor" above this RMS (mm) -- reported, never hidden.
POOR_RMS_MM = 5.0

#: Fraction of the foot length (from the heel, along +X) treated as the
#: "calcaneus / hindfoot" region of the target for landmark extraction.
HINDFOOT_FRACTION = 0.45


# --------------------------------------------------------------------------- #
#  small linear-algebra helpers
# --------------------------------------------------------------------------- #


def apply_rigid(pts: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Apply ``x -> R @ x + t`` to an ``(N, 3)`` array of points (row-wise)."""
    P = np.asarray(pts, dtype=np.float64)
    R = np.asarray(R, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64).reshape(3)
    return P @ R.T + t


def _pca(pts: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Centroid + principal axes (columns, descending variance) of a point cloud."""
    P = np.asarray(pts, dtype=np.float64)
    c = P.mean(axis=0)
    X = P - c
    cov = (X.T @ X) / max(len(P), 1)
    w, V = np.linalg.eigh(cov)            # ascending eigenvalues
    order = np.argsort(w)[::-1]
    return c, V[:, order], w[order]


def kabsch(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rigid transform mapping corresponding points ``A -> B`` (SVD, no scaling).

    Returns ``(R, t)`` such that ``R @ A_i + t ~= B_i`` (Kabsch/Umeyama, rigid
    only: det(R) = +1).
    """
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    ca = A.mean(axis=0)
    cb = B.mean(axis=0)
    H = (A - ca).T @ (B - cb)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    if d == 0:
        d = 1.0
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = cb - R @ ca
    return R, t


# --------------------------------------------------------------------------- #
#  target loading
# --------------------------------------------------------------------------- #


def load_calcn_target(vtp_path) -> tuple[np.ndarray, list[tuple[int, int, int]]]:
    """Read ``r_foot.vtp`` (calcn frame, METERS) and return millimetre points+faces.

    Returns
    -------
    points_mm : (N, 3) float64
        Vertex coordinates in the calcn frame, converted m -> mm (x1000).
    faces : list[(i, j, k)]
        Triangulated polygons (the .vtp mixes triangles and quads; quads are
        fan-triangulated here) indexing into ``points_mm``.
    """
    import vtk

    path = Path(vtp_path)
    if not path.exists():
        raise FileNotFoundError(f"calcn target geometry not found: {path}")

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    pd = reader.GetOutput()
    if pd is None or pd.GetNumberOfPoints() == 0:
        raise ValueError(f"empty target polydata: {path}")

    n = pd.GetNumberOfPoints()
    pts = np.empty((n, 3), dtype=np.float64)
    for i in range(n):
        x, y, z = pd.GetPoint(i)
        pts[i] = (x, y, z)
    pts *= 1000.0  # m -> mm

    faces: list[tuple[int, int, int]] = []
    for ci in range(pd.GetNumberOfCells()):
        ids = pd.GetCell(ci).GetPointIds()
        m = ids.GetNumberOfIds()
        poly = [int(ids.GetId(k)) for k in range(m)]
        for k in range(1, m - 1):          # fan triangulation (works for tris & quads)
            faces.append((poly[0], poly[k], poly[k + 1]))
    return pts, faces


def read_ascii_stl(path) -> tuple[np.ndarray, np.ndarray]:
    """Read an ASCII STL -> (unique_verts (V,3), triangles (F,3)).

    Vertices are de-duplicated by rounding to 1e-6 mm so the returned triangle
    indices reference the unique vertex array (this is exactly how the extracted
    ``calcaneus_r_surface.stl`` is structured).
    """
    verts: list[list[float]] = []
    with open(path, "r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("vertex"):
                parts = s.split()
                verts.append([float(parts[1]), float(parts[2]), float(parts[3])])
    V = np.asarray(verts, dtype=np.float64)
    U, inv = np.unique(np.round(V, 6), axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    tris = np.array(
        [(inv[i], inv[i + 1], inv[i + 2]) for i in range(0, len(inv), 3)],
        dtype=np.int64,
    )
    return U, tris


# --------------------------------------------------------------------------- #
#  initial alignment
# --------------------------------------------------------------------------- #


def pca_align(src_pts: np.ndarray, dst_pts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Centroid + principal-axis alignment with the calcn-axis convention.

    The source principal axes (descending variance) are mapped onto the calcn
    axes so that

        longest axis  -> X (anterior/posterior),
        second axis   -> Y (superior),
        third axis    -> Z (medio-lateral),

    and the source centroid is mapped onto the destination centroid.  The
    *sign* of each principal axis is inherently ambiguous from PCA alone; the
    deterministic sign search in :func:`register_to_calcn_frame` resolves the
    correct anatomical hemisphere (long axis -> +X anterior; superior -> +Y up).
    """
    cs, Vs, _ = _pca(src_pts)
    cd, _Vd, _ = _pca(dst_pts)
    R = np.eye(3) @ Vs.T          # maps principal axis i -> world axis i
    if np.linalg.det(R) < 0.0:    # keep a proper rotation (right-handed)
        R = np.diag([1.0, 1.0, -1.0]) @ R
    t = cd - R @ cs
    return R, t


def _proper_variants(Vs: np.ndarray) -> list[tuple[tuple, np.ndarray]]:
    """Enumerate the deterministic initial-orientation variants.

    ``Vs`` columns are the source principal axes.  We keep the long axis on X
    and try both pairings of the remaining two axes to (Y, Z) and all sign
    combinations, keeping only proper rotations.  Each variant is returned with
    a small metadata tag for the report.
    """
    variants: list[tuple[tuple, np.ndarray]] = []
    for s_long in (1, -1):
        for i2, i3 in ((1, 2), (2, 1)):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    W = np.zeros((3, 3))
                    W[0, 0] = s_long
                    W[1, i2] = s2
                    W[2, i3] = s3
                    R = W @ Vs.T
                    if np.linalg.det(R) < 0.0:
                        continue
                    variants.append(((s_long, i2, i3, s2, s3), R))
    return variants


# --------------------------------------------------------------------------- #
#  ICP refinement (rigid, vtk)
# --------------------------------------------------------------------------- #


def _polydata(points: np.ndarray, faces: Sequence[Sequence[int]] | None = None):
    """Build a vtkPolyData from points and optional triangle faces.

    When ``faces`` is None a vertex cell is created per point so the polydata is
    a valid target for vtk's cell locator.
    """
    import vtk

    vp = vtk.vtkPoints()
    for p in np.asarray(points, dtype=np.float64):
        vp.InsertNextPoint(float(p[0]), float(p[1]), float(p[2]))
    pd = vtk.vtkPolyData()
    pd.SetPoints(vp)

    cells = vtk.vtkCellArray()
    if faces is None:
        for i in range(pd.GetNumberOfPoints()):
            cells.InsertNextCell(1)
            cells.InsertCellPoint(i)
        pd.SetVerts(cells)
    else:
        for a, b, c in faces:
            cells.InsertNextCell(3)
            cells.InsertCellPoint(int(a))
            cells.InsertCellPoint(int(b))
            cells.InsertCellPoint(int(c))
        pd.SetPolys(cells)
    return pd


def icp_refine(
    src_pts: np.ndarray,
    dst_pts: np.ndarray,
    R0: np.ndarray,
    t0: np.ndarray,
    *,
    src_faces: Sequence[Sequence[int]] | None = None,
    dst_faces: Sequence[Sequence[int]] | None = None,
    max_iter: int = 200,
    tol_mm: float = 1e-5,
) -> tuple[np.ndarray, np.ndarray, float, int]:
    """Rigid ICP refinement of ``(R0, t0)`` with vtkIterativeClosestPointTransform.

    Parameters
    ----------
    src_pts, dst_pts : (N,3)/(M,3)
        Source (calcaneus) and target (r_foot) points, both in mm.
    R0, t0 :
        Initial rigid transform (``x -> R0 x + t0``) used to seed ICP.
    src_faces, dst_faces : optional triangle index arrays
        If given, the polydata carry real surfaces; otherwise a vertex cloud is
        used (still a valid vtk ICP target).
    max_iter, tol_mm :
        Iteration cap and mean-distance convergence tolerance (mm).

    Returns
    -------
    (R, t, rms_mm, n_iter)
        Final rigid transform, the point-to-point RMS residual of the source
        vertices to the target points (mm), and the vtk iteration count.
    """
    import vtk

    src_pts = np.asarray(src_pts, dtype=np.float64)
    dst_pts = np.asarray(dst_pts, dtype=np.float64)

    # seed source with the initial transform
    src_seed = apply_rigid(src_pts, R0, t0)
    src_pd = _polydata(src_seed, src_faces)
    dst_pd = _polydata(dst_pts, dst_faces)

    icp = vtk.vtkIterativeClosestPointTransform()
    icp.SetSource(src_pd)
    icp.SetTarget(dst_pd)
    icp.GetLandmarkTransform().SetModeToRigidBody()
    icp.SetMaximumNumberOfIterations(int(max_iter))
    icp.SetCheckMeanDistance(1)
    icp.SetMaximumMeanDistance(float(tol_mm))
    icp.SetStartByMatchingCentroids(0)
    icp.Modified()
    icp.Update()

    M = icp.GetMatrix()
    Rm = np.array([[M.GetElement(i, j) for j in range(3)] for i in range(3)])
    tm = np.array([M.GetElement(i, 3) for i in range(3)])
    R = Rm @ np.asarray(R0, dtype=np.float64)
    t = Rm @ np.asarray(t0, dtype=np.float64).reshape(3) + tm

    # RMS of the transformed source to the nearest target point (mm)
    from scipy.spatial import cKDTree

    P = apply_rigid(src_pts, R, t)
    d, _ = cKDTree(dst_pts).query(P)
    rms = float(np.sqrt(np.mean(d ** 2)))
    n_iter = int(icp.GetNumberOfIterations())
    return R, t, rms, n_iter


# --------------------------------------------------------------------------- #
#  landmark (Kabsch) alternative
# --------------------------------------------------------------------------- #


def _extreme(
    pts: np.ndarray, axis: np.ndarray, sign: float, top_frac: float = 0.03
) -> np.ndarray:
    """Robust extreme point along ``sign * axis`` (mean of the top ``top_frac``)."""
    proj = pts @ (sign * axis)
    k = max(1, int(round(len(pts) * top_frac)))
    idx = np.argsort(proj)[-k:]
    return pts[idx].mean(axis=0)


def _landmarks(pts: np.ndarray) -> dict[str, np.ndarray]:
    """Six anatomical extremes of a cloud already expressed in calcn axes."""
    return {
        "posterior": _extreme(pts, AXIS_ANTERIOR, -1.0),
        "anterior": _extreme(pts, AXIS_ANTERIOR, +1.0),
        "plantar": _extreme(pts, AXIS_SUPERIOR, -1.0),
        "superior": _extreme(pts, AXIS_SUPERIOR, +1.0),
        "medial": _extreme(pts, AXIS_LATERAL, -1.0),
        "lateral": _extreme(pts, AXIS_LATERAL, +1.0),
    }


_LANDMARK_ORDER = ("posterior", "anterior", "plantar", "superior", "medial", "lateral")


def landmark_fit(
    src_pts: np.ndarray,
    dst_pts: np.ndarray,
    *,
    hindfoot_fraction: float = HINDFOOT_FRACTION,
) -> dict:
    """Rigid fit from corresponding anatomical extremes (Kabsch via numpy SVD).

    The target is first cropped to its hindfoot region -- the posterior
    ``hindfoot_fraction`` of the foot length measured from the heel -- because
    the anterior-most point of the whole foot is the toe, which has no
    calcaneus counterpart.  The source's own PCA frame is then oriented with the
    same deterministic variant search and matched to those target extremes by
    axis; the variant with the smallest Kabsch residual is kept.

    Returns a dict with ``R``, ``t``, ``src_landmarks``, ``dst_landmarks``,
    ``variant`` and ``landmark_rms_mm``.
    """
    dst_pts = np.asarray(dst_pts, dtype=np.float64)
    x_min = float(dst_pts[:, 0].min())
    x_max = float(dst_pts[:, 0].max())
    x_cut = x_min + hindfoot_fraction * (x_max - x_min)
    mask = dst_pts[:, 0] <= x_cut
    dh = dst_pts[mask]
    if len(dh) < 6:
        dh = dst_pts
    dst_lm = _landmarks(dh)
    dst_arr = np.array([dst_lm[k] for k in _LANDMARK_ORDER])

    cs, Vs, _ = _pca(src_pts)
    best = None
    for meta, R in _proper_variants(Vs):
        P = apply_rigid(src_pts, R, np.zeros(3))
        src_lm = _landmarks(P)
        src_arr = np.array([src_lm[k] for k in _LANDMARK_ORDER])
        Rk, tk = kabsch(src_arr, dst_arr)
        resid = float(np.sqrt(np.mean(np.sum((apply_rigid(src_arr, Rk, tk) - dst_arr) ** 2, axis=1))))
        if best is None or resid < best["landmark_rms_mm"]:
            best = {
                "R": Rk,
                "t": tk,
                "variant": meta,
                "landmark_rms_mm": resid,
                "src_landmarks": src_arr,
                "dst_landmarks": dst_arr,
                "x_cut_mm": float(x_cut),
            }
    assert best is not None
    return best


# --------------------------------------------------------------------------- #
#  geometry validation helpers
# --------------------------------------------------------------------------- #


def mesh_volume_mm3(pts: np.ndarray, faces: np.ndarray) -> float:
    """Enclosed volume (mm^3) of a triangle surface via the divergence theorem."""
    A = pts[faces[:, 0]]
    B = pts[faces[:, 1]]
    C = pts[faces[:, 2]]
    return float(abs(np.einsum("ij,ij->i", A, np.cross(B, C)).sum() / 6.0))


def outward_face_normals(pts: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Unit face normals re-oriented to point away from the mesh centroid."""
    A, B, C = pts[faces[:, 0]], pts[faces[:, 1]], pts[faces[:, 2]]
    n = np.cross(B - A, C - A)
    ln = np.linalg.norm(n, axis=1)
    ln[ln == 0.0] = 1.0
    n = n / ln[:, None]
    fem = (A + B + C) / 3.0
    c = pts.mean(axis=0)
    s = np.sign(np.einsum("ij,ij->i", n, fem - c))
    s[s == 0.0] = 1.0
    return n * s[:, None]


def plantar_fraction(pts: np.ndarray, faces: np.ndarray) -> float:
    """Fraction of outward faces with ``normal . (-Y) > 0.5`` (plantar-facing)."""
    n = outward_face_normals(pts, faces)
    return float((n @ np.array([0.0, -1.0, 0.0]) > 0.5).mean())


def _nearest(pts: np.ndarray, q: np.ndarray) -> tuple[float, np.ndarray]:
    from scipy.spatial import cKDTree

    d, i = cKDTree(pts).query(q)
    return float(d), pts[i].copy()


def validate_candidate(pts: np.ndarray, faces: np.ndarray, rms_mm: float) -> dict:
    """Compute the physical validation block for one candidate surface."""
    tree_pts = pts
    origin_d, origin_near = _nearest(tree_pts, np.zeros(3))
    ach_d, ach_near = _nearest(tree_pts, ACHILLES_LANDMARK_MM)
    vol = mesh_volume_mm3(pts, faces)
    pf = plantar_fraction(pts, faces)

    apex_idx = int(pts[:, 0].argmin())
    sup_idx = int(pts[:, 1].argmax())
    bbox_min = pts.min(axis=0)
    bbox_max = pts.max(axis=0)
    height = float(bbox_max[1] - bbox_min[1])
    width = float(bbox_max[2] - bbox_min[2])

    valid = (
        origin_d <= 25.0
        and pf >= 0.05
        and VOLUME_RANGE_MM3[0] <= vol <= VOLUME_RANGE_MM3[1]
        and pts[apex_idx, 0] <= 0.0
        and height >= width
    )
    return {
        "rms_mm": float(rms_mm),
        "origin_distance_to_surface_mm": origin_d,
        "origin_nearest_point_mm": origin_near.tolist(),
        "achilles_nearest_surface_distance_mm": ach_d,
        "achilles_nearest_point_mm": ach_near.tolist(),
        "volume_mm3": vol,
        "volume_cm3": vol / 1000.0,
        "volume_in_plausible_range": bool(VOLUME_RANGE_MM3[0] <= vol <= VOLUME_RANGE_MM3[1]),
        "plantar_facing_face_fraction": pf,
        "posterior_apex_mm": pts[apex_idx].tolist(),
        "superior_apex_mm": pts[sup_idx].tolist(),
        "bbox_min_mm": bbox_min.tolist(),
        "bbox_max_mm": bbox_max.tolist(),
        "height_mm": height,
        "width_mm": width,
        "valid": bool(valid),
    }


# --------------------------------------------------------------------------- #
#  top-level registration
# --------------------------------------------------------------------------- #


def register_to_calcn_frame(
    src_stl,
    target_vtp_str_or_path,
    *,
    max_iter: int = 200,
    hindfoot_fraction: float = HINDFOOT_FRACTION,
) -> dict:
    """Register the THUMS calcaneus surface into the calcn frame.

    Parameters
    ----------
    src_stl :
        Path to the source ``calcaneus_r_surface.stl`` (ASCII, mm, THUMS global
        frame) -- or a pre-loaded ``(verts, tris)`` tuple.
    target_vtp_str_or_path :
        Path to ``r_foot.vtp`` (calcn frame, metres).

    Returns
    -------
    dict
        ``R`` (3x3), ``t`` (3), ``rms_mm``, ``n_iter``, ``method``,
        ``n_src_points``, ``n_dst_points``, transformed surface arrays
        (``surface_vertices_mm``, ``src_tris``) and a ``landmark_checks`` block,
        plus the full ``candidates`` list that was evaluated.
    """
    if isinstance(src_stl, (str, Path)):
        src_verts, src_tris = read_ascii_stl(src_stl)
    else:
        src_verts, src_tris = src_stl
        src_verts = np.asarray(src_verts, dtype=np.float64)
        src_tris = np.asarray(src_tris, dtype=np.int64)

    dst_pts, dst_faces = load_calcn_target(target_vtp_str_or_path)
    dst_faces_arr = np.asarray(dst_faces, dtype=np.int64)

    cs, Vs, _ = _pca(src_verts)
    cd = dst_pts.mean(axis=0)
    heel = dst_pts[int(dst_pts[:, 0].argmin())]

    candidates: list[dict] = []

    # ---- PCA-initialised candidates + ICP -------------------------------- #
    for meta, R_axes in _proper_variants(Vs):
        base = apply_rigid(src_verts, R_axes, np.zeros(3))
        translations = {
            "centroid": cd - base.mean(axis=0),
            "posterior_anchor": heel - base[int(base[:, 0].argmin())],
        }
        for tname, t0 in translations.items():
            R1, t1, rms, n_it = icp_refine(
                src_verts, dst_pts, R_axes, t0,
                src_faces=src_tris, dst_faces=dst_faces_arr, max_iter=max_iter,
            )
            P = apply_rigid(src_verts, R1, t1)
            v = validate_candidate(P, src_tris, rms)
            v.update(
                method="pca+icp", variant=meta, translation=tname,
                n_iter=n_it, R=R1, t=t1, pts=P,
            )
            candidates.append(v)

    # ---- independent landmark (Kabsch) candidate + ICP ------------------- #
    lf = landmark_fit(src_verts, dst_pts, hindfoot_fraction=hindfoot_fraction)
    R_l, t_l, rms_l, it_l = icp_refine(
        src_verts, dst_pts, lf["R"], lf["t"],
        src_faces=src_tris, dst_faces=dst_faces_arr, max_iter=max_iter,
    )
    P_l = apply_rigid(src_verts, R_l, t_l)
    v_l = validate_candidate(P_l, src_tris, rms_l)
    v_l.update(
        method="landmark+icp", variant=("kabsch",) + tuple(lf["variant"]),
        translation="kabsch", n_iter=it_l, R=R_l, t=t_l, pts=P_l,
        landmark_rms_mm=lf["landmark_rms_mm"],
    )
    candidates.append(v_l)
    lf_summary = {
        "variant": lf["variant"],
        "landmark_rms_mm": lf["landmark_rms_mm"],
        "x_cut_mm": lf["x_cut_mm"],
    }

    # ---- selection: best physical validation, tie-broken by RMS ---------- #
    valid = [c for c in candidates if c["valid"]]
    pool = valid if valid else candidates
    winner = min(pool, key=lambda c: c["rms_mm"])

    n_src = int(len(src_verts))
    n_dst = int(len(dst_pts))
    checks = {
        k: winner[k]
        for k in (
            "rms_mm", "origin_distance_to_surface_mm", "origin_nearest_point_mm",
            "achilles_nearest_surface_distance_mm", "achilles_nearest_point_mm",
            "volume_mm3", "volume_cm3", "volume_in_plausible_range",
            "plantar_facing_face_fraction", "posterior_apex_mm", "superior_apex_mm",
            "bbox_min_mm", "bbox_max_mm", "height_mm", "width_mm", "valid",
        )
    }
    checks["achilles_landmark_mm"] = ACHILLES_LANDMARK_MM.tolist()
    checks["achilles_note"] = (
        "Given landmark lies outside the r_foot.vtp bounding box in X (~100 mm "
        "posterior) and Z (~59 mm lateral), so it is NOT expressed in the "
        "r_foot.vtp calcn body frame (probably the OpenSim ground/default-stance "
        "frame). This check is informational only -- see REGISTRATION_REPORT.md."
    )

    candidates_public = [
        {
            "method": c["method"],
            "variant": list(c["variant"]),
            "translation": c["translation"],
            "rms_mm": c["rms_mm"],
            "n_iter": c["n_iter"],
            "origin_distance_to_surface_mm": c["origin_distance_to_surface_mm"],
            "plantar_facing_face_fraction": c["plantar_facing_face_fraction"],
            "volume_mm3": c["volume_mm3"],
            "posterior_apex_mm": c["posterior_apex_mm"],
            "height_mm": c["height_mm"],
            "width_mm": c["width_mm"],
            "valid": c["valid"],
        }
        for c in candidates
    ]

    return {
        "R": np.asarray(winner["R"], dtype=np.float64),
        "t": np.asarray(winner["t"], dtype=np.float64),
        "rms_mm": float(winner["rms_mm"]),
        "n_iter": int(winner["n_iter"]),
        "method": str(winner["method"]),
        "variant": list(winner["variant"]),
        "translation": str(winner["translation"]),
        "n_src_points": n_src,
        "n_dst_points": n_dst,
        "surface_vertices_mm": winner["pts"],
        "src_tris": src_tris,
        "landmark_checks": checks,
        "landmark_fit": lf_summary,
        "candidates": candidates_public,
    }


def write_ascii_stl(path, verts: np.ndarray, tris: np.ndarray, name: str = "surface") -> None:
    """Write an ASCII STL (matches the kmesh_io style). ``verts`` (V,3), ``tris`` (F,3)."""
    lines = [f"solid {name}"]
    for a, b, c in np.asarray(tris, dtype=np.int64):
        pa, pb, pc = verts[a], verts[b], verts[c]
        u = pb - pa
        v = pc - pa
        n = np.cross(u, v)
        ln = np.linalg.norm(n)
        if ln == 0.0:
            n = np.zeros(3)
        else:
            n = n / ln
        lines.append(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}")
        lines.append("    outer loop")
        for p in (pa, pb, pc):
            lines.append(f"      vertex {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
        lines.append("    endloop")
        lines.append("  endfacet")
    lines.append(f"endsolid {name}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii")


__all__ = [
    "load_calcn_target",
    "pca_align",
    "icp_refine",
    "apply_rigid",
    "register_to_calcn_frame",
    "kabsch",
    "landmark_fit",
    "read_ascii_stl",
    "mesh_volume_mm3",
    "outward_face_normals",
    "plantar_fraction",
    "validate_candidate",
    "write_ascii_stl",
    "ACHILLES_LANDMARK_MM",
]
