"""THUMS AM50 right-calcaneus -> pyfebio **two-material-domain** model builder.

Option C of the OpenSim->FEBio calcaneus study: swap the crude whole-foot
deformable geometry (``r_foot.vtp``) for the registered THUMS AM50 right
calcaneus, and split it into the two material domains the deck actually
carries:

    "calcaneus"   (CORT)  hex8   E = 15000 MPa, nu = 0.30   <- fracture-relevant
    "trabecular"  (SPON)  tet4   E = 73.4  MPa, nu = 0.45

Source mesh (mm, OpenSim ``calcn_r`` frame, origin = subtalar joint centre)::

    temp/opensim_fe/thums_calcaneus/calcaneus_r_calcnframe.npz
        node_ids (1544,)   global THUMS node ids
        nodes_xyz (1544,3) coordinates, calcn frame, mm
        conn (3957,8)      per-element connectivity (global ids)
        pid (3957,)        81001300 = CORT (634 rows, 8 distinct -> hex8)
                           81001200 = SPON (3323 rows, 4 distinct -> tet4,
                                      stored as a collapsed hex)

The SPON rows are tetrahedra written as degenerate 8-node solids: the 5th-8th
columns repeat the 4th node.  We recover the tet by taking the **first four
distinct** node ids in row order (matching ``kmesh_io``'s topology probe).

Load path (IDENTICAL to ``src/climbing/coupling/febio_model.py::
build_calcaneus_feb``)
-------------------------------------------------------------------------------
* joint load   : ghost rigid body + ``RigidForceLoad`` on ``subtalar_joint``
                 (fallback ``PressureLoad`` p = F/A);
* Achilles     : optional ``TractionLoad`` on ``achilles`` (traction vector);
* plantar      : three-way zero displacement;
* analysis     : ``STATIC``, single step, then the FEBio-4.13 ``<solver/>``
                 adaptation (``febio_model._fix_febio413``).

Boundary faces are extracted for the **mixed hex+tet** mesh.  A hex exposes
quad faces while a tet exposes triangles, so at a conformal hex/tet interface a
boundary-looking quad is fully covered by two boundary triangles; both are
removed as internal (same rule as ``kmesh_io.free_surface_triangles``), leaving
only the true outer surface.  Named surfaces reuse
``meshing._classify_surfaces`` verbatim so the thresholds are bit-identical to
``src/climbing/coupling/meshing.py``.

Units: mm-N-MPa-s throughout.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np

from climbing.coupling import meshing
from climbing.coupling.febio_model import (
    DEFAULT_JOINT_FACE,
    _fix_febio413,
    joint_area_mm2,
    joint_face_frame,
    linear_pressure_bands,
    moment_to_pressure_gradient,
)
from climbing.coupling.plantar_bc import apply_plantar_bc

__all__ = [
    "load_thums_mesh",
    "build_thums_feb",
    "hex_centroids_volumes",
    "gauge_von_mises_mixed",
    "select_fascia_attachment",
    "PID_CORT",
    "PID_SPON",
    "DEFAULT_NPZ",
    "FASCIA_X_FRAC",
    "FASCIA_Y_FRAC",
    "FASCIA_Z_FRAC",
]

# scripts/opensim_fe/thums_feb.py -> scripts/opensim_fe -> scripts -> repo root
ROOT = Path(__file__).resolve().parents[2]
THUMS_DIR = ROOT / "temp" / "opensim_fe" / "thums_calcaneus"
DEFAULT_NPZ = THUMS_DIR / "calcaneus_r_calcnframe.npz"
DEFAULT_STL = THUMS_DIR / "calcaneus_r_calcnframe.stl"

#: THUMS part ids (see scripts/opensim_fe/extract_thums_calcaneus.py)
PID_CORT = 81001300  # R_CALCANEUS_CORT  (hex8, cortical)
PID_SPON = 81001200  # R_CALCANEUS_SPON  (collapsed-hex -> tet4, trabecular)

#: element face templates.  Local ordering follows the LS-DYNA/VTK/FEBio solid
#: convention (bottom 0-3, top 4-7), so a positive-volume element yields
#: outward-facing normals once oriented away from the element centroid.
_HEX_FACES = (
    (0, 1, 2, 3),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)
_TET_FACES = ((0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3))

#: Hex volume decomposition: a fan of 6 tets about the 0-6 body diagonal.
_HEX_TET_FAN = ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6))


# ---------------------------------------------------------------------------
# mesh assembly
# ---------------------------------------------------------------------------
def _first_distinct(row: np.ndarray, k: int = 4) -> list[int]:
    """First ``k`` distinct node ids of ``row``, preserving order."""
    seen: set[int] = set()
    out: list[int] = []
    for x in row:
        xi = int(x)
        if xi not in seen:
            seen.add(xi)
            out.append(xi)
            if len(out) == k:
                break
    if len(out) != k:
        raise ValueError(f"元素行去重后不足 {k} 个节点：{row!r}")
    return out


def _mixed_boundary_faces(elements: list[np.ndarray]):
    """Outer **polygon** faces of a mixed hex+tet solid, with owner indices.

    Returns ``(polys, owners)`` where ``polys`` is a list of tuples of local
    node indices (length 4 for a hex quad face, 3 for a tet triangle) and
    ``owners[k]`` is the index into ``elements`` of the single solid owning
    face ``k``.

    A face is a boundary candidate when exactly one element references it.  A
    hex quad fully covered by two boundary triangles (the conformal hex/tet
    interface) is treated as internal together with those triangles, so the
    returned faces are the true outer shell only.  Hex faces stay quads because
    a hexahedron's boundary face *is* a quad -- classifying a quad's two
    triangles independently would split one physical face.
    """
    face_owner: dict[tuple[int, ...], list[int]] = defaultdict(list)
    face_verts: dict[tuple[int, ...], tuple[int, ...]] = {}
    for ei, e in enumerate(elements):
        tmpl = _HEX_FACES if len(e) == 8 else _TET_FACES
        for f in tmpl:
            vs = tuple(int(e[i]) for i in f)
            key = tuple(sorted(vs))
            face_owner[key].append(ei)
            face_verts.setdefault(key, vs)

    b_quads = [k for k, o in face_owner.items() if len(o) == 1 and len(k) == 4]
    b_tris = {k for k, o in face_owner.items() if len(o) == 1 and len(k) == 3}

    polys: list[tuple[int, ...]] = []
    owners: list[int] = []
    internal: set[tuple[int, ...]] = set()
    for k in b_quads:
        a, b, c, d = face_verts[k]
        covered = False
        for t1, t2 in (((a, b, c), (a, c, d)), ((a, b, d), (b, c, d))):
            k1, k2 = tuple(sorted(t1)), tuple(sorted(t2))
            if k1 in b_tris and k2 in b_tris:
                internal.update((k1, k2))
                covered = True
                break
        if not covered:
            polys.append((a, b, c, d))
            owners.append(face_owner[k][0])

    for k in b_tris - internal:
        polys.append(face_verts[k])  # type: ignore[arg-type]
        owners.append(face_owner[k][0])

    return polys, np.asarray(owners, dtype=np.int64)


def _polygon_centroids_normals(
    nodes: np.ndarray,
    polys: list[tuple[int, ...]],
    owners: np.ndarray,
    elem_centroids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Centroids + outward unit normals of a mixed tri/quad polygon list.

    Normals use Newell's method (valid for planar and warped quads) and are
    flipped to point away from the owning element's centroid.
    """
    centroids = np.empty((len(polys), 3), dtype=np.float64)
    normals = np.empty((len(polys), 3), dtype=np.float64)
    for i, poly in enumerate(polys):
        p = nodes[list(poly)]
        centroids[i] = p.mean(axis=0)
        nxt = np.roll(p, -1, axis=0)
        n = np.array(
            [
                float(np.sum((p[:, 1] - nxt[:, 1]) * (p[:, 2] + nxt[:, 2]))),
                float(np.sum((p[:, 2] - nxt[:, 2]) * (p[:, 0] + nxt[:, 0]))),
                float(np.sum((p[:, 0] - nxt[:, 0]) * (p[:, 1] + nxt[:, 1]))),
            ]
        )
        nn = float(np.linalg.norm(n))
        if nn < 1e-30:
            raise RuntimeError(f"节点 {poly} 构成的边界面法向退化")
        normals[i] = n / nn
    away = centroids - elem_centroids[owners]
    flip = np.einsum("ij,ij->i", normals, away) < 0.0
    normals[flip] = -normals[flip]
    return centroids, normals


#: Achilles region (posterior-superior tuberosity) as a FRACTION of this mesh's
#: bounding box.  The absolute ``meshing.ACHILLES_X_MAX / _Y_MIN`` constants were
#: calibrated on the IITD STL (posterior extreme ~ x = -38 mm); the registered
#: THUMS calcaneus sits ~30 mm more anterior, so that absolute box selects
#: NOTHING here.  Calibrate relative to the mesh instead.
ACHILLES_X_FRAC = 0.35  # posterior 35 % along X (X is anterior +)
ACHILLES_Y_FRAC = 0.50  # superior 50 % along Y (Y is up +)

#: Plantar-fascia (plantar aponeurosis) attachment selection.  The central band
#: of the fascia originates on the **medial process of the calcaneal
#: tuberosity** -- posterior, inferior and medial.  Fractions are of the
#: plantar-surface triangle-centroid bounding box.  ``+Z`` is **medial** for
#: this registered right calcaneus: the frame is right handed with
#: ``Z = X_anterior x Y_superior``, which points toward the body midline for a
#: right limb (verified by this bone's own anatomy: the medial tuberosity
#: process reaches ~3 mm further inferior on ``+Z``, and the medial
#: sustentaculum/articular shelf is the ``+Z`` extreme of the joint patch).
FASCIA_X_FRAC = 0.35  # posterior 35 % along X (X is anterior +)
FASCIA_Y_FRAC = 0.35  # inferior 35 % along Y (Y is up +)
FASCIA_Z_FRAC = 0.50  # medial 50 % along Z (+Z = medial for this right bone)
#: a fascia patch smaller than this is rejected (falls back to a looser rule).
FASCIA_MIN_TRIS = 3


def select_fascia_attachment(
    mesh: dict,
    *,
    x_frac: float = FASCIA_X_FRAC,
    y_frac: float = FASCIA_Y_FRAC,
    z_frac: float = FASCIA_Z_FRAC,
    min_tris: int = FASCIA_MIN_TRIS,
) -> dict:
    """Select the plantar-fascia attachment patch on the calcaneal tuberosity.

    The single-bone model has no forefoot, so the plantar fascia cannot be a
    two-ended ``calcaneus<->metatarsal`` connection; it is modelled as a
    **tension load on its calcaneal attachment region** (see
    :func:`build_thums_feb`'s ``fascia_n``).

    Rule (deterministic, geometry-only): take the plantar surface triangles and
    keep those whose centroid lies in the **posterior** ``x_frac`` band, the
    **inferior** ``y_frac`` band and on the **medial** ``z_frac`` side of the
    plantar bounding box -- i.e. the posteromedial tuberosity.  Fallbacks (in
    order, recorded in ``rule``): drop the medial restriction, then relax the
    x/y fractions to 0.5 / 0.65 / 0.8, so a patch always exists on a valid mesh.

    Returns ``tris`` (K,3 local node ids), ``nodes`` (sorted unique), the patch
    ``area_mm2`` / ``centroid_mm`` / bbox, the thresholds used and the ``rule``.
    """
    nodes = np.asarray(mesh["nodes"], dtype=np.float64)
    tris = np.asarray(mesh["surfaces"]["plantar"], dtype=np.int64).reshape(-1, 3)
    if len(tris) == 0:
        raise ValueError("跖面（plantar）无三角片，无法定位跖腱膜附着区")
    p = nodes[tris]
    cen = p.mean(axis=1)
    area = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
    if float(area.sum()) <= 0.0:
        raise ValueError("跖面三角形面积非正，无法定位跖腱膜附着区")

    lo = cen.min(axis=0)
    hi = cen.max(axis=0)
    span = np.maximum(hi - lo, 1e-12)
    thr = lo + np.asarray([x_frac, y_frac, z_frac], dtype=float) * span

    med = cen[:, 2] > thr[2]
    rules: list[tuple[str, np.ndarray]] = [
        ("posteromedial-tuberosity (x&y&medial)", (cen[:, 0] < thr[0]) & (cen[:, 1] < thr[1]) & med),
        ("posterior+inferior (medial dropped)", (cen[:, 0] < thr[0]) & (cen[:, 1] < thr[1])),
    ]
    for f in (0.50, 0.65, 0.80):
        rules.append(
            (
                f"relaxed x/y frac={f:g}",
                (cen[:, 0] < lo[0] + f * span[0]) & (cen[:, 1] < lo[1] + f * span[1]),
            )
        )

    sel: np.ndarray | None = None
    rule = ""
    for name, m in rules:
        if int(np.count_nonzero(m)) >= int(min_tris):
            sel, rule = m, name
            break
    if sel is None:
        raise RuntimeError(
            "跖腱膜附着区选取失败：所有规则都得到 < %d 个三角片" % int(min_tris)
        )

    sub = tris[sel]
    sel_cen = cen[sel]
    node_ids = np.unique(sub)
    return {
        "tris": sub,
        "nodes": node_ids,
        "area_mm2": float(area[sel].sum()),
        "n_tris": int(len(sub)),
        "n_nodes": int(len(node_ids)),
        "centroid_mm": [float(v) for v in sel_cen.mean(axis=0)],
        "bbox_min_mm": [float(v) for v in sel_cen.min(axis=0)],
        "bbox_max_mm": [float(v) for v in sel_cen.max(axis=0)],
        "thresholds_mm": [float(v) for v in thr],
        "plantar_bbox_min_mm": [float(v) for v in lo],
        "plantar_bbox_max_mm": [float(v) for v in hi],
        "medial_is_plus_z": True,
        "rule": rule,
    }


def _classify_surfaces_poly(
    centroids: np.ndarray, normals: np.ndarray
) -> dict[str, np.ndarray]:
    """Classify polygon faces: ``subtalar_joint`` / ``plantar`` / ``achilles``.

    ``joint`` and ``plantar`` use the **exact** ``meshing.py`` thresholds.
    ``achilles`` is calibrated **relative to this mesh** (posterior-superior
    tuberosity) because the absolute ``meshing.ACHILLES_*`` box was calibrated on
    the IITD STL and selects nothing on the THUMS calcaneus.  Takes precomputed
    polygon centroids/normals so a hex quad is classified as one face (not two).
    """
    dist = np.linalg.norm(centroids - np.asarray(meshing.JOINT_CENTER_MM, dtype=float), axis=1)
    joint = (dist < meshing.JOINT_RADIUS_MM) & (normals[:, 1] > meshing.JOINT_NORMAL_Y_MIN)
    plantar = normals[:, 1] < meshing.PLANTAR_NORMAL_Y_MAX
    if len(centroids):
        ax_lo, ax_hi = float(centroids[:, 0].min()), float(centroids[:, 0].max())
        ay_lo, ay_hi = float(centroids[:, 1].min()), float(centroids[:, 1].max())
        ax_thr = ax_lo + ACHILLES_X_FRAC * (ax_hi - ax_lo)
        ay_thr = ay_lo + ACHILLES_Y_FRAC * (ay_hi - ay_lo)
    else:
        ax_thr, ay_thr = meshing.ACHILLES_X_MAX, meshing.ACHILLES_Y_MIN
    achilles = (centroids[:, 0] < ax_thr) & (centroids[:, 1] > ay_thr)
    return {"subtalar_joint": joint, "plantar": plantar, "achilles": achilles}


def _triangulate_polys(
    nodes: np.ndarray, polys: list[tuple[int, ...]], mask: np.ndarray
) -> np.ndarray:
    """Triangulate the selected polygons into ``(K,3)`` triangle node indices."""
    tris: list[tuple[int, int, int]] = []
    for i in np.flatnonzero(mask):
        poly = polys[int(i)]
        a, b, c = poly[0], poly[1], poly[2]
        tris.append((a, b, c))
        if len(poly) == 4:
            tris.append((a, c, poly[3]))
    if not tris:
        return np.empty((0, 3), dtype=np.int64)
    return np.asarray(tris, dtype=np.int64)


def hex_centroids_volumes(nodes: np.ndarray, hexes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Hex8 element centroids ``(N,3)`` and absolute volumes ``(N,)`` (mm)."""
    nodes = np.asarray(nodes, dtype=np.float64)
    hexes = np.asarray(hexes, dtype=np.int64)
    centroids = nodes[hexes].mean(axis=1)
    vol = np.zeros(len(hexes), dtype=np.float64)
    for t in _HEX_TET_FAN:
        q = nodes[hexes[:, t]]
        v6 = np.einsum(
            "ij,ij->i", q[:, 1] - q[:, 0], np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0])
        )
        vol += np.abs(v6) / 6.0
    return centroids, vol


def load_thums_mesh(npz_path: str | Path = DEFAULT_NPZ) -> dict:
    """Read the registered THUMS calcaneus -> two-domain mesh dict (mm).

    Returns a dict with, at least:
        ``nodes`` (1544,3), ``elements_cort`` (634,8) hex8,
        ``elements_spon`` (3323,4) tet4, ``surfaces`` / ``node_sets`` for
        ``subtalar_joint`` / ``plantar`` / ``achilles``, ``nodes_global``,
        ``joint_center_mm``, per-domain element/node counts and volumes.
    """
    npz_path = Path(npz_path)
    if not npz_path.is_file():
        raise FileNotFoundError(f"THUMS 配准网格不存在：{npz_path}")
    d = np.load(npz_path)
    node_ids = np.asarray(d["node_ids"], dtype=np.int64)
    nodes = np.asarray(d["nodes_xyz"], dtype=np.float64)
    conn = np.asarray(d["conn"], dtype=np.int64)
    pid = np.asarray(d["pid"], dtype=np.int64)

    loc = {int(n): i for i, n in enumerate(node_ids)}
    cort_rows: list[list[int]] = []
    spon_rows: list[list[int]] = []
    for row, p in zip(conn, pid):
        pi = int(p)
        if pi == PID_CORT:
            cort_rows.append([loc[int(x)] for x in row])  # full hex8
        elif pi == PID_SPON:
            spon_rows.append([loc[x] for x in _first_distinct(row, 4)])  # tet4
        else:
            raise ValueError(f"未知 pid {pi}；只接受 CORT={PID_CORT} / SPON={PID_SPON}")
    elements_cort = np.asarray(cort_rows, dtype=np.int64)
    elements_spon = np.asarray(spon_rows, dtype=np.int64)

    # All solid elements in a single list: index space used by face owners.
    elements: list[np.ndarray] = [elements_cort[i] for i in range(len(elements_cort))]
    elements += [elements_spon[i] for i in range(len(elements_spon))]
    elem_centroids = np.asarray([nodes[e].mean(axis=0) for e in elements], dtype=np.float64)

    polys, owners = _mixed_boundary_faces(elements)
    centroids, normals = _polygon_centroids_normals(nodes, polys, owners, elem_centroids)
    masks = _classify_surfaces_poly(centroids, normals)

    # meshing.py rule: joint nodes win over plantar; achilles excludes both.
    joint_nodes = np.unique(
        np.concatenate([polys[i] for i in np.flatnonzero(masks["subtalar_joint"])])
    ) if masks["subtalar_joint"].any() else np.empty(0, dtype=np.int64)
    plantar_poly_nodes = (
        np.unique(np.concatenate([polys[i] for i in np.flatnonzero(masks["plantar"])]))
        if masks["plantar"].any()
        else np.empty(0, dtype=np.int64)
    )
    achilles_poly_nodes = (
        np.unique(np.concatenate([polys[i] for i in np.flatnonzero(masks["achilles"])]))
        if masks["achilles"].any()
        else np.empty(0, dtype=np.int64)
    )
    plantar_nodes = np.setdiff1d(plantar_poly_nodes, joint_nodes)
    achilles_nodes = np.setdiff1d(
        achilles_poly_nodes, np.union1d(joint_nodes, plantar_nodes)
    )
    node_sets = {
        "subtalar_joint": joint_nodes,
        "plantar": plantar_nodes,
        "achilles": achilles_nodes,
    }
    surfaces = {
        "subtalar_joint": _triangulate_polys(nodes, polys, masks["subtalar_joint"]),
        "plantar": _triangulate_polys(nodes, polys, masks["plantar"]),
        "achilles": _triangulate_polys(nodes, polys, masks["achilles"]),
    }

    # tet volumes for reporting
    q = nodes[elements_spon]
    tet_vol6 = np.einsum(
        "ij,ij->i", q[:, 1] - q[:, 0], np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0])
    )
    _, cort_vol = hex_centroids_volumes(nodes, elements_cort)

    # full outer shell as triangles (matches the registered STL)
    all_tris = _triangulate_polys(nodes, polys, np.ones(len(polys), dtype=bool))

    return {
        "nodes": nodes,
        "nodes_global": node_ids,
        "elements_cort": elements_cort,
        "elements_spon": elements_spon,
        "surfaces": surfaces,
        "node_sets": node_sets,
        "boundary_faces": all_tris,
        "boundary_polys": polys,
        "boundary_owners": owners,
        "boundary_normals": normals,
        "boundary_centroids": centroids,
        "element_centroids": elem_centroids,
        "joint_center_mm": tuple(float(x) for x in meshing.JOINT_CENTER_MM),
        "n_nodes": int(len(nodes)),
        "n_elements_cort": int(len(elements_cort)),
        "n_elements_spon": int(len(elements_spon)),
        "n_nodes_cort": int(len(np.unique(elements_cort))),
        "n_nodes_spon": int(len(np.unique(elements_spon))),
        "volume_cort_mm3": float(cort_vol.sum()),
        "volume_spon_mm3": float(np.abs(tet_vol6).sum() / 6.0),
        "n_boundary_faces": int(len(polys)),
        "n_boundary_tris": int(len(all_tris)),
        "n_joint_faces": int(masks["subtalar_joint"].sum()),
        "n_plantar_faces": int(masks["plantar"].sum()),
        "n_achilles_faces": int(masks["achilles"].sum()),
        "n_joint_tris": int(len(surfaces["subtalar_joint"])),
        "n_plantar_tris": int(len(surfaces["plantar"])),
        "n_achilles_tris": int(len(surfaces["achilles"])),
        "n_joint_nodes": int(len(joint_nodes)),
        "n_plantar_nodes": int(len(plantar_nodes)),
        "n_achilles_nodes": int(len(achilles_nodes)),
        "source_npz": str(npz_path),
    }


# ---------------------------------------------------------------------------
# .feb writer
# ---------------------------------------------------------------------------
def _gradient_band_loads(
    model,
    fmesh,
    mesh: dict,
    joint_face: str,
    load_n: float,
    coeff: float,
    direction,
    n_bands: int,
) -> list[tuple[str, float]]:
    """Eccentric (linearly varying) joint pressure via banded sub-surfaces.

    The joint surface triangles are projected onto ``direction`` and binned into
    ``n_bands`` bands.  Each facet gets a linear weight ``w = 1 + coeff·t`` with
    ``t ∈ [-1, 1]`` the normalised projection; facet pressures are normalised so
    that ``Σ p_f·A_f == load_n`` **exactly**, then each band gets the
    area-weighted mean pressure.  Because ``Σ_b P_b·A_b == Σ_f p_f·A_f``, the
    banded representation still transmits exactly ``load_n``.

    Returns ``[(surface_name, pressure_MPa), ...]`` and registers each banded
    surface on ``model.mesh_``.
    """
    coeff = float(coeff)
    if not (0.0 <= coeff < 1.0):
        raise ValueError(f"gradient_coeff 必须在 [0, 1)，得到 {coeff!r}")
    n_bands = int(n_bands)
    if n_bands < 2:
        raise ValueError(f"gradient_bands 至少 2，得到 {n_bands}")
    load_n = float(load_n)
    if load_n <= 0.0:
        raise ValueError(f"load_n 必须为正，得到 {load_n!r}")

    nodes = np.asarray(mesh["nodes"], dtype=np.float64)
    tris = np.asarray(mesh["surfaces"][joint_face], dtype=np.int64).reshape(-1, 3)
    if len(tris) == 0:
        raise ValueError(f"关节面 {joint_face!r} 无三角片，无法施加梯度压力")

    g = np.asarray(direction, dtype=float)
    ng = float(np.linalg.norm(g))
    if ng == 0.0:
        raise ValueError("gradient_dir 不能是零向量")
    g = g / ng

    p = nodes[tris]
    centroids = p.mean(axis=1)
    areas = 0.5 * np.linalg.norm(
        np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1
    )
    if float(areas.sum()) <= 0.0:
        raise ValueError("关节面面积非正")
    t = centroids @ g
    tmid = 0.5 * (float(t.min()) + float(t.max()))
    thalf = max(0.5 * (float(t.max()) - float(t.min())), 1e-30)
    tn = (t - tmid) / thalf  # normalised projection in [-1, 1]
    w = 1.0 + coeff * tn  # > 0 because coeff < 1
    p0 = load_n / float((w * areas).sum())
    pf = p0 * w  # facet pressure, Σ pf·A == load_n exactly

    band_idx = np.clip(((tn + 1.0) * 0.5 * n_bands).astype(int), 0, n_bands - 1)
    bands: list[tuple[int, np.ndarray, float, float]] = []
    for bi in range(n_bands):
        sel = np.flatnonzero(band_idx == bi)
        if len(sel) == 0:
            continue
        ab = float(areas[sel].sum())
        pb = float((pf[sel] * areas[sel]).sum()) / ab
        bands.append((bi, sel, ab, pb))

    total = sum(ab * pb for _, _, ab, pb in bands)
    if abs(total - load_n) > 1e-6 * max(1.0, abs(load_n)):
        raise RuntimeError(
            f"梯度面压合力未守恒：Σ P·A={total:.6f} N vs 目标 {load_n:.6f} N"
        )

    out: list[tuple[str, float]] = []
    for bi, sel, _ab, pb in bands:
        name = f"{joint_face}_g{bi}"
        surf = fmesh.Surface(name=name)
        for i, fi in enumerate(sel):
            tri = tris[int(fi)]
            surf.add_tri3(fmesh.Tri3Element(id=i + 1, text=",".join(map(str, tri + 1))))
        model.mesh_.add_surface(surf)
        out.append((name, pb))
    return out


def build_thums_feb(
    mesh: dict,
    out_feb,
    *,
    E_cort_mpa: float = 15000.0,
    nu_cort: float = 0.3,
    E_spon_mpa: float = 73.4,
    nu_spon: float = 0.45,
    rho_cort_tonne_mm3: float = 1.0e-6,
    rho_spon_tonne_mm3: float = 1.0e-6,
    load_n: float = 4000.0,
    load_dir=(0.0, -1.0, 0.0),
    joint_face: str | None = None,
    use_rigid: bool = True,
    joint_load: str = "pressure",
    gradient_coeff: float = 0.5,
    gradient_dir=(1.0, 0.0, 0.0),
    gradient_bands: int = 6,
    analysis: str = "STATIC",
    time_steps: int = 1,
    step_size: float | None = None,
    ramp_time_s: float | None = None,
    achilles_n: float = 0.0,
    achilles_dir=(0.15, 0.985, 0.0),
    achilles_surface: str = "achilles",
    plantar_bc: str = "fixed",
    spring_k: float = 1000.0,
    subtalar_force=None,
    subtalar_moment=None,
    moment_bands: int = 6,
    fascia_n: float = 0.0,
    fascia_dir=(1.0, -0.1, 0.0),
    fascia_surface: str = "fascia",
    fascia_x_frac: float = FASCIA_X_FRAC,
    fascia_y_frac: float = FASCIA_Y_FRAC,
    fascia_z_frac: float = FASCIA_Z_FRAC,
) -> Path:
    """Write the two-domain THUMS calcaneus ``.feb`` (mm-N-MPa-s).

    The load/BC conventions mirror ``febio_model.build_calcaneus_feb`` exactly
    (ghost rigid body + ``RigidForceLoad`` with ``PressureLoad`` fallback,
    optional Achilles ``TractionLoad``, single STATIC step, FEBio-4.13
    ``<solver/>`` adaptation).  The only structural difference is the two
    element domains ``"calcaneus"`` (hex8/CORT) and ``"trabecular"``
    (tet4/SPON) bound to distinct materials.

    ``plantar_bc`` selects the plantar boundary condition (default ``"fixed"``
    preserves the historical tri-axial fixity bit-for-bit):

    * ``"fixed"``  : ``BCZeroDisplacement(plantar, x=y=z=1)`` (original);
    * ``"roller"`` : normal-only zero displacement + 2 lateral pins;
    * ``"spring"`` : 3-direction Winkler foundation with stiffness ``spring_k``
                     (N/mm) per plantar node.

    See :mod:`climbing.coupling.plantar_bc`.

    Phase-S4 opt-in 3D wrench (defaults ``None`` reproduce the historical
    scalar ``load_n × load_dir`` behaviour bit-for-bit):

    * ``subtalar_force`` — full 3D force (N,3).  Rigid path applies it per DOF
      via ``RigidForceLoad``; pressure path decomposes it into a **normal**
      ``PressureLoad`` plus a **tangential** ``TractionLoad`` (so shear is
      carried faithfully), overriding ``joint_load``.
    * ``subtalar_moment`` — full 3D moment (N·m,3).  Rigid path uses
      ``RigidMomentLoad`` (values written in N·mm, i.e. ×1000); pressure path
      expresses it as a joint-face **linear pressure gradient** ``p(x)=p0+g·x``
      via :func:`climbing.coupling.febio_model.linear_pressure_bands`
      (banded constant pressure, net force conserved).  A normal pressure field
      cannot carry the moment component along the face normal (torsion); that
      component is dropped and recorded.
    * ``moment_bands`` — number of gradient bands (default 6).

    Joint-load representation (only used when ``use_rigid=False``; defaults
    reproduce the historical uniform ``PressureLoad p=F/A`` bit-for-bit):

    * ``joint_load="pressure"`` — single uniform ``PressureLoad`` (default);
    * ``joint_load="gradient"`` — the joint surface is split into
      ``gradient_bands`` bands along ``gradient_dir`` and each band gets a
      constant pressure, so the pressure varies linearly (eccentric) across the
      face; ``gradient_coeff`` in ``[0, 1)`` is the normalised tilt
      (``1 ± coeff``).  Band pressures are normalised so the **total force is
      exactly ``load_n``**;
    * ``joint_load="nodal"`` — one ``NodalForce`` per joint node (equal share
      ``load_n / N`` along ``load_dir``), total exactly ``load_n``.

    ``analysis`` selects ``"STATIC"`` (default) or ``"DYNAMIC"``.  Dynamic runs
    must carry a physical mass: pass ``rho_*_tonne_mm3`` (cortical bone
    ≈ ``1.8e-9``, trabecular ≈ ``0.6e-9``).  When ``ramp_time_s`` is given the
    joint loads reference a ``LoadCurve`` that ramps ``0→1`` linearly over
    ``ramp_time_s`` (then holds constant).  The historical default density
    (``1.0e-6``) is kept so all existing callers write identical files.
    """
    from pyfebio.model import Model
    from pyfebio import (
        boundary as fbc,
        control,
        loaddata as fld,
        loads as floads,
        material,
        mesh as fmesh,
        meshdomains,
        output as fout,
        rigid,
        step as fstep,
    )

    joint_face = joint_face or DEFAULT_JOINT_FACE
    out_feb = Path(out_feb)
    out_feb.parent.mkdir(parents=True, exist_ok=True)

    if joint_load not in ("pressure", "gradient", "nodal"):
        raise ValueError(
            f"未知 joint_load={joint_load!r}；可选 'pressure' / 'gradient' / 'nodal'"
        )
    if analysis not in ("STATIC", "DYNAMIC"):
        raise ValueError(f"未知 analysis={analysis!r}；可选 'STATIC' / 'DYNAMIC'")
    if not (float(rho_cort_tonne_mm3) > 0.0 and float(rho_spon_tonne_mm3) > 0.0):
        raise ValueError(
            f"密度必须为正 (tonne/mm³)：cort={rho_cort_tonne_mm3!r} spon={rho_spon_tonne_mm3!r}"
        )
    if analysis == "DYNAMIC" and step_size is None:
        raise ValueError("analysis='DYNAMIC' 必须显式给 step_size (s)")
    if analysis == "DYNAMIC" and (ramp_time_s is None or float(ramp_time_s) <= 0.0):
        raise ValueError("analysis='DYNAMIC' 需要正的 ramp_time_s (s)")

    nodes = np.asarray(mesh["nodes"], dtype=np.float64)
    cort = np.asarray(mesh["elements_cort"], dtype=np.int64)
    spon = np.asarray(mesh["elements_spon"], dtype=np.int64)
    if joint_face not in mesh.get("surfaces", {}):
        raise KeyError(f"mesh 里没有具名面 {joint_face!r}；有：{list(mesh.get('surfaces', {}))}")
    if "plantar" not in mesh.get("node_sets", {}):
        raise KeyError("mesh 里没有 'plantar' 节点集（跖面约束所需）")

    d = np.asarray(load_dir, dtype=float)
    nrm = float(np.linalg.norm(d))
    if nrm == 0.0:
        raise ValueError("load_dir 不能是零向量")
    d = d / nrm

    # --- Phase-S4：可选 3D wrench（默认 None → 旧标量口径逐位不变）----------
    force_vec: np.ndarray | None = None
    if subtalar_force is not None:
        force_vec = np.asarray(subtalar_force, dtype=float).ravel()
        if force_vec.shape != (3,):
            raise ValueError(f"subtalar_force 必须是 3 分量 (N)，得到 {subtalar_force!r}")
        if not np.all(np.isfinite(force_vec)):
            raise ValueError(f"subtalar_force 含 NaN/Inf：{subtalar_force!r}")
    moment_vec: np.ndarray | None = None
    if subtalar_moment is not None:
        moment_vec = np.asarray(subtalar_moment, dtype=float).ravel()
        if moment_vec.shape != (3,):
            raise ValueError(f"subtalar_moment 必须是 3 分量 (N·m)，得到 {subtalar_moment!r}")
        if not np.all(np.isfinite(moment_vec)):
            raise ValueError(f"subtalar_moment 含 NaN/Inf：{subtalar_moment!r}")
    wrench = force_vec is not None or moment_vec is not None

    model = Model()
    # --- Mesh: ONE node domain, TWO element domains ----------------------
    model.mesh_.add_node_domain(fmesh.numpy_to_nodes(nodes, name="calcaneus"))
    model.mesh_.add_element_domain(
        fmesh.numpy_to_elements(cort + 1, "hex8", name="calcaneus")
    )
    model.mesh_.add_element_domain(
        fmesh.numpy_to_elements(spon + 1, "tet4", name="trabecular")
    )

    joint_nodes = np.asarray(mesh["node_sets"][joint_face], dtype=np.int64)
    plantar_nodes = np.asarray(mesh["node_sets"]["plantar"], dtype=np.int64)
    if len(joint_nodes) == 0:
        raise ValueError("关节节点集为空")
    if len(plantar_nodes) == 0:
        raise ValueError("跖面节点集为空")
    model.mesh_.add_node_set(
        fmesh.NodeSet(name=joint_face, text=",".join(map(str, joint_nodes + 1)))
    )
    model.mesh_.add_node_set(
        fmesh.NodeSet(name="plantar", text=",".join(map(str, plantar_nodes + 1)))
    )

    # Named surfaces (tri3) for PressureLoad fallback + traceability.
    surf = fmesh.Surface(name=joint_face)
    for i, tri in enumerate(np.asarray(mesh["surfaces"][joint_face], dtype=np.int64)):
        surf.add_tri3(fmesh.Tri3Element(id=i + 1, text=",".join(map(str, tri + 1))))
    model.mesh_.add_surface(surf)

    has_achilles = (
        float(achilles_n) > 0.0
        and achilles_surface in mesh.get("surfaces", {})
        and len(np.asarray(mesh["surfaces"][achilles_surface])) > 0
    )
    # Only emit the Achilles surface/node_set when it is non-empty: an empty
    # <Surface/> makes pyfebio's xplt parser KeyError on 'faces'.
    ach_faces = np.asarray(mesh.get("surfaces", {}).get(achilles_surface, []), dtype=np.int64).reshape(-1, 3)
    if len(ach_faces) > 0:
        ach_nodes = np.asarray(mesh["node_sets"][achilles_surface], dtype=np.int64)
        if len(ach_nodes) > 0:
            model.mesh_.add_node_set(
                fmesh.NodeSet(name=achilles_surface, text=",".join(map(str, ach_nodes + 1)))
            )
        ach_surf = fmesh.Surface(name=achilles_surface)
        for i, tri in enumerate(ach_faces):
            ach_surf.add_tri3(fmesh.Tri3Element(id=i + 1, text=",".join(map(str, tri + 1))))
        model.mesh_.add_surface(ach_surf)

    # --- Materials: cortical + trabecular --------------------------------
    model.material_.add_material(
        material.IsotropicElastic(
            name="bone_cort",
            E=material.MaterialParameter(text=float(E_cort_mpa)),
            v=material.MaterialParameter(text=float(nu_cort)),
            density=material.MaterialParameter(text=float(rho_cort_tonne_mm3)),
        )
    )
    model.material_.add_material(
        material.IsotropicElastic(
            name="bone_spon",
            E=material.MaterialParameter(text=float(E_spon_mpa)),
            v=material.MaterialParameter(text=float(nu_spon)),
            density=material.MaterialParameter(text=float(rho_spon_tonne_mm3)),
        )
    )
    model.meshdomains_.add_solid_domain(
        meshdomains.SolidDomain(name="calcaneus", mat="bone_cort")
    )
    model.meshdomains_.add_solid_domain(
        meshdomains.SolidDomain(name="trabecular", mat="bone_spon")
    )

    # --- Step -------------------------------------------------------------
    if step_size is None:
        step_size = 1.0 / max(int(time_steps), 1)
    st = fstep.StepEntry(id=1, name="Step")
    st.control = control.Control(
        analysis=analysis,
        time_steps=int(time_steps),
        step_size=float(step_size),
        time_stepper=None,  # FEBio 4.13 拒绝 <time_stepper type=...>
    )

    # Dynamic runs ramp the joint load over ramp_time_s via a load controller
    # (id=1).  Load objects reference it with lc=1; static runs use lc=None so
    # the written XML is bit-identical to the historical constant-pressure form.
    lc_ref: int | None = None
    if analysis == "DYNAMIC":
        ramp = float(ramp_time_s)
        lc_ref = 1
        model.loaddata_.add_load_curve(
            fld.LoadCurve(
                id=lc_ref,
                interpolate="LINEAR",
                extend="CONSTANT",
                points=fld.CurvePoints(points=["0.0,0.0", f"{ramp:.10g},1.0"]),
            )
        )

    bcs: list = []
    plantar_meta = apply_plantar_bc(
        model,
        nodes=nodes,
        plantar_local=plantar_nodes,
        plantar_tris=np.asarray(
            mesh.get("surfaces", {}).get("plantar", []), dtype=np.int64
        ).reshape(-1, 3),
        bcs=bcs,
        kind=plantar_bc,
        spring_k=spring_k,
    )

    surface_loads: list = []
    nodal_loads: list = []
    gradient_meta: dict | None = None
    if use_rigid:
        joint_center = tuple(float(x) for x in mesh.get("joint_center_mm", (0.0, 0.0, 0.0)))
        model.add_simple_rigid_body(joint_center, "ghost")
        bcs.append(fbc.BCRigid(name="joint_to_ghost", node_set=joint_face, rb="ghost"))
        st.boundary = fbc.Boundary(all_bcs=bcs)
        # 纯力时锁转动；给力矩时放开承载该力矩的转动 DOF（否则 RigidMomentLoad 被吸收）。
        if moment_vec is None:
            rigid_bcs = [rigid.RigidFixed(rb="ghost", Ru_dof=1, Rv_dof=1, Rw_dof=1)]
        else:
            lock = {dof: 1 for i, dof in enumerate(("Ru", "Rv", "Rw")) if moment_vec[i] == 0.0}
            rigid_bcs = [rigid.RigidFixed(rb="ghost", **lock)] if lock else []
        rigid_loads = []
        fv = force_vec if force_vec is not None else load_n * d
        for comp, dof in enumerate(("Rx", "Ry", "Rz")):
            if fv[comp] != 0.0:
                rigid_loads.append(
                    rigid.RigidForceLoad(
                        rb="ghost", dof=dof, value=rigid.Value(text=float(fv[comp]))
                    )
                )
        if moment_vec is not None:
            # RigidMomentLoad.value 用 N·mm（力×长度）；输入是 N·m。
            for comp, dof in enumerate(("Ru", "Rv", "Rw")):
                if moment_vec[comp] != 0.0:
                    rigid_loads.append(
                        rigid.RigidMomentLoad(
                            rb="ghost", dof=dof,
                            value=rigid.Value(text=float(moment_vec[comp] * 1000.0)),
                        )
                    )
        st.rigid = rigid.Rigid(all_rigid_bcs=rigid_bcs, all_rigid_loads=rigid_loads)
        load_desc = (f"幽灵刚体 RigidForceLoad |F|={np.linalg.norm(fv):.1f} N × "
                     f"({fv[0]:.3f},{fv[1]:.3f},{fv[2]:.3f})")
        if moment_vec is not None:
            load_desc += (f" + RigidMomentLoad |M|={np.linalg.norm(moment_vec):.2f} N·m "
                          f"({moment_vec[0]:.3f},{moment_vec[1]:.3f},{moment_vec[2]:.3f})")
    else:
        st.boundary = fbc.Boundary(all_bcs=bcs)
        if wrench:
            frame = joint_face_frame(mesh, joint_face)
            fv = force_vec if force_vec is not None else load_n * d
            f_n = float(fv @ frame["normal_unit"])
            f_t = fv - f_n * frame["normal_unit"]
            p0 = -f_n / frame["area_mm2"]
            if moment_vec is None:
                surface_loads.append(
                    floads.PressureLoad(
                        surface=joint_face,
                        pressure=floads.Scale(lc=lc_ref, text=float(p0)),
                    )
                )
                wrench_meta = None
            else:
                slope, direction, gmeta = moment_to_pressure_gradient(
                    moment_vec, frame, base_uniform_pressure_mpa=p0
                )
                bands, wrench_meta = linear_pressure_bands(
                    model, fmesh, mesh, joint_face, frame,
                    p0_mpa=p0, slope=slope, direction=direction,
                    n_bands=moment_bands, name_prefix=f"{joint_face}_g",
                )
                for name, pb in bands:
                    surface_loads.append(
                        floads.PressureLoad(
                            surface=name,
                            pressure=floads.Scale(lc=lc_ref, text=float(pb)),
                        )
                    )
                wrench_meta["dropped_normal_axis_nm"] = gmeta["dropped_normal_axis_nm"]
                wrench_meta["moment_inplane_nm"] = gmeta["moment_inplane_nm"]
            ft_mag = float(np.linalg.norm(f_t))
            if ft_mag > 1e-12:
                ta = f_t / frame["area_mm2"]
                surface_loads.append(
                    floads.TractionLoad(
                        surface=joint_face,
                        scale=floads.Scale(text="1.0"),
                        traction=f"{ta[0]:.6f},{ta[1]:.6f},{ta[2]:.6f}",
                    )
                )
            load_desc = (
                f"3D wrench: F=({fv[0]:.1f},{fv[1]:.1f},{fv[2]:.1f}) N，法向压力 "
                f"p0={p0:.4f} MPa + 切向牵引 |F_t|={ft_mag:.1f} N"
                + (f" + 力矩梯度 {wrench_meta['n_bands']} 带" if wrench_meta else "")
            )
        elif joint_load == "pressure":
            area = joint_area_mm2(mesh, joint_face)
            pressure = float(load_n) / area  # N / mm² = MPa
            surface_loads.append(
                floads.PressureLoad(
                    surface=joint_face,
                    pressure=floads.Scale(lc=lc_ref, text=pressure),
                )
            )
            load_desc = (
                f"PressureLoad p={pressure:.4f} MPa（F={load_n:.1f} N / A={area:.2f} mm²）"
            )
        elif joint_load == "gradient":
            bands = _gradient_band_loads(
                model,
                fmesh,
                mesh,
                joint_face,
                float(load_n),
                gradient_coeff,
                gradient_dir,
                int(gradient_bands),
            )
            for name, pb in bands:
                surface_loads.append(
                    floads.PressureLoad(
                        surface=name,
                        pressure=floads.Scale(lc=lc_ref, text=float(pb)),
                    )
                )
            gradient_meta = {
                "coeff": float(gradient_coeff),
                "dir": [float(x) for x in gradient_dir],
                "n_bands": int(len(bands)),
                "bands": [{"surface": nm, "pressure_mpa": float(pb)} for nm, pb in bands],
            }
            pmin = min(pb for _, pb in bands)
            pmax = max(pb for _, pb in bands)
            load_desc = (
                f"梯度面压 {len(bands)} 带 pressure=[{pmin:.4f},{pmax:.4f}] MPa "
                f"（coeff={gradient_coeff:g}，dir=({gradient_dir[0]:g},{gradient_dir[1]:g},"
                f"{gradient_dir[2]:g})，ΣF={load_n:.1f} N 精确）"
            )
        else:  # nodal
            jn = np.asarray(mesh["node_sets"][joint_face], dtype=np.int64)
            n_j = int(len(jn))
            if n_j == 0:
                raise ValueError("关节节点集为空，无法施加节点力")
            f_node = float(load_n) / n_j
            vec = f_node * d
            for k, li in enumerate(jn):
                nm = f"{joint_face}_n{k:04d}"
                model.mesh_.add_node_set(
                    fmesh.NodeSet(name=nm, text=str(int(li) + 1))
                )
                nodal_loads.append(
                    floads.NodalForce(
                        node_set=nm,
                        value=floads.Scale(
                            lc=lc_ref,
                            text=f"{vec[0]:.10g},{vec[1]:.10g},{vec[2]:.10g}",
                        ),
                    )
                )
            load_desc = (
                f"NodalForce 逐节点力 {n_j} × {f_node:.3f} N × "
                f"({d[0]:.3f},{d[1]:.3f},{d[2]:.3f}) = {load_n:.1f} N"
            )

    if has_achilles:
        ad = np.asarray(achilles_dir, dtype=float)
        nn = float(np.linalg.norm(ad))
        if nn == 0.0:
            raise ValueError("achilles_dir 不能是零向量")
        ad = ad / nn
        a_area = joint_area_mm2(mesh, achilles_surface)
        trac = float(achilles_n) / a_area  # MPa
        vec = trac * ad
        surface_loads.append(
            floads.TractionLoad(
                surface=achilles_surface,
                scale=floads.Scale(text="1.0"),
                traction=f"{vec[0]:.6f},{vec[1]:.6f},{vec[2]:.6f}",
            )
        )
        load_desc += (
            f" + 跟腱牵引 F={achilles_n:.1f} N / A={a_area:.1f} mm² = {trac:.4f} MPa"
            f" × ({ad[0]:.2f},{ad[1]:.2f},{ad[2]:.2f})"
        )

    # 跖腱膜（plantar fascia）：单骨无前足，建模为其跟骨**跖腱膜附着区**
    # （后内侧结节）上的张力载荷 TractionLoad（traction = F/A × 单位方向，MPa）。
    # 大小 `fascia_n` 是**参数锚定**（无多体来源）；默认 0 表示不激活，
    # 此时不写任何额外 node_set/surface/load，输出与历史逐位一致。
    if float(fascia_n) > 0.0:
        fas = select_fascia_attachment(
            mesh, x_frac=fascia_x_frac, y_frac=fascia_y_frac, z_frac=fascia_z_frac
        )
        fd = np.asarray(fascia_dir, dtype=float)
        nf = float(np.linalg.norm(fd))
        if nf == 0.0:
            raise ValueError("fascia_dir 不能是零向量")
        fd = fd / nf
        model.mesh_.add_node_set(
            fmesh.NodeSet(name=fascia_surface, text=",".join(map(str, fas["nodes"] + 1)))
        )
        fas_surf = fmesh.Surface(name=fascia_surface)
        for i, tri in enumerate(np.asarray(fas["tris"], dtype=np.int64)):
            fas_surf.add_tri3(
                fmesh.Tri3Element(id=i + 1, text=",".join(map(str, tri + 1)))
            )
        model.mesh_.add_surface(fas_surf)
        f_area = float(fas["area_mm2"])
        f_trac = float(fascia_n) / f_area  # MPa (N / mm²)
        f_vec = f_trac * fd
        surface_loads.append(
            floads.TractionLoad(
                surface=fascia_surface,
                scale=floads.Scale(text="1.0"),
                traction=f"{f_vec[0]:.6f},{f_vec[1]:.6f},{f_vec[2]:.6f}",
            )
        )
        load_desc += (
            f" + 跖腱膜牵引 F={float(fascia_n):.1f} N / A={f_area:.1f} mm² = {f_trac:.4f} MPa"
            f" × ({fd[0]:.3f},{fd[1]:.3f},{fd[2]:.3f})"
            f" [附着区 {fas['n_tris']} 面/{fas['n_nodes']} 节点, {fas['rule']}]"
        )

    if surface_loads or nodal_loads:
        st.loads = floads.Loads(
            all_surface_loads=surface_loads,
            all_nodal_loads=nodal_loads,
        )

    st.control.plot_level = "PLOT_MAJOR_ITRS"
    model.step_.add_step(st)
    model.output_.add_plotfile(
        fout.OutputPlotfile(
            type="febio",
            file=out_feb.stem + ".xplt",
            all_vars=[fout.Var(type="displacement"), fout.Var(type="stress")],
        )
    )

    model.save(out_feb)
    notes = _fix_febio413(out_feb)
    for n in notes:
        print(f"[thums_feb] FEBio 4.13 适配：{n}")
    print(f"[thums_feb] 载荷：{load_desc}")
    print(f"[thums_feb] 跖面 BC：{plantar_meta}")
    return out_feb


# ---------------------------------------------------------------------------
# mixed-domain regularization (gauge) for the CORT hex cells
# ---------------------------------------------------------------------------
def gauge_von_mises_mixed(
    vm: np.ndarray,
    centroids: np.ndarray,
    volumes: np.ndarray,
    radius_mm: float,
) -> np.ndarray:
    """Per-element volume-weighted average von Mises within ``radius_mm``.

    Mirrors ``fe_post.gauge_von_mises`` but works for an arbitrary element set
    (the CORT hex8 cells), where element centroids/volumes are supplied.
    """
    from scipy.spatial import cKDTree

    vm = np.asarray(vm, dtype=np.float64)
    centroids = np.asarray(centroids, dtype=np.float64)
    volumes = np.asarray(volumes, dtype=np.float64)
    if not (len(vm) == len(centroids) == len(volumes)):
        raise ValueError(
            f"长度不匹配：vm={len(vm)} centroids={len(centroids)} volumes={len(volumes)}"
        )
    tree = cKDTree(centroids)
    nb = tree.query_ball_point(centroids, float(radius_mm))
    out = np.empty(len(vm), dtype=np.float64)
    for i, idx in enumerate(nb):
        w = volumes[idx]
        out[i] = float((vm[idx] * w).sum() / w.sum())
    return out
