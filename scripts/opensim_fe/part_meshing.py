"""part_meshing.py -- mm-native bone STL -> gmsh tet with *working* size control.

Why this exists
---------------
``meshing._gmsh_tet_remesh`` calls ``gmsh.model.mesh.classifySurfaces(angle)``
with the default ``forReparametrization=False``.  That builds a geometry whose
input STL vertices remain mandatory nodes, so ``generate(2)`` can never coarsen
the surface and ``Mesh.MeshSizeMax`` has no effect (char_len 3 -> 10 mm leaves
the tet count nearly constant).  Probe evidence::

    tibia iso-dec  baseline       cl=3 470605 -> cl=10 362545 tets   (flat)

The fix (validated by ``gmsh_size_probe.py``) is to ask gmsh for a
*reparametrizable* geometry and to switch off the point/boundary/curvature size
sources (the canonical gmsh t10 recipe)::

    classifySurfaces(angle, boundary=True, forReparametrization=True)
    Mesh.MeshSizeFromPoints        = 0
    Mesh.MeshSizeFromCurvature     = 0
    Mesh.MeshSizeExtendFromBoundary= 0
    Mesh.MeshSizeMin = Mesh.MeshSizeMax = char_len

which yields a true h-family::

    tibia iso-dec  reparam        cl=3  63115 -> cl=10  2547 tets  (monotone)

This module is **new and mm-native**: it does not touch ``meshing.py`` and does
not go through ``stl_to_tet_gmsh`` (which assumes metre input,  x1000).  It
imports only the private helpers ``_isosurface`` and ``_finalize``.

Surface routing
---------------
* ``raw``      used when the cleaned surface is a closed 2-manifold
               (0 open edges and 0 non-manifold edges) -- reparametrisation is
               only safe/robust on manifold input; a non-manifold surface can
               hang gmsh.
* ``iso-dec``  otherwise: marching-cubes isosurface -> decimate -> clean, then
               reparametrise.  This repairs the self-intersecting THUMS union
               surfaces so gmsh can build a solid at all.

CLI (repo root, PYTHONPATH=src)::

    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\part_meshing.py `
        --name tibia_r `
        --stl temp\\opensim_fe\\parts\\tibia_r\\tibia_r_surface.stl `
        --char-len 4 --ref-volume-mm3 346130.49
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from climbing.coupling import meshing as M  # noqa: E402

SDK: Path = M.SDK
LIB: Path = M.LIB
OUT_BASE = ROOT / "temp" / "opensim_fe" / "parts_gmsh"

#: classify angles tried for reparametrisation (radians).  Large angles merge
#: near-coplanar facets into reparametrizable patches; the fallback protects
#: awkward geometries.  0.001 is the legacy value kept as a last resort.
ANGLES_REPARAM = (math.pi / 2, 5 * math.pi / 12, math.pi / 3, math.pi / 4,
                  math.pi / 6, 0.1, 0.001)

__all__ = ["mesh_bone_stl_mm", "gmsh_tet_reparam", "ANGLES_REPARAM"]


# ---------------------------------------------------------------------------
# surface helpers
# ---------------------------------------------------------------------------
def load_stl_mm(stl_path: Path):
    """Read STL as-is (already mm), clean / triangulate / fill small holes."""
    import pyvista as pv

    surf = pv.read(str(stl_path))
    surf = surf.clean(point_merging=True, tolerance=1e-6, absolute=True).triangulate()
    if surf.n_open_edges:
        surf = surf.fill_holes(1e9)
        surf = surf.clean(point_merging=True, tolerance=1e-6, absolute=True)
    return surf


def manifold_stats(surf) -> tuple[int, int]:
    """(open_edges, non_manifold_edges) from the triangle-edge count table."""
    f = np.asarray(surf.faces, dtype=np.int64)
    tri = f if (f.ndim == 2 and f.shape[1] == 3) else f.reshape(-1, 4)[:, 1:]
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    e = np.sort(e, axis=1)
    _, counts = np.unique(e, axis=0, return_counts=True)
    return int((counts == 1).sum()), int((counts > 2).sum())


def repair_surface(surf, reduction: float, spacing: float):
    """marching-cubes isosurface -> decimate -> clean (mm)."""
    iso = M._isosurface(surf, spacing=max(1.0, spacing))
    dec = iso.decimate_pro(
        reduction=float(reduction), splitting=True, feature_angle=30.0,
        preserve_topology=True, boundary_vertex_deletion=False,
    ).triangulate().clean(point_merging=True, tolerance=1e-6, absolute=True)
    if dec.n_open_edges:
        dec = dec.fill_holes(1e9)
        dec = dec.clean(point_merging=True, tolerance=1e-6, absolute=True)
    return dec


# ---------------------------------------------------------------------------
# gmsh
# ---------------------------------------------------------------------------
def _import_gmsh():
    if not LIB.is_dir():
        raise RuntimeError(f"gmsh SDK lib not found: {LIB}")
    os.environ["PATH"] = f"{LIB};{SDK / 'bin'};{os.environ['PATH']}"
    if str(LIB) not in sys.path:
        sys.path.insert(0, str(LIB))
    import gmsh  # noqa: E402

    return gmsh


def gmsh_tet_reparam(surf, char_len: float, angles: tuple = ANGLES_REPARAM) -> dict:
    """STL -> reparametrised geo solid -> tet.  Returns nodes/tets + diagnostics.

    Raises ``RuntimeError`` with the verbatim gmsh error on failure.
    """
    gmsh = _import_gmsh()
    tmp_stl = Path(tempfile.gettempdir()) / f"_part_meshing_{os.getpid()}.stl"
    surf.save(str(tmp_stl))

    gmsh.initialize(str(LIB))
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.merge(str(tmp_stl))
        gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
        gmsh.model.mesh.removeDuplicateNodes()

        last: Exception | None = None
        used_angle: float | None = None
        for ang in angles:
            try:
                gmsh.model.mesh.classifySurfaces(ang, True, True)
                gmsh.model.mesh.createGeometry()
                used_angle = ang
                break
            except Exception as exc:  # gmsh uses bare Exception
                last = exc
                gmsh.model.mesh.clear()
                gmsh.merge(str(tmp_stl))
        if used_angle is None:
            raise RuntimeError(f"classifySurfaces(forReparametrization=True) failed: {last}")

        gmsh.model.geo.synchronize()
        surfs = [t for _, t in gmsh.model.getEntities(2)]
        if not surfs:
            raise RuntimeError("createGeometry produced no surfaces")
        gmsh.model.geo.addVolume([gmsh.model.geo.addSurfaceLoop(surfs)])
        gmsh.model.geo.synchronize()
        if not gmsh.model.getEntities(3):
            raise RuntimeError("surface loop did not close into a volume")

        # canonical t10: fully field-driven size -> MeshSizeMin/Max actually act
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", float(char_len))
        gmsh.option.setNumber("Mesh.MeshSizeMax", float(char_len))
        gmsh.option.setNumber("Mesh.Algorithm", 6)     # Frontal-Delaunay (surface)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)   # Delaunay (volume)

        gmsh.model.mesh.generate(2)
        n_surf_tri = len(gmsh.model.mesh.getElementsByType(2)[1])
        gmsh.model.mesh.generate(3)

        tags, coords, _ = gmsh.model.mesh.getNodes()
        tag_arr = np.asarray(tags, dtype=np.int64)
        nodes = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
        _, nidx = gmsh.model.mesh.getElementsByType(4)
        if nidx is None or len(nidx) == 0:
            raise RuntimeError("generate(3) produced no tets")
        lut = np.zeros(int(tag_arr.max()) + 1, dtype=np.int64)
        lut[tag_arr] = np.arange(len(tag_arr), dtype=np.int64)
        tets = lut[np.asarray(nidx, dtype=np.int64)].reshape(-1, 4)
        return {
            "nodes": nodes,
            "tets": tets,
            "classify_angle_deg": round(math.degrees(used_angle), 3),
            "n_surface_tri": int(n_surf_tri),
        }
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# public
# ---------------------------------------------------------------------------
def mesh_bone_stl_mm(
    stl_path,
    char_len: float,
    *,
    name: str | None = None,
    reduction: float = 0.70,
    ref_volume_mm3: float | None = None,
    out_dir=None,
) -> dict:
    """mm-native bone STL -> tet with the working reparametrised gmsh path.

    ``name`` / ``out_dir``: when ``name`` is given, writes
    ``<out_dir or OUT_BASE>/<name>/`` containing the mesh npz, the surface used
    and a summary json.  Returns the same metrics as ``meshing._finalize``.
    """
    stl_path = Path(stl_path)
    if not stl_path.is_file():
        raise FileNotFoundError(stl_path)
    raw = load_stl_mm(stl_path)
    n_open, n_nonman = manifold_stats(raw)

    t0 = time.time()
    if n_open == 0 and n_nonman == 0:
        surf = raw
        surface_used = "raw"
    else:
        surf = repair_surface(raw, reduction=reduction, spacing=max(1.0, char_len * 0.5))
        surface_used = "iso-dec"
    n_open2, n_nonman2 = manifold_stats(surf)
    print(
        f"[part_meshing] {stl_path.name}: raw open={n_open} nonmanifold={n_nonman} "
        f"-> use {surface_used} (open={n_open2} nonmanifold={n_nonman2}, "
        f"faces={surf.n_faces})",
        flush=True,
    )

    g = gmsh_tet_reparam(surf, char_len)
    nodes, tets = g["nodes"], g["tets"]

    result = M._finalize(
        surf, nodes, tets,
        method=f"gmsh-reparam({surface_used})",
        joint_center=(1e9, 1e9, 1e9), char_len=char_len,
    )
    ref = float(ref_volume_mm3) if ref_volume_mm3 else float(result["reference_volume_mm3"])
    result["reference_volume_mm3"] = ref
    result["volume_error_pct"] = abs(float(result["volume_mm3"]) - ref) / ref * 100.0 if ref else float("nan")
    result["surface_used"] = surface_used
    result["raw_open_edges"] = n_open
    result["raw_nonmanifold_edges"] = n_nonman
    result["classify_angle_deg"] = g["classify_angle_deg"]
    result["n_surface_tri"] = g["n_surface_tri"]
    result["n_faces_input"] = int(raw.n_faces)
    result["reduction"] = float(reduction) if surface_used == "iso-dec" else 0.0
    result["mesh_time_s"] = round(time.time() - t0, 2)

    if name:
        od = Path(out_dir) if out_dir else (OUT_BASE / name)
        od.mkdir(parents=True, exist_ok=True)
        np.savez(od / f"{name}_mesh.npz", nodes_xyz=nodes, tets=tets)
        try:
            surf.save(str(od / f"{name}_surface_used.stl"))
        except Exception as exc:  # non-fatal: mesh npz is the product
            print(f"[part_meshing] surface save failed ({type(exc).__name__}: {exc})")
        summary = _json_scalars(result)
        (od / f"{name}_summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        result["out_dir"] = str(od)
    return result


def _json_scalars(d: dict) -> dict:
    """Keep only JSON-serialisable scalar/tuple values (drop arrays/dicts)."""
    out = {}
    for k, v in d.items():
        if isinstance(v, np.generic):
            v = v.item()
        if isinstance(v, (str, bool, int, float)) or v is None:
            out[k] = v
        elif isinstance(v, (list, tuple)) and all(
            isinstance(x, (str, bool, int, float)) for x in v
        ):
            out[k] = list(v)
    return out


def _main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--stl", required=True)
    ap.add_argument("--char-len", type=float, default=4.0, dest="char_len")
    ap.add_argument("--reduction", type=float, default=0.70)
    ap.add_argument("--ref-volume-mm3", type=float, default=None, dest="ref_volume_mm3")
    ap.add_argument("--out-dir", default=None, dest="out_dir")
    args = ap.parse_args()

    r = mesh_bone_stl_mm(
        args.stl, args.char_len, name=args.name, reduction=args.reduction,
        ref_volume_mm3=args.ref_volume_mm3, out_dir=args.out_dir,
    )
    keys = ("method", "surface_used", "classify_angle_deg", "n_tets", "n_nodes",
            "volume_mm3", "reference_volume_mm3", "volume_error_pct",
            "min_jacobian", "jacobian_frac_lt_0_3", "n_surface_tri", "mesh_time_s")
    print("[part_meshing] " + json.dumps({k: r.get(k) for k in keys}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
