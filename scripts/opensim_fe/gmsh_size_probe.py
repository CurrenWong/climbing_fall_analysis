"""gmsh_size_probe.py -- diagnose & fix gmsh size control on bone STLs.

Motivation
----------
``meshing._gmsh_tet_remesh`` sets only ``Mesh.MeshSizeMin/Max``.  Empirically
raising those has no effect (``char_len`` 5 -> 12 mm leaves the tet count
invariant), so every real bone overshoots the 50k-tet HDF5 ceiling and the
project is forced onto the voxel fallback (``J<0.3`` ~ 17%).

This script re-implements the *same* pipeline (merge STL -> classifySurfaces ->
createGeometry -> surface loop/volume -> generate(2) -> generate(3)) but makes
the gmsh size options a first-class parameter, and sweeps ``char_len`` for each
configuration.  It is mm-native and read-only w.r.t. ``meshing``: it imports
only private helpers (``_isosurface``, ``_finalize``) and never touches the
existing artifacts.

Configurations probed
---------------------
baseline       exactly ``_gmsh_tet_remesh``'s options (MeshSizeMin=0.4*cl,
               MeshSizeMax=cl).
no-boundary    baseline + FromPoints=0, FromCurvature=0,
               ExtendFromBoundary=0, CharacteristicLengthMin/Max=cl
               (the canonical gmsh t10 recipe); also attempts the
               non-existent option ``Mesh.MeshSizeFromBoundary`` and records
               the error instead of ignoring it.
field-constant a ``Constant`` background size field = cl (+ the t10 switches).
field-matheval a ``MathEval`` background size field "F" = cl (+ t10 switches).
factor         keeps point/boundary-derived sizes, scales them by
               ``Mesh.MeshSizeFactor = cl / mean_input_edge``.

Surfaces
--------
``raw``     the input STL as-is.
``iso-dec`` marching-cubes isosurface -> decimate -> clean (repairs the
            non-manifold/self-intersecting raw bones so gmsh can build a
            solid at all).  This mirrors ``extract_part.mesh_stl`` path 2.

Usage (repo root)::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\gmsh_size_probe.py `
        --stl temp\\opensim_fe\\parts\\tibia_r\\tibia_r_surface.stl `
        --ref-volume-mm3 346130.49 --surfaces raw,iso-dec `
        --chars 3,4,5,6,8,10

Outputs ``temp/opensim_fe/parts_gmsh/gmsh_size_probe_<label>.json`` plus a
console table.  Failures/errors are recorded verbatim (never swallowed).
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

OUT_DIR = ROOT / "temp" / "opensim_fe" / "parts_gmsh"

# re-exported SDK constants from meshing for a single source of truth
SDK: Path = M.SDK
LIB: Path = M.LIB

CONFIG_NAMES = [
    "baseline",
    "reparam-baseline",
    "no-boundary",
    "field-constant",
    "field-matheval",
    "factor",
    "reparam",
]

#: classify-angle fallback lists (radians).  The legacy pipeline uses
#: ``forReparametrization=False``; the fix uses ``True``.
ANGLES_LEGACY = (0.001, 0.01, 0.1, math.pi / 2)
ANGLES_REPARAM = (math.pi / 2, 5 * math.pi / 12, math.pi / 3, math.pi / 4,
                  math.pi / 6, 0.1, 0.001)


# ---------------------------------------------------------------------------
# gmsh bootstrapping (same as meshing._gmsh_tet_remesh)
# ---------------------------------------------------------------------------
def _import_gmsh():
    if not LIB.is_dir():
        raise RuntimeError(f"gmsh SDK lib not found: {LIB}")
    os.environ["PATH"] = f"{LIB};{SDK / 'bin'};{os.environ['PATH']}"
    if str(LIB) not in sys.path:
        sys.path.insert(0, str(LIB))
    import gmsh  # noqa: E402

    return gmsh


def _safe_set(gmsh, log: list, name: str, value: float) -> None:
    """Set an option, recording (not swallowing) any gmsh error."""
    try:
        gmsh.option.setNumber(name, float(value))
    except Exception as exc:  # gmsh raises bare Exception for unknown options
        log.append(f"setNumber({name}={value:g}) -> {type(exc).__name__}: {exc}")


def _build_solid(gmsh, stl_path: Path, angles: tuple, for_reparam: bool) -> float:
    """merge STL -> classifySurfaces (angle fallback) -> geo solid volume.

    Returns the classify angle (rad) that succeeded.  ``for_reparam=True`` asks
    gmsh to build curves/surfaces it can actually reparametrize -- without it
    the input facets stay mandatory nodes and ``generate(2)`` cannot coarsen
    the surface (the root cause of "MeshSizeMax has no effect").

    Raises with the *last* gmsh error if every classify angle fails, or if the
    result cannot be closed into a volume.
    """
    stl_path = Path(stl_path)
    gmsh.merge(str(stl_path))
    gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
    gmsh.model.mesh.removeDuplicateNodes()

    last: Exception | None = None
    for ang in angles:
        try:
            gmsh.model.mesh.classifySurfaces(ang, True, for_reparam)
            gmsh.model.mesh.createGeometry()
            break
        except Exception as exc:  # gmsh uses bare Exception
            last = exc
            gmsh.model.mesh.clear()
            gmsh.merge(str(stl_path))
    else:
        raise RuntimeError(f"classifySurfaces/createGeometry failed at all angles: {last}")

    gmsh.model.geo.synchronize()
    surfs = [t for _, t in gmsh.model.getEntities(2)]
    if not surfs:
        raise RuntimeError("createGeometry produced no surfaces")
    gmsh.model.geo.addVolume([gmsh.model.geo.addSurfaceLoop(surfs)])
    gmsh.model.geo.synchronize()
    if not gmsh.model.getEntities(3):
        raise RuntimeError("surface loop did not close into a volume (open/self-intersecting)")
    return ang


def _mean_edge_mm(surf) -> float:
    """Mean unique triangle-edge length (mm) -- the size the raw points imply."""
    f = np.asarray(surf.faces, dtype=np.int64)
    if f.ndim == 2 and f.shape[1] == 3:
        tri = f
    else:  # packed VTK cell array [3, i, j, k, 3, i, j, k, ...]
        tri = f.reshape(-1, 4)[:, 1:]
    p = np.asarray(surf.points, dtype=np.float64)
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    e = np.sort(e, axis=1)
    e = np.unique(e, axis=0)
    return float(np.linalg.norm(p[e[:, 0]] - p[e[:, 1]], axis=1).mean())


def _apply_config(gmsh, kind: str, char_len: float, mean_edge: float, log: list) -> None:
    # every variant shares the meshing algorithm choices
    _safe_set(gmsh, log, "Mesh.Algorithm", 6)     # Frontal-Delaunay (surface)
    _safe_set(gmsh, log, "Mesh.Algorithm3D", 1)   # Delaunay (volume)

    if kind in ("baseline", "reparam-baseline"):
        _safe_set(gmsh, log, "Mesh.MeshSizeMin", char_len * 0.4)
        _safe_set(gmsh, log, "Mesh.MeshSizeMax", char_len)
        return

    # all non-baseline variants first suppress the point/boundary/curvature size
    # sources that pin the mesh (gmsh t10 lines 114-133).
    for name in (
        "Mesh.MeshSizeFromPoints",
        "Mesh.MeshSizeFromCurvature",
        "Mesh.MeshSizeExtendFromBoundary",
    ):
        _safe_set(gmsh, log, name, 0)
    for alias in (
        "Mesh.CharacteristicLengthFromPoints",
        "Mesh.CharacteristicLengthFromCurvature",
        "Mesh.CharacteristicLengthExtendFromBoundary",
    ):
        _safe_set(gmsh, log, alias, 0)
    # task asks for this option explicitly; it does *not* exist in gmsh 4.15.
    _safe_set(gmsh, log, "Mesh.MeshSizeFromBoundary", 0)

    if kind in ("no-boundary", "reparam"):
        _safe_set(gmsh, log, "Mesh.MeshSizeMin", char_len)
        _safe_set(gmsh, log, "Mesh.MeshSizeMax", char_len)
        _safe_set(gmsh, log, "Mesh.CharacteristicLengthMin", char_len)
        _safe_set(gmsh, log, "Mesh.CharacteristicLengthMax", char_len)
        return

    if kind == "field-constant":
        _safe_set(gmsh, log, "Mesh.MeshSizeMin", char_len)
        _safe_set(gmsh, log, "Mesh.MeshSizeMax", char_len)
        try:
            fid = gmsh.model.mesh.field.add("Constant", 1)
            gmsh.model.mesh.field.setNumber(fid, "VIn", char_len)
            gmsh.model.mesh.field.setNumber(fid, "VOut", char_len)
            gmsh.model.mesh.field.setAsBackgroundMesh(fid)
        except Exception as exc:
            log.append(f"Constant field setup -> {type(exc).__name__}: {exc}")
        return

    if kind == "field-matheval":
        _safe_set(gmsh, log, "Mesh.MeshSizeMin", char_len)
        _safe_set(gmsh, log, "Mesh.MeshSizeMax", char_len)
        try:
            fid = gmsh.model.mesh.field.add("MathEval", 1)
            gmsh.model.mesh.field.setString(fid, "F", str(float(char_len)))
            gmsh.model.mesh.field.setAsBackgroundMesh(fid)
        except Exception as exc:
            log.append(f"MathEval field setup -> {type(exc).__name__}: {exc}")
        return

    if kind == "factor":
        factor = char_len / max(mean_edge, 1e-9)
        _safe_set(gmsh, log, "Mesh.MeshSizeMin", 0.0)
        _safe_set(gmsh, log, "Mesh.MeshSizeMax", 1e22)
        _safe_set(gmsh, log, "Mesh.MeshSizeFactor", factor)
        _safe_set(gmsh, log, "Mesh.CharacteristicLengthFactor", factor)
        log.append(f"factor = char_len/mean_edge = {char_len:g}/{mean_edge:.4g} = {factor:.4g}")
        return

    raise ValueError(f"unknown config kind: {kind}")


def _extract_tets(gmsh):
    tags, coords, _ = gmsh.model.mesh.getNodes()
    tag_arr = np.asarray(tags, dtype=np.int64)
    nodes = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
    _, nidx = gmsh.model.mesh.getElementsByType(4)
    if nidx is None or len(nidx) == 0:
        raise RuntimeError("generate(3) produced no tets")
    lut = np.zeros(int(tag_arr.max()) + 1, dtype=np.int64)
    lut[tag_arr] = np.arange(len(tag_arr), dtype=np.int64)
    return nodes, lut[np.asarray(nidx, dtype=np.int64)].reshape(-1, 4)


def run_trial(surf, ref_vol: float, kind: str, char_len: float, surf_label: str) -> dict:
    rec: dict = {
        "config": kind,
        "char_len_mm": float(char_len),
        "surface": surf_label,
        "ok": False,
        "log": [],
    }
    mean_edge = _mean_edge_mm(surf)
    rec["mean_input_edge_mm"] = mean_edge

    if kind in ("reparam", "reparam-baseline"):
        angles, for_reparam = ANGLES_REPARAM, True
    else:
        angles, for_reparam = ANGLES_LEGACY, False
    rec["for_reparametrization"] = for_reparam

    # forReparametrization=True on an open/self-intersecting discrete surface can
    # hang inside gmsh's reparametrisation; require a closed manifold input.
    if for_reparam and surf.n_open_edges > 0:
        rec["time_s"] = 0.0
        rec["error"] = (
            f"skipped: surface has {surf.n_open_edges} open edges; "
            "forReparametrization=True needs a closed manifold surface"
        )
        return rec

    try:
        gmsh = _import_gmsh()
    except Exception as exc:
        rec["time_s"] = 0.0
        rec["error"] = f"import gmsh: {type(exc).__name__}: {exc}"
        return rec

    tmp_stl = Path(tempfile.gettempdir()) / f"_gmsh_probe_{os.getpid()}.stl"
    surf.save(str(tmp_stl))

    gmsh.initialize(str(LIB))
    t0 = time.time()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)
        try:
            ang = _build_solid(gmsh, tmp_stl, angles, for_reparam)
            rec["classify_angle_deg"] = round(math.degrees(ang), 3)
        except Exception as exc:
            rec["error"] = f"build_solid: {type(exc).__name__}: {exc}"
            return _finish(rec, t0, gmsh)

        _apply_config(gmsh, kind, char_len, mean_edge, rec["log"])

        try:
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.generate(3)
        except Exception as exc:
            rec["error"] = f"generate: {type(exc).__name__}: {exc}"
            return _finish(rec, t0, gmsh)

        try:
            nodes, tets = _extract_tets(gmsh)
        except Exception as exc:
            rec["error"] = f"extract: {type(exc).__name__}: {exc}"
            return _finish(rec, t0, gmsh)

        # quality via the project's own metrics; no joint faces -> centre far away
        res = M._finalize(
            surf, nodes, tets, method=f"probe-{kind}",
            joint_center=(1e9, 1e9, 1e9), char_len=char_len,
        )
        vol = float(res["volume_mm3"])
        rec.update(
            ok=True,
            n_tets=int(res["n_tets"]),
            n_nodes=int(res["n_nodes"]),
            volume_mm3=vol,
            volume_error_pct=(
                abs(vol - ref_vol) / ref_vol * 100.0 if ref_vol else float("nan")
            ),
            min_jacobian=float(res["min_jacobian"]),
            jacobian_frac_lt_0_3=float(res["jacobian_frac_lt_0_3"]),
        )
        return _finish(rec, t0, gmsh)
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass


def _finish(rec: dict, t0: float, gmsh) -> dict:
    rec["time_s"] = round(time.time() - t0, 2)
    return rec


# ---------------------------------------------------------------------------
# surfaces
# ---------------------------------------------------------------------------
def load_surfaces(stl_path: Path, modes: list[str], reduction: float) -> dict:
    import pyvista as pv

    raw = pv.read(str(stl_path))
    raw = raw.clean(point_merging=True, tolerance=1e-6, absolute=True).triangulate()
    if raw.n_open_edges:
        raw = raw.fill_holes(1e9)
        raw = raw.clean(point_merging=True, tolerance=1e-6, absolute=True)

    out = {}
    if "raw" in modes:
        out["raw"] = raw
    if "iso-dec" in modes:
        iso = M._isosurface(raw, spacing=max(1.0, 2.0))
        dec = iso.decimate_pro(
            reduction=float(reduction), splitting=True, feature_angle=30.0,
            preserve_topology=True, boundary_vertex_deletion=False,
        ).triangulate().clean(point_merging=True, tolerance=1e-6, absolute=True)
        if dec.n_open_edges:
            dec = dec.fill_holes(1e9)
            dec = dec.clean(point_merging=True, tolerance=1e-6, absolute=True)
        out["iso-dec"] = dec
    return out


def _default_ref_volume(stl_path: Path) -> float:
    sib = stl_path.parent / "summary.json"
    if sib.is_file():
        try:
            j = json.loads(sib.read_text(encoding="utf-8"))
            v = j.get("surface", {}).get("enclosed_volume_mm3")
            if v:
                return float(v)
        except Exception:
            pass
    return 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stl", default=str(
        ROOT / "temp" / "opensim_fe" / "parts" / "tibia_r" / "tibia_r_surface.stl"))
    ap.add_argument("--label", default=None, help="tag for output json (default: stl stem)")
    ap.add_argument("--ref-volume-mm3", type=float, default=None,
                    help="reference (solid-boundary) volume; else sibling summary.json")
    ap.add_argument("--surfaces", default="raw,iso-dec")
    ap.add_argument("--configs", default=",".join(CONFIG_NAMES))
    ap.add_argument("--chars", default="3,4,5,6,8,10")
    ap.add_argument("--reduction", type=float, default=0.70)
    args = ap.parse_args()

    stl_path = Path(args.stl).resolve()
    if not stl_path.is_file():
        raise FileNotFoundError(stl_path)
    label = args.label or stl_path.stem
    ref_vol = args.ref_volume_mm3 if args.ref_volume_mm3 is not None else _default_ref_volume(stl_path)
    modes = [m.strip() for m in args.surfaces.split(",") if m.strip()]
    configs = [c.strip() for c in args.configs.split(",") if c.strip()]
    chars = [float(c) for c in args.chars.split(",") if c.strip()]

    print(f"[probe] stl={stl_path}")
    print(f"[probe] ref_volume={ref_vol:.3f} mm3  surfaces={modes}  configs={configs}  chars={chars}",
          flush=True)

    surfaces = load_surfaces(stl_path, modes, args.reduction)
    for k, s in surfaces.items():
        print(f"[probe] surface {k}: faces={s.n_faces} open_edges={s.n_open_edges} "
              f"mean_edge={_mean_edge_mm(s):.3f} mm", flush=True)

    records: list[dict] = []
    for sname, surf in surfaces.items():
        for kind in configs:
            for cl in chars:
                r = run_trial(surf, ref_vol, kind, cl, sname)
                records.append(r)
                if r["ok"]:
                    print(f"[probe] {sname:8s} {kind:14s} cl={cl:<5g} "
                          f"tets={r['n_tets']:>8d} volerr={r['volume_error_pct']:6.2f}% "
                          f"J<0.3={r['jacobian_frac_lt_0_3']*100:5.2f}% "
                          f"t={r['time_s']:6.1f}s", flush=True)
                else:
                    print(f"[probe] {sname:8s} {kind:14s} cl={cl:<5g} "
                          f"FAIL {r.get('error')} t={r['time_s']:.1f}s", flush=True)
                if r["log"]:
                    for ln in r["log"]:
                        print(f"          ~ {ln}", flush=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"gmsh_size_probe_{label}.json"
    payload = {
        "stl": str(stl_path.relative_to(ROOT)) if str(stl_path).startswith(str(ROOT)) else str(stl_path),
        "label": label,
        "ref_volume_mm3": ref_vol,
        "surfaces": {k: {"n_faces": int(v.n_faces),
                         "n_open_edges": int(v.n_open_edges),
                         "mean_input_edge_mm": _mean_edge_mm(v)}
                     for k, v in surfaces.items()},
        "records": records,
    }
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[probe] wrote {out.relative_to(ROOT)}")

    print("\n=== summary (tet counts) ===")
    print(f"{'surface':8s} {'config':14s} " + " ".join(f"{c:>9g}" for c in chars))
    for sname in surfaces:
        for kind in configs:
            cells = []
            for cl in chars:
                rr = next((x for x in records
                           if x["surface"] == sname and x["config"] == kind
                           and x["char_len_mm"] == cl), None)
                if rr is None:
                    cells.append("        -")
                elif rr["ok"]:
                    cells.append(f"{rr['n_tets']:>9d}")
                else:
                    cells.append(f"{'FAIL':>9s}")
            print(f"{sname:8s} {kind:14s} " + " ".join(cells))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
