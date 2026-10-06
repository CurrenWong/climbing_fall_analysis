"""extract_part.py -- extract an arbitrary THUMS part (by PID) to STL and mesh it.

Generalises ``extract_thums_calcaneus.py`` from the right calcaneus to *any*
set of THUMS solid PIDs (a bone = one or more parts, e.g. trabecular + cortical
layers).  S3.1 scout tool: proves the OpenSim-FE meshing pipeline (calcaneus)
transfers to other bones.

Unlike the calcaneus driver this script does NOT assume the calcaneus joint
frame and does NOT touch the existing artifacts.  It reuses, read-only:

    scripts/opensim_fe/kmesh_io.py          (deck parse + surface)
    src/climbing/coupling/meshing.py        (gmsh tet remesh + quality metrics)

Coordinates are kept in the THUMS native frame (mm); nothing is re-centred or
scaled -- note ``meshing._load_surface_mm`` multiplies by 1000 (it assumes a
metre-unit VTP) so it is deliberately NOT used here.

Usage (repo root, PYTHONPATH=src)::

    python scripts/opensim_fe/extract_part.py --name tibia_r \
        --pids 81000600,81000601,81000700 --char-len 3.0

    # extract + surface only (no meshing)
    python scripts/opensim_fe/extract_part.py --name skull_pilot \
        --pids 88000001,88000002,88000003 --report-only

Outputs under ``temp/opensim_fe/parts/<name>/``:
    <name>_surface.stl      raw outer surface (mm, THUMS frame)
    <name>_volume.vtk       solid elements as unstructured grid
    <name>_elements.npz     eid / pid / connectivity
    node_ids.npy nodes_xyz.npy surface_nodes.npy
    <name>_mesh.npz         meshed (nodes_xyz, tets) if meshed
    summary.json            counts + mesh quality
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import kmesh_io as kio  # noqa: E402

ROOT = _HERE.parents[1]
MASTER = ROOT / "model" / "AM50_V71_Occupant" / "main_THUMS_AM50_V71.k"
PARTS_INDEX = ROOT / "temp" / "opensim_fe" / "parts" / "parts_index.json"
OUT_BASE = ROOT / "temp" / "opensim_fe" / "parts"

#: single-mesh upper bound before pyfebio.xplt.to_hdf5 trips its HDF5
#: attribute limit (observed ~5e4 tet in the calcaneus sweep, S3 note).
HDF5_TET_LIMIT = 50_000


def _parse_pids(text: str) -> list[int]:
    out: list[int] = []
    for tok in text.replace(" ", ",").split(","):
        tok = tok.strip()
        if tok:
            out.append(int(tok))
    return out


def _names(pids: list[int]) -> dict[int, str]:
    if not PARTS_INDEX.is_file():
        return {}
    idx = json.loads(PARTS_INDEX.read_text(encoding="utf-8"))
    by_pid = {int(r["pid"]): (r.get("cid_name") or r.get("title") or "") for r in idx}
    return {p: by_pid.get(p, "") for p in pids}


def extract(name: str, pids: list[int], out_root: Path) -> tuple[dict, dict, list]:
    """Extract raw solid mesh + outer surface. Returns (nodes, elems, tris)."""
    if not MASTER.is_file():
        raise FileNotFoundError(f"master deck not found: {MASTER}")
    tree = kio.read_include_tree(MASTER)
    print(f"[extract] *INCLUDE tree: {len(tree)} file(s)", flush=True)

    t0 = time.time()
    nodes, elems = kio.scan(tree, pids=pids)
    print(f"[extract] scan done in {time.time()-t0:.1f}s: "
          f"{len(nodes):,} nodes, {len(elems):,} solid elements", flush=True)
    return nodes, elems, tree


def surface_and_write(name: str, pids: list[int], nodes: dict, elems: list,
                      out_dir: Path) -> dict:
    if not elems:
        raise RuntimeError(f"no solid elements for pids={pids}")

    node_ids = np.array(sorted(nodes), dtype=np.int64)
    xyz = np.array([nodes[int(n)] for n in node_ids], dtype=np.float64)

    per_pid: dict[str, dict] = {}
    names = _names(pids)
    for pid in pids:
        sel = kio.pid_node_elements(pid, elems)
        n_per = np.array([len(ns) for _, _, ns in sel.elements], dtype=np.int64)
        dist = ({int(k): int(v) for k, v in zip(*np.unique(n_per, return_counts=True))}
                if len(n_per) else {})
        per_pid[str(pid)] = {
            "pid": pid, "name": names.get(pid, ""),
            "n_elements": len(sel.elements), "n_nodes": len(sel.node_ids),
            "nodes_per_element": dist,
            "distinct_nodes_per_element": (
                {int(k): int(v) for k, v in zip(*np.unique(
                    np.array([len(set(ns)) for _, _, ns in sel.elements]),
                    return_counts=True))} if sel.elements else {}),
        }

    tris = kio.free_surface_triangles(node_ids, elems)
    surf_nodes = sorted({n for t in tris for n in t})
    estats = kio.boundary_edge_stats(tris)
    signed_vol = kio.triangle_volume(tris, nodes)
    print(f"[surface] {len(tris):,} triangles, {len(surf_nodes):,} surface nodes, "
          f"boundary_edges={estats['boundary_edges']}, "
          f"non_manifold={estats['non_manifold_edges']}", flush=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    stl_path = out_dir / f"{name}_surface.stl"
    vtk_path = out_dir / f"{name}_volume.vtk"
    kio.write_stl(stl_path, tris, nodes, solid_name=name)
    vtk_summary = kio.write_vtk(vtk_path, node_ids, elems, nodes)
    np.save(out_dir / "node_ids.npy", node_ids)
    np.save(out_dir / "nodes_xyz.npy", xyz)
    np.save(out_dir / "surface_nodes.npy", np.array(surf_nodes, dtype=np.int64))
    eid = np.array([e[0] for e in elems], dtype=np.int64)
    pid_arr = np.array([e[1] for e in elems], dtype=np.int64)
    conn_w = max(len(e[2]) for e in elems)
    conn = np.full((len(elems), conn_w), -1, dtype=np.int64)
    for i, e in enumerate(elems):
        conn[i, : len(e[2])] = e[2]
    np.savez(out_dir / f"{name}_elements.npz", eid=eid, pid=pid_arr, conn=conn)

    return {
        "n_nodes": int(len(node_ids)),
        "n_elements": int(len(elems)),
        "node_id_range": [int(node_ids.min()), int(node_ids.max())],
        "per_part": per_pid,
        "surface": {
            "n_triangles": int(len(tris)),
            "n_surface_nodes": int(len(surf_nodes)),
            **estats,
            "enclosed_volume_mm3": abs(float(signed_vol)),
            "enclosed_volume_cm3": abs(float(signed_vol)) / 1000.0,
        },
        "vtk": vtk_summary,
        "stl_path": str(stl_path.relative_to(ROOT)),
    }


def _load_stl_mm(stl_path: Path):
    """Read STL as-is (already mm) -- clean / triangulate / fill small holes."""
    import pyvista as pv

    surf = pv.read(str(stl_path))
    surf = surf.clean(point_merging=True, tolerance=1e-6, absolute=True)
    surf = surf.triangulate()
    if surf.n_open_edges:
        surf = surf.fill_holes(1e9)
        surf = surf.clean(point_merging=True, tolerance=1e-6, absolute=True)
    return surf


def mesh_stl(stl_path: Path, char_len: float, reduction: float,
             max_tets: int = HDF5_TET_LIMIT, ref_volume_mm3: float | None = None) -> dict:
    """Mesh a bone STL, preferring gmsh and enforcing the HDF5 tet ceiling.

    Path strategy (first that succeeds AND stays <= ``max_tets``):

    1. ``gmsh-raw``      -- gmsh tet remesh of the raw union surface;
    2. ``iso-dec-gmsh``  -- marching-cubes isosurface -> decimate -> gmsh
                            (repairs non-manifold seams before parametrising);
    3. ``voxel``         -- ``meshing._delaunay_tet`` (voxel + scipy Delaunay),
                            the project's established bone fallback.

    Reuses meshing internals directly (no coordinate rescale; no calcaneus
    joint-face checks).  ``ref_volume_mm3`` (the divergence-theorem volume of
    the raw solid boundary) is used as the reference when given, because
    pyvista's ``.volume`` is unreliable on these non-manifold union surfaces.
    """
    from climbing.coupling import meshing as M

    raw = _load_stl_mm(stl_path)
    ref_vol = float(ref_volume_mm3) if ref_volume_mm3 else abs(float(raw.volume))
    print(f"[mesh] raw surface: faces={raw.n_faces} open_edges={raw.n_open_edges} "
          f"ref_volume={ref_vol:.0f} mm3", flush=True)

    nodes = tets = None
    method = None

    # --- 1. gmsh on the raw surface -------------------------------------
    try:
        n1, t1 = M._gmsh_tet_remesh(raw, char_len)
        if len(t1) <= max_tets:
            nodes, tets, method = n1, t1, "gmsh-raw"
        else:
            print(f"[mesh] gmsh-raw: {len(t1)} tets > {max_tets} ceiling; next path", flush=True)
    except Exception as exc:  # noqa: BLE001 - GmshSolidError / gmsh internals
        print(f"[mesh] gmsh-raw failed: {type(exc).__name__}: {exc}", flush=True)

    # --- 2. isosurface -> decimate -> gmsh ------------------------------
    if method is None:
        try:
            iso = M._isosurface(raw, spacing=max(1.0, char_len * 0.5))
            dec = iso.decimate_pro(
                reduction=float(reduction), splitting=True, feature_angle=30.0,
                preserve_topology=True, boundary_vertex_deletion=False)
            dec = dec.triangulate().clean(point_merging=True, tolerance=1e-6, absolute=True)
            if dec.n_open_edges:
                raise RuntimeError(f"isosurface still open ({dec.n_open_edges} edges)")
            n2, t2 = M._gmsh_tet_remesh(dec, char_len)
            if len(t2) <= max_tets:
                nodes, tets, method = n2, t2, f"iso-dec-gmsh(reduction={reduction:g})"
            else:
                print(f"[mesh] iso-dec-gmsh: {len(t2)} tets > {max_tets} ceiling; next path",
                      flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"[mesh] iso-dec-gmsh failed: {type(exc).__name__}: {exc}", flush=True)

    # --- 3. voxel Delaunay repair (established fallback) ----------------
    if method is None:
        n3, t3 = M._delaunay_tet(raw, char_len)
        nodes, tets, method = n3, t3, "delaunay-voxel-repair"
        if len(t3) > max_tets:
            print(f"[mesh] WARNING: voxel path still {len(t3)} tets > {max_tets}; "
                  f"increase --char-len", flush=True)

    # joint_center far away => no face is tagged as the calcaneus joint face;
    # pass the RAW surface so volume_error_pct references the true geometry.
    result = M._finalize(
        raw, nodes, tets, method=method,
        joint_center=(1e9, 1e9, 1e9), char_len=char_len,
    )
    result["n_faces_input"] = int(raw.n_faces)
    result["reduction"] = float(reduction)
    # Anchor the reference to the trusted solid-boundary volume.
    result["reference_volume_mm3"] = ref_vol
    result["volume_error_pct"] = (
        abs(float(result["volume_mm3"]) - ref_vol) / ref_vol * 100.0 if ref_vol else float("nan")
    )
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True, help="part tag, e.g. tibia_r")
    ap.add_argument("--pids", required=True, help="comma list of THUMS PIDs")
    ap.add_argument("--char-len", type=float, default=3.0, help="target tet size (mm)")
    ap.add_argument("--reduction", type=float, default=0.70, help="isosurface decimate ratio")
    ap.add_argument("--max-tets", type=int, default=HDF5_TET_LIMIT,
                    help="single-mesh tet ceiling (HDF5 attribute limit)")
    ap.add_argument("--report-only", action="store_true",
                    help="extract + surface stats only (no meshing)")
    ap.add_argument("--out", default=None, help="output dir (default parts/<name>)")
    args = ap.parse_args()

    name = args.name
    pids = _parse_pids(args.pids)
    out_dir = Path(args.out).resolve() if args.out else (OUT_BASE / name)
    names = _names(pids)
    print(f"[part] {name}: pids={pids}")
    for p in pids:
        print(f"    {p}  {names.get(p,'')}")

    nodes, elems, _tree = extract(name, pids, out_dir)
    summary = {
        "part": name,
        "pids": pids,
        "pid_names": {str(p): names.get(p, "") for p in pids},
        "source_deck": str(MASTER.relative_to(ROOT)),
        "frame": "THUMS global, mm",
    }
    summary.update(surface_and_write(name, pids, nodes, elems, out_dir))
    print(f"[extract] wrote {out_dir.relative_to(ROOT)}")

    if not args.report_only:
        t0 = time.time()
        mesh = mesh_stl(out_dir / f"{name}_surface.stl", args.char_len,
                        args.reduction, max_tets=args.max_tets,
                        ref_volume_mm3=summary["surface"]["enclosed_volume_mm3"])
        nodes_m = mesh["nodes"]
        tets = mesh["tets"]
        np.savez(out_dir / f"{name}_mesh.npz", nodes_xyz=nodes_m, tets=tets)
        mesh_summary = {
            "char_len_mm": float(args.char_len),
            "reduction": float(mesh.get("reduction", args.reduction)),
            "method": mesh.get("method"),
            "n_mesh_nodes": int(mesh["n_nodes"]),
            "n_tets": int(mesh["n_tets"]),
            "n_boundary_faces": int(mesh["n_boundary_faces"]),
            "volume_mm3": float(mesh["volume_mm3"]),
            "reference_volume_mm3": float(mesh["reference_volume_mm3"]),
            "volume_error_pct": float(mesh["volume_error_pct"]),
            "min_jacobian": float(mesh["min_jacobian"]),
            "jacobian_frac_lt_0_3": float(mesh["jacobian_frac_lt_0_3"]),
            "exceeds_hdf5_tet_limit": bool(mesh["n_tets"] > args.max_tets),
            "hdf5_tet_limit": int(args.max_tets),
            "mesh_npz": f"{name}_mesh.npz",
            "mesh_time_s": round(time.time() - t0, 1),
        }
        summary["mesh"] = mesh_summary
        print("[mesh] " + json.dumps(mesh_summary, ensure_ascii=False))

    (out_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[done] {out_dir.relative_to(ROOT)}/summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
