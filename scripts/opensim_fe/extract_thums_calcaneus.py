"""Extract the right-calcaneus solid mesh + bone material from THUMS AM50 V7.1.

Option C of the OpenSim->FEBio calcaneus study: transplant the THUMS AM50
right calcaneus geometry + material into our single-bone FEBio model.

The calcaneus in the THUMS deck is two standalone *solid* parts:

    pid 81001200  R_CALCANEUS_SPON   *MAT_ELASTIC                       (trabecular)
    pid 81001300  R_CALCANEUS_CORT   *MAT_PIECEWISE_LINEAR_PLASTICITY  (cortical)

Run from the repo root (or anywhere) with::

    D:\\Project\\climbing_fall_analysis\\.venv\\Scripts\\python.exe \
        scripts/opensim_fe/extract_thums_calcaneus.py

Outputs (all under temp/opensim_fe/thums_calcaneus/):
    calcaneus_r_surface.stl      outer surface, mm, THUMS global coords
    calcaneus_r_volume.vtk       solid elements as unstructured grid
    node_ids.npy / nodes_xyz.npy volumes nodes + coords (aligned)
    surface_nodes.npy            unique node ids on the surface
    elements.npy                 connectivity (global node ids)
    calcaneus_material.json      raw *MAT / *MAT_ADD_EROSION cards + params
    summary.json                 counts / validation numbers

Nothing is re-centred: coordinates are exactly the THUMS global frame (mm).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import kmesh_io as kio  # noqa: E402

ROOT = _HERE.parents[1]
DECK_DIR = ROOT / "model" / "AM50_V71_Occupant"
MASTER = DECK_DIR / "main_THUMS_AM50_V71.k"
OUT = ROOT / "temp" / "opensim_fe" / "thums_calcaneus"

# right-calcaneus parts (pid -> (name, kind))
CALC_PIDS = {81001200: "R_CALCANEUS_SPON", 81001300: "R_CALCANEUS_CORT"}

# best-effort field names (first line of the card) -- raw values are also kept
_MAT_ELASTIC_NAMES = ["MID", "RO", "E", "PR", "DA", "DB", "DC"]
_MAT_PLP_NAMES = ["MID", "RO", "E", "PR", "SIGY", "ETAN", "BETA", "SRC", "SRP"]


def _nums(line: str) -> list[float]:
    out: list[float] = []
    for t in line.replace(",", " ").split():
        try:
            out.append(float(t))
        except ValueError:
            pass
    return out


def _is_ascii(path: Path, limit: int = 8 << 20) -> bool:
    with open(path, "rb") as fh:
        chunk = fh.read(limit)
    return all(b < 128 for b in chunk)


def parse_material_card(keyword: str, data: list[str], mid: int) -> dict:
    """Parse one *MAT block into a raw + best-effort-named record."""
    numbers_per_line = [_nums(ln) for ln in data]
    flat = [v for row in numbers_per_line for v in row]
    names = _MAT_ELASTIC_NAMES if keyword == "MAT_ELASTIC" else _MAT_PLP_NAMES
    named: dict[str, float] = {}
    for i, v in enumerate(flat):
        named[names[i] if i < len(names) else f"field_{i}"] = v
    rec = {
        "mid": mid,
        "card": keyword,
        "raw_lines": data,
        "numbers_per_line": numbers_per_line,
        "flat_numbers": flat,
        "named_leading": named,
        "raw_card_text": f"*{keyword}\n" + "\n".join(data),
    }
    if "RO" in named:
        rec["conversions"] = {
            "RO_kg_m3": named["RO"] * 1e12,  # RO raw is tonne/mm^3 == 1e-12 kg/m^3
            "E_MPa": named.get("E"),
            "nu": named.get("PR"),
            "SIGY_MPa": named.get("SIGY"),
            "note": "raw RO in tonne/mm^3; E/SIGY already in MPa (THUMS mm-N-MPa)",
        }
    return rec


def extract_materials(blocks) -> dict:
    mats: dict[int, dict] = {}
    sections: dict[int, list[str]] = {}
    erosions: list[dict] = []
    part_titles: dict[int, str] = {}
    part_raw: dict[int, str] = {}

    mat_keywords = {k for k in blocks if k.startswith("MAT_")}
    for kw in mat_keywords:
        for data in blocks[kw]:
            if not data:
                continue
            first = _nums(data[0])
            if not first:
                continue
            mid = int(first[0])
            if mid in CALC_PIDS and kw != "MAT_ADD_EROSION" and mid not in mats:
                mats[mid] = parse_material_card(kw, data, mid)
            if kw == "MAT_ADD_EROSION":
                erosions.append({"mid": mid, "raw_lines": data,
                                 "flat_numbers": [v for ln in data for v in _nums(ln)]})

    for data in blocks.get("SECTION_SOLID", []):
        for ln in data:
            v = _nums(ln)
            if v and int(v[0]) in CALC_PIDS:
                sections[int(v[0])] = ln.split()

    for data in blocks.get("PART", []):
        if len(data) < 2:
            continue
        v = _nums(data[1]) if len(data) > 1 else []
        if v and int(v[0]) in CALC_PIDS:
            pid = int(v[0])
            part_titles[pid] = data[0]
            part_raw[pid] = data[0] + "\n" + data[1]

    return {
        "materials": mats,
        "sections": sections,
        "erosion": erosions,
        "part_titles": part_titles,
        "part_raw": part_raw,
    }


def main() -> int:
    if not MASTER.exists():
        print(f"ERROR: master deck not found: {MASTER}", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)

    tree = kio.read_include_tree(MASTER)
    print(f"*INCLUDE tree: {len(tree)} file(s)")
    for p in tree:
        print(f"    {p.relative_to(ROOT)}  ({p.stat().st_size/1e6:.1f} MB)")

    # ---- geometry (single streaming pass over the whole tree) -------------
    print("\nscanning nodes + solid elements ...")
    nodes, elems = kio.scan(tree, pids=list(CALC_PIDS))
    print(f"  {len(nodes):,} nodes, {len(elems):,} elements for "
          f"{sorted(CALC_PIDS)}")
    if not elems:
        print("ERROR: no solid elements found for the calcaneus pids", file=sys.stderr)
        return 3

    node_ids = np.array(sorted(nodes), dtype=np.int64)
    xyz = np.array([nodes[int(n)] for n in node_ids], dtype=np.float64)
    node_range = [int(node_ids.min()), int(node_ids.max())]

    # per-pid counts + node-count distribution
    per_pid: dict[str, dict] = {}
    for pid, name in CALC_PIDS.items():
        sel = kio.pid_node_elements(pid, elems)
        n_per = np.array([len(ns) for _, _, ns in sel.elements], dtype=np.int64)
        dist = {int(k): int(v) for k, v in zip(*np.unique(n_per, return_counts=True))}
        per_pid[name] = {
            "pid": pid,
            "n_elements": len(sel.elements),
            "n_nodes": len(sel.node_ids),
            "nodes_per_element": dist,
            "distinct_nodes_per_element": {
                int(k): int(v)
                for k, v in zip(*np.unique(
                    np.array([len(set(ns)) for _, _, ns in sel.elements]), return_counts=True))
            },
        }

    # ---- surface ----------------------------------------------------------
    tris = kio.free_surface_triangles(node_ids, elems)
    surf_nodes = sorted({n for t in tris for n in t})
    estats = kio.boundary_edge_stats(tris)
    signed_vol = kio.triangle_volume(tris, nodes)
    print(f"  surface: {len(tris):,} triangles, {len(surf_nodes):,} surface nodes")
    print(f"  boundary edges: {estats['boundary_edges']:,} / {estats['total_edges']:,} "
          f"({estats['boundary_ratio']*100:.2f}%), "
          f"non-manifold: {estats['non_manifold_edges']}")
    print(f"  enclosed volume (divergence thm): {abs(signed_vol)/1000.0:.1f} cm^3")

    assert len(tris) > 0, "empty surface"
    assert len(node_ids) > 0, "no nodes"
    assert estats["boundary_ratio"] < 0.25, (
        f"surface not closed-ish: boundary ratio {estats['boundary_ratio']:.3f}")

    # ---- material / part / section / erosion ------------------------------
    blocks: dict[str, list[list[str]]] = {}
    for kw, data in kio.iter_keyword_blocks(
            tree, ["PART", "SECTION_SOLID", "MAT_ADD_EROSION"]
            + [f"MAT_{k}" for k in ("ELASTIC", "PIECEWISE_LINEAR_PLASTICITY")]):
        blocks.setdefault(kw, []).append(data)
    Minfo = extract_materials(blocks)

    # also: do the mat_bone fracture/no_fracture files override calcaneus?
    override_files = [p for p in tree if "mat_bone" in p.name]
    override_mids: dict[str, list[int]] = {}
    for p in override_files:
        mids: list[int] = []
        for kw, data in kio.iter_keyword_blocks([p], ["MAT_PIECEWISE_LINEAR_PLASTICITY",
                                                      "MAT_ELASTIC", "MAT_DAMAGE_2",
                                                      "MAT_PLASTICITY_WITH_DAMAGE"]):
            v = _nums(data[0]) if data else []
            if v:
                mids.append(int(v[0]))
        override_mids[p.name] = sorted(set(mids))

    erosion_for_calc = [e for e in Minfo["erosion"] if e["mid"] in CALC_PIDS]
    deck_ascii = {p.name: _is_ascii(p) for p in tree if p.suffix == ".k"}

    # ---- write geometry artifacts ----------------------------------------
    stl_path = OUT / "calcaneus_r_surface.stl"
    vtk_path = OUT / "calcaneus_r_volume.vtk"
    kio.write_stl(stl_path, tris, nodes, solid_name="calcaneus_r")
    vtk_summary = kio.write_vtk(vtk_path, node_ids, elems, nodes)
    np.save(OUT / "node_ids.npy", node_ids)
    np.save(OUT / "nodes_xyz.npy", xyz)
    np.save(OUT / "surface_nodes.npy", np.array(surf_nodes, dtype=np.int64))
    eid_arr = np.array([e[0] for e in elems], dtype=np.int64)
    pid_arr = np.array([e[1] for e in elems], dtype=np.int64)
    conn_w = max(len(e[2]) for e in elems)
    conn = np.full((len(elems), conn_w), -1, dtype=np.int64)
    for i, e in enumerate(elems):
        conn[i, : len(e[2])] = e[2]
    np.savez(OUT / "calcaneus_elements.npz", eid=eid_arr, pid=pid_arr, conn=conn)

    print(f"\nwrote STL  {stl_path.relative_to(ROOT)}")
    print(f"wrote VTK  {vtk_path.relative_to(ROOT)}  {vtk_summary}")

    # ---- material JSON ----------------------------------------------------
    mat_json = {
        "source_deck": str(MASTER.relative_to(ROOT)),
        "include_tree": [str(p.relative_to(ROOT)) for p in tree],
        "part": {
            str(pid): {"pid": pid, "name": name,
                       "part_title": Minfo["part_titles"].get(pid, name),
                       "section_solid": " ".join(Minfo["sections"].get(pid, [])),
                       "raw_part_card": Minfo["part_raw"].get(pid, "")}
            for pid, name in CALC_PIDS.items()
        },
        "materials": {str(mid): rec for mid, rec in sorted(Minfo["materials"].items())},
        "mat_add_erosion": {
            "present_for_calcaneus": bool(erosion_for_calc),
            "cards_for_calcaneus": erosion_for_calc,
            "all_cards_in_deck": Minfo["erosion"],
            "note": (
                "No *MAT_ADD_EROSION card references the calcaneus MIDs "
                "(81001200 / 81001300). The only erosion cards in the deck "
                "attach to *MAT_SIMPLIFIED_RUBBER ligament/tissue materials "
                "(MIDs 81100000-81100400 right / 82100000-82100400 left) and "
                "all their criteria are 0.0 (disabled). Calcaneus failure is "
                "therefore governed only by *MAT_PIECEWISE_LINEAR_PLASTICITY "
                "yield, not erosion."),
        },
        "fracture_override_files": {
            "files": [p.name for p in override_files],
            "mids_defined": override_mids,
            "overrides_calcaneus": sorted(
                set(m for mids in override_mids.values() for m in mids) & set(CALC_PIDS)),
            "note": ("mat_bone_*_fracture.k re-defines femur/tibia/fibula bone "
                     "MIDs only; it does NOT re-define the calcaneus MIDs."),
        },
        "ASSUMED": {
            "note": ("No literature/assumed values are injected here. Every "
                     "number above is copied verbatim from the THUMS deck. Use "
                     "the 'conversions' blocks (pure unit scaling) to map into "
                     "FEBio mm-N-MPa."),
        },
    }
    (OUT / "calcaneus_material.json").write_text(
        json.dumps(mat_json, ensure_ascii=False, indent=2), encoding="utf-8")

    # ---- summary ----------------------------------------------------------
    summary = {
        "master_deck": str(MASTER.relative_to(ROOT)),
        "deck_is_ascii": deck_ascii,
        "calcaneus_pids": {str(k): v for k, v in CALC_PIDS.items()},
        "part_titles": {str(k): v for k, v in Minfo["part_titles"].items()},
        "n_nodes": int(len(node_ids)),
        "n_elements": int(len(elems)),
        "node_id_range": node_range,
        "per_part": per_pid,
        "surface": {
            "n_triangles": int(len(tris)),
            "n_surface_nodes": int(len(surf_nodes)),
            **estats,
            "enclosed_volume_mm3": abs(float(signed_vol)),
            "enclosed_volume_cm3": abs(float(signed_vol)) / 1000.0,
        },
        "vtk": vtk_summary,
        "mat_add_erosion_for_calcaneus": bool(erosion_for_calc),
        "artifacts": [
            "calcaneus_r_surface.stl", "calcaneus_r_volume.vtk",
            "node_ids.npy", "nodes_xyz.npy", "surface_nodes.npy",
            "calcaneus_elements.npz", "calcaneus_material.json",
        ],
    }
    (OUT / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print("\n=== SUMMARY ===")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
