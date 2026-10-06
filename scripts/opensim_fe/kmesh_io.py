"""kmesh_io.py -- Lightweight LS-DYNA keyword (.k) mesh I/O for THUMS decks.

Reusable, streaming reader that pulls *NODE / *ELEMENT_SOLID data out of an
LS-DYNA keyword file (or an *INCLUDE tree) and exports the geometry as STL /
legacy VTK.  Written for the THUMS AM50 V7.1 occupant deck but kept generic.

Deck quirks (verified empirically on THUMS AM50 V7.1):

    *NODE            whitespace-separated   I8 + 3*F16 (+ optional TC / RC)
    *ELEMENT_SOLID   8-char fixed-width, NO separators (80-char lines)
    *PART            10-char fixed-width
    files are ASCII and may be split across *INCLUDE cards.

The element reader therefore prefers fixed-width 8-char slicing for solid
element lines and only falls back to token/continuity splitting when a line
carries explicit separators (comma) or is not a multiple of 8 chars.

Public API
----------
read_include_tree(master_k)              -> list[Path]   recursively resolved
iter_nodes(paths)                        -> (nid, x, y, z)
iter_solid_elements(paths)               -> (eid, pid, [nid, ...])
pid_node_elements(pid, elems)            -> PidSelection(node_ids, elements)
free_surface_triangles(node_ids, elems)  -> [(a, b, c), ...]  (global node ids)
write_stl(path, tris, nodes)
write_vtk(path, node_ids, elems, nodes)

Extra helpers (used by the calcaneus extraction driver):
scan(paths, pids=None)                   -> (nodes, elems)   single pass
iter_keyword_blocks(paths, keywords=None)-> (keyword, [raw data lines])

Coordinates are returned in the deck's native global frame and units
(THUMS = mm).  Nothing is re-centred or scaled.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable, Iterator, Mapping, NamedTuple, Sequence

__all__ = [
    "read_include_tree",
    "iter_nodes",
    "iter_solid_elements",
    "pid_node_elements",
    "free_surface_triangles",
    "write_stl",
    "write_vtk",
    "scan",
    "iter_keyword_blocks",
    "PidSelection",
]

# --------------------------------------------------------------------------- #
#  low-level parsing helpers
# --------------------------------------------------------------------------- #

_KW_RE = re.compile(r"^\*+\s*([A-Za-z0-9_\-]+)")


def _keyword(line: str) -> str | None:
    """Return the upper-cased base keyword of a `*KEYWORD` line, else None."""
    m = _KW_RE.match(line.lstrip())
    return m.group(1).upper() if m else None


def _as_paths(paths) -> list[Path]:
    if isinstance(paths, (str, Path)):
        return [Path(paths)]
    return [Path(p) for p in paths]


def _ints_fixed(line: str, width: int) -> list[int] | None:
    """Slice a fixed-width numeric record. Returns None if any field is bad."""
    out: list[int] = []
    end = len(line) - (len(line) % width)
    for i in range(0, end, width):
        f = line[i : i + width].strip()
        if not f:
            continue
        try:
            out.append(int(f))
        except ValueError:
            return None
    return out


def _parse_int_tokens(line: str) -> list[int] | None:
    toks = [t for t in re.split(r"[,\s]+", line.strip()) if t]
    try:
        return [int(t) for t in toks]
    except ValueError:
        return None


def _parse_element_line(line: str) -> list[int] | None:
    """Parse one *ELEMENT_SOLID connectivity line.

    Handles the THUMS fixed-width 8-char format (no separators) as well as
    comma/space separated 'continuity' style records.
    """
    raw = line.rstrip("\r\n")
    st = raw.strip()
    if not st:
        return None
    if "," in st:
        return _parse_int_tokens(st)
    if len(st) % 8 == 0 and len(st) >= 24:
        vals = _ints_fixed(st, 8)
        if vals is not None:
            return vals
    toks = st.split()
    if len(toks) >= 3:
        return _parse_int_tokens(st)
    return _ints_fixed(st, 8)


def _stream_sections(paths) -> Iterator[tuple]:
    """Single-pass stream of records from a list of keyword files.

    Yields either ``("NODE", nid, x, y, z)`` or
    ``("SOLID", eid, pid, [nid, ...])``.
    """
    for path in _as_paths(paths):
        with open(path, "r", encoding="latin-1", errors="replace") as fh:
            cur: str | None = None
            for raw in fh:
                if not raw or raw[0] == "$":
                    continue
                st = raw.strip()
                if not st or st[0] == "$":
                    continue
                if st[0] == "*":
                    cur = _keyword(st)
                    continue
                if cur == "NODE":
                    toks = st.replace(",", " ").split()
                    if len(toks) >= 4:
                        try:
                            yield (
                                "NODE",
                                int(toks[0]),
                                float(toks[1]),
                                float(toks[2]),
                                float(toks[3]),
                            )
                        except ValueError:
                            pass
                elif cur is not None and cur.startswith("ELEMENT_SOLID"):
                    vals = _parse_element_line(raw)
                    if vals and len(vals) >= 3:
                        yield ("SOLID", vals[0], vals[1], vals[2:])


# --------------------------------------------------------------------------- #
#  *INCLUDE handling
# --------------------------------------------------------------------------- #


def _resolve_include(name: str, *bases: Path) -> Path | None:
    cand = Path(name)
    if cand.is_absolute():
        return cand if cand.exists() else None
    for base in bases:
        p = (base / cand).resolve()
        if p.exists():
            return p
    # last resort: relative to cwd
    p = cand.resolve()
    return p if p.exists() else None


def _file_includes(path: Path, fallback: Path) -> list[Path]:
    """Return include targets declared inside `path` (non-recursive)."""
    out: list[Path] = []
    pending = False
    try:
        with open(path, "r", encoding="latin-1", errors="replace") as fh:
            for raw in fh:
                st = raw.strip()
                if not st or st[0] == "$":
                    continue
                if st[0] == "*":
                    if _keyword(st) == "INCLUDE":
                        rest = st.split("$", 1)[0]
                        rest = rest[len("*INCLUDE") :].strip()
                        pending = not rest
                        if rest:
                            rp = _resolve_include(rest, path.parent, fallback)
                            if rp is not None:
                                out.append(rp)
                    else:
                        pending = False
                    continue
                if pending:
                    fn = st.split("$", 1)[0].strip()
                    if fn:
                        rp = _resolve_include(fn, path.parent, fallback)
                        if rp is not None:
                            out.append(rp)
                    pending = False
    except OSError:
        pass
    return out


def read_include_tree(master_k) -> list[Path]:
    """Recursively resolve every ``*INCLUDE`` under `master_k`.

    Returns all referenced files (the master itself is not returned), in a
    depth-first order, de-duplicated.  Relative paths are resolved against the
    including file's directory, then the master's directory.
    """
    master = Path(master_k).resolve()
    base = master.parent
    result: list[Path] = []
    seen: set[Path] = set()

    def visit(src: Path) -> None:
        for inc in _file_includes(src, base):
            rp = inc.resolve()
            if rp in seen:
                continue
            seen.add(rp)
            result.append(rp)
            visit(rp)

    visit(master)
    return result


# --------------------------------------------------------------------------- #
#  node / element iterators
# --------------------------------------------------------------------------- #


def iter_nodes(paths) -> Iterator[tuple[int, float, float, float]]:
    """Yield ``(nid, x, y, z)`` for every *NODE line in `paths`."""
    for rec in _stream_sections(paths):
        if rec[0] == "NODE":
            yield rec[1], rec[2], rec[3], rec[4]


def iter_solid_elements(paths) -> Iterator[tuple[int, int, list[int]]]:
    """Yield ``(eid, pid, [node ids])`` for every *ELEMENT_SOLID line."""
    for rec in _stream_sections(paths):
        if rec[0] == "SOLID":
            yield rec[1], rec[2], list(rec[3])


class PidSelection(NamedTuple):
    node_ids: list[int]
    elements: list[tuple[int, int, list[int]]]


def pid_node_elements(pid: int, elems) -> PidSelection:
    """Collect all elements of part `pid`.

    `elems` is any iterable of ``(eid, pid, [node ids])`` (e.g. the generator
    returned by :func:`iter_solid_elements`).  Returns the sorted unique node
    ids referenced by the part plus the list of its elements.
    """
    selected: list[tuple[int, int, list[int]]] = []
    nid: set[int] = set()
    for eid, p, ns in elems:
        if p == pid:
            ns = list(ns)
            selected.append((eid, p, ns))
            nid.update(ns)
    return PidSelection(sorted(nid), selected)


def scan(paths, pids: Sequence[int] | None = None):
    """Single-pass read of an include tree.

    Returns ``(nodes, elems)`` where ``nodes`` is ``{nid: (x, y, z)}`` and
    ``elems`` is ``[(eid, pid, [nid, ...]), ...]``.  When `pids` is given,
    only elements of those parts are kept and `nodes` is subset to the node
    ids they reference.
    """
    wanted = set(pids) if pids is not None else None
    all_nodes: dict[int, tuple[float, float, float]] = {}
    elems: list[tuple[int, int, list[int]]] = []
    for rec in _stream_sections(paths):
        if rec[0] == "NODE":
            all_nodes[rec[1]] = (rec[2], rec[3], rec[4])
        elif wanted is None or rec[2] in wanted:
            elems.append((rec[1], rec[2], list(rec[3])))

    if wanted is not None:
        used = {n for _, _, ns in elems for n in ns}
        nodes = {k: v for k, v in all_nodes.items() if k in used}
    else:
        nodes = all_nodes
    return nodes, elems


def iter_keyword_blocks(paths, keywords: Iterable[str] | None = None):
    """Yield ``(keyword, [data lines])`` for keyword blocks in `paths`.

    `$` comments and blank lines are dropped.  If `keywords` is given, only
    blocks whose base keyword is in that set are yielded.  Note that a block's
    first data line may be a title string (e.g. ``*PART``) or a numeric card
    (e.g. ``*MAT_*``); both are returned verbatim.
    """
    wanted = {k.upper() for k in keywords} if keywords is not None else None
    for path in _as_paths(paths):
        cur: str | None = None
        data: list[str] = []
        with open(path, "r", encoding="latin-1", errors="replace") as fh:
            for raw in fh:
                st = raw.strip()
                if not st or st[0] == "$":
                    continue
                if st[0] == "*":
                    if cur is not None and (wanted is None or cur in wanted):
                        yield cur, data
                    cur = _keyword(st)
                    data = []
                elif cur is not None:
                    data.append(st)
        if cur is not None and (wanted is None or cur in wanted):
            yield cur, data


# --------------------------------------------------------------------------- #
#  surface extraction
# --------------------------------------------------------------------------- #

_HEX_FACES = (
    (0, 1, 2, 3),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)
_TET_FACES = ((0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3))
_WEDGE_FACES = (
    (0, 1, 2),
    (3, 5, 4),
    (0, 3, 4, 1),
    (1, 4, 5, 2),
    (2, 5, 3, 0),
)


def _distinct(ns: Sequence[int]) -> list[int]:
    seen: set[int] = set()
    out: list[int] = []
    for n in ns:
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out


def _element_faces(ns: Sequence[int]) -> list[tuple[int, ...]]:
    """Local face templates for a solid element, keyed by its topology."""
    d = _distinct(ns)
    if len(d) == 4:
        return [tuple(d[i] for i in f) for f in _TET_FACES]
    if len(ns) == 8 and len(d) == 8:
        return [tuple(ns[i] for i in f) for f in _HEX_FACES]
    if len(d) == 6:
        return [tuple(d[i] for i in f) for f in _WEDGE_FACES]
    # fallback: treat as a tet on the first four distinct nodes
    if len(d) >= 4:
        return [tuple(d[i] for i in f) for f in _TET_FACES]
    return []


def _unwrap_elem(e) -> list[int]:
    if len(e) == 3 and isinstance(e[2], (list, tuple)):
        return list(e[2])
    return list(e)


def free_surface_triangles(node_ids, elems) -> list[tuple[int, int, int]]:
    """Extract the outer (boundary) surface of a solid element set as triangles.

    `elems` is an iterable of node-id lists (``[nid, ...]``) or of
    ``(eid, pid, [nid, ...])`` tuples.  A polygonal face is on the boundary
    when no other element shares it.  Because a hexahedron exposes quad faces
    while a tetrahedron exposes triangles, a quad that is fully covered by two
    boundary triangles (the typical hex/tet conformal interface) is also
    treated as internal and removed.  Returns triangles as tuples of global
    node ids.  `node_ids` is the set/sequence of ids in scope, used only to
    validate connectivity.
    """
    scope = {int(x) for x in node_ids} if node_ids is not None else None
    quad_count: dict[tuple[int, ...], int] = {}
    quad_or: dict[tuple[int, ...], tuple[int, ...]] = {}
    tri_count: dict[tuple[int, ...], int] = {}
    tri_or: dict[tuple[int, ...], tuple[int, ...]] = {}

    for e in elems:
        ns = _unwrap_elem(e)
        if scope is not None and not set(ns) <= scope:
            raise KeyError(f"element references nodes outside scope: {set(ns) - scope}")
        for f in _element_faces(ns):
            d = _distinct(f)
            if len(d) < 3:
                continue
            key = tuple(sorted(d))
            if len(d) == 4:
                quad_count[key] = quad_count.get(key, 0) + 1
                quad_or.setdefault(key, tuple(d))
            else:
                tri_count[key] = tri_count.get(key, 0) + 1
                tri_or.setdefault(key, tuple(d))

    b_quads = [k for k, c in quad_count.items() if c == 1]
    b_tris = {k for k, c in tri_count.items() if c == 1}

    internal_tris: set[tuple[int, ...]] = set()
    tris: list[tuple[int, int, int]] = []
    for k in b_quads:
        a, b, c, d = quad_or[k]
        covered = False
        for t1, t2 in (((a, b, c), (a, c, d)), ((a, b, d), (b, c, d))):
            k1, k2 = tuple(sorted(t1)), tuple(sorted(t2))
            if k1 in b_tris and k2 in b_tris:
                internal_tris.update((k1, k2))
                covered = True
                break
        if not covered:
            tris.append((a, b, c))
            tris.append((a, c, d))

    for k in b_tris - internal_tris:
        f = tri_or[k]
        tris.append((f[0], f[1], f[2]))
    return tris


def boundary_edge_stats(tris: Sequence[tuple[int, int, int]]) -> dict:
    """Count edges that are not shared by exactly two triangles."""
    edge: dict[tuple[int, int], int] = {}
    for a, b, c in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            k = (u, v) if u < v else (v, u)
            edge[k] = edge.get(k, 0) + 1
    boundary = sum(1 for n in edge.values() if n == 1)
    non_manifold = sum(1 for n in edge.values() if n > 2)
    return {
        "total_edges": len(edge),
        "boundary_edges": boundary,
        "non_manifold_edges": non_manifold,
        "boundary_ratio": boundary / len(edge) if edge else 0.0,
    }


def triangle_volume(tris: Sequence[tuple[int, int, int]], nodes: Mapping[int, Sequence[float]]) -> float:
    """Signed volume of a closed triangle surface via the divergence theorem."""
    vol = 0.0
    for a, b, c in tris:
        xa = nodes[a]
        xb = nodes[b]
        xc = nodes[c]
        vol += (
            xa[0] * (xb[1] * xc[2] - xb[2] * xc[1])
            - xa[1] * (xb[0] * xc[2] - xb[2] * xc[0])
            + xa[2] * (xb[0] * xc[1] - xb[1] * xc[0])
        )
    return vol / 6.0


# --------------------------------------------------------------------------- #
#  writers
# --------------------------------------------------------------------------- #


def _norm(v: Sequence[float]) -> tuple[float, float, float]:
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if n == 0.0:
        return (0.0, 0.0, 0.0)
    return (v[0] / n, v[1] / n, v[2] / n)


def write_stl(path, tris, nodes, *, solid_name: str = "surface") -> None:
    """Write an ASCII STL. `nodes` maps nid -> (x, y, z)."""
    lines = [f"solid {solid_name}"]
    for a, b, c in tris:
        pa, pb, pc = nodes[a], nodes[b], nodes[c]
        u = (pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2])
        v = (pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2])
        n = _norm((u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]))
        lines.append(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}")
        lines.append("    outer loop")
        for p in (pa, pb, pc):
            lines.append(f"      vertex {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
        lines.append("    endloop")
        lines.append("  endfacet")
    lines.append(f"endsolid {solid_name}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii")


def _vtk_cell(ns: Sequence[int]):
    """Return (vtk_cell_type, [global node ids]) for one element."""
    d = _distinct(ns)
    if len(d) == 4 and len(ns) in (4, 5, 6, 8):
        return 10, d  # VTK_TETRA
    if len(ns) == 8 and len(d) == 8:
        return 12, list(ns)  # VTK_HEXAHEDRON
    if len(d) == 6:
        return 13, d  # VTK_WEDGE
    if len(d) >= 4:
        return 10, d[:4]
    return None, None


def write_vtk(path, node_ids, elems, nodes) -> dict:
    """Write a legacy-ASCII VTK unstructured grid of the solid elements.

    `node_ids` defines the order of the points; `elems` may be raw node-id
    lists or ``(eid, pid, [nid, ...])`` tuples.  Degenerate 8-node records
    (collapsed hexes = tets) are written with their true topology.  Returns a
    small summary dict.
    """
    node_ids = [int(n) for n in node_ids]
    local = {n: i for i, n in enumerate(node_ids)}
    cells: list[tuple[int, list[int]]] = []
    for e in elems:
        if len(e) == 3 and isinstance(e[2], (list, tuple)):
            ns = list(e[2])
        else:
            ns = list(e)
        ctype, conn = _vtk_cell(ns)
        if ctype is None:
            continue
        cells.append((ctype, conn))

    pts = [nodes[n] for n in node_ids]
    n_cells = len(cells)
    cell_size = sum(1 + len(c) for _, c in cells)

    out = [
        "# vtk DataFile Version 3.0",
        "LS-DYNA solid mesh exported by kmesh_io",
        "ASCII",
        "DATASET UNSTRUCTURED_GRID",
        f"POINTS {len(pts)} double",
    ]
    for p in pts:
        out.append(f"{p[0]:.10g} {p[1]:.10g} {p[2]:.10g}")
    out.append(f"CELLS {n_cells} {cell_size}")
    for _, c in cells:
        out.append(f"{len(c)} " + " ".join(str(local[n]) for n in c))
    out.append(f"CELL_TYPES {n_cells}")
    for ctype, _ in cells:
        out.append(str(ctype))
    Path(path).write_text("\n".join(out) + "\n", encoding="ascii")

    counts: dict[int, int] = {}
    for ctype, _ in cells:
        counts[ctype] = counts.get(ctype, 0) + 1
    return {"n_points": len(pts), "n_cells": n_cells, "cell_types": counts}
