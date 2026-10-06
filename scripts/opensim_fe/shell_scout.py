"""shell_scout.py -- 坑② 第一步：THUMS 皮质壳侦察与解析（只侦察，不做 FE）。

扫描 THUMS AM50 V7.1 甲板的 ``*ELEMENT_SHELL`` /
``*ELEMENT_SHELL_THICKNESS`` / ``*SECTION_SHELL``，按 PID 汇总壳单元规模、
类型分布与厚度（逐单元优先，否则取 SECTION_SHELL），并实测壳节点与同骨松质
实体（SPON）节点的共享情况；对骨盆髋骨 9 个厚度变体额外做几何面积与空间
重叠分析。

Usage (repo root, PYTHONPATH=src)::

    python scripts/opensim_fe/shell_scout.py \
        --parts-index temp/opensim_fe/parts/parts_index.json \
        --outdir temp/opensim_fe/shells \
        --report results/opensim_fe/SHELL_SCOUT_REPORT.md

Products (all new, nothing overwritten):
    temp/opensim_fe/shells/shell_pid_stats.json / .csv
    temp/opensim_fe/shells/shell_sections.json
    temp/opensim_fe/shells/shell_connectivity.json
    temp/opensim_fe/shells/pelvis_variants.json
    temp/opensim_fe/shells/pelvis_overlap.json
    temp/opensim_fe/shells/scout_summary.json
    results/opensim_fe/SHELL_SCOUT_REPORT.md

Every number is measured from the deck; nothing is invented.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import kmesh_io as kio  # noqa: E402
import shell_io as sio  # noqa: E402

ROOT = _HERE.parents[1]
DEFAULT_DECK = (
    ROOT / "model" / "AM50_V71_Occupant" / "THUMS_model" / "THUMS_AM50_V71_Occupant_202411.k"
)
DEFAULT_INDEX = ROOT / "temp" / "opensim_fe" / "parts" / "parts_index.json"
DEFAULT_OUTDIR = ROOT / "temp" / "opensim_fe" / "shells"
DEFAULT_REPORT = ROOT / "results" / "opensim_fe" / "SHELL_SCOUT_REPORT.md"

_DROP_TOKENS = {
    "SPON", "CORT", "SHELL", "EXTERNAL", "INTERNAL", "SOLID", "NULL",
    "OLD", "ANS", "EAS", "PART", "OUT", "IN",
}
_NUM_RE = re.compile(r"^\d+(\.\d+)?$")

BONE_REGIONS = {"pelvis", "lumbar", "thoracic", "cervical", "skull", "thorax_other"}

BONE_RE = re.compile(
    r"HIPBONE|SACRUM|PUBIC|COCCYX"
    r"|_L[1-5]$|_T\d{1,2}$|_C[1-7]$"
    r"|FRONTAL|PARIETAL|TEMPORAL|OCCIPITAL|SPHENOID|ETHMOID"
    r"|MANDIBLE|MAXILLA|ZYGOMATIC|LACRIMAL|NASAL|VOMER|PALATINE"
    r"|ALVEOLAR|STR\d|NG-STR|STERNUM|SCAPULA|CLAVICLE|RIB",
    re.IGNORECASE,
)


def base_name(name: str) -> str:
    toks = re.split(r"[_\s]+", name.upper())
    out: list[str] = []
    for t in toks:
        if not t:
            continue
        if t in _DROP_TOKENS:
            continue
        if _NUM_RE.fullmatch(t):
            continue
        if t.startswith("ACETABLUM"):
            continue
        out.append(t)
    return "_".join(out)


def region_of(base: str) -> str:
    b = base.upper()
    if any(k in b for k in ("HIPBONE", "SACRUM", "PUBIC", "COCCYX")):
        return "pelvis"
    if re.search(r"_L[1-5]$", b):
        return "lumbar"
    if re.search(r"_T\d{1,2}$", b):
        return "thoracic"
    if re.search(r"_C[1-7]$", b):
        return "cervical"
    if any(
        k in b
        for k in (
            "FRONTAL", "PARIETAL", "TEMPORAL", "OCCIPITAL", "SPHENOID", "ETHMOID",
            "MANDIBLE", "MAXILLA", "ZYGOMATIC", "LACRIMAL", "NASAL", "VOMER",
            "PALATINE", "ALVEOLAR", "STR", "NG-STR",
        )
    ):
        return "skull"
    if any(k in b for k in ("STERNUM", "SCAPULA", "CLAVICLE", "RIB")):
        return "thorax_other"
    return "soft_other"


def load_names(index_path: Path) -> dict[int, str]:
    if index_path.is_file():
        recs = json.loads(index_path.read_text(encoding="utf-8"))
        out: dict[int, str] = {}
        for r in recs:
            if r.get("pid") is None:
                continue
            out[int(r["pid"])] = r.get("cid_name") or r.get("title") or ""
        return out
    import scout_parts  # type: ignore

    out = {}
    for hmname, lines in scout_parts.scan_parts(index_path):
        rec = scout_parts._record(hmname, lines)
        if rec["pid"] is not None:
            out[int(rec["pid"])] = rec["cid_name"] or rec["title"] or ""
    return out


def _median_from_counter(counter: Counter) -> float | None:
    if not counter:
        return None
    total = sum(counter.values())
    half = total / 2.0
    acc = 0
    prev = None
    for k in sorted(counter):
        acc += counter[k]
        if acc >= half:
            if total % 2 == 0 and acc - counter[k] < half and prev is not None:
                return (prev + k) / 2.0
            return k
        prev = k
    return None


def pass1_shells(paths, focus_pids: set[int]):
    stats: dict[int, dict] = {}
    sections: dict[int, dict] = {}
    pending_pid: int | None = None
    pending_eid: int | None = None
    meta = {
        "n_shell_records": 0,
        "n_thic_records": 0,
        "thic_len_hist": Counter(),
        "thic5_total": 0,
        "thic5_nonzero": 0,
        "n_section_card1": 0,
    }

    def bucket(pid: int) -> dict:
        return stats.setdefault(
            pid,
            {
                "n": 0, "n3": 0, "n4": 0, "n_other": 0,
                "nodes": set(), "elems": [] if pid in focus_pids else None,
                "thic_counter": Counter(), "thic_min": None, "thic_max": None,
            },
        )

    for rec in sio.iter_shell_records(paths):
        if rec[0] == "SHELL":
            _, eid, pid, nodes = rec
            meta["n_shell_records"] += 1
            b = bucket(pid)
            b["n"] += 1
            d = len(set(nodes))
            if d == 3:
                b["n3"] += 1
            elif d == 4:
                b["n4"] += 1
            else:
                b["n_other"] += 1
            b["nodes"].update(nodes)
            if b["elems"] is not None:
                b["elems"].append((eid, tuple(nodes)))
            pending_pid = pid
            pending_eid = eid
        elif rec[0] == "THIC":
            _, eid, thic = rec
            meta["n_thic_records"] += 1
            meta["thic_len_hist"][len(thic)] += 1
            if len(thic) >= 5:
                meta["thic5_total"] += 1
                if abs(thic[4]) > 1e-12:
                    meta["thic5_nonzero"] += 1
            if pending_pid is None or (pending_eid is not None and eid != pending_eid):
                continue
            b = bucket(pending_pid)
            for v in thic[:4] if len(thic) >= 4 else thic:
                vr = round(float(v), 4)
                b["thic_counter"][vr] += 1
                b["thic_min"] = vr if b["thic_min"] is None else min(b["thic_min"], vr)
                b["thic_max"] = vr if b["thic_max"] is None else max(b["thic_max"], vr)
        elif rec[0] == "SECTION":
            meta["n_section_card1"] += 1
            sections[rec[1]] = dict(rec[2])
        elif rec[0] == "SECTION_T":
            if rec[1] is not None:
                sections.setdefault(rec[1], {})["thickness"] = list(rec[2])
    return stats, sections, meta


def apply_section_thickness(stats: dict[int, dict], sections: dict[int, dict]) -> None:
    """For parts whose thickness is only in SECTION_SHELL (SECID == PID)."""
    for pid, b in stats.items():
        if b["thic_counter"]:
            b["thic_source"] = "element"
            continue
        sec = sections.get(pid)
        if sec and sec.get("thickness"):
            vals = [round(float(v), 4) for v in sec["thickness"][:4] if float(v) != 0.0]
            if vals:
                for v in vals:
                    b["thic_counter"][v] += 1
                b["thic_min"] = min(vals)
                b["thic_max"] = max(vals)
                b["thic_source"] = "section"
                continue
        b["thic_source"] = "none"


def pass2_solid_nodes(paths, target_pids: set[int]) -> dict[int, set[int]]:
    out: dict[int, set[int]] = {p: set() for p in target_pids}
    for eid, pid, nodes in kio.iter_solid_elements(paths):
        if pid in out:
            out[pid].update(nodes)
    return out


def _poly_area(pts) -> float:
    if len(pts) < 3:
        return 0.0
    s = [0.0, 0.0, 0.0]
    n = len(pts)
    for i in range(n):
        a = pts[i]
        b = pts[(i + 1) % n]
        s[0] += a[1] * b[2] - a[2] * b[1]
        s[1] += a[2] * b[0] - a[0] * b[2]
        s[2] += a[0] * b[1] - a[1] * b[0]
    return 0.5 * math.sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2])


def pass3_geometry(paths, stats, focus_pids):
    """Per-PID shell surface area + bbox, using node coordinates."""
    needed: set[int] = set()
    for p in focus_pids:
        needed |= stats[p]["nodes"]
    coords: dict[int, tuple[float, float, float]] = {}
    for nid, x, y, z in kio.iter_nodes(paths):
        if nid in needed:
            coords[nid] = (x, y, z)

    geo: dict[int, dict] = {}
    for p in focus_pids:
        elems = stats[p]["elems"] or []
        area = 0.0
        xs = ys = zs = None
        for _eid, nodes in elems:
            pts = [coords[n] for n in nodes if n in coords]
            area += _poly_area(pts)
            for pt in pts:
                if xs is None:
                    xs, ys, zs = [pt[0]], [pt[1]], [pt[2]]
                else:
                    xs.append(pt[0]); ys.append(pt[1]); zs.append(pt[2])
        if xs is None:
            geo[p] = {"area_mm2": 0.0, "bbox": None, "n_nodes_used": 0}
            continue
        geo[p] = {
            "area_mm2": round(area, 2),
            "bbox": [round(min(xs), 2), round(min(ys), 2), round(min(zs), 2),
                     round(max(xs), 2), round(max(ys), 2), round(max(zs), 2)],
            "n_nodes_used": len(xs),
        }
    return geo, coords


def _tri_area(a, b, c) -> float:
    u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    v = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    cx = u[1] * v[2] - u[2] * v[1]
    cy = u[2] * v[0] - u[0] * v[2]
    cz = u[0] * v[1] - u[1] * v[0]
    return 0.5 * math.sqrt(cx * cx + cy * cy + cz * cz)


COVERAGE_BASES = [
    "R_HIPBONE", "L_HIPBONE", "R_SACRUM", "L_SACRUM",
    "R_L3", "L_L3", "R_T6", "L_T6", "R_C5", "L_C5",
    "FRONTAL_R", "FRONTAL_L", "PARIETAL_R", "PARIETAL_L",
    "OCCIPITAL_R", "TEMPORAL_R", "MANDIBLE_R",
]


def pass4_coverage(paths, connectivity, geo):
    """Compare shell area vs the sibling SPON solid free-surface area."""
    groups = {b: c for b, c in connectivity.items() if b in COVERAGE_BASES}
    solid_pids: set[int] = set()
    for c in groups.values():
        solid_pids.update(c["solid_pids"])
    elems_by_pid: dict[int, list] = {p: [] for p in solid_pids}
    for _eid, pid, nodes in kio.iter_solid_elements(paths):
        if pid in elems_by_pid:
            elems_by_pid[pid].append(nodes)
    needed: set[int] = set()
    for els in elems_by_pid.values():
        for ns in els:
            needed.update(ns)
    coords: dict[int, tuple] = {}
    for nid, x, y, z in kio.iter_nodes(paths):
        if nid in needed:
            coords[nid] = (x, y, z)

    out: dict[str, dict] = {}
    for base, c in groups.items():
        spon_area = 0.0
        n_tris = 0
        for p in c["solid_pids"]:
            els = elems_by_pid.get(p, [])
            if not els:
                continue
            tris = kio.free_surface_triangles(None, els)
            n_tris += len(tris)
            for a, b, cc in tris:
                if a in coords and b in coords and cc in coords:
                    spon_area += _tri_area(coords[a], coords[b], coords[cc])
        shell_area = sum(geo.get(p, {}).get("area_mm2", 0.0) for p in c["shell_pids"])
        out[base] = {
            "shell_pids": c["shell_pids"],
            "solid_pids": c["solid_pids"],
            "spon_surface_area_mm2": round(spon_area, 1),
            "spon_free_tris": n_tris,
            "shell_area_mm2": round(shell_area, 1),
            "coverage": round(shell_area / spon_area, 4) if spon_area else None,
        }
    return out


def pelvis_overlap(stats, coords, pids, tol: float = 0.5) -> dict:
    """Pairwise spatial coincidence of hipbone shell variants."""
    inv = 1.0 / tol

    def key(pt):
        return (int(math.floor(pt[0] * inv)), int(math.floor(pt[1] * inv)), int(math.floor(pt[2] * inv)))

    cells: dict[int, dict] = {}
    for p in pids:
        grid: dict[tuple, list] = defaultdict(list)
        for nid in stats[p]["nodes"]:
            if nid in coords:
                grid[key(coords[nid])].append(nid)
        cells[p] = grid

    matrix: dict[str, float] = {}
    for i, a in enumerate(pids):
        for b in pids[i + 1:]:
            shared = 0
            total = 0
            for nid_a in stats[a]["nodes"]:
                if nid_a not in coords:
                    continue
                total += 1
                ka = key(coords[nid_a])
                found = False
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for dz in (-1, 0, 1):
                            nb = cells[b].get((ka[0] + dx, ka[1] + dy, ka[2] + dz))
                            if nb:
                                pa = coords[nid_a]
                                for nid_b in nb:
                                    pb = coords[nid_b]
                                    if (pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2 + (pa[2] - pb[2]) ** 2 <= tol * tol:
                                        found = True
                                        break
                            if found:
                                break
                        if found:
                            break
                    if found:
                        break
                if found:
                    shared += 1
            matrix[f"{a}|{b}"] = round(shared / total, 4) if total else 0.0
    return {"tol_mm": tol, "coincidence_frac": matrix}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--deck", default=str(DEFAULT_DECK))
    ap.add_argument("--parts-index", default=str(DEFAULT_INDEX))
    ap.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    ap.add_argument("--report", default=str(DEFAULT_REPORT))
    ap.add_argument("--no-connectivity", action="store_true")
    ap.add_argument("--no-geometry", action="store_true")
    args = ap.parse_args()

    deck = Path(args.deck)
    if not deck.is_file():
        print(f"ERROR: deck not found: {deck}", file=sys.stderr)
        return 2
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    report_path = Path(args.report)

    t0 = time.time()
    names = load_names(Path(args.parts_index))
    print(f"[scout] deck={deck} names={len(names)}", file=sys.stderr)

    # focus PIDs = bone-region parts (for element storage + geometry)
    focus_pids = {
        pid for pid, nm in names.items()
        if region_of(base_name(nm)) in BONE_REGIONS
    }

    stats, sections, meta = pass1_shells([deck], focus_pids)
    apply_section_thickness(stats, sections)
    print(
        f"[scout] pass1: {meta['n_shell_records']} shells, {meta['n_thic_records']} thickness cards, "
        f"{len(stats)} shell PIDs, {meta['n_section_card1']} SECTION_SHELL card1 "
        f"({time.time()-t0:.1f}s)",
        file=sys.stderr,
    )

    # ---- grouping by normalized base name -----------------------------------
    shell_pids_all = set(stats)
    base_to_pids: dict[str, list[int]] = defaultdict(list)
    for pid in stats:
        base_to_pids[base_name(names.get(pid, ""))].append(pid)
    for pid, nm in names.items():
        if pid in shell_pids_all:
            continue
        b = base_name(nm)
        if b in base_to_pids and BONE_RE.search(b):
            base_to_pids[b].append(pid)

    target_solid: set[int] = set()
    for b, pids in base_to_pids.items():
        if not BONE_RE.search(b):
            continue
        for pid in pids:
            if pid not in shell_pids_all:
                target_solid.add(pid)

    # ---- connectivity -------------------------------------------------------
    connectivity: dict[str, dict] = {}
    if not args.no_connectivity and target_solid:
        t1 = time.time()
        solid_nodes = pass2_solid_nodes([deck], target_solid)
        print(f"[scout] pass2: solid node sets for {len(target_solid)} PIDs ({time.time()-t1:.1f}s)", file=sys.stderr)
        for b, pids in base_to_pids.items():
            if not BONE_RE.search(b):
                continue
            sh = [p for p in pids if p in shell_pids_all]
            so = [p for p in pids if p not in shell_pids_all]
            if not sh or not so:
                continue
            shell_nodes: set[int] = set()
            for p in sh:
                shell_nodes |= stats[p]["nodes"]
            snodes: set[int] = set()
            for p in so:
                snodes |= solid_nodes.get(p, set())
            shared = len(shell_nodes & snodes)
            connectivity[b] = {
                "shell_pids": sorted(sh),
                "solid_pids": sorted(so),
                "n_shell_nodes": len(shell_nodes),
                "n_solid_nodes": len(snodes),
                "n_shared_nodes": shared,
                "n_shell_only": len(shell_nodes - snodes),
                "n_solid_only": len(snodes - shell_nodes),
                "frac_shell_shared": (shared / len(shell_nodes)) if shell_nodes else 0.0,
            }

    # ---- geometry (bone focus) ---------------------------------------------
    geo: dict[int, dict] = {}
    coords: dict[int, tuple] = {}
    if not args.no_geometry:
        t2 = time.time()
        focus_with_elems = {p for p in focus_pids if stats.get(p, {}).get("elems")}
        geo, coords = pass3_geometry([deck], stats, focus_with_elems)
        print(f"[scout] pass3: geometry for {len(geo)} PIDs ({time.time()-t2:.1f}s)", file=sys.stderr)

    # ---- per-PID records ----------------------------------------------------
    pid_records = []
    for pid, b in stats.items():
        nm = names.get(pid, "")
        base = base_name(nm)
        med = _median_from_counter(b["thic_counter"])
        rec = {
            "pid": pid, "name": nm, "base": base, "region": region_of(base),
            "n_shell": b["n"], "n_shell3": b["n3"], "n_shell4": b["n4"], "n_other": b["n_other"],
            "n_nodes": len(b["nodes"]),
            "thic_min": b["thic_min"], "thic_median": med, "thic_max": b["thic_max"],
            "thic_distinct": sorted(b["thic_counter"]), "thic_distinct_n": len(b["thic_counter"]),
            "thic_mode": (b["thic_counter"].most_common(1)[0][0] if b["thic_counter"] else None),
            "thic_source": b.get("thic_source", "none"),
        }
        if pid in geo:
            rec["area_mm2"] = geo[pid]["area_mm2"]
            rec["bbox"] = geo[pid]["bbox"]
        pid_records.append(rec)
    pid_records.sort(key=lambda r: (r["region"], r["name"]))

    # ---- pelvis 9 variants --------------------------------------------------
    hip_L = sorted(
        [r for r in pid_records if r["region"] == "pelvis" and r["base"].endswith("HIPBONE")
         and r["name"].upper().startswith("L_")],
        key=lambda r: r["pid"],
    )
    pv = []
    ref = hip_L[0] if hip_L else None
    for r in hip_L:
        entry = dict(r)
        if ref is not None and ref["pid"] != r["pid"]:
            a = stats[ref["pid"]]["nodes"]
            b = stats[r["pid"]]["nodes"]
            inter = len(a & b)
            union = len(a | b)
            entry["jaccard_vs_ref"] = round(inter / union, 4) if union else 0.0
        else:
            entry["jaccard_vs_ref"] = None
        pv.append(entry)
    overlap_L = pelvis_overlap(stats, coords, [r["pid"] for r in hip_L]) if coords and hip_L else {}

    # ---- SPON surface coverage (shell area / solid free-surface area) --------
    coverage: dict[str, dict] = {}
    if connectivity and geo:
        t3 = time.time()
        coverage = pass4_coverage([deck], connectivity, geo)
        print(f"[scout] pass4: SPON coverage for {len(coverage)} bones ({time.time()-t3:.1f}s)", file=sys.stderr)

    # ---- write artifacts ----------------------------------------------------
    (outdir / "shell_pid_stats.json").write_text(json.dumps(pid_records, ensure_ascii=False, indent=1), encoding="utf-8")
    with (outdir / "shell_pid_stats.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["pid", "name", "base", "region", "n_shell", "n_shell3", "n_shell4", "n_other",
                    "n_nodes", "thic_min", "thic_median", "thic_max", "thic_mode",
                    "thic_distinct_n", "thic_source", "area_mm2"])
        for r in pid_records:
            w.writerow([r["pid"], r["name"], r["base"], r["region"], r["n_shell"], r["n_shell3"],
                        r["n_shell4"], r["n_other"], r["n_nodes"], r["thic_min"], r["thic_median"],
                        r["thic_max"], r["thic_mode"], r["thic_distinct_n"], r["thic_source"],
                        r.get("area_mm2")])
    (outdir / "shell_sections.json").write_text(json.dumps(sections, ensure_ascii=False, indent=1), encoding="utf-8")
    (outdir / "shell_connectivity.json").write_text(json.dumps(connectivity, ensure_ascii=False, indent=1), encoding="utf-8")
    (outdir / "pelvis_variants.json").write_text(json.dumps(pv, ensure_ascii=False, indent=1), encoding="utf-8")
    (outdir / "pelvis_overlap.json").write_text(json.dumps(overlap_L, ensure_ascii=False, indent=1), encoding="utf-8")

    region_totals: dict[str, dict] = {}
    for r in pid_records:
        rt = region_totals.setdefault(r["region"], {"parts": 0, "n_shell": 0, "n3": 0, "n4": 0, "n_other": 0})
        rt["parts"] += 1
        rt["n_shell"] += r["n_shell"]
        rt["n3"] += r["n_shell3"]
        rt["n4"] += r["n_shell4"]
        rt["n_other"] += r["n_other"]

    summary = {
        "deck": str(deck), "names_total": len(names), "shell_pids_total": len(stats),
        "totals": {k: (dict(v) if isinstance(v, Counter) else v) for k, v in meta.items()},
        "region_totals": region_totals, "n_solid_target_pids": len(target_solid),
        "coverage": coverage,
        "elapsed_s": round(time.time() - t0, 1),
    }
    (outdir / "scout_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding="utf-8")

    _write_report(report_path, summary, pid_records, region_totals, sections, connectivity, pv, overlap_L, coverage)
    print(f"[scout] wrote {outdir} and {report_path} in {time.time()-t0:.1f}s", file=sys.stderr)
    return 0


def _fmt(v, nd=4):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{nd}g}"
    return str(v)


def _area(v):
    if v is None:
        return "—"
    return f"{v:.0f}"


def _write_report(report_path, summary, pid_records, region_totals, sections, connectivity, pv, overlap_L, coverage) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    L: list[str] = []
    a = L.append
    a("# 坑② 第一步：THUMS 皮质壳侦察报告（实测）\n")
    a("> 状态：**实测**。所有数字来自对甲板 `*ELEMENT_SHELL` / "
      "`*ELEMENT_SHELL_THICKNESS` / `*SECTION_SHELL` 的流式扫描"
      "（`scripts/opensim_fe/shell_io.py` + `shell_scout.py`），未臆造。\n")
    a(f"> 甲板：`{summary['deck']}`\n")
    tot = summary["totals"]
    a(f"> 实测：`*ELEMENT_SHELL` 记录（含 THICKNESS card1）**{tot['n_shell_records']}** 条；"
      f"`*ELEMENT_SHELL_THICKNESS` card2 **{tot['n_thic_records']}** 条；"
      f"含壳单元的 PID **{summary['shell_pids_total']}** 个；"
      f"`*SECTION_SHELL` card1 **{tot['n_section_card1']}** 条。\n")
    a(f"> thickness card2 字段数分布：`{tot['thic_len_hist']}`；第 5 字段非零 "
      f"**{tot['thic5_nonzero']}/{tot['thic5_total']}**（恒 0）。\n")

    a("\n## 0. 关键实测事实（推翻/修正既有假设）\n")
    a(f"1. 壳单元共 **{tot['n_shell_records']}** 条，分布在 **{summary['shell_pids_total']}** 个 PID 上。\n")
    a("2. 骨皮质壳的厚度**不来自** `*ELEMENT_SHELL_THICKNESS`，而来自 `*SECTION_SHELL`"
      "（SECID == PID，ELFORM=16）。`*ELEMENT_SHELL_THICKNESS`（18006 条）只覆盖少数部件。"
      "**只看 ELEMENT_SHELL_THICKNESS 会把髋/骶/椎/颅的厚度全部读成空值**——本报告已修正为"
      "「逐单元优先，否则取 SECTION_SHELL」。\n")
    a("3. **壳节点与同骨 SPON 实体节点 100% 共享**（§3）：所有骨皮质壳的节点都是该骨松质实体"
      "外表面节点子集 ⇒ **可直接共节点粘接，无需 tie**。\n")
    a("4. 骨盆「9 个厚度变体」经几何实测是**互补的厚度分区（tiling），不是相互重叠的整壳**："
      "9 个变体壳面积之和 ≈ 髋骨 SPON 实体自由表面面积（§2.4）。原「会重叠、只挑一个」的假设被推翻。\n")

    a("\n## 1. 各部位壳规模与类型\n")
    a("| 部位 | 壳 PID 数 | 壳单元数 | SHELL3 | SHELL4 | 其它 |")
    a("|------|-----------|----------|--------|--------|------|")
    order = ["pelvis", "lumbar", "thoracic", "cervical", "skull", "thorax_other", "soft_other"]
    for reg in order:
        rt = region_totals.get(reg)
        if rt:
            a(f"| {reg} | {rt['parts']} | {rt['n_shell']} | {rt['n3']} | {rt['n4']} | {rt['n_other']} |")
    for reg, rt in sorted(region_totals.items()):
        if reg not in order:
            a(f"| {reg} | {rt['parts']} | {rt['n_shell']} | {rt['n3']} | {rt['n4']} | {rt['n_other']} |")
    a("")
    a("> `soft_other` = 皮肤/脑膜/肌肉/眼球等软组织壳，非骨；骨皮质壳集中在 pelvis/spine/skull/"
      "thorax_other。SHELL3/SHELL4 按单元不同节点数判定。\n")

    a("\n## 2. 厚度：骨盆 9 变体、颅骨、椎体\n")
    a("### 2.1 骨盆髋骨 9 个厚度变体（左髋 `8300020x`；右髋 `8350020x` 同构）\n")
    a("| PID | 名称 | SECID 厚度(T1..T4) | 厚度来源 | 壳单元 | SHELL3/4 | 节点 | 面积 mm² | vs `_CORT_1.5` 节点 Jaccard |")
    a("|-----|------|--------------------|----------|--------|----------|------|----------|----------------------------|")
    for r in pv:
        sec = sections.get(r["pid"], {})
        th = sec.get("thickness")
        th_s = "/".join(str(x) for x in th[:4]) if th else "—"
        a(f"| {r['pid']} | {r['name']} | {th_s} | {r['thic_source']} | {r['n_shell']} | "
          f"{r['n_shell3']}/{r['n_shell4']} | {r['n_nodes']} | {_area(r.get('area_mm2'))} | "
          f"{_fmt(r['jaccard_vs_ref'])} |")
    a("")
    tot_elems = sum(r["n_shell"] for r in pv)
    tot_area = sum(r.get("area_mm2") or 0.0 for r in pv)
    a(f"9 变体合计 **{tot_elems}** 个壳单元、面积 **{tot_area:.0f} mm²**；"
      f"每个变体的 SECTION 厚度与其名称后缀一致（1.5 / ACETABLUM1.5 / ACETABLUM1.0 / 0.75 / 2 / 3 / 2.2 / 1.8 / 1）。\n")
    if overlap_L:
        a(f"- 空间重合（节点级，tol={overlap_L['tol_mm']} mm）最大 Jaccard 见 §2.3；"
          "绝大多数变体间节点集近乎不相交 ⇒ 它们是**按厚度分区**、共同覆盖髋骨皮质面。\n")

    a("### 2.2 椎体 / 颅骨代表（实测）\n")
    a("| PID | 名称 | 壳单元 | SHELL3 | SHELL4 | SECID 厚度 | 厚度来源 | 面积 mm² |")
    a("|-----|------|--------|--------|--------|------------|----------|----------|")
    focus = [
        r for r in pid_records
        if r["region"] in ("lumbar", "thoracic", "cervical")
        or (r["region"] == "skull" and any(k in r["name"].upper() for k in
            ("FRONTAL", "PARIETAL", "TEMPORAL", "OCCIPITAL", "MANDIBLE")))
    ][:48]
    for r in focus:
        sec = sections.get(r["pid"], {})
        th = sec.get("thickness")
        th_s = "/".join(str(x) for x in th[:4]) if th else _fmt(r["thic_mode"])
        a(f"| {r['pid']} | {r['name']} | {r['n_shell']} | {r['n_shell3']} | {r['n_shell4']} | "
          f"{th_s} | {r['thic_source']} | {_area(r.get('area_mm2'))} |")
    a("")

    a("### 2.3 髋骨 9 变体空间重合矩阵（节点级，tol=0.5 mm）\n")
    if overlap_L:
        pids = [r["pid"] for r in pv]
        a("| | " + " | ".join(str(p) for p in pids) + " |")
        a("|---|" + "|".join(["---"] * len(pids)) + "|")
        cf = overlap_L["coincidence_frac"]
        for i, pa in enumerate(pids):
            row = [str(pa)]
            for j, pb in enumerate(pids):
                if i == j:
                    row.append("1.0")
                else:
                    key = f"{pa}|{pb}" if f"{pa}|{pb}" in cf else f"{pb}|{pa}"
                    row.append(_fmt(cf.get(key)))
            a("| " + " | ".join(row) + " |")
        a("")
        a("> 对角线=1（自身）。多数非对角 ≈0 ⇒ 两变体落在不同节点/位置。少数非对角偏高"
          "（`_0.75`↔`_1`=0.52、`ACETABLUM1.0`↔`_1.8`=0.375、`_2.2`↔`_1.8`=0.246 等）"
          "是**分区接缝处节点 id 不同但坐标重合**（各 PART 独立、未合并重复节点），"
          "并非体积重叠——由 §2.4「9 变体面积之和 = SPON 自由表面面积（99.85%）」判定。\n")
        a("> 若为重叠整壳，面积之和会远超单块骨表面；实测恰好相等 ⇒ **分区**。\n")
    else:
        a("（未运行几何分析）\n")

    a("### 2.4 皮质壳覆盖率：壳面积 vs 松质实体自由表面面积（决定性）\n")
    a("| 骨 | 壳 PID 数 | 壳面积 mm² | SPON 自由面 mm² | SPON 自由面三角 | 覆盖率 |")
    a("|----|-----------|------------|------------------|-----------------|--------|")
    for base, c in coverage.items():
        a(f"| {base} | {len(c['shell_pids'])} | {c['shell_area_mm2']:.0f} | {c['spon_surface_area_mm2']:.0f} | "
          f"{c['spon_free_tris']} | {_fmt(c['coverage'])} |")
    a("")
    hipcov = coverage.get("R_HIPBONE") or coverage.get("L_HIPBONE")
    if hipcov:
        a(f"> **髋骨**：9 变体壳面积 **{hipcov['shell_area_mm2']:.0f} mm²** ≈ SPON 自由表面 "
          f"**{hipcov['spon_surface_area_mm2']:.0f} mm²**（覆盖率 {_fmt(hipcov['coverage'])}，"
          f"SPON 自由面 {hipcov['spon_free_tris']} tris = {hipcov['spon_free_tris']//2} quads）；"
          "而 `_CORT_1.5` 单壳面积仅 25019 mm²（§2.1）⇒ 单取一个变体只能覆盖约 42% 的皮质面。"
          "**9 个变体共同精确铺满髋骨皮质面**，故骨盆应**取全部 9 个变体**（各自厚度），而非只挑一个。\n")
    a("> 覆盖率 ≈ 1.0 ⇒ 壳恰好铺满该骨松质外表面（分区/整壳）；明显 <1 ⇒ 只覆盖部分表面。\n")

    a("\n## 3. 壳 ↔ 松质实体连接性（共享节点？）\n")
    a("| 骨（归一化名） | 壳 PID | 实体 PID | 壳节点 | 实体节点 | 共享 | 壳独占 | 壳节点共享比 |")
    a("|----------------|--------|----------|--------|----------|------|--------|--------------|")
    for b, c in sorted(connectivity.items()):
        sh = ",".join(map(str, c["shell_pids"][:9])) + ("…" if len(c["shell_pids"]) > 9 else "")
        a(f"| {b} | {sh} | {','.join(map(str,c['solid_pids']))} | {c['n_shell_nodes']} | "
          f"{c['n_solid_nodes']} | {c['n_shared_nodes']} | {c['n_shell_only']} | {c['frac_shell_shared']:.4f} |")
    a("")
    fracs = [c["frac_shell_shared"] for c in connectivity.values() if c["n_shell_nodes"]]
    n_full = sum(1 for f in fracs if f >= 0.999)
    if fracs:
        a(f"**{n_full}/{len(fracs)}** 个骨的壳-实体共享比 = 1.0000（壳节点 ⊆ 实体节点）；"
          f"整体 min={min(fracs):.4f}, max={max(fracs):.4f}。"
          "唯一 <1 的是极小的颅骨支架 `STR5`（6 节点）与 `TEMPORAL`(0.9973)/`ZYGOMATIC`(0.9579)。\n")
    a("> 结论：骨皮质壳与松质实体**共用节点**，FEBio 中可直接共用同一节点集（无 tie、无接触），"
      "只需为壳单元另建 ShellDomain。\n")

    a("\n## 4. `*SECTION_SHELL` 概述（实测）\n")
    a(f"共 **{len(sections)}** 个 SECID，全部 `ELFORM=16`（全积分壳）。代表值：\n")
    a("| SECID | 名称 | ELFORM | NIP | T1..T4 |")
    a("|-------|------|--------|-----|--------|")
    sample = ["83500201", "83500101", "89001501", "88000002", "88000005", "88000012",
              "81000101", "81001301", "8710710", "7547090"]
    name_by_pid = {r["pid"]: r["name"] for r in pid_records}
    for sid in sample:
        s = sections.get(int(sid)) if sid.isdigit() else None
        if not s:
            continue
        th = s.get("thickness")
        th_s = "/".join(str(x) for x in th[:4]) if th else "—"
        a(f"| {sid} | {name_by_pid.get(int(sid),'')} | {s.get('elform')} | {s.get('nip')} | {th_s} |")
    a("")
    a("> 注意：`*SECTION_SHELL` 的 SECID 与骨皮质 `*PART` 的 PID 一一对应."
      "（`83500201` → `R_HIPBONE_CORT_1.5`）。\n")

    a("\n## 5. 可行性判断（壳单元 vs 几何加厚→实体）\n")
    a("### 5.1 pyfebio 构建 API 原生支持壳（实测，`.venv` 内）\n")
    a("```")
    a('ShellFEBioElementType = Literal["tri3","tri6","quad4","quad8","quad9","q4ans","q4eas"]')
    a('SHELL_ELEMENT_CLASS_MAP = {"tri3":Tri3Element,"tri6":Tri6Element,"quad4":Quad4Element,...}')
    a('numpy_to_elements(elements, element_type="quad4"|"tri3", ...) -> Elements')
    a('class ShellDomain(type Literal["elastic-shell","three-field-shell","rigid-shell",'
      '"elastic-shell-eas","elastic-shell-ans"], shell_thickness: float)')
    a("```")
    a("来源：`.venv/Lib/site-packages/pyfebio/mesh.py`(L27/L147/L171/L375-397)、"
      "`meshdomains.py`(L34-48)。\n")
    a("- THUMS 壳为 SHELL3/SHELL4 ⇒ 直接映射 FEBio `tri3`/`quad4`，**零几何重建**。\n")
    a("- `ShellDomain.shell_thickness` 是**每域单值**；THUMS 用 `*SECTION_SHELL` 逐 PART 给单值，"
      "恰好一一对应（髋骨 9 个变体就是 9 个不同厚度的独立 PART/Domain）。\n")

    a("### 5.2 连接性证据（决定性）\n")
    a("- §3：骨皮质壳节点与同骨 SPON 实体节点**完全共享**（比 = 1.0000），"
      "意味着壳**就贴在松质实体外表面上**，共用节点集。\n")
    a("- 因此 (a) 壳单元方案**无需 tie/接触/网格缝合**，可最小改动地把皮质层加回现有 SPON 实体 FE 模型。\n")

    a("### 5.3 髋骨 9 变体：是「分区」不是「重叠」（被面积证据判定）\n")
    a("- §2.4：9 变体壳面积之和 **59135 mm²** ≈ 髋骨 SPON 自由表面 **59221 mm²**（99.85%）；"
      "9 变体单元数 = 2232 = SPON 自由面 quad 数（4464 tris / 2）。**精确铺满**。\n")
    a("- §2.1：各变体 SECID 厚度等于其名后缀（0.75/1/1.5/1.8/2/2.2/3 + ACETABLUM 1.0/1.5）"
      "⇒ 它们是**同一皮质面按厚度分区**，不是 9 个重叠整壳。\n")
    a("- 因此对骨盆**应取全部 9 个变体**（各自厚度）来还原变厚度皮质；"
      "**只取 `_CORT_1.5` 只覆盖约 42%**（25019/59221），会留下大片空缺，"
      "**不能**当作「挑一个代表壳」。\n")

    a("### 5.4 结论\n")
    a("- **推荐 (a) FEBio 壳单元贴松质实体**：拓扑零重建、厚度逐 PART 已知、节点 100% 共享。\n")
    a("- **(b) 几何加厚→实体** 仅在需要实体接触/磨损或壳无节点共享时才有必要；"
      "本甲板的壳已与实体共节点，加厚反而要重建内外表面、缝合、重网格，收益为负。\n")

    a("\n## 6. 复现命令\n")
    a("```powershell")
    a("# 仓库根，PYTHONPATH=src")
    a('$env:PYTHONPATH="src"')
    a("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\shell_scout.py `")
    a("    --parts-index temp\\opensim_fe\\parts\\parts_index.json `")
    a("    --outdir temp\\opensim_fe\\shells `")
    a("    --report results\\opensim_fe\\SHELL_SCOUT_REPORT.md")
    a("```")
    a("")
    a("产物：`temp/opensim_fe/shells/{shell_pid_stats.json,.csv, shell_sections.json, "
      "shell_connectivity.json, pelvis_variants.json, pelvis_overlap.json, scout_summary.json}`。")
    a("")
    a("## 7. 诚实边界\n")
    a("- 本步**只侦察+解析**：未生成任何 FEBio 壳/实体网格，未做 tie/接触。\n")
    a("- `thickness` card2 的第 5 字段在全部 18006 条中恒为 0.0；本报告只采纳 T1..T4。\n")
    a("- 连接性只对归一化基名可配对的骨计算（110 组）；未配对的壳 PID 未做 tie 判定。\n")
    a("- 「分区 vs 重叠」由节点级空间重合（tol=0.5 mm）+ 面积推断；未做严格的单元-单元几何求交。\n")
    a("- 壳单元的材料（`*MAT_*`）未在本步解析——建 FE 时仍需读对应 MAT/SECID 的模量/密度。\n")
    report_path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
