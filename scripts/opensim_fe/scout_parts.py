"""scout_parts.py -- memory-safe scanner of *PART (title, pid) in a THUMS deck.

The generic ``kmesh_io.iter_keyword_blocks`` buffers every data line of the
current block, which is fatal on a 230 MB deck whose *ELEMENT_SOLID blocks hold
millions of lines.  This scanner keeps only the first two lines after each
``*PART`` (LS-DYNA *PART card = title line + PID/SECID/MID card).

Also reads the ``$HMNAME COMPS <cid> <name>`` comment immediately preceding a
``*PART`` (HyperMesh part naming), which is the authoritative human name.

Usage (from repo root, PYTHONPATH=src)::

    python scripts/opensim_fe/scout_parts.py --out temp/opensim_fe/parts/parts_index.json
    python scripts/opensim_fe/scout_parts.py --filter "TIBIA|FIBULA|FEMUR|SKULL|..."

Output: JSON list of {cid, cid_name, pid, title, raw_pid_line}, plus a printed
filtered table.  Nothing is invented: every field is read verbatim from the deck.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DECK = ROOT / "model" / "AM50_V71_Occupant" / "THUMS_model" / "THUMS_AM50_V71_Occupant_202411.k"

# e.g. "$HMNAME COMPS81000700R_TIBIA_CORT"  or  "$HMNAME COMPS 81000700 R_TIBIA_CORT"
_HMNAME_RE = re.compile(r"^\$HMNAME\s+COMPS\s*(\d+)\s*(.*)$", re.IGNORECASE)
_FIXED = 10  # *PART PID card fixed width


def _first_int(line: str) -> int | None:
    s = line[:10].strip()
    try:
        return int(s)
    except ValueError:
        # tolerate comma / whitespace separated variants
        toks = re.split(r"[,\s]+", line.strip())
        if toks and re.fullmatch(r"\d+", toks[0]):
            return int(toks[0])
    return None


def scan_parts(deck: Path):
    """Yield dicts for every *PART block. Memory-safe single pass."""
    pending_hmname: tuple[str, str] | None = None
    cur: str | None = None
    part_lines: list[str] = []

    with open(deck, "r", encoding="latin-1", errors="replace") as fh:
        for raw in fh:
            st = raw.rstrip("\r\n")
            s = st.strip()
            if not s:
                continue
            if s[0] == "$":
                m = _HMNAME_RE.match(s)
                if m:
                    pending_hmname = (m.group(1), m.group(2).strip())
                continue
            if s[0] == "*":
                # flush previous PART
                if cur == "PART":
                    yield pending_hmname, list(part_lines)
                cur = s[1:].strip().split()[0].upper() if len(s) > 1 else ""
                part_lines = []
                if cur != "PART":
                    pending_hmname = None
                continue
            if cur == "PART" and len(part_lines) < 2:
                part_lines.append(st)
        if cur == "PART":
            yield pending_hmname, list(part_lines)


def _record(hmname, lines) -> dict:
    title = lines[0].strip() if lines else ""
    pid_line = lines[1].strip() if len(lines) > 1 else ""
    pid = _first_int(pid_line) if pid_line else None
    cid = hmname[0] if hmname else None
    cid_name = hmname[1] if hmname else None
    if pid is None and cid is not None and re.fullmatch(r"\d+", cid):
        pid = int(cid)  # HMNAME cid mirrors PID for these decks (verify in report)
    return {
        "cid": cid,
        "cid_name": cid_name,
        "pid": pid,
        "title": title,
        "pid_line": pid_line,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--deck", default=str(DECK))
    ap.add_argument("--out", default=None, help="write full index JSON here")
    ap.add_argument("--filter", default=None, help="regex over cid_name/title")
    args = ap.parse_args()

    deck = Path(args.deck)
    if not deck.is_file():
        print(f"ERROR: deck not found: {deck}", file=sys.stderr)
        return 2

    recs = []
    for hmname, lines in scan_parts(deck):
        r = _record(hmname, lines)
        if r["pid"] is not None:
            recs.append(r)
    print(f"*PART blocks parsed: {len(recs)}", file=sys.stderr)

    if args.out:
        outp = Path(args.out)
        outp.parent.mkdir(parents=True, exist_ok=True)
        outp.write_text(json.dumps(recs, ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"wrote {outp}", file=sys.stderr)

    if args.filter:
        rx = re.compile(args.filter, re.IGNORECASE)
        rows = [r for r in recs if rx.search(r["cid_name"] or "") or rx.search(r["title"] or "")]
        print(f"matched {len(rows)}:", file=sys.stderr)
        for r in rows:
            print(f"{r['pid']:>10}  {r['cid_name'] or r['title']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
