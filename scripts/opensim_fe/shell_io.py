"""shell_io.py -- streaming LS-DYNA *ELEMENT_SHELL / *SECTION_SHELL reader.

Companion to :mod:`kmesh_io` (which only understands ``*ELEMENT_SOLID``).
THUMS AM50 V7.1 stores the cortical shells of the axial skeleton and skull as
shell elements, so the solid-only reader silently loses them.  This module adds
a memory-safe, single-pass reader for:

    *ELEMENT_SHELL               I8 fixed-width  EID PID N1 N2 N3 [N4]
    *ELEMENT_SHELL_THICKNESS     card1 = the same I8 connectivity,
                                 card2 = 5 x F16 nodal thicknesses (T1..T4 + pad)
    *SECTION_SHELL               card1 = 8 x I10  SECID ELFORM SHRF NIP PROPT QR ICOMP SETYP
                                 card2 = 8 x I10  T1 T2 T3 T4 ...

Verified deck quirks (THUMS AM50 V7.1, 202411):

    * Shell connectivity lines are 8-char fixed width and, unlike the solid
      deck, are written as a single packed run of 8-char fields with no
      separators, e.g.::

        89001556 89000101 89059586 89059688 89059651 89059651
        => EID=89001556 PID=89000101 N1..N4

      Therefore the raw line MUST be sliced at width 8 *before* ``strip()``
      (stripping drops/reflows the column alignment and can concatenate
      pid+n1 into a bogus 15-digit integer).
    * ``*ELEMENT_SHELL_THICKNESS`` appears ~1172x (one per cortical shell part)
      and overrides the per-element thickness; ``*SECTION_SHELL`` appears only
      3x and does not cover the bone shells.
    * card2 carries 5 numbers (T1 T2 T3 T4 + a trailing 0.0 in every sample).

Public API
----------
iter_shell_records(paths)     -> ("SHELL", eid, pid, (nids))
                                 ("THIC", eid, (floats))        # element card2
                                 ("SECTION", secid, card1_dict)
                                 ("SECTION_T", secid, (floats)) # section card2
iter_shell_elements(paths)    -> (eid, pid, [nids])
iter_shell_thickness(paths)   -> (eid, [floats])   # from ELEMENT_SHELL_THICKNESS
iter_section_shell(paths)     -> (secid, dict)     # card1 + card2 merged

Coordinates are NOT read here; combine with :func:`kmesh_io.iter_nodes` when the
node positions are needed.  Nothing is invented: thicknesses are returned
verbatim from the deck.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterator

__all__ = [
    "iter_shell_records",
    "iter_shell_elements",
    "iter_shell_thickness",
    "iter_section_shell",
]

_KW_RE = re.compile(r"^\*+\s*([A-Za-z0-9_\-]+)")


def _keyword(line: str) -> str | None:
    m = _KW_RE.match(line.lstrip())
    return m.group(1).upper() if m else None


def _as_paths(paths) -> list[Path]:
    if isinstance(paths, (str, Path)):
        return [Path(paths)]
    return [Path(p) for p in paths]


def _fixed_ints(text: str, width: int) -> list[int] | None:
    """Slice a fixed-width numeric record; None if any non-empty field is bad."""
    out: list[int] = []
    end = len(text) - (len(text) % width)
    for i in range(0, end, width):
        f = text[i : i + width].strip()
        if not f:
            continue
        try:
            out.append(int(f))
        except ValueError:
            return None
    return out


def _parse_shell_element(line: str) -> tuple[int, int, list[int]] | None:
    """Parse one shell connectivity line -> (eid, pid, [nids])."""
    raw = line.rstrip("\r\n")
    if not raw.strip():
        return None

    # 1) fixed-width I8, preserving alignment (do NOT strip first)
    for cand in (raw, raw.rstrip()):
        if len(cand) >= 16 and len(cand) % 8 == 0:
            vals = _fixed_ints(cand, 8)
            if vals is not None and len(vals) >= 3:
                return vals[0], vals[1], vals[2:]

    # 2) comma / whitespace separated
    st = raw.strip()
    toks = [t for t in re.split(r"[,\s]+", st) if t]
    try:
        vals = [int(t) for t in toks]
    except ValueError:
        return None
    if len(vals) < 3:
        return None
    return vals[0], vals[1], vals[2:]


def _parse_thickness(line: str) -> list[float] | None:
    """Parse one thickness card (card2 of *ELEMENT_SHELL_THICKNESS / section)."""
    raw = line.rstrip("\r\n")
    if not raw.strip():
        return None

    st = raw.strip()
    if "," in st:
        toks = [t for t in re.split(r"[,\s]+", st) if t]
        try:
            return [float(t) for t in toks]
        except ValueError:
            return None

    # fixed-width float cards: F16 (element card2) or F10 (section card2)
    for w in (16, 10, 8):
        if len(raw) % w == 0 and len(raw) >= w:
            fields = [raw[i : i + w].strip() for i in range(0, len(raw), w)]
            fields = [f for f in fields if f]
            if not fields:
                continue
            try:
                return [float(f) for f in fields]
            except ValueError:
                continue

    toks = st.split()
    try:
        return [float(t) for t in toks]
    except ValueError:
        return None


def _to_int(v: str) -> int | None:
    try:
        return int(v)
    except ValueError:
        return None


def _to_float(v: str) -> float | None:
    try:
        return float(v)
    except ValueError:
        return None


def _section_slice(raw: str) -> list[str]:
    body = raw.rstrip("\r\n")
    if len(body) % 10 == 0 and len(body) >= 10:
        fields = [body[i : i + 10].strip() for i in range(0, len(body), 10)]
    else:
        fields = re.split(r"[,\s]+", body.strip())
    return [f for f in fields if f]


def iter_shell_records(paths) -> Iterator[tuple]:
    """Single-pass stream of shell / thickness / section records.

    Yields ``("SHELL", eid, pid, (nids))``, ``("THIC", eid, (thicknesses))``,
    ``("SECTION", secid, card1_dict)`` and ``("SECTION_T", secid, (floats))``,
    in deck order.
    """
    for path in _as_paths(paths):
        with open(path, "r", encoding="latin-1", errors="replace") as fh:
            cur: str | None = None
            want_elem = True  # pairing inside ELEMENT_SHELL_THICKNESS
            pending_eid: int | None = None
            pending_secid: int | None = None
            for raw in fh:
                if not raw or raw[0] == "$":
                    continue
                st = raw.strip()
                if not st or st[0] == "$":
                    continue
                if st[0] == "*":
                    cur = _keyword(st)
                    want_elem = True
                    pending_eid = None
                    pending_secid = None
                    continue
                if cur is None:
                    continue

                if cur == "ELEMENT_SHELL":
                    parsed = _parse_shell_element(raw)
                    if parsed is not None:
                        eid, pid, nodes = parsed
                        yield ("SHELL", eid, pid, tuple(nodes))

                elif cur.startswith("ELEMENT_SHELL_THICKNESS"):
                    if want_elem:
                        parsed = _parse_shell_element(raw)
                        if parsed is None:
                            continue
                        eid, pid, nodes = parsed
                        pending_eid = eid
                        yield ("SHELL", eid, pid, tuple(nodes))
                        want_elem = False
                    else:
                        thic = _parse_thickness(raw)
                        if thic is not None:
                            yield ("THIC", pending_eid, tuple(thic))
                        want_elem = True

                elif cur == "SECTION_SHELL":
                    fields = _section_slice(raw)
                    if not fields:
                        continue
                    secid = _to_int(fields[0])
                    if secid is not None:
                        pending_secid = secid
                        nums = fields[1:8]
                        yield (
                            "SECTION",
                            secid,
                            {
                                "elform": _to_int(nums[0]) if len(nums) > 0 else None,
                                "shrf": _to_float(nums[1]) if len(nums) > 1 else None,
                                "nip": _to_int(nums[2]) if len(nums) > 2 else None,
                                "propt": _to_float(nums[3]) if len(nums) > 3 else None,
                                "raw": " ".join(fields),
                            },
                        )
                    else:
                        thic = _parse_thickness(raw)
                        if thic is not None:
                            yield ("SECTION_T", pending_secid, tuple(thic))


def iter_shell_elements(paths) -> Iterator[tuple[int, int, list[int]]]:
    """Yield ``(eid, pid, [node ids])`` for every shell element (both cards)."""
    for rec in iter_shell_records(paths):
        if rec[0] == "SHELL":
            yield rec[1], rec[2], list(rec[3])


def iter_shell_thickness(paths) -> Iterator[tuple[int, list[float]]]:
    """Yield ``(eid, [thickness, ...])`` for *ELEMENT_SHELL_THICKNESS card2."""
    for rec in iter_shell_records(paths):
        if rec[0] == "THIC":
            yield rec[1], list(rec[2])


def iter_section_shell(paths) -> Iterator[tuple[int, dict]]:
    """Yield ``(secid, fields)`` for every *SECTION_SHELL (card1 + card2)."""
    cards: dict[int, dict] = {}
    order: list[int] = []
    for rec in iter_shell_records(paths):
        if rec[0] == "SECTION":
            secid, fields = rec[1], dict(rec[2])
            if secid not in cards:
                order.append(secid)
            cards.setdefault(secid, {}).update(fields)
        elif rec[0] == "SECTION_T" and rec[1] is not None:
            cards.setdefault(rec[1], {})["thickness"] = list(rec[2])
    for secid in order:
        yield secid, cards[secid]
