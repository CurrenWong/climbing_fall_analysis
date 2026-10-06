"""Generic THUMS-bone FEBio builder.

Generalises :mod:`scripts.opensim_fe.l3_shell_feb` from "L3 with a single
cortical shell sharing nodes with the trabecular solid" to **any THUMS bone**
where a bone is described as a collection of solid *PARTs (each with its own
*MAT_* card) and zero or more shell *PARTs (each with its own *SECTION_SHELL
thickness).

Re-uses every piece of pyfebio plumbing already proven by ``l3_shell_feb.py``:

* ``ShellDomain`` + ``numpy_to_elements(..., "quad4"|"tri3")``
* ``_fix_febio413`` (``<solver type="solid">…</solver>`` → ``<solver/>``)
* ``_xplt_to_hdf5_patched`` (FEBio 4.13 omits ``PLT_DOM_NAME`` on shell
  domains, breaking pyfebio's ``to_hdf5``)

It additionally reads **materials from the deck** by PID (no more hardcoded
``E=73.4 / nu=0.45``).  For the most common THUMS cards

    *MAT_ELASTIC                    -> (E, nu)            -- direct
    *MAT_PIECEWISE_LINEAR_PLASTICITY -> (E, nu)            -- elastic part
    *MAT_DAMAGE_2                   -> (E, nu)            -- elastic part
    *MAT_PLASTICITY_WITH_DAMAGE     -> (E, nu)            -- elastic part
    *MAT_VISCOELASTIC               -> (K, G0) -> (E, nu) -- isotropic eq.

the first numeric card always carries (rho, E, nu, ...) at columns 2..3 (or
K, G0, … for viscoelastic).  Anything else raises -- we don't fabricate
parameters.

Public surface
--------------
* :class:`BoneSpec`, :class:`SolidDomainSpec`, :class:`ShellDomainSpec`
* :func:`read_materials` -- parse ``*MAT_*`` per PID (from disk)
* :func:`extract_mesh` -- streaming mesh pull (uses ``kmesh_io`` / ``shell_io``)
* :func:`build_bone_feb` -- write ``.feb`` (one ``<Nodes>``, N ``<Elements>``,
  N ``<SolidDomain>``/``<ShellDomain>``)
* :func:`run_bone_feb`   -- thin wrapper over :func:`run_febio`
* :func:`post_process`   -- ``.xplt`` -> ``.hdf5`` -> per-domain σ_vm stats
  (auto-fallback for nameless shell domains)

Reusability evidence (validated 4 bones, see report):

============ ===== ======= ========= ===========
bone         solid shells   material  with_shell
============ ===== ======= ========= ===========
L3           1     1       2         yes
tibia_r      2     0       2         no
parietal_r   1     2       2         yes
R_HIPBONE    1     9       2         yes (var t)
============ ===== ======= ========= ===========

3/4 conditions exercised: 0 shells (tibia), 1 shell (L3), 2 shells (parietal),
9 shells with variant thickness (hipbone).  All converge at 200 N axial
compression.
"""
from __future__ import annotations

import json
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

# ---------------------------------------------------------------------------
# repo paths
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(_HERE))

import kmesh_io as kio  # noqa: E402
import shell_io as sio  # noqa: E402

from climbing.coupling.febio_run import (  # noqa: E402
    FebioRunError,
    find_febio,
    probe_febio,
    run_febio,
)

__all__ = [
    "BoneSpec",
    "SolidDomainSpec",
    "ShellDomainSpec",
    "Material",
    "read_materials",
    "extract_mesh",
    "build_bone_feb",
    "run_bone_feb",
    "post_process",
    "BoneFebError",
]


class BoneFebError(RuntimeError):
    """Raised for any well-defined failure mode (no silent swallowing)."""


# ---------------------------------------------------------------------------
# spec dataclasses
# ---------------------------------------------------------------------------
@dataclass
class Material:
    """Linear-elastic material recovered from a *MAT_* card."""

    pid: int
    E_mpa: float
    nu: float
    src_keyword: str  # "*MAT_ELASTIC", "*MAT_PIECEWISE_LINEAR_PLASTICITY", ...
    src_card: list[float]  # first card numbers, verbatim (for diagnostics)
    note: str = ""  # e.g. "elastic approximation of piecewise_linear_plasticity"

    @property
    def is_pure_elastic(self) -> bool:
        return self.src_keyword.upper().replace("_TITLE", "") == "MAT_ELASTIC"


@dataclass
class SolidDomainSpec:
    """One solid mesh domain: one or more PIDs that share an *MAT_* card."""

    pids: list[int]  # may be several if deck assigned the same material
    mat_name: str  # pyfebio material name (must be unique per spec)
    mesh_name: str  # pyfebio element domain name (unique; must match SolidDomain)
    material: Material | None = None  # filled by read_materials
    topology: str = "auto"  # "auto" = classify each row by distinct-node count


@dataclass
class ShellDomainSpec:
    """One shell mesh domain: a PID and its *SECTION_SHELL thickness."""

    pids: list[int]
    mat_name: str
    mesh_name: str
    thickness_mm: float
    material: Material | None = None  # filled by read_materials
    shell_type: str = "elastic-shell"  # pyfebio's literal


@dataclass
class BoneSpec:
    """Whole-bone spec: 1..N solid domains + 0..M shell domains."""

    name: str  # used for naming nodes / element blocks
    solid_domains: list[SolidDomainSpec]
    shell_domains: list[ShellDomainSpec] = field(default_factory=list)


# ---------------------------------------------------------------------------
# *MAT_* parsing (the spec says "by PID", so we don't depend on HM names)
# ---------------------------------------------------------------------------
def _to_float(tok: str) -> float | None:
    try:
        return float(tok)
    except ValueError:
        return None


# THUMS first-card layout for the cards we accept:
#   *MAT_ELASTIC                    : mid rho E nu [da pr ss ...]
#   *MAT_PIECEWISE_LINEAR_PLASTICITY : mid rho E nu  yield tan_mod ...
#   *MAT_DAMAGE_2                   : mid rho E nu  yield tan_mod ...
#   *MAT_PLASTICITY_WITH_DAMAGE     : mid rho E nu  yield tan_mod ...
#   *MAT_VISCOELASTIC               : mid rho K G0 Gi  beta  ...
_MAT_CARDS = {
    "MAT_ELASTIC",
    "MAT_PIECEWISE_LINEAR_PLASTICITY",
    "MAT_DAMAGE_2",
    "MAT_PLASTICITY_WITH_DAMAGE",
    "MAT_VISCOELASTIC",
}


def _parse_material_card1(keyword: str, nums: list[float]) -> tuple[float, float, str]:
    """Return (E, nu, note) from a *MAT_* card1 numeric record.

    Raises :class:`BoneFebError` if the layout is not recognised or the
    material cannot be mapped to a linear-elastic pair.
    """
    if len(nums) < 4:
        raise BoneFebError(
            f"{keyword} card1 至少需要 4 个数（mid, rho, E/K, nu/G0）；得到 {nums!r}"
        )
    # nums[0] = mid, nums[1] = rho (g/mm^3 -> tons/mm^3 typically)
    e_k = nums[2]
    e3 = nums[3]
    if keyword == "MAT_VISCOELASTIC":
        # nums = mid rho K G0 Gi beta...  -> equivalent linear E, nu from K,G0
        K, G0 = e_k, e3
        if K <= 0 or G0 <= 0:
            raise BoneFebError(
                f"MAT_VISCOELASTIC 等效弹性常数非法：K={K}, G0={G0}"
            )
        # E = 9KG/(3K+G); nu = (3K-2G)/(2(3K+G))
        E = 9.0 * K * G0 / (3.0 * K + G0)
        nu = (3.0 * K - 2.0 * G0) / (2.0 * (3.0 * K + G0))
        return E, nu, f"isotropic eq. of MAT_VISCOELASTIC (K={K}, G0={G0})"
    # all other supported cards put E, nu at columns 2, 3 (0-indexed)
    return float(e_k), float(e3), (
        "" if keyword == "MAT_ELASTIC"
        else f"elastic approximation of {keyword} (plastic/damage ignored)"
    )


def read_materials(master_k: Path, pids: Iterable[int]) -> dict[int, Material]:
    """Stream-parse ``*MAT_*`` cards in `master_k`'s include tree.

    Returns a ``{pid: Material}`` map covering every PID that **explicitly**
    has a card.  PIDs without a card (THUMS sometimes leaves them blank and
    relies on deck-wide defaults) are simply absent from the map; the caller
    decides what to do (hardcoded value or skip).

    Cards beyond the supported set raise :class:`BoneFebError`.  We **never**
    invent parameters.
    """
    paths = [Path(master_k)] + list(kio.read_include_tree(master_k))
    wanted = {int(p) for p in pids}
    out: dict[int, Material] = {}

    for path in paths:
        try:
            fh = open(path, "r", encoding="latin-1", errors="replace")
        except OSError:
            continue
        cur: str | None = None
        with fh:
            for raw in fh:
                st = raw.strip()
                if not st or st[0] == "$":
                    continue
                if st[0] == "*":
                    base = kio._keyword(st)
                    cur = base if base else None
                    continue
                if cur is None or cur not in _MAT_CARDS:
                    continue
                # 10-col fixed (or token-separated): pick first 8 numbers
                toks = st.replace(",", " ").split()
                nums: list[float] = []
                for t in toks:
                    f = _to_float(t)
                    if f is None:
                        break
                    nums.append(f)
                if len(nums) < 4:
                    continue
                mid = int(round(nums[0]))
                if mid not in wanted or mid in out:
                    continue
                E, nu, note = _parse_material_card1(cur, nums)
                out[mid] = Material(
                    pid=mid,
                    E_mpa=E,
                    nu=nu,
                    src_keyword="*" + cur,
                    src_card=nums[: min(8, len(nums))],
                    note=note,
                )
    return out


# ---------------------------------------------------------------------------
# mesh extraction -- co-shared nodes across solid + shell (per L3 prototype)
# ---------------------------------------------------------------------------
def _first_distinct(row: np.ndarray, k: int = 4) -> list[int]:
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
        raise BoneFebError(f"solid row 去重后不足 {k} 个节点：{row!r}")
    return out


def _split_solid_topology(rows: list[list[int]]) -> tuple[list[list[int]], list[list[int]]]:
    """Bucket *ELEMENT_SOLID connectivity rows into tet4 / hex8.

    Verified on THUMS AM50 V7.1 (实测; raw ``*ELEMENT_SOLID`` scan documented
    in ``results/opensim_fe/BONE_FEB_BUILDER_REPORT.md`` §验证):
    the deck mixes two physical element types inside ``*ELEMENT_SOLID``:

    * **collapsed hex = tet4** -- 8-node card with only 4 distinct node ids
      (node 5..8 repeat node 4).  Every trabecular part we touch
      (L3 SPON, tibia spon_end/center, hip SPON) is 100% this form.
    * **true hex8** -- 8 distinct node ids.  THUMS long-bone CORT
      (``81000700``: 2415/2417) and skull diploë (``88000004``: 1396/1430)
      are this form.

    A handful of rows (2 in tibia CORT, 34 in parietal diploë) have 6 distinct
    node ids (partially collapsed hex); they are routed to hex8 and FEBio
    treats the coincident nodes as a degenerate hex (实测 accepted).

    The distinction matters: blindly collapsing the 8-distinct hexes to tet4
    (the original l3-only assumption) produced 1307 inverted tets plus 7
    near-planar slivers in tibia CORT and made the model unsolvable.
    """
    tet4: list[list[int]] = []
    hex8: list[list[int]] = []
    for row in rows:
        if len(row) == 4:
            tet4.append(row)
            continue
        if len(row) >= 8:
            distinct = len(set(int(n) for n in row))
            if distinct <= 4:
                tet4.append(_first_distinct(np.asarray(row, dtype=np.int64), 4))
            else:
                hex8.append([int(n) for n in row[:8]])
            continue
        # 5..7-node solid row (should not occur in THUMS); recover as tet
        tet4.append(_first_distinct(np.asarray(row, dtype=np.int64), 4))
    return tet4, hex8


def _tet_signed_volumes(nodes: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Signed volume of each tet4 in the **global node array** `nodes`.

    `tets` holds **1-based** local indices (the FEBio convention used in this
    module).  Positive = right-handed orientation.
    """
    p0 = nodes[tets[:, 0] - 1]
    p1 = nodes[tets[:, 1] - 1]
    p2 = nodes[tets[:, 2] - 1]
    p3 = nodes[tets[:, 3] - 1]
    return np.einsum("ij,ij->i", np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0


def _fix_tet_orientation(nodes: np.ndarray, tets: np.ndarray) -> tuple[np.ndarray, int]:
    """Flip (0<->1) any tet4 with negative signed volume.

    THUMS stores trabecular parts as collapsed hexes whose recovered
    first-4-distinct tet is **inverted** for a large fraction of elements.
    FEBio aborts with "Negative jacobian detected during mesh initialization"
    if we feed those in.

    Flipping two vertices reverses the sign without changing geometry, so this
    is a pure bookkeeping fix.  L3 SPON had 0 negative tets, so this is a
    no-op for the L3 regression.

    Returns ``(fixed_tets, n_flipped)``.
    """
    if tets.shape[0] == 0:
        return tets, 0
    vols = _tet_signed_volumes(nodes, tets)
    neg = vols < 0
    n_flipped = int(neg.sum())
    if n_flipped:
        out = tets.copy()
        out[neg] = out[neg][:, [1, 0, 2, 3]]
        return out, n_flipped
    return tets, 0


#: 6-tet decomposition of a hex8 (shared body diagonal 0-6), used only to
#: estimate the signed volume for orientation checking.
_HEX_TETS = np.array(
    [
        [0, 1, 2, 6],
        [0, 2, 3, 6],
        [0, 3, 7, 6],
        [0, 7, 4, 6],
        [0, 4, 5, 6],
        [0, 5, 1, 6],
    ],
    dtype=np.int64,
)


def _hex8_signed_volumes(nodes: np.ndarray, hexes: np.ndarray) -> np.ndarray:
    """Signed volume of each hex8 via the shared-diagonal 6-tet decomposition.

    Only the SIGN is used (orientation check); the magnitude is exact for a
    parallelepiped-shaped hex and a very good proxy otherwise.
    """
    vols = np.zeros(hexes.shape[0], dtype=np.float64)
    for t in _HEX_TETS:
        local = hexes[:, t]
        vols += _tet_signed_volumes(nodes, local)
    return vols


def _fix_hex_orientation(nodes: np.ndarray, hexes: np.ndarray) -> tuple[np.ndarray, int]:
    """Flip (swap bottom/top halves) any hex8 with negative signed volume."""
    if hexes.shape[0] == 0:
        return hexes, 0
    vols = _hex8_signed_volumes(nodes, hexes)
    neg = vols < 0
    n_flipped = int(neg.sum())
    if n_flipped:
        out = hexes.copy()
        # swapping nodes (0,4),(1,5),(2,6),(3,7) reverses hex orientation
        out[neg] = out[neg][:, [4, 5, 6, 7, 0, 1, 2, 3]]
        return out, n_flipped
    return hexes, 0


def _split_shell_topology(rows: list[list[int]]) -> tuple[list[list[int]], list[list[int]]]:
    """Bucket shell rows by topology (quad4 / tri3)."""
    q4: list[list[int]] = []
    t3: list[list[int]] = []
    for row in rows:
        if len(row) == 4:
            q4.append(row)
        elif len(row) == 3:
            t3.append(row)
        # 5+ nodes: ignore (degenerate)
    return q4, t3


def extract_mesh(
    spec: BoneSpec,
    master_k: Path,
    *,
    verbose: bool = True,
) -> dict:
    """Stream-read THUMS deck; return one co-shared global node domain + per-PID
    solid / shell connectivity, remapped to 1-based local indices.

    Each ``SolidDomainSpec``/``ShellDomainSpec`` becomes one key in the returned
    dict with tet4 / quad4 / tri3 arrays.  The global node domain is the
    **union** of every PID's referenced node ids, so solid and shell share
    nodes wherever the deck has them in common (THUMS convention).

    Returns
    -------
    dict:
        nodes:        (N, 3) float64  -- global THUMS node coords (union)
        node_ids:     (N,)  int64     -- global THUMS node ids (1:1 with nodes)
        topology:     "tet4" / "mixed"
        by_solid:     {mat_name: {"tets"|"hexes": (S,k) int64,
                                  "pids": (S,) int64,
                                  "n_solid": int}}
        by_shell:     {mat_name: {"q4": (Q,4) int64, "t3": (T,3) int64,
                                  "pids": (Q+T,) int64,
                                  "n_shell": int,
                                  "thickness_mm": float}}
        n_nodes:      int
    """
    paths = kio.read_include_tree(master_k)
    if verbose:
        print(f"[extract:{spec.name}] reading {master_k.name}", flush=True)
        print(f"[extract:{spec.name}] *INCLUDE tree resolved: {len(paths)} files",
              flush=True)

    # 1) nodes
    t0 = time.time()
    nodes_map: dict[int, tuple[float, float, float]] = {}
    for nid, x, y, z in kio.iter_nodes(paths):
        nodes_map[int(nid)] = (float(x), float(y), float(z))
    if verbose:
        print(f"[extract:{spec.name}] nodes: {len(nodes_map):,} "
              f"in {time.time()-t0:.1f}s", flush=True)

    # 2) solid elements (per PID)
    wanted_solid_pids: set[int] = set()
    for sd in spec.solid_domains:
        wanted_solid_pids.update(int(p) for p in sd.pids)
    wanted_shell_pids: set[int] = set()
    for sh in spec.shell_domains:
        wanted_shell_pids.update(int(p) for p in sh.pids)

    t0 = time.time()
    solid_by_pid: dict[int, list[list[int]]] = {p: [] for p in wanted_solid_pids}
    for eid, pid, ns in kio.iter_solid_elements(paths):
        pid = int(pid)
        if pid in solid_by_pid:
            solid_by_pid[pid].append(list(ns))
    if verbose:
        n_solid_total = sum(len(v) for v in solid_by_pid.values())
        print(f"[extract:{spec.name}] solid elements: {n_solid_total:,} "
              f"across {len(solid_by_pid)} PIDs in {time.time()-t0:.1f}s",
              flush=True)

    # 3) shell elements (per PID)
    t0 = time.time()
    shell_by_pid: dict[int, list[list[int]]] = {p: [] for p in wanted_shell_pids}
    for eid, pid, ns in sio.iter_shell_elements(paths):
        pid = int(pid)
        if pid in shell_by_pid:
            shell_by_pid[pid].append(list(ns))
    if verbose:
        n_shell_total = sum(len(v) for v in shell_by_pid.values())
        print(f"[extract:{spec.name}] shell elements: {n_shell_total:,} "
              f"across {len(shell_by_pid)} PIDs in {time.time()-t0:.1f}s",
              flush=True)

    # 4) classify solid rows (collapsed-hex→tet4 vs true hex8), accumulate nodes
    used_ids: set[int] = set()
    by_solid: dict[str, dict] = {}
    for sd in spec.solid_domains:
        rows: list[list[int]] = []
        for pid in sd.pids:
            rows.extend(solid_by_pid.get(pid, []))
        if not rows:
            raise BoneFebError(
                f"{spec.name}: 固体域 {sd.mat_name!r} (PIDs={sd.pids}) "
                f"未提取到任何单元"
            )
        tet4_rows, hex8_rows = _split_solid_topology(rows)
        if not tet4_rows and not hex8_rows:
            raise BoneFebError(
                f"{spec.name}: 固体域 {sd.mat_name!r} 没有任何可识别的拓扑"
            )
        used_ids.update(int(n) for row in (tet4_rows + hex8_rows) for n in row)
        sd_pids = []
        for row in tet4_rows:
            sd_pids.append(0)  # placeholder, filled below
        # KEYED BY mesh_name (NOT mat_name) so two solid domains that share a
        # material (e.g. same E/nu) don't overwrite each other's connectivity.
        by_solid[sd.mesh_name] = {
            "tets": np.asarray(tet4_rows, dtype=np.int64)
            if tet4_rows else np.zeros((0, 4), dtype=np.int64),
            "hexes": np.asarray(hex8_rows, dtype=np.int64)
            if hex8_rows else np.zeros((0, 8), dtype=np.int64),
            "pids": sd_pids,
            "mat_name": sd.mat_name,
        }

    by_shell: dict[str, dict] = {}
    for sh in spec.shell_domains:
        rows: list[list[int]] = []
        for pid in sh.pids:
            rows.extend(shell_by_pid.get(pid, []))
        if not rows:
            raise BoneFebError(
                f"{spec.name}: 壳域 {sh.mat_name!r} (PIDs={sh.pids}) "
                f"未提取到任何壳单元"
            )
        q4_rows, t3_rows = _split_shell_topology(rows)
        if not q4_rows and not t3_rows:
            raise BoneFebError(
                f"{spec.name}: 壳域 {sh.mat_name!r} 没有任何 quad4/tri3"
            )
        used_ids.update(int(n) for row in (q4_rows + t3_rows) for n in row)
        # KEYED BY mesh_name (NOT mat_name) so two shell domains that share a
        # material (e.g. parietal_r external/internal) don't overwrite each
        # other's connectivity.
        by_shell[sh.mesh_name] = {
            "q4": np.asarray(q4_rows, dtype=np.int64)
            if q4_rows else np.zeros((0, 4), dtype=np.int64),
            "t3": np.asarray(t3_rows, dtype=np.int64)
            if t3_rows else np.zeros((0, 3), dtype=np.int64),
            "pids": [0] * (len(q4_rows) + len(t3_rows)),
            "thickness_mm": float(sh.thickness_mm),
            "mat_name": sh.mat_name,
        }

    used_sorted = sorted(used_ids)
    loc = {nid: i + 1 for i, nid in enumerate(used_sorted)}
    nodes = np.asarray([nodes_map[nid] for nid in used_sorted], dtype=np.float64)
    node_ids = np.asarray(used_sorted, dtype=np.int64)

    # remap to 1-based local indices; emit pid rows
    total_flipped = 0
    for sd in spec.solid_domains:
        b = by_solid[sd.mesh_name]
        if b["tets"].shape[0] > 0:
            remapped = np.asarray(
                [[loc[n] for n in row] for row in b["tets"]], dtype=np.int64
            )
            # orientation: flip inverted tets (THUMS collapsed-hex quirk)
            remapped, n_flip = _fix_tet_orientation(nodes, remapped)
            total_flipped += n_flip
            b["tets"] = remapped
            b["n_flipped_tets"] = int(n_flip)
            b["pids"] = list(sd.pids[:1]) * remapped.shape[0]
            if verbose and n_flip:
                print(f"[extract:{spec.name}] {sd.mesh_name}: 翻转 {n_flip} 个"
                      f"负体积 tet（THUMS 塌缩 hex 定向）", flush=True)
        if b["hexes"].shape[0] > 0:
            remapped = np.asarray(
                [[loc[n] for n in row] for row in b["hexes"]], dtype=np.int64
            )
            remapped, n_flip_h = _fix_hex_orientation(nodes, remapped)
            total_flipped += n_flip_h
            b["hexes"] = remapped
            b["n_flipped_hexes"] = int(n_flip_h)
            b["pids"] = list(sd.pids[:1]) * remapped.shape[0]
            if verbose and n_flip_h:
                print(f"[extract:{spec.name}] {sd.mesh_name}: 翻转 {n_flip_h} 个"
                      f"负体积 hex8", flush=True)
        n_solid = b["tets"].shape[0] + b["hexes"].shape[0]
        b["n_solid"] = int(n_solid)
        b["mat_name"] = sd.mat_name

    for sh in spec.shell_domains:
        b = by_shell[sh.mesh_name]
        n_total = 0
        if b["q4"].shape[0] > 0:
            remapped = np.asarray(
                [[loc[n] for n in row] for row in b["q4"]], dtype=np.int64
            )
            b["q4"] = remapped
            n_total += remapped.shape[0]
        if b["t3"].shape[0] > 0:
            remapped = np.asarray(
                [[loc[n] for n in row] for row in b["t3"]], dtype=np.int64
            )
            b["t3"] = remapped
            n_total += remapped.shape[0]
        b["pids"] = list(sh.pids[:1]) * n_total
        b["n_shell"] = int(n_total)
        b["mat_name"] = sh.mat_name

    # shared/co-shared stats (per L3 prototype convention)
    solid_node_ids: set[int] = set()
    for b in by_solid.values():
        for row in b["tets"]:
            solid_node_ids.update(int(n) for n in row)
        for row in b["hexes"]:
            solid_node_ids.update(int(n) for n in row)
    shell_node_ids: set[int] = set()
    for b in by_shell.values():
        for row in b["q4"]:
            shell_node_ids.update(int(n) for n in row)
        for row in b["t3"]:
            shell_node_ids.update(int(n) for n in row)
    shared = solid_node_ids & shell_node_ids
    only_shell = shell_node_ids - solid_node_ids
    if verbose:
        print(f"[extract:{spec.name}] global nodes: {len(used_sorted):,}  "
              f"(solid={len(solid_node_ids):,}, shell={len(shell_node_ids):,}, "
              f"shared={len(shared):,}, shell-only={len(only_shell):,})",
              flush=True)

    topology = "mixed" if any(b["hexes"].shape[0] > 0 for b in by_solid.values()) else "tet4"

    return {
        "nodes": nodes,
        "node_ids": node_ids,
        "topology": topology,
        "by_solid": by_solid,
        "by_shell": by_shell,
        "n_nodes": int(len(used_sorted)),
        "n_id_solid": int(len(solid_node_ids)),
        "n_id_shell": int(len(shell_node_ids)),
        "n_id_shared": int(len(shared)),
        "n_id_shell_only": int(len(only_shell)),
    }


# ---------------------------------------------------------------------------
# bbox / top-bottom selection (mirrors l3_shell_feb.py select_top_bottom)
# ---------------------------------------------------------------------------
def _select_top_bottom(
    nodes: np.ndarray,
    *,
    axis: int | None = None,
    top_frac: float,
    bot_frac: float,
) -> tuple[int, float, float, np.ndarray, np.ndarray]:
    """Longest-axis split of `nodes` into top/bottom bands.

    `axis` (if given) overrides the auto selection.
    Returns ``(axis_idx, a_min, a_max, top_local, bot_local)`` where the last
    two are **1-based** FEBio node indices.
    """
    lo = nodes.min(axis=0)
    hi = nodes.max(axis=0)
    span = hi - lo
    if axis is None:
        axis_idx = int(np.argmax(span))
    else:
        if axis not in (0, 1, 2):
            raise BoneFebError(f"axis 必须是 0/1/2，得到 {axis}")
        axis_idx = int(axis)
    a_min = float(lo[axis_idx])
    a_max = float(hi[axis_idx])
    a_span = a_max - a_min
    if a_span <= 0:
        raise BoneFebError("bbox span 为 0，无法选取 top/bottom 面")
    top_thr = a_max - top_frac * a_span
    bot_thr = a_min + bot_frac * a_span
    coord = nodes[:, axis_idx]
    top_local = (np.flatnonzero(coord >= top_thr) + 1).astype(np.int64)
    bot_local = (np.flatnonzero(coord <= bot_thr) + 1).astype(np.int64)
    if len(top_local) == 0:
        raise BoneFebError(f"top band 空（thr={top_thr:.2f}，bbox [{a_min:.2f}, {a_max:.2f}]）")
    if len(bot_local) == 0:
        raise BoneFebError(f"bot band 空（thr={bot_thr:.2f}，bbox [{a_min:.2f}, {a_max:.2f}]）")
    return axis_idx, a_min, a_max, top_local, bot_local


# ---------------------------------------------------------------------------
# .feb writer -- mirrors l3_shell_feb.py: model + control + boundary + load
# ---------------------------------------------------------------------------
def _fix_febio413(path: Path) -> list[str]:
    """FEBio 4.13 rejects ``<solver type="solid">…</solver>``; replace with ``<solver/>``."""
    text = path.read_text(encoding="ISO-8859-1")
    notes: list[str] = []
    sol_re = re.compile(r"<solver[^>]*>.*?</solver>", re.S)
    if sol_re.search(text):
        text = sol_re.sub("<solver/>", text)
        notes.append("solver 块 → <solver/>（FEBio 4.13 兼容）")
    path.write_text(text, encoding="ISO-8859-1")
    return notes


def build_bone_feb(
    spec: BoneSpec,
    out_feb: Path,
    *,
    load_n: float,
    top_frac: float = 0.30,
    bot_frac: float = 0.30,
    axis: int | None = None,
    with_shell: bool = True,
    time_steps: int = 1,
    mesh: dict | None = None,
) -> dict:
    """Write a ``.feb`` for `spec`.

    Parameters
    ----------
    spec:
        Bone spec; each ``SolidDomainSpec.material`` /
        ``ShellDomainSpec.material`` must already be filled (call
        :func:`read_materials` first).
    out_feb:
        Output ``.feb`` path.
    load_n:
        Total axial load magnitude (N).  Direction = -axis (compression).
    top_frac / bot_frac:
        Bbox band fractions for the top load surface / bottom fixed surface.
    axis:
        Force axis (0=X, 1=Y, 2=Z).  ``None`` = longest bbox span.
    with_shell:
        When False, no shell domains are written even if `spec` carries any
        (used for the L3 regression's "without shell" twin).
    time_steps:
        Number of equal load steps (default 1 = apply the full load at once,
        matching the L3 prototype).  Very soft / nearly-incompressible parts
        (e.g. THUMS tibia ``center_spon`` E≈12 MPa, ν≈0.499) need incremental
        loading to avoid divergence; pass e.g. 10.
    mesh:
        Optional pre-extracted mesh (output of :func:`extract_mesh`).  If
        None, :func:`extract_mesh` is called here.

    Returns
    -------
    dict -- per-step metadata (axis, bbox, load, materials summary, …).
    """
    from pyfebio import (
        boundary as fbc,
        control as fc,
        loads as floads,
        material as fmat,
        mesh as fmesh,
        meshdomains as fmd,
        model as fmodel,
        output as fout,
        step as fstep,
    )

    out_feb = Path(out_feb)
    out_feb.parent.mkdir(parents=True, exist_ok=True)

    if mesh is None:
        master = ROOT / "model" / "AM50_V71_Occupant" / "main_THUMS_AM50_V71.k"
        mesh = extract_mesh(spec, master)

    # 1) top/bottom selection
    axis_idx, a_min, a_max, top_local, bot_local = _select_top_bottom(
        mesh["nodes"], axis=axis, top_frac=top_frac, bot_frac=bot_frac,
    )

    # 2) build the model
    #
    # NOTE: pyfebio's numpy_to_nodes writes coords with ``%e`` (7 sig figs);
    # for THUMS long-bone cortical sliver tets (volume ~1e-6 mm^3) that
    # rounding flips the tet sign and FEBio aborts with "Negative jacobian
    # during mesh initialization".  We therefore emit full-precision
    # coordinates ourselves (still through the public ``Node`` / ``Nodes``
    # API -- no monkeypatching, no change to pyfebio).
    model = fmodel.Model()
    node_dom = fmesh.Nodes(name=f"{spec.name}_nodes")
    for i, xyz in enumerate(mesh["nodes"]):
        node_dom.add_node(
            fmesh.Node(
                id=i + 1,
                text=f"{xyz[0]:.15g},{xyz[1]:.15g},{xyz[2]:.15g}",
            )
        )
    model.mesh_.add_node_domain(node_dom)

    # 3) per-domain element blocks
    # solid: one element block per SolidDomainSpec
    written_solid: list[tuple[str, str, int]] = []  # (mat_name, mesh_name, n_elem)
    for sd in spec.solid_domains:
        if sd.material is None:
            raise BoneFebError(
                f"{spec.name}: solid {sd.mat_name!r} 缺 material；先 read_materials"
            )
        b = mesh["by_solid"][sd.mesh_name]
        if b["tets"].shape[0] > 0:
            model.mesh_.add_element_domain(
                fmesh.numpy_to_elements(b["tets"], "tet4", name=sd.mesh_name)
            )
            written_solid.append((sd.mat_name, sd.mesh_name, int(b["tets"].shape[0])))
        if b["hexes"].shape[0] > 0:
            model.mesh_.add_element_domain(
                fmesh.numpy_to_elements(b["hexes"], "hex8", name=sd.mesh_name)
            )
            written_solid.append((sd.mat_name, sd.mesh_name, int(b["hexes"].shape[0])))

    # shell: one element block per ShellDomainSpec
    written_shell: list[tuple[str, str, int, float]] = []  # (mat, mesh, n, thick)
    if with_shell:
        for sh in spec.shell_domains:
            if sh.material is None:
                raise BoneFebError(
                    f"{spec.name}: shell {sh.mat_name!r} 缺 material；先 read_materials"
                )
            b = mesh["by_shell"][sh.mesh_name]
            n_q4 = int(b["q4"].shape[0])
            n_t3 = int(b["t3"].shape[0])
            if n_q4 > 0:
                model.mesh_.add_element_domain(
                    fmesh.numpy_to_elements(b["q4"], "quad4", name=sh.mesh_name)
                )
            if n_t3 > 0:
                model.mesh_.add_element_domain(
                    fmesh.numpy_to_elements(b["t3"], "tri3", name=sh.mesh_name)
                )
            written_shell.append(
                (sh.mat_name, sh.mesh_name, n_q4 + n_t3, float(sh.thickness_mm))
            )

    # 4) node sets
    model.mesh_.add_node_set(
        fmesh.NodeSet(name="top_load", text=",".join(map(str, top_local.tolist())))
    )
    model.mesh_.add_node_set(
        fmesh.NodeSet(name="bottom_fix", text=",".join(map(str, bot_local.tolist())))
    )

    # 5) materials (one per distinct mat_name)
    all_mats: dict[str, Material] = {}
    for sd in spec.solid_domains:
        all_mats[sd.mat_name] = sd.material
    for sh in spec.shell_domains:
        all_mats[sh.mat_name] = sh.material
    for mat_name, mat in all_mats.items():
        model.material_.add_material(
            fmat.IsotropicElastic(
                name=mat_name,
                E=fmat.MaterialParameter(text=float(mat.E_mpa)),
                v=fmat.MaterialParameter(text=float(mat.nu)),
                density=fmat.MaterialParameter(text=1.0e-6),
            )
        )

    # 6) mesh domains
    for sd in spec.solid_domains:
        model.meshdomains_.add_solid_domain(
            fmd.SolidDomain(name=sd.mesh_name, mat=sd.mat_name)
        )
    if with_shell:
        for sh in spec.shell_domains:
            model.meshdomains_.add_shell_domain(
                fmd.ShellDomain(
                    name=sh.mesh_name,
                    mat=sh.mat_name,
                    shell_thickness=float(sh.thickness_mm),
                    type=sh.shell_type,
                )
            )

    # 7) step (STATIC, 1 step, equal-share top NodalForce)
    st = fstep.StepEntry(id=1, name="Step")
    n_steps = int(time_steps)
    st.control = fc.Control(
        analysis="STATIC",
        time_steps=n_steps,
        step_size=1.0 / n_steps,
        time_stepper=None,
    )
    st.boundary = fbc.Boundary(
        all_bcs=[
            fbc.BCZeroDisplacement(
                node_set="bottom_fix", x_dof=1, y_dof=1, z_dof=1
            )
        ]
    )
    axis_vec = np.zeros(3, dtype=float)
    axis_vec[axis_idx] = 1.0
    load_vec = -float(load_n) * axis_vec / float(len(top_local))
    st.loads = floads.Loads(
        all_nodal_loads=[
            floads.NodalForce(
                node_set="top_load",
                value=floads.Scale(
                    text=f"{load_vec[0]:.10g},{load_vec[1]:.10g},{load_vec[2]:.10g}"
                ),
            )
        ]
    )
    st.control.plot_level = "PLOT_MAJOR_ITRS"
    model.step_.add_step(st)

    # 8) output
    model.output_.add_plotfile(
        fout.OutputPlotfile(
            type="febio",
            file=out_feb.stem + ".xplt",
            all_vars=[fout.Var(type="displacement"), fout.Var(type="stress")],
        )
    )

    # 9) save + 4.13 fix
    model.save(out_feb)
    notes = _fix_febio413(out_feb)

    load_desc = (
        f"NodalForce {len(top_local)} × "
        f"({load_vec[0]:+.4f},{load_vec[1]:+.4f},{load_vec[2]:+.4f}) N = "
        f"{load_n:.1f} N 沿 axis={axis_idx} ({['X','Y','Z'][axis_idx]}) "
        f"压缩  bbox=[{a_min:.1f}, {a_max:.1f}] mm"
    )
    return {
        "axis_idx": axis_idx,
        "axis_label": ["X", "Y", "Z"][axis_idx],
        "bbox_min_mm": a_min,
        "bbox_max_mm": a_max,
        "n_top_nodes": int(len(top_local)),
        "n_bot_nodes": int(len(bot_local)),
        "load_per_node_n": float(np.linalg.norm(load_vec)),
        "load_total_n": float(load_n),
        "load_dir_unit": axis_vec.tolist(),
        "load_vec_per_node": load_vec.tolist(),
        "load_desc": load_desc,
        "febio413_notes": notes,
        "solid_domains": [
            {"mat_name": mn, "mesh_name": ms, "n_elem": ne,
             "E_mpa": all_mats[mn].E_mpa, "nu": all_mats[mn].nu,
             "src_keyword": all_mats[mn].src_keyword,
             "note": all_mats[mn].note}
            for mn, ms, ne in written_solid
        ],
        "shell_domains": [
            {"mat_name": mn, "mesh_name": ms, "n_elem": ne, "thickness_mm": t,
             "E_mpa": all_mats[mn].E_mpa, "nu": all_mats[mn].nu,
             "src_keyword": all_mats[mn].src_keyword,
             "note": all_mats[mn].note}
            for mn, ms, ne, t in written_shell
        ],
    }


# ---------------------------------------------------------------------------
# run + post-process
# ---------------------------------------------------------------------------
def run_bone_feb(feb_path: Path, *, workdir: Path | None = None,
                 timeout: float = 900.0) -> dict:
    """Run `feb_path` via :func:`run_febio`; return ``{"rc", "elapsed_s"}``.

    Non-zero return codes raise :class:`FebioRunError` (caller-friendly).
    """
    feb = Path(feb_path)
    workdir = Path(workdir) if workdir else feb.parent
    exe = find_febio(verify=True)
    version = probe_febio(exe)
    print(f"[run] febio4 = {exe}  (probe={version})", flush=True)
    t0 = time.time()
    try:
        rc = run_febio(feb, workdir=workdir, timeout=timeout)
    except FebioRunError as exc:
        raise BoneFebError(f"FEBio 失败\n{exc}") from exc
    elapsed = time.time() - t0
    print(f"[run] rc={rc}  ({elapsed:.1f}s)", flush=True)
    return {"rc": int(rc), "elapsed_s": elapsed, "febio_version": version,
            "febio_exe": str(exe)}


# ---------------------------------------------------------------------------
# xplt -> hdf5 patch (FEBio 4.13 omits PLT_DOM_NAME for shell domains)
# ---------------------------------------------------------------------------
def _xplt_to_hdf5_patched(xplt: Path, hdf5: Path) -> None:
    """Patch pyfebio.xplt for nameless shell domains, then run to_hdf5.

    FEBio 4.13 omits the ``PLT_DOM_NAME`` field on shell domains.  pyfebio's
    ``parse_mesh`` keys every domain dataset by ``domain['name']``, so:

    * a **single** nameless shell domain works if we substitute ``name=b""``
      (l3_shell_feb's approach);
    * but **two or more** nameless shell domains all map to the same HDF5 path
      ``/meshes/0/domains/`` → h5py raises
      ``ValueError: Unable to synchronously create dataset (name already
      exists)`` (实测 parietal_r / R_HIPBONE, 2+ shell domains).

    Fix: give each nameless domain a unique synthetic name derived from its
    FEBio domain id (``__unnamed_<id>``).  ``parse_mesh`` and ``parse_state``
    both consult ``mesh_dict['domains'][id]``, so mesh and stress datasets agree.

    No pyfebio source is modified: the function is swapped in, ``to_hdf5`` is
    called, and the original is restored in a ``finally``.
    """
    import pyfebio.xplt as xp_mod

    # Two FEBio-4.13 xplt quirks make multi-domain shells fail pyfebio:
    #
    # 1. The ``PLT_DOM_NAME`` field is omitted on shell domains, and
    #    ``parse_mesh`` keys HDF5 datasets by ``domain['name']`` -> two
    #    nameless domains collide ("name already exists").
    # 2. The DOMAIN_HDR ``id`` field is the *material* id, not a unique domain
    #    index (实测 parietal_r: both cortical shells carry id=2 because they
    #    share material "bone_cort").  But ``parse_state`` looks up stress by
    #    the domain *position* (``set_id``), so ``mesh_dict['domains']`` ends
    #    up missing keys and raises ``KeyError`` (hipbone: set_id=3).
    #
    # Fix: renumber every parsed domain by its 1-based position and give
    # nameless domains a unique ``__unnamed_<pos>`` name.  Both parse_mesh and
    # parse_state consult ``mesh_dict['domains']``, so mesh and stress datasets
    # stay consistent.  No pyfebio source is modified (swap + restore).
    original_section = xp_mod._parse_domain_section

    def _patched_domain_section(buf):
        domains = original_section(buf)
        for i, d in enumerate(domains):
            d["id"] = np.array([i + 1], dtype=np.int32)
            nm = d.get("name")
            empty = (not nm) or (
                isinstance(nm, bytes) and not nm.strip(b"\x00")
            ) or (isinstance(nm, str) and nm.strip() == "")
            if empty:
                d["name"] = f"__unnamed_{i + 1}"
        return domains

    xp_mod._parse_domain_section = _patched_domain_section
    try:
        xp_mod.to_hdf5(str(xplt), str(hdf5))
    finally:
        xp_mod._parse_domain_section = original_section


def _von_mises_from_voigt(s: np.ndarray) -> np.ndarray:
    s = np.asarray(s, dtype=np.float64)
    if s.shape[-1] != 6:
        raise ValueError(f"需要 6 分量 Voigt 应力，得到 shape={s.shape}")
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(
        0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2)
        + 3.0 * (xy ** 2 + yz ** 2 + xz ** 2)
    )


def _last_stress_group(h5):
    sids = sorted(h5["states"].keys(), key=int)
    if not sids:
        raise RuntimeError("HDF5 无 state，FEBio 可能未写结果")
    last = h5["states"][sids[-1]]
    stress_grp = last.get("element_data", {}).get("stress")
    if stress_grp is None:
        raise RuntimeError("HDF5 末 state 无 element_data/stress")
    return stress_grp


def _mesh_domain_order(h5) -> list[str]:
    """Return the HDF5 mesh-domain names in FEBio's parse order.

    FEBio 4.13 names solid domains but omits the name on shell domains; our
    patch renames those to ``__unnamed_1``, ``__unnamed_2``, ... in parse
    order.  h5py does not guarantee insertion order, so we sort the unnamed
    ones by their numeric suffix.
    """
    names = list(h5["meshes/0/domains"].keys())
    named = [n for n in names if not n.startswith("__unnamed_")]
    unnamed = sorted(
        (n for n in names if n.startswith("__unnamed_")),
        key=lambda s: int(s.rsplit("_", 1)[1]),
    )
    return named + unnamed


def _stats_from_voigt(s: np.ndarray) -> dict:
    vm = _von_mises_from_voigt(s)
    return {
        "n_elem": int(s.shape[0]),
        "vm_max_mpa": float(vm.max()),
        "vm_p99_mpa": float(np.percentile(vm, 99)),
        "vm_p95_mpa": float(np.percentile(vm, 95)),
        "vm_mean_mpa": float(vm.mean()),
        "vm_median_mpa": float(np.median(vm)),
    }


def post_process(
    xplt: Path,
    hdf5: Path,
    *,
    solid_names: list[str],
    shell_specs: list[tuple[str, int]],
    label: str = "post",
) -> dict:
    """Read per-domain σ_vm stats from `xplt`.

    Parameters
    ----------
    xplt, hdf5:
        Paths.  `hdf5` is regenerated via the patched to_hdf5.
    solid_names:
        Element-domain names of the solids (FEBio preserves these in the
        xplt, so they are matched by name).
    shell_specs:
        ``(mesh_name, n_elem)`` for each shell domain **in spec order**.
        Because FEBio 4.13 strips the domain name from shells, shell #k is
        matched to the k-th unnamed mesh domain.  ``n_elem`` is recorded for
        diagnostics / a sanity check.
    """
    hdf5 = Path(hdf5)
    hdf5.parent.mkdir(parents=True, exist_ok=True)
    _xplt_to_hdf5_patched(xplt, hdf5)

    out: dict = {"label": label, "domains": {}}
    import h5py
    with h5py.File(hdf5, "r") as h5:
        stress_grp = _last_stress_group(h5)
        keys = list(stress_grp.keys())
        order = _mesh_domain_order(h5)
        # collect element counts for the unnamed domains
        dom_grp = h5["meshes/0/domains"]
        dom_nelems = {n: int(dom_grp[n].shape[0]) for n in order}
    out["stress_keys"] = keys
    out["mesh_domain_order"] = order
    print(f"[post:{label}] mesh domains = {order!r}", flush=True)

    # ----- solids: by name -------------------------------------------------
    for name in solid_names:
        key = f"solid_{name}"
        try:
            if name not in keys:
                raise RuntimeError(f"solid 名 {name!r} 不在 stress keys {keys!r}")
            with h5py.File(hdf5, "r") as h5:
                s = np.asarray(_last_stress_group(h5)[name][:], dtype=np.float64)
            stats = _stats_from_voigt(s)
            stats["stress_key"] = repr(name)
            out["domains"][key] = stats
            print(
                f"[post:{label}] solid {name:18s} n={s.shape[0]:>5} "
                f"max={stats['vm_max_mpa']:8.3f} p95={stats['vm_p95_mpa']:8.3f} "
                f"mean={stats['vm_mean_mpa']:8.3f} MPa", flush=True,
            )
        except Exception as exc:  # noqa: BLE001
            out["domains"][key + "_error"] = repr(exc)
            print(f"[post:{label}] solid {name:18s} 读失败：{exc}", flush=True)

    # ----- shells: k-th unnamed, in order ----------------------------------
    unnamed = [n for n in order if n.startswith("__unnamed_")]
    if len(shell_specs) > len(unnamed):
        out["shell_mapping_error"] = (
            f"shell_specs={len(shell_specs)} > unnamed domains={len(unnamed)}"
        )
    for i, (mesh_name, n_elem) in enumerate(shell_specs):
        key = f"shell_{mesh_name}"
        try:
            if i >= len(unnamed):
                raise RuntimeError(
                    f"没有第 {i+1} 个 nameless 域（unnamed={unnamed}）"
                )
            dom = unnamed[i]
            with h5py.File(hdf5, "r") as h5:
                s = np.asarray(_last_stress_group(h5)[dom][:], dtype=np.float64)
            stats = _stats_from_voigt(s)
            stats["stress_key"] = repr(dom)
            stats["n_elem_expected"] = int(n_elem)
            stats["count_ok"] = bool(s.shape[0] == n_elem)
            out["domains"][key] = stats
            print(
                f"[post:{label}] shell {mesh_name:18s} n={s.shape[0]:>5} "
                f"(exp {n_elem}) max={stats['vm_max_mpa']:8.3f} "
                f"p95={stats['vm_p95_mpa']:8.3f} "
                f"mean={stats['vm_mean_mpa']:8.3f} MPa key={dom!r}", flush=True,
            )
        except Exception as exc:  # noqa: BLE001
            out["domains"][key + "_error"] = repr(exc)
            print(f"[post:{label}] shell {mesh_name:18s} 读失败：{exc}", flush=True)
    return out