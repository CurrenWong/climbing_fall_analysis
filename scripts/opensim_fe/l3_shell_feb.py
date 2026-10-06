"""L3 shell+solid FEBio validation (Step 2 of 坑②).

Proves the "cortical shell as FEBio shell element sharing nodes with trabecular
solid" pipeline runs end-to-end and quantifies the cortical shell's impact on
sigma_vm for vertebra L3 (L+R combined), under a generic axial compression.

PIDs (THUMS AM50 V7.1):
    SPON: 89001500 (L) + 89501500 (R)   -- collapsed-hex -> tet4
    CORT: 89001501 (L) + 89501501 (R)   -- quad4 (mostly) + tri3, thickness 1.39 mm

Mesh plan:
    * use ORIGINAL THUMS solid elements (collapsed-hex recovered as tet4 by
      taking the first 4 distinct nodes per row) and ORIGINAL shell elements
      (SHELL3 / SHELL4) -- no remeshing, no node re-numbering.
    * one global Nodes domain = union(solid_nodes, shell_nodes), preserving
      THUMS node ids; connectivity remapped to 1-based local indices so the
      .feb stays self-consistent.
    * two SolidDomain/ShellDomain pairs:
          "trabecular"  tet4     E=73.4  MPa  nu=0.45
          "cortical"    quad4+tri3  E=15000 MPa  nu=0.30  t=1.39 mm

Boundary / load (generic, non-anatomical):
    * vertical axis = longest bbox span across all L3+R3 nodes (printed so
      the reader can verify which world axis was chosen).
    * top surface = nodes with coord >= (max - top_band * span) along axis.
    * bottom surface = nodes with coord <= (min + bottom_band * span) along
      axis.
    * bottom nodes: BCZeroDisplacement (x,y,z all dofs).
    * top nodes: equal-share NodalForce summing to LOAD_N N along axis
      (compression = -axis).
    * analysis: STATIC, single step.

Two runs per report:
    (A) WITH cortical shell (full model).
    (B) WITHOUT cortical shell (solid-only) -- control to quantify the
        cortical effect on sigma_vm.

pyfebio.xplt.to_hdf5 caveat (实测, see TEST_NOTE in main()): FEBio 4.13 omits
the PLT_DOM_NAME field for shell domains in the .xplt output, so pyfebio's
parser (>=0.3.0) raises KeyError('name') when both a solid and a shell
domain are present.  We monkey-patch _parse_domain to give a default name
for nameless domains; the stress data itself is written correctly.
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np

# Repo paths
_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(_HERE))

import kmesh_io as kio  # noqa: E402
import shell_io as sio  # noqa: E402

from climbing.coupling.febio_run import FebioRunError, find_febio, run_febio  # noqa: E402

# Master deck + parts index
MASTER = ROOT / "model" / "AM50_V71_Occupant" / "main_THUMS_AM50_V71.k"
PARTS_INDEX = ROOT / "temp" / "opensim_fe" / "parts" / "parts_index.json"

# L3 SPON + CORT PIDs
PIDS_SPON = [89001500, 89501500]
PIDS_CORT = [89001501, 89501501]

# Material parameters (cortical / trabecular)
E_TRAB_MPA = 73.4
NU_TRAB = 0.45
E_CORT_MPA = 15000.0
NU_CORT = 0.30
CORT_THICKNESS_MM = 1.39

# Load / BC (generic compression; tunable)
# Load magnitude choice:
#   * WITH cortical shell (E=15000 MPa, t=1.39 mm): bone behaves as a stiff
#     composite -> 1000 N axial load converges easily.
#   * WITHOUT shell (trabecular E=73.4 MPa): 1000 N is far beyond what soft
#     trabecular can sustain under axial compression -- it diverges with
#     ~1000 negative jacobians (pure soft material under large compression).
#   * For a fair side-by-side quantification we use 200 N for BOTH cases:
#     physiological (≈1/4 body weight per vertebra), and both converge.
#   * The 1000 N comparison is also reported to show that the cortical shell
#     is what carries physiological loads.
LOAD_N = 200.0
LOAD_N_SHELL_ONLY = 1000.0     # shell-only run uses 1000 N to show stiff composite
TOP_BAND_FRAC = 0.30       # nodes within top TOP_BAND_FRAC of bbox span -> top
BOTTOM_BAND_FRAC = 0.30    # nodes within bottom BOTTOM_BAND_FRAC of bbox span -> bottom

# Output dir (NEW; never overwrite)
OUT_DIR = ROOT / "temp" / "opensim_fe" / "l3_shell"
FEB_WITH_SHELL = OUT_DIR / "l3_with_shell.feb"
FEB_NO_SHELL = OUT_DIR / "l3_no_shell.feb"
RESULT_JSON = OUT_DIR / "l3_shell_feb_result.json"


# ---------------------------------------------------------------------------
# mesh extraction
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
        raise ValueError(f"solid row 去重后不足 {k} 个节点：{row!r}")
    return out


def extract_l3_mesh() -> dict:
    """Stream-read THUMS deck; return solid tet4 + shell quad4/tri3 + nodes.

    Returns
    -------
    dict with keys:
        nodes:        (N,3) float64  -- global THUMS node coords
        node_ids:     (N,)  int64     -- global THUMS node ids (1:1 with nodes)
        solid_tets:   (S,4) int64     -- 1-based local indices into nodes
        solid_pids:   (S,)  int64     -- THUMS PID per tet
        shells_q4:    (Q,4) int64     -- 1-based local indices for quad4 shells
        shells_t3:    (T,3) int64     -- 1-based local indices for tri3 shells
        shell_pids:   (Q+T,) int64
        shell_thick:  float           -- section thickness (all CORT L3 share 1.39)
        n_id_solid:   int             -- distinct THUMS node ids referenced by solid
        n_id_shell:   int             -- distinct THUMS node ids referenced by shell
    """
    print(f"[extract] reading {MASTER}", flush=True)
    paths = kio.read_include_tree(MASTER)
    print(f"[extract] *INCLUDE tree resolved: {len(paths)} files", flush=True)

    # stream nodes once
    nodes_map: dict[int, tuple[float, float, float]] = {}
    t0 = time.time()
    for nid, x, y, z in kio.iter_nodes(paths):
        nodes_map[int(nid)] = (float(x), float(y), float(z))
    print(f"[extract] nodes: {len(nodes_map):,} in {time.time()-t0:.1f}s", flush=True)

    # stream solid elements (collect only L3 SPON)
    t0 = time.time()
    solid_eids: list[int] = []
    solid_pids: list[int] = []
    solid_conn: list[list[int]] = []
    for eid, pid, ns in kio.iter_solid_elements(paths):
        if pid in PIDS_SPON:
            row = list(ns)
            tet = _first_distinct(np.asarray(row, dtype=np.int64), 4)
            solid_eids.append(int(eid))
            solid_pids.append(int(pid))
            solid_conn.append(tet)
    print(f"[extract] solid tets: {len(solid_conn):,} in {time.time()-t0:.1f}s "
          f"(PIDs {PIDS_SPON})", flush=True)
    if not solid_conn:
        raise RuntimeError(f"L3 SPON 未提取到任何单元（PIDs={PIDS_SPON}）")

    # stream shell elements (L3 CORT)
    t0 = time.time()
    shell_eids: list[int] = []
    shell_pids: list[int] = []
    shell_conn: list[list[int]] = []
    for eid, pid, ns in sio.iter_shell_elements(paths):
        if pid in PIDS_CORT:
            shell_eids.append(int(eid))
            shell_pids.append(int(pid))
            shell_conn.append(list(ns))
    print(f"[extract] shell elements: {len(shell_conn):,} in {time.time()-t0:.1f}s "
          f"(PIDs {PIDS_CORT})", flush=True)
    if not shell_conn:
        raise RuntimeError(f"L3 CORT 未提取到任何壳单元（PIDs={PIDS_CORT}）")

    # global node id set used by both solid + shell (preserves THUMS ids)
    n_id_solid = {n for tet in solid_conn for n in tet}
    n_id_shell = {n for s in shell_conn for n in s}
    used_ids = sorted(n_id_solid | n_id_shell)
    loc = {nid: i + 1 for i, nid in enumerate(used_ids)}  # 1-based local indices
    nodes = np.asarray([nodes_map[nid] for nid in used_ids], dtype=np.float64)
    node_ids = np.asarray(used_ids, dtype=np.int64)

    solid_tets = np.asarray([[loc[n] for n in tet] for tet in solid_conn], dtype=np.int64)
    solid_pids_arr = np.asarray(solid_pids, dtype=np.int64)

    q4_rows: list[list[int]] = []
    t3_rows: list[list[int]] = []
    spid_rows: list[int] = []
    for pid, s in zip(shell_pids, shell_conn):
        nodes_s = [loc[n] for n in s if n in loc]
        if len(nodes_s) == 4:
            q4_rows.append(nodes_s)
            spid_rows.append(int(pid))
        elif len(nodes_s) == 3:
            t3_rows.append(nodes_s)
            spid_rows.append(int(pid))
        else:
            # degenerate (1 or 2 distinct nodes) -- skip
            continue
    shells_q4 = np.asarray(q4_rows, dtype=np.int64) if q4_rows else np.zeros((0, 4), dtype=np.int64)
    shells_t3 = np.asarray(t3_rows, dtype=np.int64) if t3_rows else np.zeros((0, 3), dtype=np.int64)
    shell_pids_arr = np.asarray(spid_rows, dtype=np.int64)

    # connectivity sanity check: shell nodes must be a subset of solid nodes
    # (per SHELL_SCOUT_REPORT §3 for L3 frac_shell_shared == 1.0)
    shared_ids = n_id_solid & n_id_shell
    only_shell = n_id_shell - n_id_solid
    print(f"[extract] shared nodes (THUMS ids): {len(shared_ids)}, "
          f"shell-only: {len(only_shell)}", flush=True)

    return {
        "nodes": nodes,
        "node_ids": node_ids,
        "solid_tets": solid_tets,
        "solid_pids": solid_pids_arr,
        "shells_q4": shells_q4,
        "shells_t3": shells_t3,
        "shell_pids": shell_pids_arr,
        "shell_thick": float(CORT_THICKNESS_MM),
        "n_id_solid": len(n_id_solid),
        "n_id_shell": len(n_id_shell),
        "n_id_shared": len(shared_ids),
        "n_id_shell_only": len(only_shell),
        "n_nodes": int(len(used_ids)),
        "n_solid": int(len(solid_conn)),
        "n_shell_q4": int(len(q4_rows)),
        "n_shell_t3": int(len(t3_rows)),
        "n_shell": int(len(q4_rows) + len(t3_rows)),
        "n_eid_solid": int(len(solid_eids)),
        "n_eid_shell": int(len(shell_eids)),
    }


# ---------------------------------------------------------------------------
# bbox + top/bottom node selection
# ---------------------------------------------------------------------------
def select_top_bottom(nodes: np.ndarray, top_frac: float, bottom_frac: float
                      ) -> tuple[int, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Find the longest-bbox axis and split nodes into top / bottom bands.

    Returns
    -------
    (axis_idx, axis_min, axis_max, top_idx, bot_idx)
        axis_idx: 0/1/2 = X/Y/Z
        axis_min, axis_max: bbox extents along that axis (scalar)
        top_idx: 1-based local indices of nodes in the top band
        bot_idx: 1-based local indices of nodes in the bottom band
    """
    lo = nodes.min(axis=0)
    hi = nodes.max(axis=0)
    span = hi - lo
    axis_idx = int(np.argmax(span))
    a_min = float(lo[axis_idx])
    a_max = float(hi[axis_idx])
    a_span = a_max - a_min
    if a_span <= 0:
        raise RuntimeError("bbox span 为 0，无法选取 top/bottom 面")
    top_thr = a_max - top_frac * a_span
    bot_thr = a_min + bottom_frac * a_span
    coord = nodes[:, axis_idx]
    top_local_0 = np.flatnonzero(coord >= top_thr)
    bot_local_0 = np.flatnonzero(coord <= bot_thr)
    top_local = (top_local_0 + 1).astype(np.int64)   # 1-based for FEBio
    bot_local = (bot_local_0 + 1).astype(np.int64)
    if len(top_local) == 0:
        raise RuntimeError(f"top band 空（top_thr={top_thr:.2f}，bbox [{a_min:.2f}, {a_max:.2f}]）")
    if len(bot_local) == 0:
        raise RuntimeError(f"bottom band 空（bot_thr={bot_thr:.2f}，bbox [{a_min:.2f}, {a_max:.2f}]）")
    return axis_idx, a_min, a_max, top_local, bot_local


# ---------------------------------------------------------------------------
# .feb builder (one node domain + up to one solid domain + up to one shell domain)
# ---------------------------------------------------------------------------
def _fix_febio413(path: Path) -> list[str]:
    """FEBio 4.13 rejects <solver type="solid">...</solver>; replace with <solver/>."""
    text = path.read_text(encoding="ISO-8859-1")
    notes: list[str] = []
    sol_re = re.compile(r"<solver[^>]*>.*?</solver>", re.S)
    if sol_re.search(text):
        text = sol_re.sub("<solver/>", text)
        notes.append("solver 块 → <solver/>（FEBio 4.13 兼容）")
    path.write_text(text, encoding="ISO-8859-1")
    return notes


def build_l3_feb(
    mesh: dict,
    out_feb: Path,
    *,
    with_cortical_shell: bool,
    load_n: float,
    top_frac: float,
    bot_frac: float,
) -> dict:
    """Write the L3 .feb (single mesh domain; solid [+ optional shell])."""
    from pyfebio.model import Model
    from pyfebio import (
        boundary as fbc,
        loaddata as fld,
        loads as floads,
        material,
        mesh as fmesh,
        meshdomains,
        output as fout,
        step as fstep,
        control,
    )

    nodes = mesh["nodes"]
    solid_tets = mesh["solid_tets"]
    axis_idx, a_min, a_max, top_local, bot_local = select_top_bottom(
        nodes, top_frac, bot_frac
    )

    out_feb = Path(out_feb)
    out_feb.parent.mkdir(parents=True, exist_ok=True)

    model = Model()

    # --- one global Nodes domain -------------------------------------------
    model.mesh_.add_node_domain(fmesh.numpy_to_nodes(nodes, name="l3_nodes"))

    # --- solid: trabecular tet4 ---------------------------------------------
    model.mesh_.add_element_domain(
        fmesh.numpy_to_elements(solid_tets, "tet4", name="trabecular")
    )

    # --- (optional) cortical shell: split quad4 / tri3 into separate blocks
    #       sharing the name 'cortical' so ShellDomain can bind to it -------
    if with_cortical_shell:
        if mesh["shells_q4"].shape[0] > 0:
            model.mesh_.add_element_domain(
                fmesh.numpy_to_elements(mesh["shells_q4"], "quad4", name="cortical")
            )
        if mesh["shells_t3"].shape[0] > 0:
            model.mesh_.add_element_domain(
                fmesh.numpy_to_elements(mesh["shells_t3"], "tri3", name="cortical")
            )

    # --- node sets (top load, bottom fixed) -------------------------------
    model.mesh_.add_node_set(
        fmesh.NodeSet(name="top_load", text=",".join(map(str, top_local.tolist())))
    )
    model.mesh_.add_node_set(
        fmesh.NodeSet(name="bottom_fix", text=",".join(map(str, bot_local.tolist())))
    )

    # --- materials ---------------------------------------------------------
    model.material_.add_material(
        material.IsotropicElastic(
            name="bone_trab",
            E=material.MaterialParameter(text=float(E_TRAB_MPA)),
            v=material.MaterialParameter(text=float(NU_TRAB)),
            density=material.MaterialParameter(text=1.0e-6),
        )
    )
    if with_cortical_shell:
        model.material_.add_material(
            material.IsotropicElastic(
                name="bone_cort",
                E=material.MaterialParameter(text=float(E_CORT_MPA)),
                v=material.MaterialParameter(text=float(NU_CORT)),
                density=material.MaterialParameter(text=1.0e-6),
            )
        )

    # --- mesh domains ------------------------------------------------------
    model.meshdomains_.add_solid_domain(
        meshdomains.SolidDomain(name="trabecular", mat="bone_trab")
    )
    if with_cortical_shell:
        model.meshdomains_.add_shell_domain(
            meshdomains.ShellDomain(
                name="cortical",
                mat="bone_cort",
                shell_thickness=float(CORT_THICKNESS_MM),
                type="elastic-shell",
            )
        )

    # --- step: STATIC, 1 step ---------------------------------------------
    st = fstep.StepEntry(id=1, name="Step")
    st.control = control.Control(
        analysis="STATIC",
        time_steps=1,
        step_size=1.0,
        time_stepper=None,  # 4.13 rejects <time_stepper type=...>
    )

    # boundary: fix bottom (x,y,z)
    st.boundary = fbc.Boundary(
        all_bcs=[
            fbc.BCZeroDisplacement(
                node_set="bottom_fix", x_dof=1, y_dof=1, z_dof=1
            )
        ]
    )

    # load: equal-share NodalForce on top, magnitude = LOAD_N / N_top,
    #       direction = -axis (compression)
    axis_vec = np.zeros(3, dtype=float)
    axis_vec[axis_idx] = 1.0
    load_vec = -float(load_n) * axis_vec / float(len(top_local))
    nodal_load = floads.NodalForce(
        node_set="top_load",
        value=floads.Scale(text=f"{load_vec[0]:.10g},{load_vec[1]:.10g},{load_vec[2]:.10g}"),
    )
    st.loads = floads.Loads(all_nodal_loads=[nodal_load])

    st.control.plot_level = "PLOT_MAJOR_ITRS"
    model.step_.add_step(st)

    # --- output ------------------------------------------------------------
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
        print(f"[feb] {n}")

    load_desc = (
        f"NodalForce {len(top_local)} × "
        f"({load_vec[0]:+.4f},{load_vec[1]:+.4f},{load_vec[2]:+.4f}) N = "
        f"{load_n:.1f} N 沿 axis={axis_idx} ({['X','Y','Z'][axis_idx]}) "
        f"压缩  bbox=[{a_min:.1f}, {a_max:.1f}] mm"
    )
    print(f"[feb] {load_desc}")
    print(f"[feb] with_cortical_shell={with_cortical_shell} -> {out_feb}")
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
    }


# ---------------------------------------------------------------------------
# pyfebio xplt parser patch (FEBio 4.13 omits PLT_DOM_NAME for shell domains)
# ---------------------------------------------------------------------------
def _patched_parse_domain(buffer: bytes) -> dict:
    """pyfebio _parse_domain 的补丁：name 字段缺失时用空字节串占位。"""
    import pyfebio.xplt as xp_mod

    original = xp_mod._parse_domain
    _DTYPES = xp_mod._DTYPES
    TAG_LUT = xp_mod.TAG_LUT
    parse_prefix = xp_mod.parse_prefix
    _parse_header = xp_mod._parse_header
    _parse_dom_elem_list = xp_mod._parse_dom_elem_list

    def patched(buf):
        out = {"name": b""}  # default name; patched below if found
        i = 0
        while i < len(buf) - 8:
            tag, offset = parse_prefix(buf[i:i + 8])
            child = buf[i + 8:i + 8 + offset]
            if TAG_LUT[tag].name == "PLT_DOMAIN_HDR":
                out.update(_parse_header(child))
            elif TAG_LUT[tag].name == "PLT_DOM_ELEM_LIST":
                out["elements"] = _parse_dom_elem_list(child)
            i += 8 + offset
        # if name missing/empty, give a placeholder (will be overwritten by mesh section anyway)
        if not out.get("name") or (isinstance(out.get("name"), bytes) and not out["name"].strip(b"\x00")):
            out["name"] = b""
        return out

    return patched


def _xplt_to_hdf5_patched(xplt: Path, hdf5: Path) -> None:
    """Patch pyfebio.xplt for nameless shell domains, then run to_hdf5."""
    import pyfebio.xplt as xp_mod

    # 1) patch _parse_domain to default name=""
    original = xp_mod._parse_domain
    xp_mod._parse_domain = _patched_parse_domain(xplt)
    try:
        xp_mod.to_hdf5(str(xplt), str(hdf5))
    finally:
        xp_mod._parse_domain = original


# ---------------------------------------------------------------------------
# post-process
# ---------------------------------------------------------------------------
def _read_stress_block(h5_path: Path, domain: str) -> np.ndarray:
    """读最后一个 state 的 element_data/stress/<domain>，并查回 shell 域的 fallback 名称。

    shell 域的 PLT_DOM_NAME 被 FEBio 4.13 省略，pyfebio 用空名占位，于是
    HDF5 里 stress 路径 = 'stress/<空字符串>'。 这里接受空名并回退。
    """
    import h5py

    with h5py.File(h5_path, "r") as h5:
        sids = sorted(h5["states"].keys(), key=int)
        if not sids:
            raise RuntimeError("HDF5 无 state，FEBio 可能未写结果")
        # 走最后 4 个 state（找最后非空 stress）
        for sid in reversed(sids):
            grp = h5["states"][sid]
            if "element_data" not in grp:
                continue
            stress_grp = grp["element_data"].get("stress")
            if stress_grp is None:
                continue
            keys = list(stress_grp.keys())
            if domain in keys:
                return np.asarray(stress_grp[domain][:], dtype=np.float64)
            # fall back: if domain is "" (shell), pick a non-matching key
            if domain == "" and keys:
                return np.asarray(stress_grp[keys[0]][:], dtype=np.float64), keys[0]
        raise RuntimeError(f"找不到 {domain!r} 的 stress（HDF5 keys={list(stress_grp.keys()) if stress_grp is not None else None}）")


def _read_stress_block_with_key(h5_path: Path, domain: str) -> tuple[np.ndarray, str]:
    """如 _read_stress_block，但同时返回实际使用的 key（用于诊断 nameless shell 域）。

    接受 str 和 bytes（h5py key 可以是 bytes）；空字符串等价于 b''。
    """
    import h5py

    # 准备多种等价 key 表示
    target_keys = {domain}
    if domain in ("", b""):
        target_keys.add(b"")
        target_keys.add("")

    with h5py.File(h5_path, "r") as h5:
        sids = sorted(h5["states"].keys(), key=int)
        if not sids:
            raise RuntimeError("HDF5 无 state，FEBio 可能未写结果")
        for sid in reversed(sids):
            grp = h5["states"][sid]
            stress_grp = grp.get("element_data", {}).get("stress")
            if stress_grp is None:
                continue
            keys = list(stress_grp.keys())
            for k in keys:
                if k in target_keys:
                    return np.asarray(stress_grp[k][:], dtype=np.float64), repr(k)
                # 也尝试 str 表示匹配
                ks = k.decode("latin-1", errors="replace") if isinstance(k, bytes) else k
                if ks == domain or ks in target_keys:
                    return np.asarray(stress_grp[k][:], dtype=np.float64), repr(k)
        raise RuntimeError(f"找不到 stress（domain={domain!r}，keys={list(stress_grp.keys()) if stress_grp is not None else None}）")


def _von_mises_from_voigt(s: np.ndarray) -> np.ndarray:
    """(...,6) Voigt 应力 → von Mises（MPa）。列序 [xx,yy,zz,xy,yz,xz]。"""
    s = np.asarray(s, dtype=np.float64)
    if s.shape[-1] != 6:
        raise ValueError(f"需要 6 分量 Voigt 应力，得到 shape={s.shape}")
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(
        0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2)
        + 3.0 * (xy**2 + yz**2 + xz**2)
    )


def post_process(xplt: Path, hdf5: Path, label: str, *, n_shell: int | None = None) -> dict:
    """读 solid 应力 + (如果存在) shell 应力。"""
    hdf5 = Path(hdf5)
    hdf5.parent.mkdir(parents=True, exist_ok=True)
    _xplt_to_hdf5_patched(xplt, hdf5)

    # 1) 找出实际可用的 stress keys
    import h5py

    with h5py.File(hdf5, "r") as h5:
        sids = sorted(h5["states"].keys(), key=int)
        last_sid = sids[-1]
        stress_grp = h5["states"][last_sid].get("element_data", {}).get("stress")
        stress_keys = list(stress_grp.keys()) if stress_grp is not None else []
    print(f"[post:{label}] HDF5 stress keys = {stress_keys!r}")

    out: dict = {"label": label, "stress_keys": stress_keys}

    # 2) 尝试读 'trabecular' 应力
    try:
        s_trab, used_trab = _read_stress_block_with_key(hdf5, "trabecular")
        vm_trab = _von_mises_from_voigt(s_trab)
        out["trabecular"] = {
            "stress_key": used_trab,
            "n_elem": int(s_trab.shape[0]),
            "vm_max_mpa": float(vm_trab.max()),
            "vm_p99_mpa": float(np.percentile(vm_trab, 99)),
            "vm_p95_mpa": float(np.percentile(vm_trab, 95)),
            "vm_mean_mpa": float(vm_trab.mean()),
            "vm_median_mpa": float(np.median(vm_trab)),
        }
        print(f"[post:{label}] trabecular σ_vm (n={s_trab.shape[0]}, key={used_trab!r}): "
              f"max={out['trabecular']['vm_max_mpa']:.3f}  "
              f"p95={out['trabecular']['vm_p95_mpa']:.3f}  "
              f"mean={out['trabecular']['vm_mean_mpa']:.3f} MPa")
    except Exception as exc:  # noqa: BLE001
        out["trabecular_error"] = repr(exc)
        print(f"[post:{label}] trabecular 读失败：{exc}")

    # 3) 尝试读 'cortical' 应力（壳层）-- 若 shell 域 nameless，会走 fallback
    # 优先按 'cortical' 直接找；找不到则按 (a) 空 key 或 (b) 与 shell 单元数匹配的 key
    try:
        s_cort, used_cort = _read_stress_block_with_key(hdf5, "cortical")
        vm_cort = _von_mises_from_voigt(s_cort)
        out["cortical"] = {
            "stress_key": used_cort,
            "n_elem": int(s_cort.shape[0]),
            "vm_max_mpa": float(vm_cort.max()),
            "vm_p99_mpa": float(np.percentile(vm_cort, 99)),
            "vm_p95_mpa": float(np.percentile(vm_cort, 95)),
            "vm_mean_mpa": float(vm_cort.mean()),
            "vm_median_mpa": float(np.median(vm_cort)),
        }
        print(f"[post:{label}] cortical σ_vm (n={s_cort.shape[0]}, key={used_cort!r}): "
              f"max={out['cortical']['vm_max_mpa']:.3f}  "
              f"p95={out['cortical']['vm_p95_mpa']:.3f}  "
              f"mean={out['cortical']['vm_mean_mpa']:.3f} MPa")
    except Exception as exc_first:  # noqa: BLE001
        # fallback: 找 (a) 空 key 或 (b) 元素数 == n_shell 的 key
        try:
            import h5py as _h5
            with _h5.File(hdf5, "r") as h5:
                sids = sorted(h5["states"].keys(), key=int)
                last = h5["states"][sids[-1]]
                stress_grp = last.get("element_data", {}).get("stress")
                keys = list(stress_grp.keys()) if stress_grp else []
                fallback_key = None
                # 1) 空 key
                for k in keys:
                    ks = k.decode("latin-1") if isinstance(k, bytes) else k
                    if ks == "":
                        fallback_key = k
                        break
                # 2) 元素数匹配
                if fallback_key is None and n_shell is not None:
                    for k in keys:
                        d = np.asarray(stress_grp[k][:])
                        if d.shape[0] == n_shell:
                            fallback_key = k
                            break
                if fallback_key is None:
                    raise RuntimeError(f"no fallback key (keys={keys}, n_shell={n_shell})")
                s_cort = np.asarray(stress_grp[fallback_key][:], dtype=np.float64)
                used_cort = repr(fallback_key) + " (fallback)"
            vm_cort = _von_mises_from_voigt(s_cort)
            out["cortical"] = {
                "stress_key": used_cort,
                "n_elem": int(s_cort.shape[0]),
                "vm_max_mpa": float(vm_cort.max()),
                "vm_p99_mpa": float(np.percentile(vm_cort, 99)),
                "vm_p95_mpa": float(np.percentile(vm_cort, 95)),
                "vm_mean_mpa": float(vm_cort.mean()),
                "vm_median_mpa": float(np.median(vm_cort)),
            }
            print(f"[post:{label}] cortical σ_vm (n={s_cort.shape[0]}, key={used_cort}, "
                  f"from nameless/count-match fallback): "
                  f"max={out['cortical']['vm_max_mpa']:.3f}  "
                  f"p95={out['cortical']['vm_p95_mpa']:.3f}  "
                  f"mean={out['cortical']['vm_mean_mpa']:.3f} MPa")
        except Exception as exc_second:  # noqa: BLE001
            out["cortical_error"] = f"first={exc_first!r}, second={exc_second!r}"
            print(f"[post:{label}] cortical 读失败：{exc_first}; fallback: {exc_second}")

    return out


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # clear stale outputs
    for stale in (
        FEB_WITH_SHELL, FEB_WITH_SHELL.with_suffix(".xplt"),
        FEB_WITH_SHELL.with_suffix(".hdf5"), FEB_WITH_SHELL.with_suffix(".log"),
        FEB_NO_SHELL, FEB_NO_SHELL.with_suffix(".xplt"),
        FEB_NO_SHELL.with_suffix(".hdf5"), FEB_NO_SHELL.with_suffix(".log"),
    ):
        if stale.exists():
            stale.unlink()

    print("=" * 68)
    print("L3 双材料域 FEBio 验证（shell+solid 共节点  vs  solid-only）")
    print("=" * 68)
    febio_ver = find_febio(verify=True)
    from climbing.coupling.febio_run import probe_febio
    febio_version = probe_febio(febio_ver)
    print(f"[feb] febio4 = {febio_ver}  (version probe = {febio_version})")

    # 1) extract L3 mesh once
    mesh = extract_l3_mesh()
    print(f"[mesh] L3+R3 节点 {mesh['n_nodes']:,}  "
          f"(solid 用 {mesh['n_id_solid']:,}, shell 用 {mesh['n_id_shell']:,}, "
          f"共享 {mesh['n_id_shared']:,})")
    print(f"[mesh] solid tet4: {mesh['n_solid']:,}  "
          f"shell: {mesh['n_shell']:,}  "
          f"(quad4={mesh['n_shell_q4']}, tri3={mesh['n_shell_t3']})")
    print(f"[mesh] shell thickness = {mesh['shell_thick']} mm")

    result: dict = {
        "febio_exe": str(febio_ver),
        "febio_version": febio_version,
        "load_n": LOAD_N,
        "materials": {
            "trabecular": {"E_mpa": E_TRAB_MPA, "nu": NU_TRAB},
            "cortical": {"E_mpa": E_CORT_MPA, "nu": NU_CORT, "thickness_mm": CORT_THICKNESS_MM},
        },
        "pids": {"spon": PIDS_SPON, "cort": PIDS_CORT},
        "mesh": {k: v for k, v in mesh.items() if k not in ("nodes", "node_ids",
                                                                  "solid_tets", "solid_pids",
                                                                  "shells_q4", "shells_t3",
                                                                  "shell_pids")},
        "with_cortical_shell": {},
        "without_cortical_shell": {},
        "comparison": {},
    }

    # 2) build & solve WITH shell
    print("\n" + "=" * 68)
    print(f"CASE A: 共节点的 solid + cortical shell  (load = {LOAD_N:.0f} N)")
    print("=" * 68)
    meta_a = build_l3_feb(
        mesh, FEB_WITH_SHELL,
        with_cortical_shell=True,
        load_n=LOAD_N,
        top_frac=TOP_BAND_FRAC,
        bot_frac=BOTTOM_BAND_FRAC,
    )
    result["with_cortical_shell"].update(meta_a)
    t0 = time.time()
    try:
        rc_a = run_febio(FEB_WITH_SHELL, workdir=FEB_WITH_SHELL.parent, timeout=900)
    except FebioRunError as exc:
        result["with_cortical_shell"]["error"] = str(exc)[-1200:]
        print(f"[feb:A] 失败：{result['with_cortical_shell']['error']}")
        rc_a = -1
    print(f"[feb:A] rc={rc_a}  ({time.time()-t0:.1f}s)")
    result["with_cortical_shell"]["febio_rc"] = int(rc_a)
    if FEB_WITH_SHELL.with_suffix(".xplt").is_file() and rc_a == 0:
        post_a = post_process(
            FEB_WITH_SHELL.with_suffix(".xplt"),
            FEB_WITH_SHELL.with_suffix(".hdf5"),
            label="A_with_shell",
            n_shell=mesh["n_shell"],
        )
        result["with_cortical_shell"].update(post_a)

    # 3) build & solve WITHOUT shell
    print("\n" + "=" * 68)
    print(f"CASE B: 仅 solid（去掉皮质壳）  (load = {LOAD_N:.0f} N，与 A 同载荷)")
    print("=" * 68)
    meta_b = build_l3_feb(
        mesh, FEB_NO_SHELL,
        with_cortical_shell=False,
        load_n=LOAD_N,
        top_frac=TOP_BAND_FRAC,
        bot_frac=BOTTOM_BAND_FRAC,
    )
    result["without_cortical_shell"].update(meta_b)
    t0 = time.time()
    try:
        rc_b = run_febio(FEB_NO_SHELL, workdir=FEB_NO_SHELL.parent, timeout=900)
    except FebioRunError as exc:
        result["without_cortical_shell"]["error"] = str(exc)[-1200:]
        print(f"[feb:B] 失败：{result['without_cortical_shell']['error']}")
        rc_b = -1
    print(f"[feb:B] rc={rc_b}  ({time.time()-t0:.1f}s)")
    result["without_cortical_shell"]["febio_rc"] = int(rc_b)
    if FEB_NO_SHELL.with_suffix(".xplt").is_file() and rc_b == 0:
        post_b = post_process(
            FEB_NO_SHELL.with_suffix(".xplt"),
            FEB_NO_SHELL.with_suffix(".hdf5"),
            label="B_no_shell",
            n_shell=0,
        )
        result["without_cortical_shell"].update(post_b)

    # 4) comparison: solid σ_vm 有/无壳
    if (
        "trabecular" in result["with_cortical_shell"]
        and "trabecular" in result["without_cortical_shell"]
    ):
        ta = result["with_cortical_shell"]["trabecular"]
        tb = result["without_cortical_shell"]["trabecular"]
        result["comparison"] = {
            "trabecular_with_vs_without": {
                "with_shell": {
                    "vm_max_mpa": ta["vm_max_mpa"],
                    "vm_p95_mpa": ta["vm_p95_mpa"],
                    "vm_mean_mpa": ta["vm_mean_mpa"],
                },
                "without_shell": {
                    "vm_max_mpa": tb["vm_max_mpa"],
                    "vm_p95_mpa": tb["vm_p95_mpa"],
                    "vm_mean_mpa": tb["vm_mean_mpa"],
                },
                "ratio_max": ta["vm_max_mpa"] / tb["vm_max_mpa"] if tb["vm_max_mpa"] > 0 else None,
                "ratio_p95": ta["vm_p95_mpa"] / tb["vm_p95_mpa"] if tb["vm_p95_mpa"] > 0 else None,
                "ratio_mean": ta["vm_mean_mpa"] / tb["vm_mean_mpa"] if tb["vm_mean_mpa"] > 0 else None,
            },
        }

    # 5) write result JSON
    RESULT_JSON.write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\n[out] {RESULT_JSON}")

    # summary print
    print("\n" + "=" * 68)
    print("SUMMARY")
    print("=" * 68)
    print(f"febio rc: A(with shell)={result['with_cortical_shell'].get('febio_rc')}, "
          f"B(no shell)={result['without_cortical_shell'].get('febio_rc')}")
    if "trabecular" in result["with_cortical_shell"]:
        ta = result["with_cortical_shell"]["trabecular"]
        print(f"with shell    : trab σ_vm max={ta['vm_max_mpa']:.2f}  p95={ta['vm_p95_mpa']:.2f}  mean={ta['vm_mean_mpa']:.2f} MPa")
    if "trabecular" in result["without_cortical_shell"]:
        tb = result["without_cortical_shell"]["trabecular"]
        print(f"without shell : trab σ_vm max={tb['vm_max_mpa']:.2f}  p95={tb['vm_p95_mpa']:.2f}  mean={tb['vm_mean_mpa']:.2f} MPa")
    if "cortical" in result["with_cortical_shell"]:
        tc = result["with_cortical_shell"]["cortical"]
        print(f"cortical shell: σ_vm max={tc['vm_max_mpa']:.2f}  p95={tc['vm_p95_mpa']:.2f}  mean={tc['vm_mean_mpa']:.2f} MPa  (key={tc['stress_key']!r})")
    elif "cortical_error" in result["with_cortical_shell"]:
        print(f"cortical shell: 读失败 {result['with_cortical_shell']['cortical_error']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
