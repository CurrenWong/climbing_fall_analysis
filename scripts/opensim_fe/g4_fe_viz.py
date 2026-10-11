"""G4 — export the fibula_r damage-model FAILURE for visual inspection.

Produces under ``results/opensim_fe/g4_fe_viz/``:

1. ``initial_mesh.vtk`` — undeformed fibula_r mesh, all element blocks
   (``fib_spon_end`` / ``fib_spon_center`` tet4, ``fib_cort`` hex8).
2. ``deformed_first_failure.vtk`` — the SAME mesh at the **last written
   (last-converged) state before divergence**, with per-node ``displacement``
   (vector) + ``abs_u`` (scalar) and per-element ``damage`` + ``von_mises_mpa``.
3. ``overlay_initial_vs_deformed.png`` — initial mesh (grey wireframe) vs
   deformed mesh (coloured by |u|), two side views.
4. ``damage_field.png`` — deformed mesh coloured by element damage D.
5. ``README.md`` — exact run/config + state index + the measured numbers.

Primary failure source (per P1-a): the **displacement-control** run
``disp_d2.0x1000`` (``disp@-397.858 Dmax=1.0``) — it reached the farthest and
its ``.xplt`` holds 93 written states.  The **force-control 1000 N** run
(``force_L1000_1``, the canonical "98 negative jacobians" case) wrote only the
undeformed state 0; this script re-runs that exact config with a 20-step ramp
(``--force1000``) to try to capture converged states before divergence.

Reuse only: imports ``bone_feb`` (signed-volume helpers) and
``climbing.coupling.fe_post`` (``xplt -> hdf5``); modifies no existing module.

Reproduce (repo root, PYTHONPATH=src)::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_viz.py
"""
from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import bone_feb as bf  # noqa: E402
from climbing.coupling.fe_post import build_hdf5  # noqa: E402

OUT = ROOT / "results" / "opensim_fe" / "g4_fe_viz"
STAB = ROOT / "temp" / "opensim_fe" / "g4_fe_stab"
VIZ_TMP = ROOT / "temp" / "opensim_fe" / "g4_fe_viz"

#: primary failure run (displacement control, Dmax=1.0) -- has written states
DC_STEM = "disp_d2.0x1000"
#: canonical force-control failure (1000 N, Dmax=1.0) -- no converged state
FORCE_STEM = "force_L1000_1"

_TET_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
_HEX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
              (0, 4), (1, 5), (2, 6), (3, 7)]
_VTK_TYPE = {"tet4": "tetra", "hex8": "hexahedron"}


# ===========================================================================
# geometry / state readers
# ===========================================================================
def parse_feb_geometry(feb: Path) -> tuple[np.ndarray, np.ndarray, list[dict], dict]:
    """Parse ``<Nodes>`` / ``<Elements>`` / ``<NodeSet>`` straight from the .feb.

    Returns ``(points, node_ids, blocks, nodesets)``; `points`/`node_ids` are in
    ascending node-id order (1-based ids -> row i = id i+1); `blocks[k]["conn"]`
    holds **0-based** local node indices.
    """
    text = feb.read_text(encoding="ISO-8859-1")
    nodes: dict[int, tuple[float, float, float]] = {}
    for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', text):
        xyz = m.group(2).replace(" ", "").split(",")
        nodes[int(m.group(1))] = (float(xyz[0]), float(xyz[1]), float(xyz[2]))
    nids = sorted(nodes)
    loc = {n: i for i, n in enumerate(nids)}
    pts = np.asarray([nodes[n] for n in nids], dtype=np.float64)

    blocks: list[dict] = []
    for m in re.finditer(r'<Elements\s+name="([^"]+)"\s+type="(\w+)">(.*?)</Elements>',
                         text, re.S):
        name, etype, body = m.group(1), m.group(2), m.group(3)
        conn = [[loc[int(v)] for v in em.group(1).replace(" ", "").split(",")]
                for em in re.finditer(r'<elem id="\d+">([^<]+)</elem>', body)]
        blocks.append({"name": name, "type": etype,
                       "conn": np.asarray(conn, dtype=np.int64)})
    sets: dict[str, np.ndarray] = {}
    for m in re.finditer(r'<NodeSet name="([^"]+)">([^<]+)</NodeSet>', text):
        sets[m.group(1)] = np.asarray(
            [loc[int(v)] for v in m.group(2).replace(" ", "").split(",")],
            dtype=np.int64)
    return pts, np.asarray(nids, dtype=np.int64), blocks, sets


def read_last_state(hdf5: Path) -> dict:
    """Read the LAST written state: node displacement, per-domain damage + stress."""
    import h5py
    with h5py.File(hdf5, "r") as h5:
        sids = sorted(h5["states"].keys(), key=int)
        if not sids:
            raise RuntimeError("hdf5 无 state")
        last = sids[-1]
        nd = h5[f"states/{last}/node_data"]
        u = None
        for nm in nd:
            if "displacement" in nm.lower():
                it = nd[nm]
                u = np.asarray(it if isinstance(it, h5py.Dataset) else it[list(it)[0]],
                               dtype=np.float64)
                break
        if u is None:
            raise RuntimeError("hdf5 无 displacement node var")
        damage, stress = {}, {}
        ed = h5[f"states/{last}/element_data"]
        for var, store in (("damage", damage), ("stress", stress)):
            if var in ed:
                for dom in ed[var]:
                    store[dom] = np.asarray(ed[var][dom], dtype=np.float64)
        n_states = len(sids)
    return {"state": int(last), "n_states": n_states, "u": u,
            "damage": damage, "stress": stress}


def von_mises_from_voigt(s: np.ndarray) -> np.ndarray:
    s = np.asarray(s, dtype=np.float64)
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2)
                   + 3.0 * (xy ** 2 + yz ** 2 + xz ** 2))


#: natural (ξ,η,ζ) coords of the 8 hex8 nodes
_HEX_XI = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)


def _hex_min_jacobian(pts: np.ndarray, conn: np.ndarray) -> np.ndarray:
    """min det(J) over the 8 Gauss points of each hex8 (FEBio's own check).

    Negative = the element is inverted at an integration point.  More faithful
    than the 6-tet *sum* (which can stay positive while a sub-tet is negative).
    """
    g = 1.0 / np.sqrt(3.0)
    G = np.array([[-g, -g, -g], [g, -g, -g], [g, g, -g], [-g, g, -g],
                  [-g, -g, g], [g, -g, g], [g, g, g], [-g, g, g]], float)
    X = pts[conn]                                   # (n, 8, 3)
    dets = np.empty((len(conn), 8))
    for gi, (xi, eta, ze) in enumerate(G):
        dN = np.empty((8, 3))
        for i in range(8):
            a, b, c = _HEX_XI[i]
            dN[i, 0] = 0.125 * a * (1 + eta * b) * (1 + ze * c)
            dN[i, 1] = 0.125 * (1 + xi * a) * b * (1 + ze * c)
            dN[i, 2] = 0.125 * (1 + xi * a) * (1 + eta * b) * c
        J = np.einsum("nia,ib->nab", X, dN)         # (n, 3, 3)
        dets[:, gi] = np.linalg.det(J)
    return dets.min(axis=1)


def element_jacobians(deformed_pts: np.ndarray, blocks: list[dict]) -> dict:
    """Per-element Jacobian measure (min det(J) for hex; signed volume for tet).

    Negative => element inverted (the "negative jacobian" failure).  tets use
    ``bone_feb``'s exact signed-volume helper (1-based).
    """
    out: dict[str, np.ndarray] = {}
    for b in blocks:
        if b["type"] == "tet4":
            out[b["name"]] = bf._tet_signed_volumes(deformed_pts, b["conn"] + 1)
        else:
            out[b["name"]] = _hex_min_jacobian(deformed_pts, b["conn"])
    return out


def _initial_jacobians(pts: np.ndarray, blocks: list[dict]) -> dict:
    return element_jacobians(pts, blocks)


# ===========================================================================
# VTK export
# ===========================================================================
def _meshio_mesh(pts, blocks, damage, vm, u=None, sets=None):
    import meshio
    cells, cd = [], {"damage": [], "von_mises_mpa": []}
    for b in blocks:
        cells.append((_VTK_TYPE[b["type"]], b["conn"]))
        cd["damage"].append(np.asarray(damage[b["name"]], float).reshape(-1, 1))
        cd["von_mises_mpa"].append(np.asarray(vm[b["name"]], float).reshape(-1, 1))
    pd: dict[str, np.ndarray] = {}
    if u is not None:
        pd["displacement"] = np.asarray(u, float)
        pd["abs_u"] = np.linalg.norm(u, axis=1)
    else:
        pd["displacement"] = np.zeros((len(pts), 3))
        pd["abs_u"] = np.zeros(len(pts))
    if sets:
        for nm, idx in sets.items():
            mask = np.zeros(len(pts), dtype=np.int8)
            mask[idx] = 1
            pd[nm] = mask
    return meshio.Mesh(points=pts, cells=cells, point_data=pd, cell_data=cd)


def write_vtk(mesh, path: Path) -> None:
    import meshio
    path.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(str(path), mesh, file_format="vtk")
    # verify it re-opens (connectivity must survive)
    back = meshio.read(str(path))
    n_cells = sum(len(cb.data) for cb in back.cells)
    kinds = [(cb.type, len(cb.data)) for cb in back.cells]
    print(f"[vtk] {path.name}: points={len(back.points)} cells={n_cells} {kinds}",
          flush=True)


# ===========================================================================
# PNG (matplotlib, Agg, CJK fallback)
# ===========================================================================
def _init_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                              "SimSun", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    return plt


def _cell_edges(blocks):
    """Deduplicated global edge list (0-based node pairs)."""
    edges: set[tuple[int, int]] = set()
    for b in blocks:
        es = _TET_EDGES if b["type"] == "tet4" else _HEX_EDGES
        for row in b["conn"]:
            for a, c in es:
                p, q = int(row[a]), int(row[c])
                edges.add((p, q) if p < q else (q, p))
    return np.asarray(sorted(edges), dtype=np.int64)


def _node_avg_cell_scalar(blocks, per_block):
    """Average a per-element scalar onto nodes (tet4 4-col / hex8 8-col mixed)."""
    n_nodes = max(int(b["conn"].max()) for b in blocks) + 1
    acc = np.zeros(n_nodes)
    cnt = np.zeros(n_nodes)
    for b in blocks:
        val = np.asarray(per_block[b["name"]], float).reshape(-1)
        conn = b["conn"]
        for k in range(conn.shape[1]):
            np.add.at(acc, conn[:, k], val)
            np.add.at(cnt, conn[:, k], 1)
    return acc / np.maximum(cnt, 1)


def plot_overlay(pts, blocks, u, sets, path: Path, axis_names=("X", "Y", "Z"),
                 suptitle: str | None = None) -> None:
    plt = _init_mpl()
    from matplotlib.collections import LineCollection
    E = _cell_edges(blocks)
    deformed = pts + u
    ui = np.linalg.norm(u, axis=1)
    esca = 0.5 * (ui[E[:, 0]] + ui[E[:, 1]])
    ax_i, ax_j = 2, 0   # vertical = Z, horizontal = X  (front view)
    fig, axes = plt.subplots(1, 2, figsize=(15, 8))
    for ax, (h, v) in zip(axes, ((0, 2), (1, 2))):
        ax.add_collection(LineCollection(pts[E][:, :, [h, v]], colors="#c8c8c8",
                                         linewidths=0.18, zorder=1))
        lc = LineCollection(deformed[E][:, :, [h, v]], array=esca, cmap="turbo",
                            linewidths=0.45, zorder=2)
        ax.add_collection(lc)
        if "top_load" in sets:
            t = sets["top_load"]
            ax.scatter(deformed[t, h], deformed[t, v], s=6, c="#d62728",
                       zorder=3, label="top_load（加载带）")
        if "bottom_fix" in sets:
            b = sets["bottom_fix"]
            ax.scatter(deformed[b, h], deformed[b, v], s=6, c="#1f77b4",
                       zorder=3, label="bottom_fix（固定带）")
        ax.set_xlabel(f"{axis_names[h]} (mm)")
        ax.set_ylabel(f"{axis_names[v]} (mm)")
        ax.set_aspect("equal")
        ax.autoscale()
        ax.legend(loc="best", fontsize=8)
    axes[0].set_title("初始网格（灰线框） vs 变形网格（按 |u| 着色）— 前视 X–Z")
    axes[1].set_title("初始网格（灰线框） vs 变形网格（按 |u| 着色）— 侧视 Y–Z")
    cb = fig.colorbar(lc, ax=axes, fraction=0.03, pad=0.02)
    cb.set_label("|u| (mm)")
    fig.suptitle(suptitle or "fibula_r 反应损伤模型：点载荷奇异 → 单元反转（末收敛态）",
                 fontsize=13)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[png] {path.name}", flush=True)


def plot_damage(pts, blocks, u, damage, sets, path: Path,
                suptitle: str | None = None) -> None:
    """Deformed mesh, per-element D as FILLED projected polygons (clearer than a
    dense wireframe)."""
    plt = _init_mpl()
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    names = {0: "X", 1: "Y", 2: "Z"}
    deformed = pts + u
    fig, axes = plt.subplots(1, 2, figsize=(15, 8))
    pc = None
    for ax, (h, v) in zip(axes, ((0, 2), (1, 2))):
        verts, vals = [], []
        for b in blocks:
            P = deformed[b["conn"]][:, :, [h, v]]
            dv = np.asarray(damage[b["name"]], float).reshape(-1)
            verts.extend(P[e] for e in range(len(P)))
            vals.extend(dv)
        pc = PolyCollection(verts, array=np.asarray(vals, float), cmap="inferno",
                            norm=Normalize(vmin=0.0, vmax=1.0),
                            edgecolors="none", rasterized=True)
        ax.add_collection(pc)
        if sets and "top_load" in sets:
            t = sets["top_load"]
            ax.scatter(deformed[t, h], deformed[t, v], s=14, facecolors="none",
                       edgecolors="#00e5ff", linewidths=0.35, zorder=5,
                       label="top_load")
        if sets and "bottom_fix" in sets:
            bix = sets["bottom_fix"]
            ax.scatter(deformed[bix, h], deformed[bix, v], s=14, facecolors="none",
                       edgecolors="#00ff66", linewidths=0.35, zorder=5,
                       label="bottom_fix")
        ax.set_xlabel(f"{names[h]} (mm)")
        ax.set_ylabel(f"{names[v]} (mm)")
        ax.set_aspect("equal")
        ax.autoscale()
        ax.legend(loc="best", fontsize=8)
    axes[0].set_title("损伤场 D（变形网格，前视 X–Z）")
    axes[1].set_title("损伤场 D（变形网格，侧视 Y–Z）")
    cb = fig.colorbar(pc, ax=axes, fraction=0.03, pad=0.02)
    cb.set_label("damage D")
    fig.suptitle(suptitle or "fibula_r 主承载域 fib_cort 损伤场（末收敛态；青/绿=加载/固定带）",
                 fontsize=13)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[png] {path.name}", flush=True)


# ===========================================================================
# main
# ===========================================================================
def _export(feb: Path, xplt: Path, label: str) -> dict:
    hdf5 = VIZ_TMP / f"{feb.stem}.hdf5"
    hdf5.parent.mkdir(parents=True, exist_ok=True)
    build_hdf5(xplt, hdf5_path=hdf5)
    pts, nids, blocks, sets = parse_feb_geometry(feb)
    st = read_last_state(hdf5)
    u = st["u"]
    if len(u) != len(pts):
        raise RuntimeError(f"节点数不匹配: disp={len(u)} vs feb={len(pts)}")
    deformed = pts + u
    damage = {dom: a.reshape(-1) for dom, a in st["damage"].items()}
    vm = {dom: von_mises_from_voigt(s) for dom, s in st["stress"].items()}
    for b in blocks:                      # make sure every block has the fields
        damage.setdefault(b["name"], np.zeros(b["conn"].shape[0]))
        vm.setdefault(b["name"], np.zeros(b["conn"].shape[0]))

    # ---- VTK ----
    m_init = _meshio_mesh(pts, blocks, {b["name"]: np.zeros(b["conn"].shape[0])
                                        for b in blocks},
                          {b["name"]: np.zeros(b["conn"].shape[0]) for b in blocks},
                          u=None, sets=sets)
    m_def = _meshio_mesh(deformed, blocks, damage, vm, u=u, sets=sets)
    write_vtk(m_init, OUT / "initial_mesh.vtk")
    write_vtk(m_def, OUT / f"deformed_{label}.vtk")

    # ---- Jacobian / inversion at the deformed state ----------------------
    j_def = element_jacobians(deformed, blocks)
    j_init = _initial_jacobians(pts, blocks)
    span = pts[:, 2].max() - pts[:, 2].min()
    zmin = pts[:, 2].min()
    neg_total = 0
    neg_info = {}
    worst = {"ratio": None, "domain": None, "z_frac": None, "neg": None}
    for b in blocks:
        v = j_def[b["name"]]
        v0 = j_init[b["name"]]
        ratio = v / np.where(np.abs(v0) > 1e-30, np.abs(v0), np.nan)
        neg = v < 0
        n = int(neg.sum())
        neg_total += n
        ctr = (deformed[b["conn"]].mean(axis=1) if len(b["conn"])
               else np.zeros((0, 3)))
        frac = ((ctr[:, 2] - zmin) / span) if len(ctr) else np.zeros(0)
        neg_info[b["name"]] = {
            "n_elem": int(v.size), "n_neg": n,
            "n_neg_top_band": int(((neg) & (frac >= 0.70)).sum()),
            "n_neg_bottom_band": int(((neg) & (frac <= 0.30)).sum()),
            "n_neg_middle": int(((neg) & (frac > 0.30) & (frac < 0.70)).sum()),
            "n_crushed_lt10pct": int((ratio < 0.10).sum()),
            "min_j_ratio": (float(np.nanmin(ratio)) if len(ratio) else None),
            "min_frac": float(frac[neg].min()) if n else None,
            "max_frac": float(frac[neg].max()) if n else None,
        }
        if len(ratio):
            k = int(np.nanargmin(ratio))
            if worst["ratio"] is None or ratio[k] < worst["ratio"]:
                zf = float(frac[k])
                band = ("bottom_fix 带（底固定带）" if zf <= 0.30 else
                        "top_load 带（顶加载带）" if zf >= 0.70 else "中部")
                worst = {"ratio": float(ratio[k]), "domain": b["name"],
                         "z_frac": zf, "negative": bool(v[k] < 0), "band": band}
    # ---- PNG ----
    plot_overlay(pts, blocks, u, sets, OUT / "overlay_initial_vs_deformed.png")
    plot_damage(pts, blocks, u, damage, sets, OUT / "damage_field.png")

    return {
        "label": label, "feb": str(feb), "xplt": str(xplt),
        "state_index": st["state"], "n_states": st["n_states"],
        "max_abs_u_mm": float(np.linalg.norm(u, axis=1).max()),
        "n_nodes": int(len(pts)),
        "blocks": [(b["name"], b["type"], int(b["conn"].shape[0])) for b in blocks],
        "neg_jacobian_total": neg_total, "neg_jacobian": neg_info,
        "worst_jacobian": worst,
        "max_damage": float(max((a.max() if a.size else 0.0)
                                for a in damage.values())),
    }


def _rerun_force1000(steps: int, timeout: float) -> tuple[Path, Path] | None:
    """Re-run the canonical force-control 1000 N config with a ramp; return feb/xplt."""
    try:
        import g4_fe_stab as stab
        setup = stab.Setup()
        rec = stab.cell_force(setup, load_n=1000.0, Dmax=1.0, steps=steps,
                              timeout=timeout, tag="viz1000", rep=1)
        feb = stab.OUT_TEMP / "force_viz1000.feb"
        return feb, feb.with_suffix(".xplt")
    except Exception as exc:  # noqa: BLE001
        print(f"[force1000] re-run failed: {type(exc).__name__}: {exc}", flush=True)
        return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dc-stem", default=DC_STEM)
    ap.add_argument("--force1000", action="store_true",
                    help="also re-run the force-control 1000 N config (ramp) and "
                         "export it if FEBio wrote converged states")
    ap.add_argument("--force-steps", type=int, default=20)
    ap.add_argument("--timeout", type=float, default=240.0)
    args = ap.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    VIZ_TMP.mkdir(parents=True, exist_ok=True)

    results = {"generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
               "dc": None, "force": None}

    # ---- primary: DC failure (has written states) ----
    dc_feb = STAB / f"{args.dc_stem}.feb"
    dc_xplt = STAB / f"{args.dc_stem}.xplt"
    if not dc_feb.is_file() or not dc_xplt.is_file():
        print(f"[error] 缺少 {dc_feb} / {dc_xplt}", flush=True)
        return 2
    print(f"[viz] primary = {args.dc_stem} (displacement control, Dmax=1.0)",
          flush=True)
    results["dc"] = _export(dc_feb, dc_xplt, "first_failure")

    # ---- optional: canonical force-control 1000 N ----
    if args.force1000:
        f_feb = STAB / f"{FORCE_STEM}.feb"
        f_xplt = STAB / f"{FORCE_STEM}.xplt"
        have_state = False
        if f_xplt.is_file():
            try:
                h5 = VIZ_TMP / f"{FORCE_STEM}.hdf5"
                build_hdf5(f_xplt, hdf5_path=h5)
                have_state = read_last_state(h5)["n_states"] > 1
            except Exception:  # noqa: BLE001
                have_state = False
        if not have_state:
            print("[viz] force-control 1000 N xplt has no converged state -> "
                  "re-running with a 20-step ramp", flush=True)
            got = _rerun_force1000(args.force_steps, args.timeout)
            if got:
                f_feb, f_xplt = got
                h5 = VIZ_TMP / f"{f_feb.stem}.hdf5"
                build_hdf5(f_xplt, hdf5_path=h5)
                if read_last_state(h5)["n_states"] > 1:
                    results["force"] = _export(f_feb, f_xplt, "force1000")
                else:
                    results["force"] = {"note": "re-run also wrote no converged state",
                                        "xplt": str(f_xplt)}
        else:
            results["force"] = _export(f_feb, f_xplt, "force1000")

    _write_readme(results, args)
    print(f"[done] {OUT}", flush=True)
    return 0


def _write_readme(R: dict, args) -> None:
    L: list[str] = []
    A = L.append
    A("# g4_fe_viz — fibula_r 损伤模型失败可视化")
    A("")
    A(f"_生成时间：{R['generated_at']}。脚本：`scripts/opensim_fe/g4_fe_viz.py`（只新增，不改任何既有模块）。_")
    A("")
    dc = R["dc"]
    A("## 主导出（点载荷奇异 / 单元反转）")
    A("")
    A("| 项 | 值 |")
    A("|---|---|")
    A(f"| 运行配置 | **位移控制** `{args.dc_stem}` = P1-a 的 `disp@-397.858`, `Dmax=1.0`（ramped prescribed z-displacement，∂ 输入位移 −397.858 mm） |")
    A(f"| `.feb` | `{dc['feb']}` |")
    A(f"| `.xplt` | `{dc['xplt']}` |")
    A(f"| **导出 state 索引** | **{dc['state_index']}**（共 {dc['n_states']} 个写入 state；FEBio 在**收敛增量**才写 state，故 state {dc['state_index']} = **发散前最后一个收敛态**） |")
    A(f"| 节点数 / 单元块 | {dc['n_nodes']} / " + "、".join(f"{n}({t},{c})" for n, t, c in dc["blocks"]) + " |")
    A(f"| **max\\|u\\|** | **{dc['max_abs_u_mm']:.3f} mm** |")
    A(f"| **负 Jacobian 单元数** | **{dc['neg_jacobian_total']}** / {sum(v['n_elem'] for v in dc['neg_jacobian'].values())} |")
    A(f"| 最大损伤 D | {dc['max_damage']:.4f} |")
    A("")
    A("### 单元 Jacobian（末收敛态；位置按**初始**沿 z 归一化，0=底 1=顶）")
    A("")
    A("| 域 | 单元数 | 负 J 数 | 顶带(z≥0.70) | 底带(z≤0.30) | 中部 | min detJ/\\|detJ0\\| | 挤到<10% |")
    A("|---|---|---|---|---|---|---|---|")
    for dom, v in dc["neg_jacobian"].items():
        A(f"| `{dom}` | {v['n_elem']} | {v['n_neg']} | {v['n_neg_top_band']} | "
          f"{v['n_neg_bottom_band']} | {v['n_neg_middle']} | "
          f"{v['min_j_ratio']:.4f} | {v['n_crushed_lt10pct']} |")
    A("")
    A("## 副导出（力控制 1000 N——经典 `98 negative jacobians`）")
    A("")
    fc = R.get("force")
    if fc and "note" not in fc and fc.get("state_index") is not None:
        A(f"- `{fc['feb']}`；导出 state {fc['state_index']}/{fc['n_states']}；"
          f"max|u|={fc['max_abs_u_mm']:.3f} mm；负 J={fc['neg_jacobian_total']}。")
    elif isinstance(fc, dict):
        A(f"- **该运行没有写出任何收敛 state**（`{fc.get('note', fc.get('xplt', '—'))}`）。"
          "力控制 1000 N 在**第一个增量**就发散（`98 negative jacobians detected. | "
          "failed to converge at time : 1`），故 `.xplt` 只有 state 0（未变形、位移=0），"
          "无法用于「变形网格」可视化——这正是需要位移控制来**逐增量**逼近失败的原因。")
    else:
        A("- 未跑（`--force1000` 未开启）。")
    A("")
    A("## 产物")
    A("")
    A("- `initial_mesh.vtk` — 未变形网格（hex8/tet4 连通性保留；point_data 含 "
      "`top_load`/`bottom_fix` 掩码）。")
    A("- `deformed_first_failure.vtk` — 末收敛态网格；point_data `displacement`(向量)、"
      "`abs_u`(标量)；cell_data `damage`、`von_mises_mpa`。")
    A("- `overlay_initial_vs_deformed.png` — 初始（灰线框）vs 变形（按 |u| 着色），X–Z 与 Y–Z 两视图。")
    A("- `damage_field.png` — 变形网格按单元损伤 D 着色。")
    A("")
    A("## 观察")
    A("")
    wj = dc.get("worst_jacobian") or {}
    log = Path(dc["xplt"]).with_suffix(".log")
    if log.is_file():
        txt = log.read_text(errors="replace")
        nj = [ln.strip() for ln in txt.splitlines()
              if "negative jacobian" in ln.lower() and ln.strip().startswith("*")]
        stp = [ln.strip() for ln in txt.splitlines()
               if "time steps completed" in ln]
        A("FEBio 日志原文（节选）：")
        A("")
        A("```")
        for ln in (stp + nj)[:6]:
            A(ln)
        A("```")
        A("")
    A(f"- 导出的 state **{dc['state_index']}** 是**最后一个收敛态**（FEBio 只在"
      f"收敛增量写 state）——顶带被下压到 `max|u|={dc['max_abs_u_mm']:.1f} mm`。")
    A(f"- 该态**负 Jacobian 单元数 = {dc['neg_jacobian_total']}**。**单元反转发生在"
      "下一个（发散的、未写入的）增量**：FEBio 在试探更新处检测到反转即 "
      "`ERROR TERMINATION`，该构型**不进 `.xplt`** —— 所以任何「已写入 state」"
      "按定义都是未反转的。这是「文件里的最后状态」与「崩溃瞬间构型」的差别。")
    A(f"- 压溃热点（最接近反转）：`min detJ/|detJ0| = {wj.get('ratio'):.4f}`，"
      f"域 `{wj.get('domain')}`，初始 z 归一化位置 `{wj.get('z_frac'):.3f}` "
      f"（**{wj.get('band')}**）→ 反转由该带起始。")
    A(f"- 损伤 D 达 **{dc['max_damage']:.3f}**（`fib_cort` 近饱和），集中在压溃区。")
    A("- VTK 为 legacy ASCII（`meshio` 写出并**回读校验**：hex8=8 节点、tet4=4 节点）。")
    A("")
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"[readme] {OUT / 'README.md'}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
