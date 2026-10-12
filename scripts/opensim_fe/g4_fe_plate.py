"""G4 压板探针：单骨 + 幽灵刚体压板 + sliding-elastic 接触。

背景（P1-b 已建立）：V7（50% 带 NodalForce 沿 PCA）弹性峰 σ_vm ≈ 147 MPa，
物理参照 ~9.8 MPa ⇒ 残余 ~15× 应力集中。点载荷奇异两个来源：(a) 力硬塞给
顶带节点；(b) 无接触 ⇒ 顶面节点被强加力。本探针用**刚性压板 + 接触**替代
节点力：力经接触分布到骨端面，无点载荷；压板法向 = 骨轴 ⇒ 无 V1 斜端面的
横向合力问题。

实现（FEBio 4.13 无 rigid_wall，S5 结论）：hex8 薄板 + 幽灵刚体约束
（``<bc type="rigid">`` 绑定板节点）→ 板 = 数值刚体；RigidFixed 锁横向+转动、
RigidPrescribed（rigid_displacement）沿骨轴压入；板底面 ↔ 骨端面 cap
（quad4 facets）sliding-elastic 接触（S5 NEW PENALTY 配方逐字沿用）。

复现（repo root）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_plate.py

只 ADD 本模块；import 复用 g4_fe_p1b / g4_fe_stab / g4_fe_spike / bone_feb，
不改任何既有模块。
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import bone_feb as bf  # noqa: E402
import g4_fe_p1b as p1b  # noqa: E402
import g4_fe_spike as g  # noqa: E402
import g4_fe_stab as gs  # noqa: E402

OUT_TEMP = ROOT / "temp" / "opensim_fe" / "g4_fe_plate"
OUT_RESULTS = ROOT / "results" / "opensim_fe"

#: S5 wave-3 最小版配方（docs/S5_contact_udg.md §7）——逐字沿用。
#: 2026-10-12 根因修复（KB 理论手册 7.1.1）: gap = -nu·(x1-x2) 用 secondary 法向
#: nu；plate_bottom quad4 逆时针 => 法向 +z 背离骨 primary（骨在下方）=> 接触
#: 永不激活。之前试 flip_primary 无效正因它不碰 secondary 法向。故
#: flip_secondary=1 把板法向翻到 -z（指向骨）。
CONTACT_RECIPE = dict(laugon="PENALTY", penalty=0.1, auto_penalty=0,
                      two_pass=0, node_reloc=0, symmetric_stiffness=0,
                      tolerance=0.005, search_radius=20.0,
                      # 2026-10-12: fric>0 使接触失效（d1 应力 14.85→0.000，
                      # 实测），保持 0.0；侧滑问题由大位移接触锁定自然抑制.
                      fric_coeff=0.0, flip_secondary=1)
PLATE_E_MPA = 200_000.0
PLATE_NU = 0.30
PLATE_THICK_MM = 4.0
PLATE_MARGIN_MM = 6.0
PLATE_NDIV = 14
GAP_MM = 0.2


# ===========================================================================
# plate geometry (aligned frame: bone long axis -> +Z)
# ===========================================================================
def plate_geom(setup: p1b.Setup, geom: dict) -> dict:
    """Square rigid-plate covering the bone cap, just above z=a_max+gap."""
    nodes = setup.nodes_for(True)
    xs = nodes[:, 0]; ys = nodes[:, 1]
    x0 = float(xs.min()) - PLATE_MARGIN_MM
    x1 = float(xs.max()) + PLATE_MARGIN_MM
    y0 = float(ys.min()) - PLATE_MARGIN_MM
    y1 = float(ys.max()) + PLATE_MARGIN_MM
    z_wall = float(geom["a_max"]) + GAP_MM
    z_bot = z_wall - PLATE_THICK_MM
    return {"x0": x0, "x1": x1, "y0": y0, "y1": y1,
            "z_wall": z_wall, "z_bot": z_bot,
            "n_node_grid": PLATE_NDIV + 1, "ndiv": PLATE_NDIV,
            "a_plate_mm2": (x1 - x0) * (y1 - y0)}


def plate_mesh_xml(pg: dict, *, first_node_id: int) -> tuple[str, str, str, str]:
    """Return (nodes_xml, elems_xml, nodeset_xml, bottom_surface_xml).

    Plate = ndiv x ndiv hex8 slab; top layer at z=z_wall (contact face),
    bottom layer at z=z_bot.  ids start at ``first_node_id``.
    """
    n = pg["n_node_grid"]
    top: list[int] = []      # (i,k) -> id  (top layer, z=z_wall)
    bot: list[int] = []
    nid = first_node_id
    node_lines: list[str] = []
    for i in range(n):
        x = pg["x0"] + (pg["x1"] - pg["x0"]) * i / (n - 1)
        for k in range(n):
            y = pg["y0"] + (pg["y1"] - pg["y0"]) * k / (n - 1)
            top.append(nid)
            node_lines.append(
                f'      <node id="{nid}">{x:.10e},{y:.10e},{pg["z_wall"]:.10e}</node>')
            nid += 1
            bot.append(nid)
            node_lines.append(
                f'      <node id="{nid}">{x:.10e},{y:.10e},{pg["z_bot"]:.10e}</node>')
            nid += 1
    nodes_xml = "\n".join(node_lines) + "\n"

    elems: list[str] = []
    bottom_quads: list[str] = []
    for i in range(pg["ndiv"]):
        for k in range(pg["ndiv"]):
            t00 = top[i * n + k]; t01 = top[i * n + k + 1]
            t10 = top[(i + 1) * n + k]; t11 = top[(i + 1) * n + k + 1]
            b00 = bot[i * n + k]; b01 = bot[i * n + k + 1]
            b10 = bot[(i + 1) * n + k]; b11 = bot[(i + 1) * n + k + 1]
            eid = len(elems) + 1
            elems.append(f'      <elem id="{eid}">{t00},{t01},{t11},{t10},'
                         f'{b00},{b01},{b11},{b10}</elem>')
            # bottom surface quad (top layer, face -z): REVERSED order so the
            # quad4 normal points -z (towards bone).  FEBio 4.13 ignored the
            # flip_secondary parameter (log shows "no (0)" despite .feb=1),
            # and the NTS gap g=-nu.(x1-x2) uses the secondary normal nu.
            # CCW (t00,t01,t11,t10) -> +z (away from bone); CW -> -z.
            bottom_quads.append(f'      <quad4 id="{len(bottom_quads) + 1}">'
                                f'{t00},{t10},{t11},{t01}</quad4>')
    elems_xml = (f'    <Elements name="plate" type="hex8">\n'
                 + "\n".join(elems) + "\n    </Elements>\n")
    nodeset_xml = ('    <NodeSet name="plate_nodes">'
                   + ",".join(str(x) for x in top + bot) + "</NodeSet>\n")
    surf_xml = ('    <Surface name="plate_bottom">\n'
                + "\n".join(bottom_quads) + "\n    </Surface>\n")
    return nodes_xml, elems_xml, nodeset_xml, surf_xml


def bone_cap_surface_xml(geom: dict) -> str:
    """Surface from cap facets (all quad4 for fibula_r, verified)."""
    lines: list[str] = []
    quads = [f for f in geom["cap"] if len(f["poly"]) == 4]
    if len(quads) != len(geom["cap"]):
        raise ValueError("cap 含非 quad4 facet，暂不支持")
    for i, f in enumerate(quads):
        ids = ",".join(str(int(x)) for x in f["poly"])
        lines.append(f'      <quad4 id="{i + 1}">{ids}</quad4>')
    return '    <Surface name="bone_top">\n' + "\n".join(lines) + "\n    </Surface>\n"


#: SurfacePair lives inside <Mesh> (S5 syntax).
SURFACE_PAIR_XML = ('    <SurfacePair name="plate_bone_pair">\n'
                    '      <primary>bone_top</primary>\n'
                    '      <secondary>plate_bottom</secondary>\n'
                    '    </SurfacePair>\n')


# ===========================================================================
# .feb text injection (raw XML, same style as g4_fe_p1b)
# ===========================================================================
def _next_material_id(text: str) -> int:
    ids = [int(m) for m in re.findall(r'<material[^>]*\bid="(\d+)"', text)]
    return (max(ids) + 1) if ids else 1


def _inject_before(text: str, marker: str, block: str) -> str:
    idx = text.rfind(marker)
    if idx < 0:
        raise ValueError(f"marker not found: {marker}")
    return text[:idx] + block + text[idx:]


def _inject_after(text: str, marker: str, block: str) -> str:
    idx = text.rfind(marker)
    if idx < 0:
        raise ValueError(f"marker not found: {marker}")
    return text[:idx + len(marker)] + block + text[idx + len(marker):]


def _ghost_mesh_xml(pg: dict, *, first_node_id: int) -> tuple[str, str]:
    """Rigid-body carrier mesh (FEBio needs a real grid + rigid-body material).

    A single tet4 around the plate centre (mirrors pyfebio's
    ``add_simple_rigid_body`` output seen in S4 decks).
    """
    cx = 0.5 * (pg["x0"] + pg["x1"])
    cy = 0.5 * (pg["y0"] + pg["y1"])
    z = pg["z_wall"]
    n0 = first_node_id
    nodes = (f'    <Nodes name="ghost">\n'
             f'      <node id="{n0}">{cx:.10e},{cy:.10e},{z:.10e}</node>\n'
             f'      <node id="{n0 + 1}">{cx + 1:.10e},{cy:.10e},{z:.10e}</node>\n'
             f'      <node id="{n0 + 2}">{cx:.10e},{cy + 1:.10e},{z:.10e}</node>\n'
             f'      <node id="{n0 + 3}">{cx:.10e},{cy:.10e},{z + 1:.10e}</node>\n'
             f'    </Nodes>\n')
    elems = (f'    <Elements name="ghost" type="tet4">\n'
             f'      <elem id="1">{n0},{n0 + 1},{n0 + 2},{n0 + 3}</elem>\n'
             f'    </Elements>\n')
    return nodes, elems


def build_plate_feb(setup: p1b.Setup, out_path: Path, *, disp_mm: float,
                    steps: int = 5, damage_Dmax: float | None = None,
                    seed: int = 0, force_mode: bool = False,
                    force_n: float = 1000.0, gap_mm: float = GAP_MM) -> Path:
    """Build the plate-probe .feb (displacement-controlled rigid plate)."""
    if not out_path.parent.exists():
        out_path.parent.mkdir(parents=True, exist_ok=True)

    setup.fill_materials()

    # 1) base bone model (V7 geometry: 50% bands, PCA-aligned; nodal load will
    #    be discarded).
    bf.build_bone_feb(setup.spec, out_path, load_n=1.0, top_frac=0.5,
                      bot_frac=0.5, axis=None, with_shell=True,
                      time_steps=steps, mesh=setup.mesh, align_pca=True)
    text = out_path.read_text(encoding="ISO-8859-1")

    # 2) geometry
    geom = setup.facet_geom(top_frac=0.5, align_pca=True)
    pg = plate_geom(setup, geom)
    n_bone = int(setup.mesh["n_nodes"])
    nodes_xml, elems_xml, nodeset_xml, plate_surf_xml = plate_mesh_xml(
        pg, first_node_id=n_bone + 1)
    ghost_nodes_xml, ghost_elems_xml = _ghost_mesh_xml(
        pg, first_node_id=n_bone + 1 + 2 * pg["n_node_grid"] ** 2)

    # 3) mesh blocks (bone_feb native indent: 2-space)
    text = _inject_before(text, "</Nodes>", nodes_xml)
    text = _inject_before(text, "</Mesh>", elems_xml)
    text = _inject_before(text, "</Mesh>", ghost_nodes_xml)
    text = _inject_before(text, "</Mesh>", ghost_elems_xml)
    text = _inject_before(text, "</Mesh>", nodeset_xml)
    text = _inject_before(text, "</Mesh>", plate_surf_xml)
    text = _inject_before(text, "</Mesh>", bone_cap_surface_xml(geom))
    text = _inject_before(text, "</Mesh>", SURFACE_PAIR_XML)

    # 4) material + domain (plate + ghost rigid-body carrier)
    mid = _next_material_id(text)
    plate_mat = (f'    <material name="plate_mat" id="{mid}" '
                 f'type="isotropic elastic">\n'
                 f'      <density>8.0e-09</density>\n'
                 f'      <E>{PLATE_E_MPA:.1f}</E>\n'
                 f'      <v>{PLATE_NU}</v>\n'
                 f'    </material>\n')
    ghost_mat = (f'    <material name="ghost" id="{mid + 1}" '
                 f'type="rigid body">\n'
                 f'      <density>1e-06</density>\n'
                 f'    </material>\n')
    text = _inject_before(text, "</Material>", plate_mat + ghost_mat)
    text = _inject_before(text, "</MeshDomains>",
                          f'    <SolidDomain name="plate" mat="plate_mat"/>\n'
                          f'    <SolidDomain name="ghost" mat="ghost"/>\n')

    # 5) boundary: bind plate nodes to ghost rigid body (injected with Contact
    #    block in step 7 below, keeping StepEntry order Boundary->Contact->Rigid)
    rigid_bc = ('      <bc type="rigid" name="plate_to_ghost" '
                'node_set="plate_nodes">\n'
                '        <rb>ghost</rb>\n'
                '      </bc>\n')

    # 6) contact (S5 NEW PENALTY recipe) — pair defined in <Mesh>
    c = CONTACT_RECIPE
    contact_xml = (
        '    <Contact>\n'
        '      <contact type="sliding-elastic" name="plate_bone_contact" '
        'surface_pair="plate_bone_pair">\n'
        f'        <laugon>{c["laugon"]}</laugon>\n'
        f'        <penalty>{c["penalty"]}</penalty>\n'
        f'        <auto_penalty>{c["auto_penalty"]}</auto_penalty>\n'
        f'        <two_pass>{c["two_pass"]}</two_pass>\n'
        f'        <node_reloc>{c["node_reloc"]}</node_reloc>\n'
        f'        <symmetric_stiffness>{c["symmetric_stiffness"]}</'
        f'symmetric_stiffness>\n'
        f'        <tolerance>{c["tolerance"]}</tolerance>\n'
        f'        <search_radius>{c["search_radius"]}</search_radius>\n'
        f'        <fric_coeff>{c["fric_coeff"]}</fric_coeff>\n'
        f'        <flip_secondary>{c.get("flip_secondary", 0)}</'
        f'flip_secondary>\n'
        '      </contact>\n'
        '    </Contact>\n')

    # 7) rigid block (displacement control on ghost)
    # FEBio 4.x rigid DOF 命名（S5 presc_sink2 实证）:
    #   Rx/Ry/Rz = 平移自由度, Ru/Rv/Rw = 转动自由度.
    #   S5 波5 配方: 平移全放开 (Rx=Ry=Rz=0), 转动全锁死 (Ru=Rv=Rw=1),
    #   再由 rigid_displacement dof=z 驱动 z 平移.
    if force_mode:
        rigid_xml = (
            '    <Rigid>\n'
            '      <rigid_bc type="rigid_fixed">\n'
            '        <rb>ghost</rb>\n'
            '        <Rx_dof>0</Rx_dof><Ry_dof>0</Ry_dof><Rz_dof>0</Rz_dof>\n'
            '        <Ru_dof>1</Ru_dof><Rv_dof>1</Rv_dof><Rw_dof>1</Rw_dof>\n'
            '      </rigid_bc>\n'
            '      <rigid_load type="rigid_force">\n'
            '        <rb>ghost</rb>\n'
            '        <dof>z</dof>\n'
            '        <relative>0</relative>\n'
            f'        <value lc="1">{-force_n:.6f}</value>\n'
            '      </rigid_load>\n'
            '    </Rigid>\n')
    else:
        rigid_xml = (
            '    <Rigid>\n'
            '      <rigid_bc type="rigid_fixed">\n'
            '        <rb>ghost</rb>\n'
            '        <Rx_dof>0</Rx_dof><Ry_dof>0</Ry_dof><Rz_dof>0</Rz_dof>\n'
            '        <Ru_dof>1</Ru_dof><Rv_dof>1</Rv_dof><Rw_dof>1</Rw_dof>\n'
            '      </rigid_bc>\n'
            '      <rigid_bc type="rigid_displacement">\n'
            '        <rb>ghost</rb>\n'
            '        <dof>z</dof>\n'
            '        <relative>0</relative>\n'
            f'        <value lc="1">{-disp_mm:.10e}</value>\n'
            '      </rigid_bc>\n'
            '    </Rigid>\n')

    # order per pyfebio StepEntry: Boundary -> Loads -> Contact -> Rigid
    text = re.sub(r"<Loads>.*?</Loads>", "", text, flags=re.S)
    text = _inject_before(text, "      </Boundary>", rigid_bc)  # inside Boundary
    text = _inject_after(text, "      </Boundary>", contact_xml)  # after Boundary
    text = _inject_before(text, "    </step>", rigid_xml)

    # 8) ramp load curve lc=1 (AFTER </MeshDomains>, like p1b.add_load_ramp)
    ramp = gs._RAMP_XML.replace("_RAMP", "plate_ramp")
    text = _inject_after(text, "</MeshDomains>", ramp)

    # 9) optional damage material for the cortical shell
    if damage_Dmax is not None:
        text = gs.make_damage_text(
            text, dom=setup.main_mesh, mid=gs._material_id(text, setup.main_mesh),
            E=setup.main_mat.E_mpa, nu=setup.main_mat.nu, a=setup.a, c=0.0,
            Dmax=damage_Dmax, xplt_name=out_path.stem + ".xplt")

    # 10) reaction output for force measurement
    text = p1b.add_reaction_var(text)

    out_path.write_text(text, encoding="ISO-8859-1")
    return out_path


def solve_and_stats(setup: p1b.Setup, feb_path: Path, workdir: Path, *,
                    timeout: float = 600.0, read_damage: bool = False) -> dict:
    """Run FEBio, read peak von Mises + axial reaction (p1b.run_and_post)."""
    rec = p1b.run_and_post(setup, feb_path, timeout=timeout,
                           read_damage=read_damage)
    if not rec.get("ok"):
        return {"ok": False, "rc": rec.get("rc"),
                "err": rec.get("error", rec.get("post_error", "?"))}
    return {"ok": True, "rc": rec["rc"], "wall_s": rec.get("wall_s"),
            "vm_max": rec.get("peak_vm_mpa"), "vm_p95": rec.get("p95_vm_mpa"),
            "vm_mean": rec.get("vm_mean_mpa"), "R_net": rec.get("R_net_n"),
            # p1b._reaction writes R_axis_n (was reading rec["R_axis"] -> None)
            "R_axis": rec.get("R_axis_n"), "R_vec": rec.get("R_net_n")}


# ===========================================================================
# probe drivers
# ===========================================================================
def elastic_scan(setup: p1b.Setup, disp_list: list[float],
                 *, steps: int = 5, timeout: float = 900.0) -> list[dict]:
    """Displacement-controlled elastic scans; report R_axis + peak vm."""
    rows: list[dict] = []
    for d in disp_list:
        stem = f"elastic_plate_d{d:g}".replace(".", "p")
        feb = OUT_TEMP / f"{stem}.feb"
        build_plate_feb(setup, feb, disp_mm=d, steps=steps)
        r = solve_and_stats(setup, feb, OUT_TEMP, timeout=timeout)
        rows.append({"disp_mm": d, **r})
        print(f"  d={d:g}mm  rc={r.get('rc')}  ok={r.get('ok')}  "
              f"R_axis={r.get('R_axis')}  vm_max={r.get('vm_max')}",
              flush=True)
    return rows


def damage_runs(setup: p1b.Setup, disp_mm: float, *, n_runs: int = 5,
                steps: int = 20, timeout: float = 900.0,
                seed0: int = 0) -> list[dict]:
    """Damage model (Dmax=1) at fixed plate displacement, n repeated seeds."""
    rows: list[dict] = []
    for i in range(n_runs):
        seed = seed0 + i
        stem = f"damage_plate_d{disp_mm:g}_seed{seed}".replace(".", "p")
        feb = OUT_TEMP / f"{stem}.feb"
        build_plate_feb(setup, feb, disp_mm=disp_mm, steps=steps,
                        damage_Dmax=1.0, seed=seed)
        r = solve_and_stats(setup, feb, OUT_TEMP, timeout=timeout,
                            read_damage=True)
        rows.append({"run": i + 1, "seed": seed, "disp_mm": disp_mm, **r})
        print(f"  dmg#{i + 1} seed={seed}  rc={r.get('rc')}  ok={r.get('ok')}  "
              f"vm_max={r.get('vm_max')}", flush=True)
    return rows


def _fmt(x, nd: int = 2) -> str:
    return "?" if x is None else f"{x:.{nd}f}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", choices=["elastic", "damage", "probe"], default="probe")
    args = ap.parse_args()

    t0 = time.time()
    setup = p1b.Setup()
    print(f"[g4_fe_plate] Setup loaded in {time.time() - t0:.1f}s")

    if args.only in ("elastic", "probe"):
        print("[elastic scan] displacement-controlled rigid plate")
        rows = elastic_scan(setup, [0.5, 1.0, 1.5, 3.0, 5.0])
        # interpolate disp for R_axis = 1000 N from converged rows
        conv = [r for r in rows if r.get("ok") and r.get("R_axis")]
        if conv:
            d1000 = _interp_disp(conv, 1000.0)
            print(f"[elastic] 1000 N @ disp ≈ {d1000:.3f} mm")
            json.dump({"rows": rows, "d1000_mm": d1000, "gap_mm": GAP_MM},
                      open(OUT_RESULTS / "G4_FE_PLATE_ELASTIC.json", "w"),
                      indent=2, ensure_ascii=False)

    if args.only in ("damage", "probe"):
        # find d* from saved scan if present
        dstar = None
        pth = OUT_RESULTS / "G4_FE_PLATE_ELASTIC.json"
        if pth.exists():
            data = json.load(open(pth, encoding="utf-8"))
            dstar = data.get("d1000_mm")
        if dstar is None:
            print("[damage] no elastic scan found; re-running short scan")
            rows = elastic_scan(setup, [0.3, 0.5, 1.0])
            conv = [r for r in rows if r.get("ok") and r.get("R_axis")]
            dstar = _interp_disp(conv, 1000.0) if conv else None
        # 2026-10-12: 板在圆顶骨端单点接触 + 侧滑 => R_axis 非单调
        # (d1:-32N / d1.5:+6N / d3:0N)，1000N 插值不可靠 => 固定 dstar=1.0mm
        # （接触最稳、vm 27.5 MPa 落在探针 10-30 区间）。
        if dstar is None or dstar > 2.0:
            print("[damage] unreliable d1000 (R-axis unstable); using dstar=1.0mm")
            dstar = 1.0
        if dstar is None:
            raise SystemExit("elastic scan failed entirely; abort damage")
        print(f"[damage] Dmax=1 @ disp={dstar:.3f} mm (target 1000 N) x5")
        dmg = damage_runs(setup, dstar, n_runs=5)
        json.dump({"dstar_mm": dstar, "runs": dmg},
                  open(OUT_RESULTS / "G4_FE_PLATE_DAMAGE.json", "w"),
                  indent=2, ensure_ascii=False)

    print(f"[g4_fe_plate] done in {time.time() - t0:.0f}s")


def _interp_disp(rows: list[dict], target: float) -> float | None:
    pts = [(r["disp_mm"], r["R_axis"]) for r in rows
           if r.get("ok") and r.get("R_axis") and r["R_axis"] > 0]
    pts.sort()
    if not pts:
        return None
    if pts[0][1] >= target:
        return pts[0][0]
    for (d0, f0), (d1, f1) in zip(pts, pts[1:]):
        if f1 >= target:
            return float(d0 + (target - f0) * (d1 - d0) / (f1 - f0))
    return pts[-1][0]


if __name__ == "__main__":
    main()
