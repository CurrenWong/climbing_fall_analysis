# -*- coding: utf-8 -*-
"""G4 软适配层探针：单骨 + 贴合圆顶的软垫 + 刚体压头（ghost 刚性驱动）。

背景（plan v2.6 §9.11，PASS·有条件）：平板压圆顶 = 单点接触（反力仅 6-43 N、
板侧滑、d1000 插值不可靠）。本探针在骨端圆顶加**贴合软垫**：垫底面 = 圆顶面片
顶点 + gap（与骨 bone_top 共形，Aproj=230 mm²，是单点接触的 ~100×），垫顶面 =
+thick，hex8 每面片一个单元；垫节点 rigid 绑定 ghost。接触 bone_top(primary)
↔ pad_bottom(secondary)，pad_bottom quad4 **手工反转**（法向指向骨；FEBio 4.13
静默忽略 flip_secondary，须直接改元素顺序 —— §9.11 根因）。

== 结论（2026-10-12，探针判定 FAIL：软适配层方案 (a) 被证伪）==

1. **几何/装配全部打通**：接触激活、骨应力出现（vm 5.2-5.7 MPa）、FEBio 接受网格。
   但 1000 N 目标达不到：位移控 d=2mm 即失败；力控 250N 就上百负 Jacobian。
2. **物理根因**：贴合垫把载荷摊到整个骨端（230 mm²），细长腓骨于是**弯曲让位**
   而非局部受压 —— 1000 N 需要厘米级挠度。垫本身（3mm 薄层）也随之翻面。
   ⇒ §9.11 建议的"软适配层贴合圆顶"不能作为 P1-b 立项路径。
3. **横向空空间（重大，既往结论需加限定）**：S5 平板配方平移全放开
   （Rx=Ry=Rz=0），fric=0 时横向 DOF 是零刚度空空间，求解器往里填垃圾
   （实测：平板 1e9 mm、垫 19-40 mm）。小垫直接滑出圆顶 ⇒ 零接触零应力；
   轴向导向（Rx=Ry=1）后接触正常。但**平板探针在导向后每步都不收敛** ——
   §9.11 的 27.5 MPa 是在横向漂移状态下测得的，保留该结论必须带上此限定。

== 过程中踩到并已修正的 4 个坑（详见 .learnings）==

- ``facet_geom()['cap'][i]['poly']`` 是 **1-based** 网格索引：取坐标必须
  ``nodes[poly-1]``（poly 直取面积误差 1313%，poly-1 误差 0.6%）。
- ``bone_feb._hex8_signed_volumes`` 收 **1-based** 局部索引：传 0-based 会让
  idx 0 回绕读末行，体积值变成垃圾（同一网格跨脚本算出 ±2 数量级差异）。
- 仓库 6-tet 体积检查只是**符号代理**，对强翘曲面片不足以判 FEBio 接受性；
  真判据是三线性 hex8 在 8 角点 + 2×2×2 Gauss 点的 detJ > 0。
- FEBio ``rigid_load`` 的平移 DOF 拼写是 ``Rz``，``dof=z`` 只对
  ``rigid_displacement`` 有效（``g4_fe_plate.py`` force_mode 里同一个 bug）。

复现（repo root）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_pad.py --only elastic --d 1.0

只 ADD 本模块；import 复用 g4_fe_plate / g4_fe_p1b / g4_fe_stab / bone_feb，
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

import g4_fe_plate as gp  # noqa: E402  (reuse CONTACT_RECIPE / injections / stats)
import g4_fe_p1b as p1b  # noqa: E402
import g4_fe_stab as gs  # noqa: E402
import bone_feb as bf  # noqa: E402

OUT_TEMP = ROOT / "temp" / "opensim_fe" / "g4_fe_pad"
OUT_RESULTS = ROOT / "results" / "opensim_fe"

#: 圆顶面片截止：所有顶点 z >= PAD_Z0 的面片构成贴合区。poly 为 **1-based**
#: 网格索引（2026-10-12 实测：poly-1 面积误差 0.6%，0-based 误差 1313%）。
#: 实测 fibula_r：Z0=92 → 24 面片、A=348.6 mm²、Aproj=230.4 mm²，
#: 且三线性 hex8 全部正 detJ（minCorner 0.42 / minGauss 0.92）；Z0=90 亦可
#: (34 面片, Aproj 316.7) 但 minGauss 0.62 偏临界；Z0≤88 含近乎竖直的
#: 肩部面片（axdiv≈0）⇒ 拉伸单元翻面。
PAD_Z0 = 92.0
PAD_GAP_MM = 0.05      # 垫底与骨面间隙（共形，小间隙快速建立接触）
PAD_THICK_MM = 3.0     # 垫厚（沿骨轴 z 均匀偏移）
#: 2026-10-12 实测：E=20 MPa 时垫被接触罚压穿（24 个垫单元翻转，日志
#: "24 negative jacobians"）⇒ 贴合由**几何共形**实现、不靠软化；垫刚度只需
#: 远低于骨（18500）以摊开载荷即可。默认 2000 MPa（≈骨/9）。
PAD_E_MPA = 2000.0
PAD_NU = 0.30

#: .feb stem prefix (CLI --tag changes it so parameter sweeps don't overwrite)
STEM = "pad"

#: True  -> bind **all** pad nodes to the ghost (pad = rigid conforming platen,
#:          exactly the working plate-probe recipe); False -> bind only the top
#:          face (pad deforms).
PAD_BIND_ALL = False


# ===========================================================================
# dome facets + pad geometry (aligned frame: bone long axis -> +Z)
# ===========================================================================
def dome_facets(setup: p1b.Setup, geom: dict, *, z0: float = PAD_Z0) -> list[dict]:
    """Cap facets whose ALL vertices are at z >= z0 (the proximal dome).

    NOTE: facet_geom()['cap'][i]['poly'] holds **1-based** mesh node indices
    (verified 2026-10-12: mean rel area err 0.006 vs 13.1 for 0-based), so
    coordinates come from nodes[poly - 1].  Writing poly straight into the
    .feb is correct (FEBio node ids are 1-based too).
    """
    nodes = setup.nodes_for(True)
    sel = []
    for f in geom["cap"]:
        zs = np.array([nodes[int(i) - 1, 2] for i in f["poly"]])
        if float(zs.min()) >= z0:
            sel.append(f)
    return sel


def pad_mesh_xml(setup: p1b.Setup, geom: dict, *, first_node_id: int
                 ) -> tuple[str, str, str, str, dict]:
    """Return (nodes_xml, elems_xml, pad_top_nodeset_xml, pad_bottom_surf_xml,
    meta).  One hex8 per dome facet: bottom face = dome+gap, top face = +thick.
    Bottom surface quads REVERSED (normal -z pointing to bone)."""
    facets = dome_facets(setup, geom)
    if not facets:
        raise ValueError(f"no dome facets at z0={PAD_Z0}")
    nodes = setup.nodes_for(True)

    # coordinate keys (dedup across shared facet vertices)
    bkey: dict[int, tuple] = {}   # bone mesh node id -> bottom key
    tkey: dict[int, tuple] = {}   # bone mesh node id -> top key
    for f in facets:
        for i in f["poly"]:
            p = nodes[int(i) - 1]
            bkey[int(i)] = (float(p[0]), float(p[1]), float(p[2]) + PAD_GAP_MM)
            tkey[int(i)] = (float(p[0]), float(p[1]),
                            float(p[2]) + PAD_GAP_MM + PAD_THICK_MM)

    uniq_b = sorted({k for k in bkey.values()},
                    key=lambda k: (round(k[0], 6), round(k[1], 6),
                                   round(k[2], 6)))
    uniq_t = sorted({k for k in tkey.values()},
                    key=lambda k: (round(k[0], 6), round(k[1], 6),
                                   round(k[2], 6)))
    all_coords = uniq_b + uniq_t
    # NOTE: bone_feb._tet_signed_volumes takes **1-based** local indices
    bidx = {k: j + 1 for j, k in enumerate(uniq_b)}        # 1-based local idx
    tidx = {k: len(uniq_b) + j + 1 for j, k in enumerate(uniq_t)}

    # local hexes (1-based), fix orientation with the module's own validator
    hex_local: list[list[int]] = []
    for f in facets:
        poly = [int(x) for x in f["poly"]]
        hex_local.append([bidx[bkey[i]] for i in poly]
                         + [tidx[tkey[i]] for i in poly])
    vols = bf._hex8_signed_volumes(np.array(all_coords, dtype=float),
                                   np.array(hex_local, dtype=np.int64))
    n_flip = int((vols < 0).sum())
    if n_flip:
        fixed = []
        for k, h in enumerate(hex_local):
            if vols[k] < 0:
                h = h[4:8] + h[0:4]          # swap bottom/top halves
            fixed.append(h)
        hex_local = fixed

    # allocate node ids (bottom then top, matching all_coords order)
    id_by_local = {j + 1: first_node_id + j for j in range(len(all_coords))}
    bnode = {k: id_by_local[bidx[k]] for k in bkey.values()}
    tnode = {k: id_by_local[tidx[k]] for k in tkey.values()}
    node_lines = [f'      <node id="{first_node_id + j}">'
                  f'{all_coords[j][0]:.10e},{all_coords[j][1]:.10e},'
                  f'{all_coords[j][2]:.10e}</node>'
                  for j in range(len(all_coords))]
    nodes_xml = "\n".join(node_lines) + "\n"

    elems: list[str] = []
    bottom_quads: list[str] = []
    for f, h in zip(facets, hex_local):
        poly = [int(x) for x in f["poly"]]
        b = [id_by_local[j] for j in h[0:4]]     # element bottom (after fix)
        t = [id_by_local[j] for j in h[4:8]]     # element top
        eid = len(elems) + 1
        elems.append(f'      <elem id="{eid}">{b[0]},{b[1]},{b[2]},{b[3]},'
                     f'{t[0]},{t[1]},{t[2]},{t[3]}</elem>')
        # secondary surface: pad BOTTOM geometry (bnode), reversed winding
        bb = [bnode[bkey[i]] for i in poly]
        bottom_quads.append(f'      <quad4 id="{len(bottom_quads) + 1}">'
                            f'{bb[0]},{bb[3]},{bb[2]},{bb[1]}</quad4>')

    elems_xml = ('    <Elements name="pad" type="hex8">\n'
                 + "\n".join(elems) + "\n    </Elements>\n")
    bind_ids = (sorted(set(bnode.values()) | set(tnode.values())) if PAD_BIND_ALL
                else sorted(tnode.values()))
    nodeset_xml = ('    <NodeSet name="pad_bind">'
                   + ",".join(str(v) for v in bind_ids) + "</NodeSet>\n")
    surf_xml = ('    <Surface name="pad_bottom">\n'
                + "\n".join(bottom_quads) + "\n    </Surface>\n")
    meta = {"n_pad_nodes": len(all_coords), "n_pad_elems": len(elems),
            "n_dome_facets": len(facets), "n_hex_flipped": n_flip,
            "a_pad_mm2": float(sum(f["area"] for f in facets)),
            "a_pad_proj_mm2": float(sum(f["area"] * f["axdiv"] for f in facets))}
    return nodes_xml, elems_xml, nodeset_xml, surf_xml, meta


def bone_dome_surface_xml(setup: p1b.Setup, geom: dict) -> str:
    """Primary surface = the same dome facets (original winding)."""
    facets = dome_facets(setup, geom)
    lines = []
    for i, f in enumerate(facets):
        ids = ",".join(str(int(x)) for x in f["poly"])
        lines.append(f'      <quad4 id="{i + 1}">{ids}</quad4>')
    return '    <Surface name="bone_top">\n' + "\n".join(lines) + "\n    </Surface>\n"


SURFACE_PAIR_XML = ('    <SurfacePair name="bone_pad_pair">\n'
                    '      <primary>bone_top</primary>\n'
                    '      <secondary>pad_bottom</secondary>\n'
                    '    </SurfacePair>\n')


def _dome_pg(setup: p1b.Setup, geom: dict) -> dict:
    """Ghost-carrier geometry dict (centered on the dome top)."""
    facets = dome_facets(setup, geom)
    nodes = setup.nodes_for(True)
    xs = [float(nodes[int(i) - 1, 0]) for f in facets for i in f["poly"]]
    ys = [float(nodes[int(i) - 1, 1]) for f in facets for i in f["poly"]]
    zs = [float(nodes[int(i) - 1, 2]) for f in facets for i in f["poly"]]
    return {"x0": min(xs) - 1.0, "x1": max(xs) + 1.0,
            "y0": min(ys) - 1.0, "y1": max(ys) + 1.0,
            "z_wall": max(zs) + PAD_GAP_MM + PAD_THICK_MM}


def build_pad_feb(setup: p1b.Setup, out_path: Path, *, disp_mm: float | None = None,
                  steps: int = 5, damage_Dmax: float | None = None,
                  seed: int = 0, force_n: float | None = None) -> Path:
    """Build the pad-probe .feb.

    Displacement control (``disp_mm``) drives the ghost's z; force control
    (``force_n``) instead applies that axial force to the ghost.  Force control
    is the physically usable one here: a conforming pad spreads the load over
    ~230 mm^2 of dome, so the slender fibula *bends away* instead of compressing
    locally — at 1000 N that needs centimetres of deflection, which
    displacement control cannot reach (d=2 mm already failed).
    """
    if not out_path.parent.exists():
        out_path.parent.mkdir(parents=True, exist_ok=True)

    setup.fill_materials()

    # 1) base bone model (V7 geometry; nodal load discarded)
    bf.build_bone_feb(setup.spec, out_path, load_n=1.0, top_frac=0.5,
                      bot_frac=0.5, axis=None, with_shell=True,
                      time_steps=steps, mesh=setup.mesh, align_pca=True)
    text = out_path.read_text(encoding="ISO-8859-1")

    # 2) pad geometry
    geom = setup.facet_geom(top_frac=0.5, align_pca=True)
    n_bone = int(setup.mesh["n_nodes"])
    nodes_xml, elems_xml, nodeset_xml, pad_surf_xml, meta = pad_mesh_xml(
        setup, geom, first_node_id=n_bone + 1)
    ghost_pg = _dome_pg(setup, geom)
    ghost_nodes_xml, ghost_elems_xml = gp._ghost_mesh_xml(
        ghost_pg, first_node_id=n_bone + 1 + meta["n_pad_nodes"])

    # 3) mesh blocks
    text = gp._inject_before(text, "</Nodes>", nodes_xml)
    text = gp._inject_before(text, "</Mesh>", elems_xml)
    text = gp._inject_before(text, "</Mesh>", ghost_nodes_xml)
    text = gp._inject_before(text, "</Mesh>", ghost_elems_xml)
    text = gp._inject_before(text, "</Mesh>", nodeset_xml)
    text = gp._inject_before(text, "</Mesh>", pad_surf_xml)
    text = gp._inject_before(text, "</Mesh>", bone_dome_surface_xml(setup, geom))
    text = gp._inject_before(text, "</Mesh>", SURFACE_PAIR_XML)

    # 4) material + domain (pad + ghost rigid-body carrier)
    mid = gp._next_material_id(text)
    pad_mat = (f'    <material name="pad_mat" id="{mid}" '
               f'type="isotropic elastic">\n'
               f'      <density>8.0e-09</density>\n'
               f'      <E>{PAD_E_MPA:.1f}</E>\n'
               f'      <v>{PAD_NU}</v>\n'
               f'    </material>\n')
    ghost_mat = (f'    <material name="ghost" id="{mid + 1}" '
                 f'type="rigid body">\n'
                 f'      <density>1e-06</density>\n'
                 f'    </material>\n')
    text = gp._inject_before(text, "</Material>", pad_mat + ghost_mat)
    text = gp._inject_before(text, "</MeshDomains>",
                             f'    <SolidDomain name="pad" mat="pad_mat"/>\n'
                             f'    <SolidDomain name="ghost" mat="ghost"/>\n')

    # 5) boundary: bind pad top nodes to ghost rigid body
    rigid_bc = ('      <bc type="rigid" name="pad_to_ghost" '
                'node_set="pad_bind">\n'
                '        <rb>ghost</rb>\n'
                '      </bc>\n')

    # 6) contact (S5 NEW PENALTY recipe) — pair defined in <Mesh>
    c = gp.CONTACT_RECIPE
    contact_xml = (
        '    <Contact>\n'
        '      <contact type="sliding-elastic" name="bone_pad_contact" '
        'surface_pair="bone_pad_pair">\n'
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

    # 7) rigid block: axial-guided platen.  2026-10-12: the S5 plate recipe
    # frees ALL translations (Rx=Ry=Rz=0); with fric=0 the lateral DOFs are a
    # zero-stiffness null space and the solver fills them with garbage
    # (measured: plate 1e9 mm, pad 19-40 mm).  The small conforming pad slides
    # clean off the dome => zero contact => zero bone stress.  Fix: guide the
    # platen axially (x,y fixed; z free and driven; rotations fixed).
    rigid_fixed = (
        '      <rigid_bc type="rigid_fixed">\n'
        '        <rb>ghost</rb>\n'
        '        <Rx_dof>1</Rx_dof><Ry_dof>1</Ry_dof><Rz_dof>0</Rz_dof>\n'
        '        <Ru_dof>1</Ru_dof><Rv_dof>1</Rv_dof><Rw_dof>1</Rw_dof>\n'
        '      </rigid_bc>\n')
    if force_n is not None:
        # rigid_load translation DOF uses Rx/Ry/Rz (NOT single-letter z — that
        # is the rigid_displacement spelling; FEBio rejects it here with
        # 'tag "dof": invalid value: z'.  Confirmed against thums_feb.py:1033.)
        drive = ('      <rigid_load type="rigid_force">\n'
                 '        <rb>ghost</rb>\n'
                 '        <dof>Rz</dof>\n'
                 '        <relative>0</relative>\n'
                 f'        <value lc="1">{-force_n:.6f}</value>\n'
                 '      </rigid_load>\n')
    else:
        drive = ('      <rigid_bc type="rigid_displacement">\n'
                 '        <rb>ghost</rb>\n'
                 '        <dof>z</dof>\n'
                 '        <relative>0</relative>\n'
                 f'        <value lc="1">{-disp_mm:.10e}</value>\n'
                 '      </rigid_bc>\n')
    rigid_xml = '    <Rigid>\n' + rigid_fixed + drive + '    </Rigid>\n'

    # order per pyfebio StepEntry: Boundary -> Loads -> Contact -> Rigid
    text = re.sub(r"<Loads>.*?</Loads>", "", text, flags=re.S)
    text = gp._inject_before(text, "      </Boundary>", rigid_bc)
    text = gp._inject_after(text, "      </Boundary>", contact_xml)
    text = gp._inject_before(text, "    </step>", rigid_xml)

    # 8) ramp load curve lc=1 (AFTER </MeshDomains>, like p1b.add_load_ramp)
    ramp = gs._RAMP_XML.replace("_RAMP", "pad_ramp")
    text = gp._inject_after(text, "</MeshDomains>", ramp)

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
    """Run FEBio, read peak von Mises + axial reaction (plate recipe)."""
    return gp.solve_and_stats(setup, feb_path, workdir, timeout=timeout,
                              read_damage=read_damage)


# ===========================================================================
# probe drivers
# ===========================================================================
def elastic_scan(setup: p1b.Setup, disp_list: list[float],
                 *, steps: int = 5, timeout: float = 900.0) -> list[dict]:
    """Displacement-controlled elastic scans; report R_axis + peak vm."""
    rows: list[dict] = []
    for d in disp_list:
        stem = f"{STEM}_d{d:g}".replace(".", "p")
        feb = OUT_TEMP / f"{stem}.feb"
        build_pad_feb(setup, feb, disp_mm=d, steps=steps)
        r = solve_and_stats(setup, feb, OUT_TEMP, timeout=timeout)
        rows.append({"disp_mm": d, **r})
        print(f"  d={d:g}mm  rc={r.get('rc')}  ok={r.get('ok')}  "
              f"R_axis={r.get('R_axis')}  vm_max={r.get('vm_max')}", flush=True)
    return rows


def damage_runs(setup: p1b.Setup, disp_mm: float, *, n_runs: int = 5,
                steps: int = 20, timeout: float = 900.0,
                seed0: int = 0) -> list[dict]:
    """Damage model (Dmax=1) at fixed pad displacement, n repeated seeds."""
    rows: list[dict] = []
    for i in range(n_runs):
        seed = seed0 + i
        stem = f"{STEM}_damage_d{disp_mm:g}_seed{seed}".replace(".", "p")
        feb = OUT_TEMP / f"{stem}.feb"
        build_pad_feb(setup, feb, disp_mm=disp_mm, steps=steps,
                      damage_Dmax=1.0, seed=seed)
        r = solve_and_stats(setup, feb, OUT_TEMP, timeout=timeout,
                            read_damage=True)
        rows.append({"run": i + 1, "seed": seed, "disp_mm": disp_mm, **r})
        print(f"  dmg#{i + 1} seed={seed}  rc={r.get('rc')}  ok={r.get('ok')}  "
              f"vm_max={r.get('vm_max')}", flush=True)
    return rows


def main() -> None:
    global PAD_E_MPA, PAD_THICK_MM, PAD_GAP_MM, STEM, PAD_BIND_ALL
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", choices=["elastic", "damage", "probe"], default="probe")
    ap.add_argument("--e-pad", type=float, default=None,
                    help=f"pad Young modulus override (default {PAD_E_MPA})")
    ap.add_argument("--thick", type=float, default=None,
                    help=f"pad thickness override (default {PAD_THICK_MM})")
    ap.add_argument("--gap", type=float, default=None,
                    help=f"pad bottom gap override (default {PAD_GAP_MM})")
    ap.add_argument("--d", type=float, nargs="*", default=None,
                    help="elastic scan displacement list (mm)")
    ap.add_argument("--bind-all", action="store_true",
                    help="bind all pad nodes to the ghost (rigid conforming platen)")
    ap.add_argument("--force", type=float, default=None,
                    help="force mode: axial force on the ghost (N), instead of --d")
    ap.add_argument("--tag", default="", help="suffix for result filenames")
    args = ap.parse_args()

    if args.e_pad is not None:
        PAD_E_MPA = args.e_pad
    if args.thick is not None:
        PAD_THICK_MM = args.thick
    if args.gap is not None:
        PAD_GAP_MM = args.gap
    STEM = f"pad{args.tag}" if args.tag else "pad"
    if args.bind_all:
        PAD_BIND_ALL = True
    sfx = args.tag
    out_elastic = OUT_RESULTS / f"G4_FE_PAD_ELASTIC{sfx}.json"
    out_damage = OUT_RESULTS / f"G4_FE_PAD_DAMAGE{sfx}.json"

    t0 = time.time()
    setup = p1b.Setup()
    print(f"[g4_fe_pad] Setup loaded in {time.time() - t0:.1f}s  "
          f"E_pad={PAD_E_MPA} thick={PAD_THICK_MM} gap={PAD_GAP_MM}")

    if args.force is not None:
        print(f"[force mode] rigid_force Rz = {args.force:.1f} N on the ghost")
        feb = OUT_TEMP / f"{STEM}_force{args.force:g}".replace(".", "p")
        build_pad_feb(setup, feb, force_n=args.force, steps=5)
        r = solve_and_stats(setup, feb, OUT_TEMP, timeout=900)
        print(f"  F={args.force:.1f} N  rc={r.get('rc')} ok={r.get('ok')} "
              f"R_axis={r.get('R_axis')} vm_max={r.get('vm_max')}")
        json.dump({"force_n": args.force, "row": r, "bind_all": PAD_BIND_ALL,
                   "E_pad_mpa": PAD_E_MPA, "thick_mm": PAD_THICK_MM,
                   "gap_mm": PAD_GAP_MM},
                  open(out_elastic, "w"), indent=2, ensure_ascii=False)
    elif args.only in ("elastic", "probe"):
        disp_list = args.d if args.d else [0.2, 0.5, 1.0, 1.5, 2.0, 3.0]
        print(f"[elastic scan] displacement-controlled rigid hat on conforming pad")
        rows = elastic_scan(setup, disp_list)
        conv = [r for r in rows if r.get("ok") and r.get("R_axis") is not None]
        r_max = max((abs(r["R_axis"]) for r in conv), default=0.0)
        # only report d1000 when the curve actually REACHES 1000 N; otherwise
        # _interp_disp would silently clamp to the last displacement
        d1000 = (gp._interp_disp(conv, 1000.0)
                 if conv and r_max >= 1000.0 else None)
        print(f"[elastic] peak |R_axis| = {r_max:.3f} N over scan "
              f"(1000 N reached: {r_max >= 1000.0})")
        print(f"[elastic] 1000 N @ disp ≈ {d1000:.3f} mm" if d1000
              else "[elastic] 载荷未达 1000 N ⇒ 无 d1000（方案判据未达成）")
        json.dump({"rows": rows, "d1000_mm": d1000, "r_max_n": r_max,
                   "bind_all": PAD_BIND_ALL, "gap_mm": PAD_GAP_MM,
                   "E_pad_mpa": PAD_E_MPA, "thick_mm": PAD_THICK_MM},
                  open(out_elastic, "w"), indent=2, ensure_ascii=False)

    if args.only in ("damage", "probe"):
        dstar = None
        if out_elastic.exists():
            dstar = json.load(open(out_elastic, encoding="utf-8")).get("d1000_mm")
        if dstar is None:
            print("[damage] no reliable d1000; short scan fallback")
            rows = elastic_scan(setup, [0.5, 1.0, 1.5])
            conv = [r for r in rows if r.get("ok") and r.get("R_axis")]
            dstar = gp._interp_disp(conv, 1000.0) if conv else None
        if dstar is None:
            raise SystemExit("elastic scan failed entirely; abort damage")
        print(f"[damage] Dmax=1 @ disp={dstar:.3f} mm (target 1000 N) x5")
        dmg = damage_runs(setup, dstar, n_runs=5)
        json.dump({"dstar_mm": dstar, "runs": dmg}, open(out_damage, "w"),
                  indent=2, ensure_ascii=False)

    print(f"[g4_fe_pad] done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()

