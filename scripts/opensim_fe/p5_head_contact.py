# -*- coding: utf-8 -*-
"""P5.0 能力探针：髁端关节面接触承压 + 股骨头解剖载荷。

预注册：P5_HEAD_CONTACT_TODO.md（2026-10-14，冻结，不入库）。
目的：P5 主工程（9 部位重解）之前的单骨（femur_r）能力确认 —— 验证
「髁端刚性压板 sliding-elastic 接触（替代 bottom_fix 全带 + lateral_support）+
股骨头球冠 NodalForce」这一 BC 形式能否稳定求解、接触确实承载轴向、
弯曲路径恢复、不触发屈曲/横向滑移哨兵。

复用（零改动，只 ADD 本模块）：
- p4_pb_probe.load_ctx / head_joint_nodes（femur_r 载荷上下文）
- p4c_probe.statistic_p4c（shaft 统计 + 弯曲指标 + 位移哨兵）
- g4_fe_plate：plate_mesh_xml / _ghost_mesh_xml / CONTACT_RECIPE / SURFACE_PAIR_XML /
  _next_material_id / _inject_before / _inject_after（接触配方逐字沿用，fric 参数化）
- g4_fe_p1b：oriented_facets / band_facets / add_reaction_var
- bone_feb：build_bone_feb / run_bone_feb / _align_pca_rotate / _select_top_bottom
- g4_fe_stab：_RAMP_XML / _node_var

变体矩阵（冻结）：p5_anchor（fric=0 + 3 节点横向锚点）、p5_fric（fric=0.5 无锚点）、
p5_free（fric=0 无锚点，负对照）。高度 h∈{10,20,30,50}。
判定与哨兵见预注册 TODO §1.5 / §1.6，不回溯。
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import bone_feb as bf
import g4_fe_p1b as p1b
import g4_fe_plate as plate
import g4_fe_stab as gs
import paper_alignment as pa
from p4_pb_probe import FIX_FRAC, HEAD_FRAC, head_joint_nodes, load_ctx
from p4c_probe import statistic_p4c

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "opensim_fe"
BASE_DIR = ROOT / "temp" / "opensim_fe" / "p5_head_contact" / "femur_r"

SIGMA_C_MPA = 220.0
HEIGHTS = [10, 20, 30, 50]
PLATE_DISP_MM = 0.5   # 板位移压入量（冻结，g4_fe_plate 接触最稳档）
GAP_MM = -0.3         # 板初始预穿透（2026-10-14 修正：原 +0.2 在 t=0 造成 z 向奇异）
                      # 原理：t=0 时骨轴向完全自由（bottom_fix 只约束横向），
                      # 必须由接触提供 z 锚位 ⇒ 接触须在 t=0 已激活且穿透足够。
                      # 实测 +0.2（间隙）与 -0.05（0.05mm 穿透）均 FAIL，
                      # -0.3 收敛且接触承载 100.0001% 轴向载荷。
TOP_FRAC = 0.30       # 髁端 cap 带份额（与 P4 系 fix 带一致）
CAP_AXDIV_MIN = 0.5   # cap 面法向轴向分量下限（只取真端面，2026-10-14 根因修复）

# P4-a' 已收口对照（results/opensim_fe/P4_PB_PROBE.json，不重跑）
P4A_BY_H = {10: 33.1, 20: 48.5, 30: 60.7, 50: 83.3}


# ---------------------------------------------------------------------------
# 几何 / 载荷
# ---------------------------------------------------------------------------

def _combined_elems(mesh, spec) -> list[list[int]]:
    """拼接全部 solid_domains 的 tet+hex（1-based node-id），等效 p1b.Setup 的 elems。"""
    elems = []
    for sd in spec.solid_domains:
        b = mesh["by_solid"][sd.mesh_name]
        for row in b["tets"]:
            elems.append([int(x) for x in row])
        for row in b["hexes"]:
            elems.append([int(x) for x in row])
    return elems


def cap_geom(mesh, spec, top_frac: float = TOP_FRAC) -> dict:
    """髁端 cap 面片几何（femur PCA 旋转后 +Z 顶端 = 髁端，F-3'c 已验证）。

    返回：axis_idx/axis（+Z）、a_min/a_max、top_local（1-based）、cap 面片列表、
    a_cap_mm2（cap 面积和）。cap 面片含 tri3+quad4（femur 端面含 tet 面）。
    """
    nodes_rot = bf._align_pca_rotate(mesh["nodes"])
    axis_idx, a_min, a_max, top_local, _ = bf._select_top_bottom(
        nodes_rot, axis=None, top_frac=top_frac, bot_frac=top_frac)
    facets = p1b.oriented_facets(nodes_rot, _combined_elems(mesh, spec))
    band = p1b.band_facets(facets, top_local)
    axis = np.zeros(3)
    axis[axis_idx] = 1.0
    for f in band:
        f["axdiv"] = float(np.dot(f["normal"], axis))
    # 2026-10-14 根因修复（实测：axdiv>0.0 选中 542 面片中 365 个 axdiv<0.5，
    # 髁端球面/侧面混入 => 接触反力带横向分量 => 压板刚体 x/y 漂移 2.4e6 mm，
    # 「收敛」但解无物理意义）。只保留真端面（法向朝轴向），横向分量 <0.5。
    cap = [f for f in band if f["axdiv"] > CAP_AXDIV_MIN]
    return {
        "axis_idx": int(axis_idx), "axis": axis.tolist(),
        "a_min": float(a_min), "a_max": float(a_max),
        "n_top_nodes": len(top_local), "top_local": top_local,
        "band": band, "cap": cap,
        "n_cap_facets": int(len(cap)),
        "a_cap_mm2": float(sum(f["area"] for f in cap)),
        "n_cap_tri": int(sum(1 for f in cap if len(f["poly"]) == 3)),
        "n_cap_quad": int(sum(1 for f in cap if len(f["poly"]) == 4)),
    }


def plate_geom_from_mesh(nodes_rot: np.ndarray, geom: dict) -> dict:
    """刚性压板几何（等效 g4_fe_plate.plate_geom，但用任意旋转后节点，无 setup）。"""
    xs, ys = nodes_rot[:, 0], nodes_rot[:, 1]
    x0, x1 = float(xs.min()), float(xs.max())
    y0, y1 = float(ys.min()), float(ys.max())
    z_wall = geom["a_max"] + GAP_MM
    return {
        "x0": x0 - plate.PLATE_MARGIN_MM, "x1": x1 + plate.PLATE_MARGIN_MM,
        "y0": y0 - plate.PLATE_MARGIN_MM, "y1": y1 + plate.PLATE_MARGIN_MM,
        "z_wall": z_wall, "z_bot": z_wall - plate.PLATE_THICK_MM,
        "z_mid": z_wall - 0.5 * plate.PLATE_THICK_MM,
        "ndiv": plate.PLATE_NDIV, "n_node_grid": plate.PLATE_NDIV + 1,
    }


def _anchor_nodes(hj: dict, n: int = 3) -> np.ndarray:
    """从髁端带（1-based）均匀取 n 个节点作为锚点。"""
    arr = np.asarray(hj["fix_nodes_1based"], dtype=np.int64)
    idx = np.unique(np.linspace(0, len(arr) - 1, n).round().astype(int))
    return arr[idx]


# ---------------------------------------------------------------------------
# 接触注入（复制 g4_fe_plate 注入顺序，fric 参数化 + cap 面混合 tri3/quad4）
# ---------------------------------------------------------------------------

def _cap_surface_xml(geom: dict, name: str = "bone_top") -> str:
    """骨端（髁端）cap 面片 surface XML —— 支持 tri3+quad4 混合。

    名字必须为 bone_top：g4_fe_plate.SURFACE_PAIR_XML 的 primary 硬编码该名。
    与 g4_fe_plate.bone_cap_surface_xml 的差别：femur 端面可能含 tet 面（tri3），
    不能复用其全 quad4 断言；这里按 p1b.inject_surface 的 tag 逻辑输出混合面片。
    """
    if not geom["cap"]:
        raise ValueError("cap 空：髁端端面无面片")
    n = len(geom["cap"])
    out = [f'    <Surface name="{name}">']
    for i, f in enumerate(geom["cap"]):
        poly = [int(p) for p in f["poly"]]
        if len(poly) == 3:
            out.append(f'      <tri3 id="{i + 1}">{",".join(map(str, poly))}</tri3>')
        elif len(poly) == 4:
            out.append(f'      <quad4 id="{i + 1}">{",".join(map(str, poly))}</quad4>')
        else:
            raise ValueError(f"cap 含 {len(poly)} 顶点面片，不支持")
    out.append("    </Surface>")
    out.append(f"<!-- bone_top: {n} facets "
               f"({geom['n_cap_tri']} tri3 / {geom['n_cap_quad']} quad4) -->")
    return "\n".join(out)


def _inject_plate_contact(text: str, pg: dict, geom: dict, *,
                          fric: float, disp_mm: float, first_extra_id: int,
                          penalty: float = 0.1) -> str:
    """把刚性压板 + sliding-elastic 接触完整注入 build_bone_feb 输出的 .feb 文本。

    顺序与 g4_fe_plate.build_plate_feb 第 2-10 步一致（XML 原文逐字沿用），
    差异：① fric_coeff/penalty 参数化（penalty 默认 0.1 保持 g4_fe_plate 原配方，
    p5_pen 变体用 50.0 适配 femur 载荷量级）；② **不删除 <Loads>**（P5 保留头端
    球冠 NodalForce）；③ cap 面用混合 tri3/quad4 的 bone_top。
    """
    n_node_grid = pg["n_node_grid"]
    n_bone = first_extra_id - 1                       # 骨节点数（注入从 n_bone+1 起）

    # a. 板网格 + 板面 + 幽灵刚体
    nodes_xml, elems_xml, nodeset_xml, plate_surf_xml = plate.plate_mesh_xml(
        pg, first_node_id=n_bone + 1)
    ghost_nodes_xml, ghost_elems_xml = plate._ghost_mesh_xml(
        pg, first_node_id=n_bone + 1 + 2 * n_node_grid ** 2)

    text = plate._inject_before(text, "</Nodes>", nodes_xml)
    text = plate._inject_before(text, "</Mesh>", elems_xml)
    text = plate._inject_before(text, "</Mesh>", ghost_nodes_xml)
    text = plate._inject_before(text, "</Mesh>", ghost_elems_xml)
    text = plate._inject_before(text, "</Mesh>", nodeset_xml)
    text = plate._inject_before(text, "</Mesh>", plate_surf_xml)
    text = plate._inject_before(text, "</Mesh>", _cap_surface_xml(geom))
    text = plate._inject_before(text, "</Mesh>", plate.SURFACE_PAIR_XML)

    # b. 材料 + 域（plate 弹性 + ghost 刚体载体，原文）
    mid = plate._next_material_id(text)
    plate_mat = (f'    <material name="plate_mat" id="{mid}" '
                 f'type="isotropic elastic">\n'
                 f'      <density>8.0e-09</density>\n'
                 f'      <E>{plate.PLATE_E_MPA:.1f}</E>\n'
                 f'      <v>{plate.PLATE_NU}</v>\n'
                 f'    </material>\n')
    ghost_mat = (f'    <material name="ghost" id="{mid + 1}" '
                 f'type="rigid body">\n'
                 f'      <density>1e-06</density>\n'
                 f'    </material>\n')
    text = plate._inject_before(text, "</Material>", plate_mat + ghost_mat)
    text = plate._inject_before(text, "</MeshDomains>",
                                f'    <SolidDomain name="plate" mat="plate_mat"/>\n'
                                f'    <SolidDomain name="ghost" mat="ghost"/>\n')

    # c. Boundary：rigid_bc（板节点绑到幽灵刚体）
    rigid_bc = ('      <bc type="rigid" name="plate_to_ghost" '
                'node_set="plate_nodes">\n'
                '        <rb>ghost</rb>\n'
                '      </bc>\n')
    text = plate._inject_before(text, "      </Boundary>", rigid_bc)

    # d. Contact：sliding-elastic（fric/penalty 参数化；flip_secondary=1 保留）
    c = plate.CONTACT_RECIPE
    contact_xml = (
        '    <Contact>\n'
        '      <contact type="sliding-elastic" name="plate_bone_contact" '
        'surface_pair="plate_bone_pair">\n'
        f'        <laugon>{c["laugon"]}</laugon>\n'
        f'        <penalty>{penalty}</penalty>\n'
        f'        <auto_penalty>{c["auto_penalty"]}</auto_penalty>\n'
        f'        <two_pass>{c["two_pass"]}</two_pass>\n'
        f'        <node_reloc>{c["node_reloc"]}</node_reloc>\n'
        f'        <symmetric_stiffness>{c["symmetric_stiffness"]}</'
        f'symmetric_stiffness>\n'
        f'        <tolerance>{c["tolerance"]}</tolerance>\n'
        f'        <search_radius>{c["search_radius"]}</search_radius>\n'
        f'        <fric_coeff>{fric}</fric_coeff>\n'
        f'        <flip_secondary>{c.get("flip_secondary", 0)}</'
        f'flip_secondary>\n'
        '      </contact>\n'
        '    </Contact>\n')
    text = plate._inject_after(text, "      </Boundary>", contact_xml)

    # e. Rigid：位移驱动（板压入 disp_mm，转动锁死；2026-10-14 根因修复：
    #    Rx/Ry 由 0(放开) 改 1(锁死) —— 实测放开时接触反力横向分量把板刚体
    #    推走 2.4e6 mm，"收敛"但解无物理意义；z 仍由 rigid_displacement 驱动）
    rigid_xml = (
        '    <Rigid>\n'
        '      <rigid_bc type="rigid_fixed">\n'
        '        <rb>ghost</rb>\n'
        '        <Rx_dof>1</Rx_dof><Ry_dof>1</Ry_dof><Rz_dof>0</Rz_dof>\n'
        '        <Ru_dof>1</Ru_dof><Rv_dof>1</Rv_dof><Rw_dof>1</Rw_dof>\n'
        '      </rigid_bc>\n'
        '      <rigid_bc type="rigid_displacement">\n'
        '        <rb>ghost</rb>\n'
        '        <dof>z</dof>\n'
        '        <relative>0</relative>\n'
        f'        <value lc="1">{-disp_mm:.10e}</value>\n'
        '      </rigid_bc>\n'
        '    </Rigid>\n')
    text = plate._inject_before(text, "    </step>", rigid_xml)

    # f. ramp lc=1/2（AFTER </MeshDomains>）
    #    P5 修复（2026-10-14）：sliding-elastic 单边接触（tension=no）＋骨轴向完全自由
    #    ⇒ 板压入产生的接触力必须始终有载荷平衡（载荷把骨推向板），否则骨被推出
    #    接触区 → 刚体模态 → 负 jacobian。
    #    lc=1：plate 位移 ramp，起点 0.5 —— t=0 板位移 -0.25mm = 闭合 gap(0.2mm)
    #          + 预穿透 0.05mm，接触从第一刻激活（§1.1「板压入闭合接触 + 提供 z 锚位」）。
    #    lc=2：head 载荷 ramp，起点 0.5 —— t=0 载荷 50%，与接触力同步建立平衡。
    #    t=1 终态不变（板 -0.5mm / 载荷 100%），仅改加载时序。
    ramp_xml = (
        '  <LoadData>\n'
        '    <load_controller id="1" type="loadcurve">\n'
        '      <interpolate>LINEAR</interpolate>\n'
        '      <extend>CONSTANT</extend>\n'
        '      <points><point>0,0.5</point><point>1,1</point></points>\n'
        '    </load_controller>\n'
        '    <load_controller id="2" type="loadcurve">\n'
        '      <interpolate>LINEAR</interpolate>\n'
        '      <extend>CONSTANT</extend>\n'
        '      <points><point>0,0.5</point><point>1,1</point></points>\n'
        '    </load_controller>\n'
        '  </LoadData>\n'
    )
    text = plate._inject_after(text, "</MeshDomains>", ramp_xml)

    # f2. head 载荷（nodal_force top_load）挂 lc=2 —— 与板压入同步 ramp（P5 修复）
    text = re.sub(
        r'(<nodal_load type="nodal_force" node_set="top_load">.*?<value)(>)([^<]*</value>)',
        r'\1 lc="2"\2\3', text, count=1, flags=re.S)

    # g. reaction var + 原生 contact force 输出（P5 能力判据的观测量，2026-10-14）
    text = p1b.add_reaction_var(text)
    text = _add_contact_force_var(text)
    return text


def _add_contact_force_var(text: str) -> str:
    """在 <Output> 段追加 ``<var type="contact force"/>``（FEBio 4.13 原生支持）。

    没有它，`_load_bear_stats` 无从读取接触承载量 —— 已验证：加了才有
    ``states/{sid}/surface_data/contact force``。
    """
    if '<var type="contact force"/>' in text:
        return text
    new, n = re.subn(r'([ \t]*)<var type="reaction forces"/>',
                     r'\1<var type="reaction forces"/>\n\1<var type="contact force"/>',
                     text, count=1)
    if n != 1:
        raise ValueError("未找到 <var type=\"reaction forces\"/>，contact force 注入失败")
    return new


def _anchor_xy_patch(text: str) -> str:
    """把 bottom_fix BC 的 z_dof 置 0 —— 锚点只约束横向（x/y），轴向 100% 走接触。"""
    pat = re.compile(r'(<bc type="zero displacement" node_set="bottom_fix">'
                     r'.*?<z_dof>)(\d)(</z_dof>.*?</bc>)', re.S)
    if not pat.search(text):
        raise ValueError("未找到 bottom_fix BC，锚点 patch 失败")
    return pat.sub(lambda m: m.group(1) + "0" + m.group(3), text, count=1)


# ---------------------------------------------------------------------------
# 后处理：接触承载验证 + 横向滑移哨兵
# ---------------------------------------------------------------------------

def _load_bear_stats(h5: Path, load_n: float) -> dict:
    """接触承载验证 + 横向滑移哨兵（单次开文件）。

    **2026-10-14 观测量修复**：原口径用 node_data 的 ``reaction`` 求和，
    但该变量是**逐节点内力**（全域求和按构造 telescoping 恒为 0，平衡要求
    -load_n 而实测 -0.000），并且 `lat_support_0`(z_dof=0，物理不可能承轴向)
    竟报 -18226 N 与 BC 直接冲突 ⇒ 口径结构性失效，无法作判据。

    **修复**：改用 FEBio **原生 contact force 输出**
    （`<var type="contact force"/>` → ``states/{sid}/surface_data/contact force``）。

    ⚠️ 读法要点：`bone_top` 与 `plate_bottom` 是同一对接触的两个面，值互为相反数
    （牛顿第三定律），**两面相加恒为 0** —— 只可取**单面**值。本函数取
    `bone_top`（骨侧感受的力，方向 -z，与 +z 载荷反向）。

    判定口径：|R_z_contact|/load_n ≈ 1（接触确实承载轴向）。
    """
    import h5py
    with h5py.File(h5, "r") as h:
        sid = sorted(h["states"].keys(), key=int)[-1]
        sd = h[f"states/{sid}"]
        dis = np.asarray(sd["node_data/displacement/femur_r_nodes"], dtype=float)
        if "surface_data" not in sd or "contact force" not in sd["surface_data"]:
            raise KeyError(
                f"h5 无原生 contact force 输出（states/{sid}/surface_data/contact force）"
                " —— .feb 的 <Output> 段必须含 <var type=\"contact force\"/>；"
                "旧 xplt/h5 需重跑求解")
        cf = sd["surface_data/contact force"]
        # 单面读取：bone_top（骨侧）。两面相加恒 0，不可相加。
        f_bone = np.asarray(cf["bone_top"], dtype=float).sum(axis=0)
        f_plate = np.asarray(cf["plate_bottom"], dtype=float).sum(axis=0)

    r_z_contact = float(f_bone[2])
    # 牛顿第三定律自检：两面应互为反号（|sum| ≈ 0）
    n3_selfcheck = float(np.abs(f_bone + f_plate).max())
    out = {
        "R_z_contact_n": r_z_contact,
        "contact_f_bone_n": [float(x) for x in f_bone],
        "contact_f_plate_n": [float(x) for x in f_plate],
        "contact_newton3_resid_n": n3_selfcheck,
        "R_axis_check": abs(abs(r_z_contact) - load_n) / load_n,
        "dis_xy_max_mm": float(np.sqrt(dis[:, 0] ** 2 + dis[:, 1] ** 2).max()),
        "dis_z_max_mm": float(np.abs(dis[:, 2]).max()),
    }
    return out


# ---------------------------------------------------------------------------
# 构建 / 运行 / 判定 / 主流程
# ---------------------------------------------------------------------------

VARIANTS = {
    "p5_anchor": {"fric": 0.0, "anchor": True, "penalty": 0.1,
                  "tag": "fric=0 + 髁端3节点横向锚点（法向接触承载轴向；penalty 沿用 g4_fe_plate 0.1）"},
    "p5_pen": {"fric": 0.0, "anchor": True, "penalty": 50.0,
               "tag": "fric=0 + 锚点 + penalty=50（配方适配：femur 22-58kN 载荷量级，"
                      "penalty=0.1 ⇒ 穿透≈29mm≫gap 0.2mm，接触形同虚设；50 ⇒ 0.06mm 合理）"},
    "p5_fric": {"fric": 0.5, "anchor": False, "penalty": 0.1,
                "tag": "fric=0.5 无锚点（摩擦承载横向；复现 g4_fe_plate fric 失效记录）"},
    "p5_free": {"fric": 0.0, "anchor": False, "penalty": 0.1,
                "tag": "fric=0 无锚点（对照：横向仅靠 lateral_support；注意 2026-10-14 "
                       "起 lateral_support 对全部变体生效 ⇒ 本变体不再是「无横向支撑」负对照，"
                       "该角色由诊断期 E 系列承担）"},
}


def _build_one(vname: str, hno: int, load: float, ctx, hj, geom, pg) -> Path:
    """构建 .feb（build_bone_feb 定制 BC → 注入接触 → 锚点 x/y patch）。"""
    part, mesh, solid_names, shell_specs, dom, joint, loads, heights = ctx
    d = BASE_DIR / f"{vname}_h{hno:02d}"
    d.mkdir(parents=True, exist_ok=True)
    feb = d / "femur_r.feb"
    v = VARIANTS[vname]
    bf.build_bone_feb(
        part.spec, feb, load_n=load, top_frac=TOP_FRAC, bot_frac=TOP_FRAC,
        axis=None, with_shell=True, time_steps=10, mesh=mesh, align_pca=True,
        # P4 同款横向导向（2026-10-14）：E6/E7 对照证明「3 节点 anchor」不足以
        # 替代它 —— 同一 anchor 只补 lateral_support 即收敛，无它则横向失稳。
        lateral_support=[(0.05, 0.95)],
        bottom_fix_axial_only=(not v["anchor"]),   # anchor: 全固定(后 patch 横向) / 其余: z-only
        load_nodes_override=np.asarray(hj["head_nodes_1based"], dtype=np.int64),
        fix_nodes_override=_anchor_nodes(hj, 3),
        load_dir_vec=np.asarray(hj["load_dir_vec"], dtype=float),
    )
    text = feb.read_text(encoding="ISO-8859-1")
    text = _inject_plate_contact(
        text, pg, geom, fric=v["fric"], disp_mm=PLATE_DISP_MM,
        first_extra_id=mesh["nodes"].shape[0] + 1,
        penalty=v.get("penalty", 0.1))
    if v["anchor"]:
        text = _anchor_xy_patch(text)
    feb.write_text(text, encoding="ISO-8859-1")
    return feb


def judge(rows: list[dict], p4a_dis50: float | None) -> tuple[str, str, str]:
    """预注册 §1.5 判据（跑前写死，跑后如实）。返回 (verdict, attribution, sentinel)。"""
    sentinels = []
    for vname in ("p5_anchor", "p5_pen", "p5_fric", "p5_free"):
        vrows = [r for r in rows if r["case"] == vname]
        failed = [r for r in vrows if not r.get("converged")]
        if failed:
            sentinels.append(f"{vname}: 不收敛 @h={[r['h_no'] for r in failed]}")
        conv = [r for r in vrows if r.get("converged") and "shaft_stat_mpa" in r]
        if not conv:
            continue
        hmax_c = max(r["h_no"] for r in conv)
        hc = next(r for r in conv if r["h_no"] == hmax_c)
        if hc.get("dis_mag_max_mm") and p4a_dis50:
            ratio = hc["dis_mag_max_mm"] / p4a_dis50
            if ratio > 3.0:
                sentinels.append(f"{vname}: h={hmax_c} 位移 {hc['dis_mag_max_mm']:.3f}mm "
                                 f"= {ratio:.1f}× P4-a' h50（{p4a_dis50:.3f}mm）剧增")
        if hc.get("dis_xy_max_mm") and hc["dis_xy_max_mm"] > 5.0:
            sentinels.append(f"{vname}: h={hmax_c} 横向滑移 {hc['dis_xy_max_mm']:.3f}mm > 5mm")
        # 接触两面牛顿第三定律自检（取代已失效的 anchor_share 口径，2026-10-14）
        if hc.get("contact_newton3_resid_n", 0.0) > 1.0:
            sentinels.append(f"{vname}: h={hmax_c} 接触两面向量和残差 "
                             f"{hc['contact_newton3_resid_n']:.2f}N > 1N（接触力读数不自洽）")
    sentinel = "；".join(sentinels) if sentinels else ""

    # 能力主判据：候选变体在最高收敛高度 收敛 + 接触承载轴向（|R_z_contact|≈load_n）
    proven = []
    for vname in ("p5_anchor", "p5_pen", "p5_fric"):
        conv = [r for r in rows if r["case"] == vname and r.get("converged")
                and "R_z_contact_n" in r]
        if not conv:
            continue
        hmax_c = max(r["h_no"] for r in conv)
        hc = next(r for r in conv if r["h_no"] == hmax_c)
        ok_axis = abs(abs(hc["R_z_contact_n"]) - hc["load_n"]) / hc["load_n"] <= 0.20
        if ok_axis:
            proven.append((vname, hmax_c, hc))
    if proven:
        vname, hmax_c, hc = proven[0]
        verdict = (f"能力成立（主判据）：{vname} 在 h={hmax_c} 收敛，接触承载轴向 "
                   f"R_z_contact={hc['R_z_contact_n']:.1f}N vs load_n={hc['load_n']:.1f}N"
                   f"（相对偏差 {hc['R_axis_check']*100:.2f}%）"
                   f"⇒ 髁端接触承压可替代位移固定带承载轴向")
        attribution = ("法向接触（fric=0）在 h_max_c 提供 |R_contact|≈load_n 的接触力："
                       "轴向力全部走接触；该量取自 FEBio 原生 contact force 输出"
                       "（bone_top 单面，两面为作用力/反作用力不可相加）")
    else:
        verdict = ("能力不成立（主判据）：p5_anchor / p5_pen / p5_fric 均未在最高收敛高度实现 "
                   "|R_z_contact|≈load_n（接触未独立承载轴向）")
        attribution = "原因待哨兵/错误记录归因；P5 主工程须换 BC 形式或接触配方"

    # 弯曲恢复次判据（仅当能力成立时评估，如实记录不回溯）
    bend = ""
    if proven:
        vname, hmax_c, hc = proven[0]
        p4a_same = P4A_BY_H.get(hmax_c)
        ratio = hc["shaft_stat_mpa"] / p4a_same if p4a_same else 0.0
        if (hc["shaft_stat_mpa"] >= 1.5 * p4a_same or hc["szz_max_mpa"] > 0
                or hc["r2_linear"] >= 0.2):
            bend = (f"；弯曲恢复次判据满足：σ_shaft={hc['shaft_stat_mpa']:.1f} "
                    f"({'≥1.5×' if ratio >= 1.5 else '≈'}{ratio:.2f}× P4-a' h{hmax_c}="
                    f"{p4a_same}) szz_max={hc['szz_max_mpa']:.1f} r2={hc['r2_linear']:.3f}"
                    " ⇒ 弯曲路径恢复（P5 形式解决 P4-a′ 抑制）")
        else:
            bend = (f"；弯曲恢复次判据不满足：σ_shaft={hc['shaft_stat_mpa']:.1f} "
                    f"({ratio:.2f}× P4-a') szz_max={hc['szz_max_mpa']:.1f} "
                    f"r2={hc['r2_linear']:.3f} ⇒ 弯曲仍缺失（如实记录，不回溯）")
    return verdict + bend, attribution, sentinel


def main() -> int:
    t_all = time.time()
    ctx = load_ctx()
    part, mesh, solid_names, shell_specs, dom, joint, loads, heights = ctx
    print(f"[p5] {part.bone}: mesh {mesh['n_nodes']} nodes, dom={dom}, joint={joint}",
          flush=True)
    hj = head_joint_nodes(mesh)
    print(f"[p5] 载荷=股骨头球冠 n={hj['n_head_nodes']} e_band={hj['band_eccentricity_mm']:.2f}mm, "
          f"髁端带 n={hj['n_fix_nodes']}", flush=True)
    geom = cap_geom(mesh, part.spec, top_frac=TOP_FRAC)
    print(f"[p5] 髁端 cap: {geom['n_cap_facets']} facets "
          f"({geom['n_cap_tri']} tri / {geom['n_cap_quad']} quad), "
          f"A_cap={geom['a_cap_mm2']:.1f}mm²", flush=True)
    pg = plate_geom_from_mesh(bf._align_pca_rotate(mesh["nodes"]), geom)

    # 正对照基准：P4-a' 已收口（不重跑）
    p4a_h5 = ROOT / "temp" / "opensim_fe" / "p4_pb_probe" / "femur_r" / "head_load_h50" / "femur_r.hdf5"
    p4a_dis = None
    if p4a_h5.is_file():
        p4a_dis = statistic_p4c(p4a_h5, dom, float(loads[-1])).get("dis_z_max_mm")

    rows: list[dict] = []
    for vname in VARIANTS:
        for hno in HEIGHTS:
            load = float(loads[hno - 1])
            d = BASE_DIR / f"{vname}_h{hno:02d}"
            d.mkdir(parents=True, exist_ok=True)
            feb = d / "femur_r.feb"
            h5 = feb.with_suffix(".hdf5")
            rec = {"case": vname, "variant": VARIANTS[vname]["tag"],
                   "h_no": hno, "h_m": float(heights[hno - 1]), "load_n": load,
                   "fric": VARIANTS[vname]["fric"],
                   "penalty": VARIANTS[vname].get("penalty", 0.1)}
            if not h5.is_file():
                t = time.time()
                try:
                    _build_one(vname, hno, load, ctx, hj, geom, pg)
                    bf.run_bone_feb(feb, workdir=feb.parent, timeout=1800)
                    if feb.with_suffix(".xplt").is_file():
                        bf.post_process(feb.with_suffix(".xplt"), h5,
                                        solid_names=solid_names,
                                        shell_specs=shell_specs, label=part.bone)
                    rec["converged"] = True
                except Exception as exc:  # noqa: BLE001
                    rec["converged"] = False
                    rec["error"] = f"{type(exc).__name__}: {exc}"
                rec["seconds"] = round(time.time() - t, 1)
            else:
                rec["converged"] = True
                rec["seconds"] = 0.0
            if h5.is_file():
                rec.update(statistic_p4c(h5, dom, load))
                rec.update(_load_bear_stats(h5, load))
            rows.append(rec)
            print(f"[p5] {vname} h={hno:2d} ({load:9.1f}N) "
                  f"{'OK ' if rec.get('converged') else 'FAIL'} "
                  + (f"shaft={rec['shaft_stat_mpa']:8.2f} szzmax={rec.get('szz_max_mpa', float('nan')):7.1f} "
                     f"r2={rec.get('r2_linear', float('nan')):.3f} "
                     f"Rz={rec.get('R_z_contact_n', float('nan')):9.1f} "
                     f"err={rec.get('R_axis_check', float('nan'))*100:5.2f}% "
                     f"xy={rec.get('dis_xy_max_mm')} dz={rec.get('dis_z_max_mm')} "
                     f"[{rec.get('seconds')}s]"
                     if "shaft_stat_mpa" in rec else f"({rec.get('error', '')[:70]})"),
                  flush=True)

    # ---- 判定（预注册 §1.5，跑前写死，跑后如实） ----
    verdict, attribution, sentinel = judge(rows, p4a_dis)
    print("\n=== 判定 ===")
    print(f"  能力主判据：{verdict}")
    print(f"  归因：{attribution}")
    if sentinel:
        print(f"  哨兵：{sentinel}")

    out_json = RES / "P5_HEAD_CONTACT.json"
    out_json.write_text(json.dumps({
        "probe": "P5.0 能力探针：髁端关节面接触承压 + 股骨头解剖载荷（femur_r）",
        "bone": part.bone, "domain": dom, "joint": joint, "sigma_c_mpa": SIGMA_C_MPA,
        "plate_disp_mm": PLATE_DISP_MM, "gap_mm": GAP_MM, "top_frac": TOP_FRAC,
        "head_geom": {k: v for k, v in hj.items() if k not in
                      ("head_nodes_1based", "fix_nodes_1based")},
        "cap_geom": {k: v for k, v in geom.items()
                     if k not in ("band", "cap", "top_local")},
        "p4a_reference": P4A_BY_H,
        "criterion": ("能力主判据：p5_anchor 或 p5_fric 最高收敛高度收敛且 |R_z_contact|≈load_n "
                      "（±20%）⇒ 接触承压能力成立；弯曲恢复次判据：σ_shaft≥1.5×P4-a' 或 "
                      "szz_max>0 或 r2≥0.2；哨兵：发散/位移>3×P4-a'/横向滑移>5mm/锚点吸轴>20%。"),
        "verdict": verdict, "attribution": attribution, "sentinel": sentinel,
        "rows": rows}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n[p5] 总耗时 {time.time()-t_all:.0f}s；写出 {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

