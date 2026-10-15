# -*- coding: utf-8 -*-
"""``joint_bc="contact"`` —— 距下关节的**刚性距骨 + 真实接触**载荷路径。

**为什么需要它**（`results/opensim_fe/V1_LOAD_INTRO_SINGULARITY.md` §3）

V1 现行口径把 `PressureLoad p = F/A` 施加在距下关节面的一个**子集**上。
该子集的边界是「牵引从 `p` 突降到 0」的不连续面 ⇒ 线性弹性下产生**载荷引入奇异**
（实测 `SC = σ_max/p_avg` 达 **3.96 / 20.39**，而骨的体力 SC 应 2–4）。
零调参 A/B 已证明峰值**不由名义压力决定**：补片面积 222→2495 mm² 时
`p_avg` 降 **11.24×** 而 `σ_max` 只降 **2.18×** ⇒ 峰值由**边界几何**决定。

**本模块做什么**

在 `joint_face`（距下关节面，法向朝上）**上方**造一块**薄 hex8 板**（"距骨帽"，
`E = 200 GPa` 数值刚性），其**底面**与该面构成 FEBio `sliding-elastic` 接触对；
把总力改为 `PressureLoad` 施加在**帽的顶面**（`p = F/A_cap_top`），
帽只允许沿 **y** 移动（`x = z = 0` 约束）。

**为什么这能消除奇异**：接触牵引场由 **KKT 条件自洽解出**，在接触区边界
`t_n → 0` **连续过渡**，不存在人为的牵引突降。这与 [Y25] 原文 §1.1 的
**面面接触（μ = 0.5）** 载荷路径一致，也与本项目跖面已验收的
``plantar_bc="contact"``（`climbing.coupling.plantar_bc._build_pad_floor`）同构。

**recipe**：逐字沿用 `plantar_bc._apply_contact` 已验证的 NEW PENALTY 配方
（`laugon="PENALTY"`、`penalty=0.1`、`search_radius=20`、`node_reloc=0`、
`symmetric_stiffness=0`、`tolerance=0.005`）。

**默认关闭**：`joint_bc="pressure"` 走原路径，deck 逐位不变。

用法::

    from climbing.coupling.joint_contact import apply_joint_contact
    meta = apply_joint_contact(model, nodes=nodes, joint_local=joint_nodes,
                               bcs=bcs, load_n=F, mu=0.5)
    if meta.get("contact") is not None:
        st.contact = meta["contact"]
"""
from __future__ import annotations

import numpy as np

DEFAULT_MU = 0.5                 # [Y25] 原文 §1.1 面面接触 μ = 0.5
DEFAULT_PENALTY = 0.1            # 与 plantar_bc 已验证配方一致
DEFAULT_SEARCH_RADIUS = 20.0
DEFAULT_TOLERANCE = 0.005
DEFAULT_GAP_MM = 0.5             # 初始间隙（帽底 vs 关节面）
DEFAULT_CAP_E_MPA = 200_000.0    # 数值刚性（远高于皮质 15000）
DEFAULT_CAP_NU = 0.30
DEFAULT_CAP_THICK_MM = 5.0
DEFAULT_CAP_MARGIN_MM = 20.0
DEFAULT_CAP_NDIV = 4


def _local_to_fe(local) -> list[int]:
    """0-based 局部索引 → 1-based FEBio 节点 id。"""
    return [int(i) + 1 for i in local]


def _nodal_normals(nodes: np.ndarray, quads: np.ndarray) -> np.ndarray:
    """由 quad4 集算**逐节点**面积加权法向（单位化）。

    关节面是**凹面**（实测纵向跨 34 mm），平面帽只能单点接触 ⇒ 必须让帽底
    贴合各节点法向。法向取自 quad 的右手法则（`_build_thums_joint_quads`
    已保证朝外）。
    """
    n_acc = np.zeros((len(nodes), 3), dtype=np.float64)
    for q in np.asarray(quads, dtype=np.int64):
        p = nodes[q]
        nv = np.cross(p[1] - p[0], p[2] - p[0])
        for v in q:
            n_acc[int(v)] += nv
    norm = np.linalg.norm(n_acc, axis=1, keepdims=True)
    norm[norm == 0.0] = 1.0
    return n_acc / norm


def _gid_to_xyz(model, gid: int) -> np.ndarray:
    """按 FEBio 全局节点 id 从 model 的网格里取坐标（共形帽算面积用）。"""
    for nd in model.mesh_.nodes[0].all_nodes:
        if int(nd.id) == int(gid):
            xyz = [float(x) for x in str(nd.text).split(',')]
            return np.asarray(xyz, dtype=np.float64)
    raise KeyError(f"节点 id {gid} 不在 model 网格中")


def _hex8_signed_volume(p: np.ndarray) -> float:
    """hex8 的**有符号**体积（中心点 8 个四面体之和）。负号 ⇒ 节点序需反转。"""
    c = p.mean(axis=0)
    faces = ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))
    vol = 0.0
    for idx in faces:
        a, b, d, e = (p[i] for i in idx)
        for tri in ((a, b, d), (a, d, e)):
            vol += float(np.dot(tri[0] - c, np.cross(tri[1] - c, tri[2] - c))) / 6.0
    return vol


def _build_conforming_cap(
    model,
    nodes: np.ndarray,
    joint_quads: np.ndarray,
    *,
    gap_mm: float,
    thickness_mm: float,
    mat_id: int,
    mat_name: str,
    e_mpa: float,
    nu: float,
) -> dict:
    """**共形**距骨帽：底面逐节点沿关节面法向偏移 ``gap_mm``，顶面再偏 ``thickness_mm``。

    与平面版（`_build_talus_cap`）的区别：平面版在凹关节面上只有**一个最高点**
    接触 ⇒ 接触只传 ~5% 载荷；共形版贴合整个面 ⇒ 接触区延展整个关节面。

    每个 facet quad 生成 1 个 hex8（``bot`` 4 节点 + ``top`` 4 节点）；**节点序按
    有符号体积自动定向**（负体积即反转），不靠人工推导。
    """
    from pyfebio import material as fmat
    from pyfebio import mesh as fmesh
    from pyfebio import meshdomains as fdom

    q4 = np.asarray(joint_quads, dtype=np.int64).reshape(-1, 4)
    nrm = _nodal_normals(nodes, q4)

    used = np.unique(q4)
    bot_id: dict[int, int] = {}
    top_id: dict[int, int] = {}
    node_domain = model.mesh_.nodes[0]
    next_id = int(node_domain.all_nodes[-1].id) + 1
    for li in used:
        li = int(li)
        p = np.asarray(nodes, dtype=np.float64)[li]
        n = nrm[li]
        for store, off in ((bot_id, float(gap_mm)),
                           (top_id, float(gap_mm) + float(thickness_mm))):
            q = p + off * n
            store[li] = next_id
            node_domain.add_node(fmesh.Node(
                id=next_id, text=f"{q[0]:.10e},{q[1]:.10e},{q[2]:.10e}"))
            next_id += 1

    elems: list[list[int]] = []
    bot_quads: list[list[int]] = []
    top_quads: list[list[int]] = []
    flipped = 0
    for q in q4:
        b = [bot_id[int(v)] for v in q]
        t = [top_id[int(v)] for v in q]
        row = b + t
        pts = np.array([np.asarray(nodes, dtype=np.float64)[int(v)] for v in q])
        # 定向：用 offset 后的真实节点坐标算有符号体积
        allp = np.vstack([pts + float(gap_mm) * nrm[q],
                          pts + (float(gap_mm) + float(thickness_mm)) * nrm[q]])
        if _hex8_signed_volume(allp) < 0.0:
            row = t[::-1] + b[::-1]
            flipped += 1
        elems.append(row)
        # 接触面 = bot 面；法向须朝 **−n**（朝骨）⇒ 相对 hex 的 b 序反一下
        bot_quads.append(list(reversed(b)))
        top_quads.append(list(b))   # 顶面法向朝 +n（承载面）

    elems_arr = np.asarray(elems, dtype=np.int64)
    last_elem_id = int(model.mesh_.elements[-1].all_elements[-1].id)
    model.mesh_.add_element_domain(
        fmesh.numpy_to_elements(elems_arr, "hex8", name="talus_cap",
                                offset=last_elem_id))
    model.material_.add_material(fmat.IsotropicElastic(
        name=mat_name, id=int(mat_id),
        density=fmat.MaterialParameter(text=1.0e-9),
        E=fmat.MaterialParameter(text=float(e_mpa)),
        v=fmat.MaterialParameter(text=float(nu))))
    model.meshdomains_.add_solid_domain(
        fdom.SolidDomain(name="talus_cap", mat=mat_name))

    cap_ids = sorted(bot_id.values()) + sorted(top_id.values())
    model.mesh_.add_node_set(fmesh.NodeSet(
        name="talus_cap_set", text=",".join(str(i) for i in cap_ids)))
    bot_nodes = [(int(bot_id[int(li)]),
                  tuple(np.asarray(nodes, dtype=np.float64)[int(li)]
                        + float(gap_mm) * nrm[int(li)]))
                 for li in used]

    return {
        "bot_surface": bot_quads,
        "top_surface": top_quads,
        "cap_set": cap_ids,
        "bot_nodes": bot_nodes,
        "n_elems": int(len(elems)),
        "n_cap_nodes": int(len(cap_ids)),
        "n_quads": int(len(q4)),
        "n_flipped": int(flipped),
        "gap_mm": float(gap_mm),
        "thickness_mm": float(thickness_mm),
        "conforming": True,
    }


def _build_talus_cap(
    model,
    nodes: np.ndarray,
    joint_local: np.ndarray,
    *,
    gap_mm: float,
    thickness_mm: float,
    margin_mm: float,
    ndiv: int,
    mat_id: int,
    mat_name: str,
    e_mpa: float,
    nu: float,
) -> dict:
    """在关节面**上方**建一块固定 x/z 的薄 hex8 板（"距骨帽"）。

    ``_build_pad_floor`` 的 y 镜像版：
    * 板**底**面 ``y = pts.y.max() + gap`` —— 这是与关节面的接触面，**法向朝下**
    * 板**顶**面 ``y = bot + thickness`` —— 承载面（PressureLoad 施加于此）
    * 底面 quad 的节点序**反序**于 `_build_pad_floor.top_quads`，
      以便 FEBio 从**节点序**算出的法向朝下（朝骨）。

    返回 ``dict``：``bot_surface``（底面 quad，1-based 全局 id）、``top_local``
    （顶面节点的 0-based 局部索引序，供 PressureLoad 用）、``cap_set``（帽全部节点
    的 1-based id）、``bot_y_mm``、``n_elems``。
    """
    from pyfebio import mesh as fmesh

    pts = np.asarray(nodes, dtype=np.float64)[np.asarray(joint_local, dtype=np.int64)]
    x0, x1 = float(pts[:, 0].min()), float(pts[:, 0].max())
    z0, z1 = float(pts[:, 2].min()), float(pts[:, 2].max())
    x0 -= float(margin_mm)
    x1 += float(margin_mm)
    z0 -= float(margin_mm)
    z1 += float(margin_mm)
    bot_y = float(pts[:, 1].max()) + float(gap_mm)     # ← 镜像：在**上方**
    top_y = bot_y + float(thickness_mm)

    ndiv = max(int(ndiv), 1)
    node_domain = model.mesh_.nodes[0]
    next_id = int(node_domain.all_nodes[-1].id) + 1

    bot_ids: dict[tuple[int, int], int] = {}
    top_ids: dict[tuple[int, int], int] = {}
    for i in range(ndiv + 1):
        for k in range(ndiv + 1):
            x = x0 + (x1 - x0) * i / ndiv
            z = z0 + (z1 - z0) * k / ndiv
            bot_ids[(i, k)] = next_id
            node_domain.add_node(fmesh.Node(
                id=next_id, text=f"{x:.10e},{bot_y:.10e},{z:.10e}"))
            next_id += 1
            top_ids[(i, k)] = next_id
            node_domain.add_node(fmesh.Node(
                id=next_id, text=f"{x:.10e},{top_y:.10e},{z:.10e}"))
            next_id += 1

    elems: list[list[int]] = []
    bot_quads: list[list[int]] = []
    for i in range(ndiv):
        for k in range(ndiv):
            # hex8 节点序：**前 4 = y 较大的一侧**，后 4 = y 较小的一侧。
            # 与 `plantar_bc._build_pad_floor` 同款（那里前 4 是 ``top_ids``）。
            # 反过来会得到负 Jacobian（FEBio: "Negative jacobian detected
            # during mesh initialization"）。
            elems.append([
                top_ids[(i, k)], top_ids[(i + 1, k)],
                top_ids[(i + 1, k + 1)], top_ids[(i, k + 1)],
                bot_ids[(i, k)], bot_ids[(i + 1, k)],
                bot_ids[(i + 1, k + 1)], bot_ids[(i, k + 1)],
            ])
            # 底面（与骨相触）：法向必须朝 **−Y**（朝骨）。
            # 叉积复算：v1=(x_i,y,z_k) v2=(x_{i+1},y,z_k) v3=(x_{i+1},y,z_{k+1})
            #   v2−v1 = (+dx, 0, 0)；v3−v1 = (+dx, 0, +dz)
            #   (v2−v1)×(v3−v1) = (0·dz − 0·0, 0·dx − dx·dz, 0) = (0, **−dx·dz**, 0)
            #   ⇒ −Y ✓（若写成 [i,k]→[i,k+1]→… 会得到 +Y，接触把帽推向错误方向，
            #      表现为求解中 16 个帽单元全部 negative jacobian）。
            bot_quads.append([
                bot_ids[(i, k)], bot_ids[(i + 1, k)],
                bot_ids[(i + 1, k + 1)], bot_ids[(i, k + 1)],
            ])

    elems_arr = np.asarray(elems, dtype=np.int64)
    last_elem_id = int(model.mesh_.elements[-1].all_elements[-1].id)
    model.mesh_.add_element_domain(
        fmesh.numpy_to_elements(elems_arr, "hex8", name="talus_cap",
                                offset=last_elem_id)
    )

    from pyfebio import material as fmat
    from pyfebio import meshdomains as fdom

    model.material_.add_material(
        fmat.IsotropicElastic(
            name=mat_name,
            id=int(mat_id),
            density=fmat.MaterialParameter(text=1.0e-9),
            E=fmat.MaterialParameter(text=float(e_mpa)),
            v=fmat.MaterialParameter(text=float(nu)),
        )
    )
    model.meshdomains_.add_solid_domain(
        fdom.SolidDomain(name="talus_cap", mat=mat_name))

    cap_ids = sorted(bot_ids.values()) + sorted(top_ids.values())
    model.mesh_.add_node_set(fmesh.NodeSet(
        name="talus_cap_set", text=",".join(str(i) for i in cap_ids)))

    # 顶面节点按 **(i,k) 网格序**返回（供按同一分格重建 quad4）；
    # 另给一个排序版作 node_set 用。
    top_id_list = [top_ids[(i, k)] for i in range(ndiv + 1) for k in range(ndiv + 1)]
    model.mesh_.add_node_set(fmesh.NodeSet(
        name="talus_cap_top_set",
        text=",".join(str(i) for i in sorted(top_id_list))))

    return {
        "bot_surface": bot_quads,
        "cap_set": cap_ids,
        "top_ids": top_id_list,          # (i,k) 网格序
        "ndiv": ndiv,
        "bot_y_mm": bot_y,
        "top_y_mm": top_y,
        "n_elems": int(len(elems)),
        "n_cap_nodes": int(len(cap_ids)),
        # 供 hold 弹簧用：底面每个节点的 (id, 坐标)
        "bot_nodes": [(int(bot_ids[(i, k)]),
                       (x0 + (x1 - x0) * i / ndiv, bot_y,
                        z0 + (z1 - z0) * k / ndiv))
                      for i in range(ndiv + 1) for k in range(ndiv + 1)],
    }


def _apply_cap_hold_spring(
    model,
    bcs: list,
    bot_nodes,
    *,
    k_n_per_mm: float,
    offset_mm: float = 1.0,
    set_name: str = "talus_cap_hold",
) -> dict:
    """帽底节点沿 **y** 各连一根**很弱**的线性弹簧到固定锚点。

    纯**数值正则化**：单边接触在闭合前 y 方向零刚度 ⇒ 刚度矩阵奇异
    （实测帽位移爆到 **7.18e12**，8 个单元倒置）。弹簧刚度取 ≪ 骨刚度
    （默认 1 N/mm/节点 × 50 节点 = 50 N/mm，对 1e5 N/mm 量级的骨可忽略），
    仅用于消掉刚体模态；接触承载后其贡献可忽略 —— 交付时必须**报告弹簧总力
    与总载荷之比**作为「可忽略」的证据。

    与 `climbing.coupling.plantar_bc._apply_spring`（跖面同款做法，
    `contact_hold_spring_k=100`）同构。
    """
    from pyfebio import boundary as fbc
    from pyfebio import discrete
    from pyfebio import mesh as fmesh

    k = float(k_n_per_mm)
    if not (np.isfinite(k) and k > 0.0):
        raise ValueError(f"k_n_per_mm 必须为正有限值，得到 {k_n_per_mm!r}")
    node_domain = model.mesh_.nodes[0]
    next_id = int(node_domain.all_nodes[-1].id) + 1
    elems: list = []
    anchors: list[int] = []
    for gid, (x, y, z) in bot_nodes:
        anchor = next_id
        anchors.append(anchor)
        node_domain.add_node(fmesh.Node(
            id=next_id, text=f"{x:.10e},{y + float(offset_mm):.10e},{z:.10e}"))
        elems.append(fmesh.DiscreteElement(text=f"{int(gid)},{next_id}"))
        next_id += 1
    model.mesh_.add_discrete_set(fmesh.DiscreteSet(name=set_name, elements=elems))
    mat_id = int(len(model.discrete_.discrete_materials)) + 1
    model.discrete_.add_discrete_material(
        discrete.Spring(id=mat_id, name=f"{set_name}_mat", E=k))
    model.discrete_.add_discrete_element(
        discrete.DiscreteEntry(dmat=mat_id, discrete_set=set_name))
    anchor_set = f"{set_name}_anchor"
    model.mesh_.add_node_set(fmesh.NodeSet(
        name=anchor_set, text=",".join(str(i) for i in anchors)))
    bcs.append(fbc.BCZeroDisplacement(
        node_set=anchor_set, x_dof=1, y_dof=1, z_dof=1))
    return {"k_n_per_mm": k, "n_springs": int(len(elems)),
            "n_anchors": int(len(anchors)), "offset_mm": float(offset_mm)}


def apply_joint_contact(
    model,
    nodes: np.ndarray,
    joint_local: np.ndarray,
    bcs: list,
    *,
    joint_surface_name: str = "subtalar_joint",
    load_n: float,
    joint_tris: np.ndarray | None = None,
    mu: float = DEFAULT_MU,
    penalty: float = DEFAULT_PENALTY,
    search_radius: float = DEFAULT_SEARCH_RADIUS,
    tolerance: float = DEFAULT_TOLERANCE,
    gap_mm: float = DEFAULT_GAP_MM,
    cap_e_mpa: float = DEFAULT_CAP_E_MPA,
    cap_nu: float = DEFAULT_CAP_NU,
    cap_thickness_mm: float = DEFAULT_CAP_THICK_MM,
    cap_margin_mm: float = DEFAULT_CAP_MARGIN_MM,
    cap_ndiv: int = DEFAULT_CAP_NDIV,
    laugon: str = "PENALTY",
    auto_penalty: int = 0,
    two_pass: int = 0,
    node_reloc: int = 0,
    symmetric_stiffness: int = 0,
    swap_pair: bool = False,
    create_bone_surface: bool = True,
    joint_quads: np.ndarray | None = None,
    joint_surface_name_q4: str = "subtalar_joint_q4",
    cap_hold_spring_k: float | None = None,
    cap_hold_offset_mm: float = 1.0,
) -> dict:
    """把距下关节的载荷改成「刚性距骨帽 + 真实接触」。

    ``create_bone_surface``：``thums_feb.build_thums_feb`` 已经建过
    ``<Surface name=joint_face>``（FEBio 拒绝重名，报
    ``A surface with name "..." was already defined``）⇒ 从 builder 调用时传
    ``False``。

    返回 ``dict``，含 ``"contact"``（挂到 step 的 ``<Contact>``）、
    ``"cap_top_surface"``（PressureLoad 的 surface 名）、``"pressure_mpa"``。
    """
    from pyfebio import contact as fcontact
    from pyfebio import mesh as fmesh
    from pyfebio import boundary as fbc

    joint_local = np.asarray(joint_local, dtype=np.int64)
    # **共形帽**：给了 joint_quads 就走共形（逐节点沿法向偏移），
    # 否则回落平面帽（合成网格 / 单域测试用）。
    if joint_quads is not None:
        cap = _build_conforming_cap(
            model, nodes, joint_quads,
            gap_mm=float(gap_mm), thickness_mm=float(cap_thickness_mm),
            mat_id=int(len(model.material_.all_materials)) + 1,
            mat_name="talus_cap_mat", e_mpa=float(cap_e_mpa), nu=float(cap_nu),
        )
    else:
        cap = _build_talus_cap(
            model, nodes, joint_local,
            gap_mm=float(gap_mm), thickness_mm=float(cap_thickness_mm),
            margin_mm=float(cap_margin_mm), ndiv=int(cap_ndiv),
            mat_id=int(len(model.material_.all_materials)) + 1,
            mat_name="talus_cap_mat", e_mpa=float(cap_e_mpa), nu=float(cap_nu),
        )

    # --- 帽的 BC：只放开 y（上/下压），x = z = 0 ------------------------------
    bcs.append(fbc.BCZeroDisplacement(node_set="talus_cap_set", x_dof=1, z_dof=1))

    # --- opt-in：帽的弱 y 支撑弹簧（消单边接触的刚体模态） --------------------
    hold_meta: dict | None = None
    if cap_hold_spring_k is not None and float(cap_hold_spring_k) > 0.0:
        hold_meta = _apply_cap_hold_spring(
            model, bcs, cap["bot_nodes"],
            k_n_per_mm=float(cap_hold_spring_k),
            offset_mm=float(cap_hold_offset_mm))

    # --- 1. 关节面 Surface（骨侧） ------------------------------------------
    # **必须 quad4**：距下关节的节点全在 CORT hex8 域上，tri3 facet 会被 FEBio
    # 判 ``invalid facets`` ⇒ 接触不承载（实测帽的 y 位移爆到 7.18e12）。
    # 与跖面同一坑同一解（`docs/S5_contact_facet_fix.md` §1.3）。
    surf_name = joint_surface_name_q4 if joint_quads is not None else joint_surface_name
    if joint_quads is not None:
        q4 = np.asarray(joint_quads, dtype=np.int64).reshape(-1, 4)
        if len(q4) == 0:
            raise ValueError("joint_quads 非空但为 0 行 —— 拒绝静默退化到 tri3。")
        surf = fmesh.Surface(name=surf_name)
        for i, quad in enumerate(q4):
            surf.add_quad4(fmesh.Quad4Element(
                id=i + 1, text=",".join(map(str, _local_to_fe(quad)))))
        model.mesh_.add_surface(surf)
    elif create_bone_surface:
        surf = fmesh.Surface(name=joint_surface_name)
        if joint_tris is not None and len(joint_tris) > 0:
            tris = np.asarray(joint_tris, dtype=np.int64).reshape(-1, 3)
            for i, tri in enumerate(tris):
                surf.add_tri3(fmesh.Tri3Element(
                    id=i + 1, text=",".join(map(str, _local_to_fe(tri)))))
        else:
            raise ValueError(
                "apply_joint_contact 需要 joint_tris（距下关节面的 tri3 连接）；"
                "THUMS 网格请传 mesh['surfaces']['subtalar_joint']"
            )
        model.mesh_.add_surface(surf)

    # --- 2. 帽底 Surface ---------------------------------------------------
    cap_surf = fmesh.Surface(name="talus_cap_bottom")
    for i, quad in enumerate(cap["bot_surface"]):
        cap_surf.add_quad4(fmesh.Quad4Element(
            id=i + 1, text=",".join(str(j) for j in quad)))
    model.mesh_.add_surface(cap_surf)

    # --- 3. 帽顶 Surface（PressureLoad 施加于此） ---------------------------
    # 共形帽：顶面 quad 直接由 `cap["top_surface"]` 给出（与 bot 同拓扑）。
    # 平面帽：按 (ndiv+1)² 分格重建。
    if cap.get("conforming"):
        cap_top = fmesh.Surface(name="talus_cap_top")
        for i, quad in enumerate(cap["top_surface"]):
            cap_top.add_quad4(fmesh.Quad4Element(
                id=i + 1, text=",".join(str(j) for j in quad)))
        model.mesh_.add_surface(cap_top)
    else:
        nd = int(cap["ndiv"])
        top_ids = cap["top_ids"]
        top_quads: list[list[int]] = []
        for i in range(nd):
            for k in range(nd):
                a = i * (nd + 1) + k
                top_quads.append([top_ids[a], top_ids[a + 1],
                                  top_ids[a + nd + 2], top_ids[a + nd + 1]])
        cap_top = fmesh.Surface(name="talus_cap_top")
        for i, quad in enumerate(top_quads):
            cap_top.add_quad4(fmesh.Quad4Element(
                id=i + 1, text=",".join(str(j) for j in quad)))
        model.mesh_.add_surface(cap_top)

    # 帽顶面积（mm²）⇒ p = F/A_cap_top
    # 共形帽：直接对各 top quad 求面积（贴合关节面，不是包围盒）；
    # 平面帽：包围盒 + margin（保持原口径）。
    if cap.get("conforming"):
        a_cap = 0.0
        for quad in cap["top_surface"]:
            pts = np.array([_gid_to_xyz(model, g) for g in quad])
            a_cap += 0.5 * float(np.linalg.norm(
                np.cross(pts[1] - pts[0], pts[2] - pts[0])))
            a_cap += 0.5 * float(np.linalg.norm(
                np.cross(pts[2] - pts[0], pts[3] - pts[0])))
    else:
        pts_bot = np.asarray(nodes, dtype=np.float64)[joint_local]
        dx = (float(pts_bot[:, 0].max()) - float(pts_bot[:, 0].min())
              + 2.0 * float(cap_margin_mm))
        dz = (float(pts_bot[:, 2].max()) - float(pts_bot[:, 2].min())
              + 2.0 * float(cap_margin_mm))
        a_cap = dx * dz
    a_cap = max(float(a_cap), 1e-12)
    p_cap = float(load_n) / a_cap

    # --- 4. SurfacePair + Contact -----------------------------------------
    pair_primary, pair_secondary = (
        (surf_name, "talus_cap_bottom") if swap_pair
        else ("talus_cap_bottom", surf_name)
    )
    model.mesh_.add_surface_pair(fmesh.SurfacePair(
        name="subtalar_talus_pair",
        primary=pair_primary, secondary=pair_secondary))

    contact_obj = fcontact.Contact()
    contact_obj.add_contact(fcontact.SlidingElastic(
        name="subtalar_talus_contact",
        surface_pair="subtalar_talus_pair",
        type="sliding-elastic",
        laugon=str(laugon).upper(),
        penalty=float(penalty),
        auto_penalty=int(auto_penalty),
        two_pass=int(two_pass),
        node_reloc=int(node_reloc),
        symmetric_stiffness=int(symmetric_stiffness),
        tolerance=float(tolerance),
        search_radius=float(search_radius),
        fric_coeff=float(mu),
    ))

    return {
        "kind": "contact",
        "contact": contact_obj,
        "cap_top_surface": "talus_cap_top",
        "pressure_mpa": float(p_cap),
        "a_cap_mm2": float(a_cap),
        "gap_mm": float(gap_mm),
        "bot_y_mm": float(cap.get("bot_y_mm", float("nan"))),
        "top_y_mm": float(cap.get("top_y_mm", float("nan"))),
        "conforming": bool(cap.get("conforming", False)),
        "n_quads": int(cap.get("n_quads", 0)),
        "n_flipped": int(cap.get("n_flipped", 0)),
        "n_cap_nodes": int(cap["n_cap_nodes"]),
        "n_cap_elems": int(cap["n_elems"]),
        "mu": float(mu),
        "pair_primary": str(pair_primary),
        "pair_secondary": str(pair_secondary),
        "load_n": float(load_n),
        "hold_spring_k_n_per_mm": (
            None if hold_meta is None else float(hold_meta["k_n_per_mm"])),
        "hold_meta": hold_meta,
    }
