"""可切换的跖面（plantar）边界条件：``fixed`` / ``roller`` / ``spring`` / ``contact``。

背景
----
现有跟骨 FE（``febio_model.build_calcaneus_feb`` / ``scripts/opensim_fe/thums_feb.py``）
把跖面**三向全固定**（``BCZeroDisplacement(x=y=z=1)``）。地面上骨与地面之间其实
只传**压力 + 摩擦/滑移**，全固定会在约束面棱边人为制造**固定端奇异**，使
``gauge`` 虚高、首次骨折高度偏低。本模块把跖面 BC 参数化：

``fixed``
    三向零位移（原行为，**默认**，保证现有回归逐位不变）。
``roller``
    只约束**法向**（``BCNormalDisplacement``，切向自由）→ 去掉切向剪力/固定端
    弯矩；再对 2 个跖面节点补 x/z 固定，消除"绕竖轴转动 + 水平平移"的刚体奇异。
``roller_free``
    诊断用变体：全局 Y 零位移 + 极软切向弹簧（放开切向）。
``spring``
    Winkler 弹性地基：每个跖面节点向**三个正交锚点**（−Y 法向 + X/Z 切向）各连
    一根线性弹簧（``k`` N/mm），锚点三向固定。这样三向刚体模态都被弹性约束住，
    无需附加刚性销；``k → ∞`` 时代数地逼近 ``fixed`` 基线（这是正确性自检）。
``contact``
    **波3 最小版**（``docs/S5接触方案.md`` §3.1）：跖面 ⇄ 对偶面（一**全节点固定**
    的 ``hex8`` 薄板 = 数值刚性壁）+ Coulomb 摩擦（``μ``）。**不建泡沫实体**，故
    不依赖泡沫本构标定。recipe = ``laugon=PENALTY`` + ``penalty=0.1`` +
    ``search_radius=20``，取自 :mod:`climbing.coupling.pad_contact` 的 udg 验证配方
    （NEW PENALTY 配方；见 ``docs/S5_contact_udg.md`` §7 与 ``docs/S5_QUOTABLE.md``）。
    **默认不启用**（``plantar_bc`` 默认仍是 ``"fixed"``），opt-in。
    FEBio 4.0 的 ``<Mesh><rigid_wall>`` 与 ``<contact type="rigid_wall">`` 均被本机
    4.13 拒绝（见 :func:`_build_pad_floor`），故用固定薄板代替解析壁。

单位 mm–N–MPa：弹簧 ``E`` = 刚度 k（N/mm）；接触 ``penalty`` = 无量纲缩放因子，
``μ`` = 无量纲；``search_radius`` = mm。

设计约束
--------
* 顶层不 import pyfebio（``febio_model`` 要在无 pyfebio 环境可 import）。
* ``fixed`` 路径必须与旧 XML **逐位一致**：不额外添加 surface / node_set。
* ``contact`` 路径新增 1 个 ``<Surface name="plantar">``（即便已有同名
  ``NodeSet`` 也复用——不影响现有 bc）；1 个 ``<SurfacePair>``；1 个
  ``<Contact>`` 块（写入 builder 的 step 内）。``<rigid_wall>`` 由 builder
  后处理注入，因 pyfebio 0.3.0 没有 rigidwall 元素类。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "BC_KINDS",
    "apply_plantar_bc",
    "DEFAULT_CONTACT_MU",
    "DEFAULT_CONTACT_PENALTY",
    "DEFAULT_CONTACT_SEARCH_RADIUS",
    "DEFAULT_CONTACT_TOLERANCE",
    "DEFAULT_CONTACT_GAP_MM",
]

#: 支持的跖面 BC 类型（顺序即文档顺序）。
#: ``roller_free`` 是诊断用变体：全局 Y 固定 + 极软切向弹簧，真正放开切向。
#: ``contact`` 是波3 最小版：刚性壁 + Coulomb 摩擦（opt-in）。
BC_KINDS = ("fixed", "roller", "roller_free", "spring", "contact")

#: ``plantar_bc="contact"`` 的默认接触配方参数（opt-in）。
#: recipe 逐字取自 ``pad_contact.py`` 的 udg NEW PENALTY 配方：
#: ``laugon=PENALTY`` + ``penalty=0.1`` + ``search_radius=20`` +
#: ``node_reloc=0`` + ``tolerance=0.005`` + ``fric_coeff=mu``。
#: 见 ``docs/S5_contact_udg.md`` §7、``docs/S5_QUOTABLE.md``。
DEFAULT_CONTACT_MU: float = 0.6
DEFAULT_CONTACT_PENALTY: float = 0.1
DEFAULT_CONTACT_SEARCH_RADIUS: float = 20.0
DEFAULT_CONTACT_TOLERANCE: float = 0.005
#: 对偶面（固定薄板）顶面到跖面最低点的初始间隙（mm）。FEBio 接触搜索要求严格
#: 正间隙才能播种（重合面不被检测到）；与 ``pad_contact.DEFAULT_GAP_MM`` 一致。
DEFAULT_CONTACT_GAP_MM: float = 0.5

#: 方向 key -> 锚点偏移方向（单位向量，取一个正交方向即可）。
_DIRECTION_OFFSET = {
    "normal": (0.0, -1.0, 0.0),  # 跖面朝 −Y
    "x": (1.0, 0.0, 0.0),
    "z": (0.0, 0.0, 1.0),
}


def _local_to_fe(local: np.ndarray) -> list[int]:
    """0-based 局部节点索引 -> 1-based FEBio 节点 id（list）。"""
    return [int(i) + 1 for i in np.asarray(local, dtype=np.int64).ravel()]


def _add_node_set(model, fmesh, name: str, local: np.ndarray) -> None:
    model.mesh_.add_node_set(
        fmesh.NodeSet(name=name, text=",".join(map(str, _local_to_fe(local))))
    )


def _pick_lateral_pins(nodes: np.ndarray, plantar_local: np.ndarray) -> tuple[int, int]:
    """挑 2 个跖面节点做切向销：最靠后(min X) + 离它最远的节点。

    这两个节点在 x/z 上各固定一次，即可消掉"水平平移 + 绕竖直轴转动"这 3 个
    仅靠法向约束无法约束的刚体模态。
    """
    local = np.asarray(plantar_local, dtype=np.int64)
    pts = np.asarray(nodes, dtype=np.float64)[local]
    i0 = int(np.argmin(pts[:, 0]))
    dist = np.linalg.norm(pts - pts[i0], axis=1)
    i1 = int(np.argmax(dist))
    if i1 == i0:  # 退化保护（不该发生：plantar 是多节点面）
        i1 = 0 if i0 != 0 else min(1, len(local) - 1)
    return int(local[i0]), int(local[i1])


def _apply_roller(
    model, nodes: np.ndarray, plantar_local: np.ndarray, plantar_tris: np.ndarray, bcs: list
) -> dict:
    """法向零位移 + 2 个切向销。"""
    from pyfebio import boundary as fbc
    from pyfebio import loaddata as fld
    from pyfebio import mesh as fmesh

    if plantar_tris is None or len(plantar_tris) == 0:
        raise ValueError("roller BC 需要非空的 plantar 三角面（mesh['surfaces']['plantar']）")

    surf = fmesh.Surface(name="plantar")
    for i, tri in enumerate(np.asarray(plantar_tris, dtype=np.int64)):
        surf.add_tri3(fmesh.Tri3Element(id=i + 1, text=",".join(map(str, _local_to_fe(tri)))))
    model.mesh_.add_surface(surf)

    # FEBio 要求 <scale> 的 lc 指向已定义的 load controller（lc="0" 非法）。
    # 定义一个恒 =1 的 loadcurve（id=1），scale text=0.0 => 法向位移 = 0.0×1。
    model.loaddata_.add_load_curve(
        fld.LoadCurve(
            id=1,
            points=fld.CurvePoints(points=["0.0,1.0", "1.0,1.0"]),
        )
    )
    bcs.append(
        fbc.BCNormalDisplacement(surface="plantar", scale=fbc.Value(lc=1, text=0.0))
    )

    a, b = _pick_lateral_pins(nodes, plantar_local)
    _add_node_set(model, fmesh, "plantar_pin_a", np.array([a]))
    _add_node_set(model, fmesh, "plantar_pin_b", np.array([b]))
    bcs.append(fbc.BCZeroDisplacement(node_set="plantar_pin_a", x_dof=1, z_dof=1))
    bcs.append(fbc.BCZeroDisplacement(node_set="plantar_pin_b", x_dof=1, z_dof=1))
    return {"pins": [a, b], "n_surface_tris": int(len(plantar_tris))}


def _apply_spring(
    model,
    nodes: np.ndarray,
    plantar_local: np.ndarray,
    bcs: list,
    spring_k: float,
    anchor_offset_mm: float,
    directions: tuple[str, ...] = ("normal", "x", "z"),
    set_name: str = "plantar_spring",
    mat_id: int = 1,
) -> dict:
    """Winkler 弹性地基：每个跖面节点沿 ``directions`` 各连一根线性弹簧到固定锚点。"""
    from pyfebio import boundary as fbc
    from pyfebio import discrete
    from pyfebio import mesh as fmesh

    k = float(spring_k)
    if not np.isfinite(k) or k <= 0.0:
        raise ValueError(f"spring_k 必须为正有限值 (N/mm)，得到 {spring_k!r}")
    d = float(anchor_offset_mm)
    if not np.isfinite(d) or d <= 0.0:
        raise ValueError(f"anchor_offset_mm 必须为正有限值，得到 {anchor_offset_mm!r}")
    dirs = tuple(str(x) for x in directions)
    if not dirs or any(x not in _DIRECTION_OFFSET for x in dirs):
        raise ValueError(f"directions 非法：{directions!r}；可选 {tuple(_DIRECTION_OFFSET)}")

    node_domain = model.mesh_.nodes[0]
    local = np.asarray(plantar_local, dtype=np.int64)

    elems: list = []
    anchor_ids: list[int] = []
    # 用**实际最后一个节点 id** 续号（而非 ``len(nodes)``）：接触路径会在
    # ``_apply_contact`` 里先追加一块固定薄板节点，此时 ``len(nodes)`` 与真实
    # 最大 id 不符，沿用旧口径会与薄板节点 id 冲突（FEBio: invalid id）。
    # 默认 spring 路径下两者相等 ⇒ 逐位不变。
    next_id = int(node_domain.all_nodes[-1].id) + 1
    for li in local:
        p = np.asarray(nodes, dtype=np.float64)[int(li)]
        for dk in dirs:
            off = d * np.asarray(_DIRECTION_OFFSET[dk], dtype=float)
            anchor_ids.append(next_id)
            coord = p + off
            node_domain.add_node(
                fmesh.Node(id=next_id, text=f"{coord[0]:.10e},{coord[1]:.10e},{coord[2]:.10e}")
            )
            elems.append(fmesh.DiscreteElement(text=f"{int(li) + 1},{next_id}"))
            next_id += 1

    model.mesh_.add_discrete_set(fmesh.DiscreteSet(name=set_name, elements=elems))
    model.discrete_.add_discrete_material(
        discrete.Spring(id=int(mat_id), name=f"{set_name}_mat", E=k)
    )
    model.discrete_.add_discrete_element(
        discrete.DiscreteEntry(dmat=int(mat_id), discrete_set=set_name)
    )

    anchor_set = f"{set_name}_anchor"
    model.mesh_.add_node_set(
        fmesh.NodeSet(name=anchor_set, text=",".join(map(str, anchor_ids)))
    )
    bcs.append(
        fbc.BCZeroDisplacement(node_set=anchor_set, x_dof=1, y_dof=1, z_dof=1)
    )
    return {
        "spring_k_n_per_mm": k,
        "directions": list(dirs),
        "n_springs": int(len(elems)),
        "n_anchors": int(len(anchor_ids)),
        "anchor_offset_mm": d,
    }


def _apply_roller_free(
    model,
    nodes: np.ndarray,
    plantar_local: np.ndarray,
    bcs: list,
    tangential_k: float,
    anchor_offset_mm: float,
) -> dict:
    """真正的"放开切向"滚子：全局 Y（竖直法向）零位移 + 极软切向弹簧。

    ``BCNormalDisplacement`` 在**曲面**跖面上会对共享节点施加多个面法向约束，
    合起来等价于全固定（本仓库实测与 fixed 逐位一致）。因此这里改用**全局竖直
    方向**的法向约束（水平地面接触的法向），切向用**很软**的 Winkler 弹簧吸收
    三个刚体模态（无硬销 → 不引入人工应力集中）。
    """
    from pyfebio import boundary as fbc

    bcs.append(fbc.BCZeroDisplacement(node_set="plantar", y_dof=1))
    meta = _apply_spring(
        model, nodes, plantar_local, bcs, tangential_k, anchor_offset_mm,
        directions=("x", "z"), set_name="plantar_spring_tan",
    )
    meta["normal_bc"] = "zero displacement y (global)"
    return meta


def _build_pad_floor(
    model,
    nodes: np.ndarray,
    plantar_local: np.ndarray,
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
    """在跖面下方建一块**固定薄板**（``hex8``）作为接触对偶面（数值刚性壁）。

    FEBio 4.0 格式的 ``<contact type="rigid_wall">`` / ``<Mesh><rigid_wall>`` 均
    被本机 4.13 拒绝（实测 ``unrecognized tag`` / ``invalid value for attribute
    "type"``），而 pyfebio 0.3.0 也没有 rigidwall 元素类。因此最小版用一块**全节点
    固定**的薄板代替解析刚性壁：其顶面是接触对偶面，刚度不在物理上起作用（全部
    节点被 BC 锁死），数值上等价于刚性壁。这与 ``pad_contact.py`` 模块文档里
    "stiff plate 与 rigid body 数值上不可区分" 的既有结论一致。

    返回 ``dict``：``top_surface``（``(M,4)`` 0-based 局部→全局的 quad 连接）、
    ``node_ids``（1-based 全局 id 列表）、``wall_y_mm``、``n_elems``。
    """
    from pyfebio import mesh as fmesh

    pts = np.asarray(nodes, dtype=np.float64)[np.asarray(plantar_local, dtype=np.int64)]
    x0, x1 = float(pts[:, 0].min()), float(pts[:, 0].max())
    z0, z1 = float(pts[:, 2].min()), float(pts[:, 2].max())
    x0 -= float(margin_mm)
    x1 += float(margin_mm)
    z0 -= float(margin_mm)
    z1 += float(margin_mm)
    wall_y = float(pts[:, 1].min()) - float(gap_mm)
    bot_y = wall_y - float(thickness_mm)

    ndiv = max(int(ndiv), 1)
    node_domain = model.mesh_.nodes[0]
    # 全局 id 续号（与 ``Model.add_simple_rigid_body`` 同款：last id + 1），
    # 对单域 / 多域（THUMS CORT+SPON）都安全。
    next_id = int(node_domain.all_nodes[-1].id) + 1

    top_ids: dict[tuple[int, int], int] = {}
    bot_ids: dict[tuple[int, int], int] = {}
    for i in range(ndiv + 1):
        for k in range(ndiv + 1):
            x = x0 + (x1 - x0) * i / ndiv
            z = z0 + (z1 - z0) * k / ndiv
            top_ids[(i, k)] = next_id
            node_domain.add_node(fmesh.Node(id=next_id, text=f"{x:.10e},{wall_y:.10e},{z:.10e}"))
            next_id += 1
            bot_ids[(i, k)] = next_id
            node_domain.add_node(fmesh.Node(id=next_id, text=f"{x:.10e},{bot_y:.10e},{z:.10e}"))
            next_id += 1

    elems: list[list[int]] = []
    top_quads: list[list[int]] = []
    for i in range(ndiv):
        for k in range(ndiv):
            elems.append([
                top_ids[(i, k)], top_ids[(i + 1, k)], top_ids[(i + 1, k + 1)], top_ids[(i, k + 1)],
                bot_ids[(i, k)], bot_ids[(i + 1, k)], bot_ids[(i + 1, k + 1)], bot_ids[(i, k + 1)],
            ])
            top_quads.append([
                top_ids[(i, k)], top_ids[(i, k + 1)], top_ids[(i + 1, k + 1)], top_ids[(i + 1, k)],
            ])

    elems_arr = np.asarray(elems, dtype=np.int64)
    # 续号接在**最后一个**单元域之后（THUMS 有 calcaneus + trabecular 两个域），
    # 与 ``add_simple_rigid_body`` 的 ``elements[-1].all_elements[-1].id + 1`` 一致。
    last_elem_id = int(model.mesh_.elements[-1].all_elements[-1].id)
    model.mesh_.add_element_domain(
        fmesh.numpy_to_elements(elems_arr, "hex8", name="pad_floor", offset=last_elem_id)
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
    model.meshdomains_.add_solid_domain(fdom.SolidDomain(name="pad_floor", mat=mat_name))

    node_ids = sorted(top_ids.values()) + sorted(bot_ids.values())
    model.mesh_.add_node_set(
        fmesh.NodeSet(name="pad_floor_fixed", text=",".join(str(i) for i in node_ids))
    )
    return {
        "top_surface": top_quads,
        "node_ids": node_ids,
        "wall_y_mm": wall_y,
        "n_elems": int(len(elems)),
    }


def _apply_contact(
    model,
    nodes: np.ndarray,
    plantar_local: np.ndarray,
    plantar_tris: np.ndarray,
    bcs: list,
    *,
    mu: float = DEFAULT_CONTACT_MU,
    penalty: float = DEFAULT_CONTACT_PENALTY,
    laugon: str = "PENALTY",
    search_radius: float = DEFAULT_CONTACT_SEARCH_RADIUS,
    tolerance: float = DEFAULT_CONTACT_TOLERANCE,
    gap_mm: float = DEFAULT_CONTACT_GAP_MM,
    contact_type: str = "sliding-elastic",
    auto_penalty: int = 0,
    two_pass: int = 0,
    node_reloc: int = 0,
    symmetric_stiffness: int = 0,
    fric_coeff: float | None = None,
    pad_e_mpa: float = 200_000.0,
    pad_nu: float = 0.30,
    pad_thickness_mm: float = 5.0,
    pad_margin_mm: float = 20.0,
    pad_ndiv: int = 4,
    plantar_quads: np.ndarray | None = None,
    swap_pair: bool = False,
    hold_spring_k: float | None = None,
    hold_anchor_offset_mm: float = 1.0,
) -> dict:
    """波3 最小版跖面接触：跖面（骨 facet）⇄ **固定薄板**（数值刚性壁）+ Coulomb μ。

    recipe 逐字沿用 :mod:`climbing.coupling.pad_contact` 已验证的 NEW PENALTY 配方
    （``laugon="PENALTY"`` + ``penalty=0.1`` + ``search_radius=20`` + ``node_reloc=0``
    + ``symmetric_stiffness=0`` + ``tolerance=0.005`` + ``fric_coeff=μ``，见
    ``docs/S5_contact_udg.md`` §7、``docs/S5_QUOTABLE.md``）。

    全部参数是 **opt-in kwarg**；顶层 :func:`apply_plantar_bc` 只在 ``kind="contact"``
    时路由到此，故默认行为逐位不变。返回的 ``dict`` 含 ``"contact"`` 键
    （pyfebio :class:`pyfebio.contact.Contact` 对象）；builder 需把它挂到 step 的
    ``<Contact>``（``st.contact = meta["contact"]``）——因为 ``apply_plantar_bc``
    只拿到 ``model``，拿不到 ``step``。

    Plantar facet type（**新：5 选项**）

    - **未传** ``plantar_quads``（默认 ``None``）→ 走 tri3 路径（合成网格 +
      旧版测试保持兼容；这是波3 落地的语义）
    - **传** ``plantar_quads``（``(K, 4)`` 0-based 局部节点索引）→ 走 quad4 路径
      （**真实 THUMS 网格的跖面属于 CORT hex8 域**，因此正确 facet 类型是 quad4。
      tri3 facet 挂在 hex8 域上 FEBio 报 ``invalid facets`` ⇒ 接触不承载，详见
      ``docs/S5_contact_facet_fix.md`` §1）。

    防「静默退化」：``plantar_quads`` 非空但元素行数 0 → 显式抛 ``ValueError``
    （绝不回退到 tri3）。
    """
    from pyfebio import boundary as fbc
    from pyfebio import contact as fcontact
    from pyfebio import mesh as fmesh

    tris = np.asarray(plantar_tris, dtype=np.int64)
    if plantar_quads is None:
        # tri3 路径（合成 / 旧版兼容路径）
        if tris.size == 0 or tris.ndim != 2 or tris.shape[1] != 3:
            raise ValueError(
                "contact BC 需要非空的 plantar 三角面 (K,3) "
                f"（mesh['surfaces']['plantar']），得到 shape={tris.shape}"
            )
        quads_for_surface: np.ndarray = np.empty((0, 4), dtype=np.int64)
    else:
        quads_for_surface = np.asarray(plantar_quads, dtype=np.int64)
        if quads_for_surface.ndim != 2 or quads_for_surface.shape[1] != 4:
            raise ValueError(
                "plantar_quads 必须是 (K,4) 的 0-based 局部节点索引，"
                f"得到 shape={quads_for_surface.shape}"
            )
        if len(quads_for_surface) == 0:
            # 防静默退化：调用方明确要走 quad4，但提供的 quad4 集为空 ⇒ 抛错
            raise ValueError(
                "plantar_quads 为空：调用方要求走 quad4 路径但 quad4 集为空。"
                "请检查父网格的 boundary_polys / 跖面分类；绝不静默回退到 tri3。"
            )
    if not (np.isfinite(mu) and mu >= 0.0):
        raise ValueError(f"mu 必须有限且 >= 0，得到 {mu!r}")
    if not (np.isfinite(penalty) and penalty > 0.0):
        raise ValueError(f"penalty 必须有限且 > 0，得到 {penalty!r}")
    if not (np.isfinite(search_radius) and search_radius > 0.0):
        raise ValueError(f"search_radius 必须有限且 > 0，得到 {search_radius!r}")
    if not (np.isfinite(tolerance) and tolerance > 0.0):
        raise ValueError(f"tolerance 必须有限且 > 0，得到 {tolerance!r}")
    if not (np.isfinite(gap_mm) and gap_mm > 0.0):
        raise ValueError(f"gap_mm 必须有限且 > 0（接触搜索需严格正间隙），得到 {gap_mm!r}")
    if str(contact_type) != "sliding-elastic":
        raise ValueError(f"contact 只写 sliding-elastic 对偶面，得到 {contact_type!r}")
    fric = float(mu if fric_coeff is None else fric_coeff)
    if not (np.isfinite(fric) and fric >= 0.0):
        raise ValueError(f"fric_coeff 必须有限且 >= 0，得到 {fric!r}")

    # --- 1. 跖面 Surface ---------------------------------------------------
    # 优先 quad4（hex8 域的正确 facet 类型）；未给 plantar_quads 时走 tri3。
    surf = fmesh.Surface(name="plantar")
    if len(quads_for_surface) > 0:
        for i, quad in enumerate(quads_for_surface):
            surf.add_quad4(fmesh.Quad4Element(id=i + 1, text=",".join(map(str, _local_to_fe(quad)))))
    else:
        for i, tri in enumerate(tris):
            surf.add_tri3(fmesh.Tri3Element(id=i + 1, text=",".join(map(str, _local_to_fe(tri)))))
    model.mesh_.add_surface(surf)

    # --- 2. 固定薄板（数值刚性壁） ---------------------------------------
    n_materials = int(len(model.material_.all_materials))
    pad_meta = _build_pad_floor(
        model, nodes, plantar_local,
        gap_mm=float(gap_mm), thickness_mm=float(pad_thickness_mm),
        margin_mm=float(pad_margin_mm), ndiv=int(pad_ndiv),
        mat_id=n_materials + 1, mat_name="pad_floor_mat",
        e_mpa=float(pad_e_mpa), nu=float(pad_nu),
    )
    bcs.append(
        fbc.BCZeroDisplacement(node_set="pad_floor_fixed", x_dof=1, y_dof=1, z_dof=1)
    )

    pad_surf = fmesh.Surface(name="pad_top")
    for i, quad in enumerate(pad_meta["top_surface"]):
        pad_surf.add_quad4(fmesh.Quad4Element(id=i + 1, text=",".join(map(str, quad))))
    model.mesh_.add_surface(pad_surf)

    # FEBio sliding-elastic：primary = master（通常为**刚性/更刚**面），secondary =
    # slave（被投影的节点）。``pad_contact.py`` 的验证配方用 primary=indenter（刚），
    # secondary=pad（软）。wave-3 用 primary=plantar（软），secondary=pad_top（刚）——
    # 这是反向的；opt-in ``swap_pair=True`` 可切到 primary=pad_top / secondary=plantar。
    pair_primary, pair_secondary = (
        ("pad_top", "plantar") if swap_pair else ("plantar", "pad_top")
    )
    model.mesh_.add_surface_pair(
        fmesh.SurfacePair(
            name="plantar_pad_pair", primary=pair_primary, secondary=pair_secondary
        )
    )

    # --- 3. contact 块（挂到 step 由 builder 完成） ----------------------
    contact_obj = fcontact.Contact()
    contact_obj.add_contact(
        fcontact.SlidingElastic(
            name="plantar_pad_contact",
            surface_pair="plantar_pad_pair",
            type="sliding-elastic",
            laugon=str(laugon).upper(),
            penalty=float(penalty),
            auto_penalty=int(auto_penalty),
            two_pass=int(two_pass),
            node_reloc=int(node_reloc),
            symmetric_stiffness=int(symmetric_stiffness),
            tolerance=float(tolerance),
            search_radius=float(search_radius),
            fric_coeff=float(fric),
        )
    )

    # --- 4. 可选：弱 hold 弹簧（消掉接触未闭合前的刚体模态） -------------
    # 力学机理：跖面接触是**单边**约束（只在压入时提供刚度），静止时切向摩擦
    # 也为零 ⇒ 在 gap 未闭合前，骨只有接触支撑 ⇒ 切向刚度矩阵奇异 ⇒ Newton
    # 解出 1e23 级位移而发散。加一组**很弱**的弹簧（k ≪ 骨刚度，默认 None 不启用）
    # 可消掉刚体模态，让接触正常闭合；接触承载后弹簧力可忽略。
    hold_meta: dict | None = None
    if hold_spring_k is not None and float(hold_spring_k) > 0.0:
        n_discrete = int(len(model.discrete_.discrete_materials))
        hold_meta = _apply_spring(
            model, nodes, plantar_local, bcs, float(hold_spring_k),
            float(hold_anchor_offset_mm),
            directions=("normal", "x", "z"),
            set_name="plantar_contact_hold", mat_id=n_discrete + 1,
        )
    return {
        "hold_spring_k_n_per_mm": (
            None if hold_meta is None else float(hold_spring_k)
        ),
        "kind": "contact",
        "contact": contact_obj,
        "n_surface_tris": int(len(tris)),
        "n_surface_quads": int(len(quads_for_surface)),
        "plantar_facet_kind": ("quad4" if len(quads_for_surface) > 0 else "tri3"),
        "n_pad_nodes": int(len(pad_meta["node_ids"])),
        "n_pad_elems": int(pad_meta["n_elems"]),
        "wall_y_mm": float(pad_meta["wall_y_mm"]),
        "mu": float(mu),
        "penalty": float(penalty),
        "search_radius": float(search_radius),
        "tolerance": float(tolerance),
        "gap_mm": float(gap_mm),
        "laugon": str(laugon).upper(),
        "pair_primary": str(pair_primary),
        "pair_secondary": str(pair_secondary),
    }


def apply_plantar_bc(
    model,
    *,
    nodes: np.ndarray,
    plantar_local: np.ndarray,
    plantar_tris: np.ndarray | None,
    bcs: list,
    kind: str = "fixed",
    spring_k: float = 1000.0,
    anchor_offset_mm: float = 1.0,
    contact_mu: float = DEFAULT_CONTACT_MU,
    contact_penalty: float = DEFAULT_CONTACT_PENALTY,
    contact_search_radius: float = DEFAULT_CONTACT_SEARCH_RADIUS,
    contact_tolerance: float = DEFAULT_CONTACT_TOLERANCE,
    contact_gap_mm: float = DEFAULT_CONTACT_GAP_MM,
    contact_laugon: str = "PENALTY",
    contact_pad_e_mpa: float = 200_000.0,
    contact_plantar_quads: np.ndarray | None = None,
    contact_node_reloc: int = 0,
    contact_swap_pair: bool = False,
    contact_hold_spring_k: float | None = None,
    contact_hold_anchor_offset_mm: float = 1.0,
) -> dict:
    """把跖面 BC 施加到 ``model``，BC 对象追加进 ``bcs``；返回元数据。

    Parameters
    ----------
    model:
        已建好 mesh（节点域 / 节点集 / 关节面）的 pyfebio ``Model``。
    nodes:
        **原始形变网格**节点坐标 ``(N,3)`` mm（不含锚点）。
    plantar_local:
        跖面节点在 ``nodes`` 里的 0-based 局部索引。
    plantar_tris:
        ``(K,3)`` 跖面三角（0-based 局部索引）；``fixed``/``spring`` 可传 None，
        ``contact``/``roller`` 必需非空。
    bcs:
        BC 列表，本函数向其追加。
    kind:
        ``"fixed"`` / ``"roller"`` / ``"roller_free"`` / ``"spring"`` / ``"contact"``。
    spring_k:
        Winkler 刚度 k（N/mm），仅 ``spring`` 用。
    anchor_offset_mm:
        锚点偏移距离（mm），仅 ``spring`` 用（不改变 k，只决定弹簧初始长度）。
    contact_mu, contact_penalty, contact_search_radius, contact_tolerance, contact_gap_mm, contact_laugon:
        **仅 ``kind="contact"`` 用**（opt-in）。默认取 udg NEW PENALTY 配方
        （``μ=0.6`` / ``penalty=0.1`` / ``search_radius=20`` / ``tolerance=0.005`` /
        ``gap=0.5 mm`` / ``laugon="PENALTY"``）。见 :func:`_apply_contact`。
    contact_pad_e_mpa:
        仅 ``kind="contact"`` 用：对偶固定薄板的杨氏模量（MPa），默认 200 000
        （数值刚性壁）。
    contact_plantar_quads:
        仅 ``kind="contact"`` 用：跖面 quad4 facet ``(K,4)`` 0-based 局部节点索引
        （opt-in，``None`` = 走原 tri3 路径）。**真实 THUMS 网格的跖面属于
        CORT hex8 域**，因此正确 facet 类型是 ``quad4``（tri3 facet 挂在 hex8 域上
        FEBio 报 ``invalid facets`` ⇒ 接触不承载，详见
        ``docs/S5_contact_facet_fix.md`` §1）。由 ``thums_feb.build_thums_feb``
        在检测到 ``elements_cort`` + ``boundary_polys`` 时自动填入。

    Returns
    -------
    dict
        元数据（kind、节点数、pins / 弹簧数 / 刚度 / 接触配方等），供报告与自检。
        当 ``kind="contact"`` 时，额外含 ``"contact"`` 键（pyfebio ``Contact`` 对象）
        ——builder 必须把它挂到 step：``st.contact = meta["contact"]``（``transient
        contact`` 段在 FEBio 4.0 里属于 step）。
    """
    from pyfebio import boundary as fbc

    kind = str(kind)
    if kind not in BC_KINDS:
        raise ValueError(f"未知 plantar_bc={kind!r}；可选 {BC_KINDS}")
    local = np.asarray(plantar_local, dtype=np.int64).ravel()
    if len(local) == 0:
        raise ValueError("跖面节点集为空，无法施加跖面 BC")

    if kind == "fixed":
        # 与旧行为逐位一致：只加一个三向零位移 BC，不加任何 surface / 额外节点集。
        bcs.append(fbc.BCZeroDisplacement(node_set="plantar", x_dof=1, y_dof=1, z_dof=1))
        return {"kind": "fixed", "n_plantar_nodes": int(len(local))}

    if kind == "roller":
        meta = _apply_roller(model, nodes, local, np.asarray(plantar_tris) if plantar_tris is not None else np.empty((0, 3), dtype=np.int64), bcs)
        meta.update({"kind": "roller", "n_plantar_nodes": int(len(local))})
        return meta

    if kind == "roller_free":
        meta = _apply_roller_free(model, nodes, local, bcs, spring_k, anchor_offset_mm)
        meta.update({"kind": "roller_free", "n_plantar_nodes": int(len(local))})
        return meta

    if kind == "contact":
        meta = _apply_contact(
            model, nodes, local,
            np.asarray(plantar_tris, dtype=np.int64) if plantar_tris is not None
            else np.empty((0, 3), dtype=np.int64),
            bcs,
            mu=contact_mu, penalty=contact_penalty,
            laugon=contact_laugon, search_radius=contact_search_radius,
            tolerance=contact_tolerance, gap_mm=contact_gap_mm,
            pad_e_mpa=contact_pad_e_mpa,
            plantar_quads=contact_plantar_quads,
            node_reloc=int(contact_node_reloc),
            swap_pair=bool(contact_swap_pair),
            hold_spring_k=contact_hold_spring_k,
            hold_anchor_offset_mm=float(contact_hold_anchor_offset_mm),
        )
        meta["n_plantar_nodes"] = int(len(local))
        return meta

    # spring
    meta = _apply_spring(model, nodes, local, bcs, spring_k, anchor_offset_mm)
    meta.update({"kind": "spring", "n_plantar_nodes": int(len(local))})
    return meta
