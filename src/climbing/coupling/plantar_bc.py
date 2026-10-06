"""可切换的跖面（plantar）边界条件：``fixed`` / ``roller`` / ``spring``。

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
``spring``
    Winkler 弹性地基：每个跖面节点向**三个正交锚点**（−Y 法向 + X/Z 切向）各连
    一根线性弹簧（``k`` N/mm），锚点三向固定。这样三向刚体模态都被弹性约束住，
    无需附加刚性销；``k → ∞`` 时代数地逼近 ``fixed`` 基线（这是正确性自检）。

单位 mm–N–MPa：弹簧 ``E`` = 刚度 k（N/mm）。

设计约束
--------
* 顶层不 import pyfebio（``febio_model`` 要在无 pyfebio 环境可 import）。
* ``fixed`` 路径必须与旧 XML **逐位一致**：不额外添加 surface / node_set。
"""

from __future__ import annotations

import numpy as np

__all__ = ["BC_KINDS", "apply_plantar_bc"]

#: 支持的跖面 BC 类型（顺序即文档顺序）。
#: ``roller_free`` 是诊断用变体：全局 Y 固定 + 极软切向弹簧，真正放开切向。
BC_KINDS = ("fixed", "roller", "roller_free", "spring")

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
    n0 = int(len(nodes))
    local = np.asarray(plantar_local, dtype=np.int64)

    elems: list = []
    anchor_ids: list[int] = []
    next_id = n0 + 1
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
        ``(K,3)`` 跖面三角（0-based 局部索引）；``fixed``/``spring`` 可传 None。
    bcs:
        BC 列表，本函数向其 append。
    kind:
        ``"fixed"`` / ``"roller"`` / ``"spring"``。
    spring_k:
        Winkler 刚度 k（N/mm），仅 ``spring`` 用。
    anchor_offset_mm:
        锚点偏移距离（mm），仅 ``spring`` 用（不改变 k，只决定弹簧初始长度）。

    Returns
    -------
    dict
        元数据（kind、节点数、pins / 弹簧数 / 刚度等），供报告与自检。
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

    # spring
    meta = _apply_spring(model, nodes, local, bcs, spring_k, anchor_offset_mm)
    meta.update({"kind": "spring", "n_plantar_nodes": int(len(local))})
    return meta
