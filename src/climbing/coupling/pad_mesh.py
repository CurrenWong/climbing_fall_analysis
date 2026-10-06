"""矩形坍塌垫的体网格生成器（S5 接触阶段脚-垫体的提供者）。

职责
----
S5-contact 的 §3.2 入口：返回一个矩形 crash pad 的四面体网格 + 具名面
（``top`` = +Z 接触面，``bottom`` = −Z 固定面，``sides`` = 四个侧面合并），
供后续 ``fe_post.read_displacement``、FEBio ``<contact>`` 等流程直接消费。

坐标 / 单位
-----------
* 全部 mm；原点放在 ``(0, 0, 0)``，盒子沿 +X/+Y/+Z 张到 ``size_mm``。
* 坐标系约定：+X 长度方向、+Y 宽度方向、+Z 高度方向。
* "foot contact face" = ``top``（+Z 面）；"fixed face" = ``bottom``（−Z 面）。

为什么走"结构化六面体 → 每 cell 拆 6 个 tet"（而不是 gmsh SDK）
------------------------------------------------------------
本模块的输入是**矩形**，没有任何几何自由度需要 gmsh 去判断（没有曲面、没有
拓扑、没有尺寸场），用 gmsh 反而引入了以下三类风险：

1. SDK 路径（``gmsh.initialize(LIB)``）只在 ``D:\\Program\\gmsh-sdk\\...`` 存在时
   才能用，把模块与一台机器绑死；
2. gmsh 对一个长方体仍可能输出 6-tet 以外的非一致拆分（如 5-tet split），导致
   Jacobian 出现 0.3 以下，违反 §"裸方块 ``jacobian_frac_lt_0_3 == 0``"；
3. gmsh 的节点编号对外暴露（``gmsh.model.mesh.getNodes`` 用的是 tag），
   在 SDK 重启后顺序未必确定 → 同一网格两次跑可能节点编号漂移。

结构化拆分的优势：

* 完全确定性：节点和 tet 编号由 ``(i, j, k)`` 推导；
* 全部 tet 的 Jacobian **可解析对照**：cell 边长 ``(dx, dy, dz)`` 时，6 个 tet
  分成两类（3 个"标准"tet 棱长互垂、Jmin=1；3 个"对角线"tet 棱含两条对角线、
  Jmin = 0.5）。最坏情况 ``min_jacobian = 0.5`` → ``jacobian_frac_lt_0_3 == 0``；
* 只依赖 ``numpy``，与 ``meshing.py`` 的 gmsh/Delaunay 双路径完全解耦；
* 体积、面积、Jacobian 可解析对照公式值，测试无须任何 mesh-only 库。

权威参考
--------
* meshio 的 hex→tet 拆分（``src/meshio/_common.py`` 中的 ``hex_to_tets``）—— 本
  模块采用完全相同的 6-tet 模式。
* Jacobian / 体积定义与 ``climbing.coupling.meshing`` 一致（同一份公式），
  保证两个模块之间数值可直接比较。
"""
from __future__ import annotations

from pathlib import Path as _Path

import numpy as np

# 复用 meshing.py 的网格后处理工具，确保本模块与 bone FE 网格的
# "Jacobian/体积"定义逐位相同（§0 一致性）。
from climbing.coupling.meshing import (  # noqa: E402  (复用即合约)
    _orient_and_clean,
    _scaled_jacobian,
    _tet_volumes,
)

__all__ = ["pad_block_tet", "pad_block_hex"]


# ---------------------------------------------------------------------------
# 内部辅助
# ---------------------------------------------------------------------------
def _validate_inputs(
    nx: int,
    ny: int,
    nz: int,
    size_mm: tuple[float, float, float],
) -> tuple[int, int, float, float, float]:
    """参数基本校验 + 单位强转 mm。"""
    if not isinstance(nx, int) or nx <= 0:
        raise ValueError(f"nx 必须是正整数，得到 {nx!r}")
    if not isinstance(ny, int) or ny <= 0:
        raise ValueError(f"ny 必须是正整数，得到 {ny!r}")
    if not isinstance(nz, int) or nz <= 0:
        raise ValueError(f"nz 必须是正整数，得到 {nz!r}")
    if len(size_mm) != 3:
        raise ValueError(f"size_mm 必须是 3 元组，得到 {size_mm!r}")
    sx, sy, sz = (float(v) for v in size_mm)
    for name, v in (("sx", sx), ("sy", sy), ("sz", sz)):
        if not (np.isfinite(v) and v > 0):
            raise ValueError(f"{name} 必须是正有限数（mm），得到 {v!r}")
    return nx, ny, sx, sy, sz


def _grid_node_id(i: int, j: int, k: int, nx: int, ny: int) -> int:
    """结构化节点 (i, j, k) 的线性 ID（i ∈ [0, nx], j ∈ [0, ny], k ∈ [0, nz]）。

    行列优先顺序：x 最快、z 最慢。== ``i + (nx+1) * (j + (ny+1) * k)``。
    """
    return i + (nx + 1) * (j + (ny + 1) * k)


def _hex_to_six_tets(
    i: int, j: int, k: int, nx: int, ny: int, nz: int
) -> np.ndarray:
    """一个结构化 hex cell → 6 个 tet，编号顺序与 meshio 的 ``hex_to_tets`` 一致。

    八顶点编号（与 §"为什么..." 同图）::

        7----6
        /|   /|
      4----5 |     z
      | 3---|-2    |
      |/    |/  y  +---> y
      0----1   ↗
              x

    6-tet 拆分（3 个"标准"tet + 3 个"对角线"tet）：::

        tet1: (0, 1, 3, 4)   # 标准：bot-front-right
        tet2: (1, 3, 4, 5)   # 对角：bot-front-top
        tet3: (1, 2, 3, 6)   # 标准：bot-back-right
        tet4: (1, 3, 5, 6)   # 对角：bot-mid-top
        tet5: (3, 4, 5, 7)   # 标准：bot-top-left
        tet6: (3, 5, 6, 7)   # 对角：bot-mid-top-left

    全部 tet 在 (dx, dy, dz) > 0 时取向一致（正 Jacobian）；
    三个"标准"tet 的 scaled-Jacobian = 1.0；三个"对角线"tet = 0.5。
    """
    v0 = _grid_node_id(i,     j,     k,     nx, ny)
    v1 = _grid_node_id(i + 1, j,     k,     nx, ny)
    v2 = _grid_node_id(i + 1, j + 1, k,     nx, ny)
    v3 = _grid_node_id(i,     j + 1, k,     nx, ny)
    v4 = _grid_node_id(i,     j,     k + 1, nx, ny)
    v5 = _grid_node_id(i + 1, j,     k + 1, nx, ny)
    v6 = _grid_node_id(i + 1, j + 1, k + 1, nx, ny)
    v7 = _grid_node_id(i,     j + 1, k + 1, nx, ny)
    return np.array(
        [
            [v0, v1, v3, v4],
            [v1, v3, v4, v5],
            [v1, v2, v3, v6],
            [v1, v3, v5, v6],
            [v3, v4, v5, v7],
            [v3, v5, v6, v7],
        ],
        dtype=np.int64,
    )


def _structured_pad_mesh(
    nx: int, ny: int, nz: int, sx: float, sy: float, sz: float
) -> tuple[np.ndarray, np.ndarray]:
    """生成节点（(nx+1)*(ny+1)*(nz+1), 3）和 tet（6*nx*ny*nz, 4），全部 mm。

    节点编号约定：``_grid_node_id(i, j, k) = i + (nx+1) * (j + (ny+1) * k)``，即
    **i 最快、k 最慢**（C-order）。meshgrid 用 ``(kz, jz, iz_)`` 调入以匹配此线。
    """
    dx, dy, dz = sx / nx, sy / ny, sz / nz
    is_ = np.arange(nx + 1, dtype=np.float64)
    js = np.arange(ny + 1, dtype=np.float64)
    ks = np.arange(nz + 1, dtype=np.float64)
    # meshgrid(ks, js, is_, indexing='ij') → shape (nz+1, ny+1, nx+1)
    # C-order ravel：i 在最右 = 最快轴 → 线性索引 == _grid_node_id(i, j, k) ✓
    K, J, I = np.meshgrid(ks, js, is_, indexing="ij")
    nodes = np.stack([I.ravel() * dx, J.ravel() * dy, K.ravel() * dz], axis=1)
    assert nodes.shape == ((nx + 1) * (ny + 1) * (nz + 1), 3)

    # tet：每个 hex cell → 6 tet
    tets = np.empty((6 * nx * ny * nz, 4), dtype=np.int64)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                idx = 6 * (i + nx * (j + ny * k))
                tets[idx:idx + 6] = _hex_to_six_tets(i, j, k, nx, ny, nz)
    assert tets.shape == (6 * nx * ny * nz, 4)
    return nodes, tets


def _box_boundary_faces(
    nx: int,
    ny: int,
    nz: int,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """直接按盒子拓扑构造外壳三角面（保持原始顶点顺序），避免 ``_boundary_faces``
    丢方向的问题（它把 (a, b, c) ``sorted`` 后当 key，对方向不敏感）。

    Returns
    -------
    surfaces, normals
        ``surfaces`` ``{"top", "bottom", "sides"}`` → (K, 3) int64。
        ``normals`` (sum(K), 3) float64 外法向单位向量；**顺序与 surfaces 字典拼
        起来的 ``top → sides → bottom`` 一致**。

    顶点顺序保证：``(v1-v0) × (v2-v0)`` 给出**外法向**（CCW 从外看）。
    """
    top_tris: list[np.ndarray] = []
    bot_tris: list[np.ndarray] = []
    sides_tris: list[np.ndarray] = []
    normals_top: list[np.ndarray] = []
    normals_bot: list[np.ndarray] = []
    normals_sides: list[np.ndarray] = []

    # ----- +Z 面（k=nz-1 这层 cell 的顶面）：外法向 +Z --------------------
    # cell (i, j, k=nz-1) 的顶面四顶点：v4 (i,j,nz), v5 (i+1,j,nz), v6 (i+1,j+1,nz), v7 (i,j+1,nz)
    # 切成 2 三角：(v4, v5, v7) 与 (v5, v6, v7) — 都 CCW 从 +Z 看
    for j in range(ny):
        for i in range(nx):
            v4 = _grid_node_id(i,     j,     nz, nx, ny)
            v5 = _grid_node_id(i + 1, j,     nz, nx, ny)
            v6 = _grid_node_id(i + 1, j + 1, nz, nx, ny)
            v7 = _grid_node_id(i,     j + 1, nz, nx, ny)
            top_tris.append(np.array([v4, v5, v7], dtype=np.int64))
            top_tris.append(np.array([v5, v6, v7], dtype=np.int64))
            normals_top.append(np.array([0.0, 0.0, 1.0]))
            normals_top.append(np.array([0.0, 0.0, 1.0]))

    # ----- -Z 面（k=0 这层 cell 的底面）：外法向 -Z --------------------
    # cell (i, j, k=0) 的底面四顶点：v0 (i,j,0), v1 (i+1,j,0), v2 (i+1,j+1,0), v3 (i,j+1,0)
    # 切成 2 三角 CCW 从 -Z 看：(v0, v3, v1) 与 (v1, v3, v2)
    for j in range(ny):
        for i in range(nx):
            v0 = _grid_node_id(i,     j,     0, nx, ny)
            v1 = _grid_node_id(i + 1, j,     0, nx, ny)
            v2 = _grid_node_id(i + 1, j + 1, 0, nx, ny)
            v3 = _grid_node_id(i,     j + 1, 0, nx, ny)
            bot_tris.append(np.array([v0, v3, v1], dtype=np.int64))
            bot_tris.append(np.array([v1, v3, v2], dtype=np.int64))
            normals_bot.append(np.array([0.0, 0.0, -1.0]))
            normals_bot.append(np.array([0.0, 0.0, -1.0]))

    # ----- +X 面（i=nx-1 这列 cell 的右面）：外法向 +X --------------------
    # cell (i=nx-1, j, k) 的右面四顶点：v1, v2, v6, v5
    # CCW 从 +X 看（x 朝外，y 右，z 上）：(v1, v2, v5) 与 (v2, v6, v5)
    # （(v2-v1)×(v5-v1) = +X；原 (v1,v5,v2) 给出 -X，是内翻，已修正。）
    for k in range(nz):
        for j in range(ny):
            v1 = _grid_node_id(nx, j,     k,     nx, ny)
            v2 = _grid_node_id(nx, j + 1, k,     nx, ny)
            v5 = _grid_node_id(nx, j,     k + 1, nx, ny)
            v6 = _grid_node_id(nx, j + 1, k + 1, nx, ny)
            sides_tris.append(np.array([v1, v2, v5], dtype=np.int64))
            sides_tris.append(np.array([v2, v6, v5], dtype=np.int64))
            normals_sides.append(np.array([1.0, 0.0, 0.0]))
            normals_sides.append(np.array([1.0, 0.0, 0.0]))

    # ----- -X 面（i=0 这列 cell 的左面）：外法向 -X --------------------
    # cell (i=0, j, k) 的左面四顶点：v0, v4, v3, v7
    # CCW 从 -X 看（x 朝外，y 左，z 上）：(v0, v4, v3) 与 (v4, v7, v3)
    # （(v4-v0)×(v3-v0) = -X；原 (v0,v3,v4) 给出 +X，是内翻，已修正。）
    for k in range(nz):
        for j in range(ny):
            v0 = _grid_node_id(0, j,     k,     nx, ny)
            v3 = _grid_node_id(0, j + 1, k,     nx, ny)
            v4 = _grid_node_id(0, j,     k + 1, nx, ny)
            v7 = _grid_node_id(0, j + 1, k + 1, nx, ny)
            sides_tris.append(np.array([v0, v4, v3], dtype=np.int64))
            sides_tris.append(np.array([v4, v7, v3], dtype=np.int64))
            normals_sides.append(np.array([-1.0, 0.0, 0.0]))
            normals_sides.append(np.array([-1.0, 0.0, 0.0]))

    # ----- +Y 面（j=ny-1 这行 cell 的前/远面）：外法向 +Y ------------------
    # cell (i, j=ny-1, k) 的顶面前面四顶点：v3, v7, v2, v6
    # CCW 从 +Y 看（y 朝外，x 左，z 上）：(v3, v7, v2) 与 (v2, v7, v6)
    # （(v7-v3)×(v2-v3) = +Y；原 (v3,v2,v7) 给出 -Y，是内翻，已修正。）
    for k in range(nz):
        for i in range(nx):
            v2 = _grid_node_id(i + 1, ny, k,     nx, ny)
            v3 = _grid_node_id(i,     ny, k,     nx, ny)
            v6 = _grid_node_id(i + 1, ny, k + 1, nx, ny)
            v7 = _grid_node_id(i,     ny, k + 1, nx, ny)
            sides_tris.append(np.array([v3, v7, v2], dtype=np.int64))
            sides_tris.append(np.array([v2, v7, v6], dtype=np.int64))
            normals_sides.append(np.array([0.0, 1.0, 0.0]))
            normals_sides.append(np.array([0.0, 1.0, 0.0]))

    # ----- -Y 面（j=0 这行 cell 的近面）：外法向 -Y ------------------
    # cell (i, j=0, k) 的近面四顶点：v0, v1, v5, v4
    # CCW 从 -Y 看（y 朝外，x 右，z 上）：(v0, v1, v4) 与 (v1, v5, v4)
    # （(v1-v0)×(v4-v0) = -Y；原 (v0,v4,v1) 给出 +Y，是内翻，已修正。）
    for k in range(nz):
        for i in range(nx):
            v0 = _grid_node_id(i,     0, k,     nx, ny)
            v1 = _grid_node_id(i + 1, 0, k,     nx, ny)
            v4 = _grid_node_id(i,     0, k + 1, nx, ny)
            v5 = _grid_node_id(i + 1, 0, k + 1, nx, ny)
            sides_tris.append(np.array([v0, v1, v4], dtype=np.int64))
            sides_tris.append(np.array([v1, v5, v4], dtype=np.int64))
            normals_sides.append(np.array([0.0, -1.0, 0.0]))
            normals_sides.append(np.array([0.0, -1.0, 0.0]))

    surfaces = {
        "top": np.stack(top_tris, axis=0).astype(np.int64),
        "bottom": np.stack(bot_tris, axis=0).astype(np.int64),
        "sides": np.stack(sides_tris, axis=0).astype(np.int64),
    }
    normals = np.stack(normals_top + normals_sides + normals_bot, axis=0).astype(np.float64)
    return surfaces, normals


def _write_msh(
    result: dict,
    out_msh: _Path,
) -> None:
    """写一份 ``meshio`` ``.msh 2.2``，含具名 ``Physical Surface``：top/sides/bottom。

    复用 ``meshing.py._write_msh`` 的格式约定，但 ``field_data`` 改用
    pad 三面（3 个 surface tag）。给"人工检查 / 其它 FE 工具"复用。
    """
    import meshio  # 局部 import：meshio 不在本模块经常做严格依赖；meshing.py 同款延迟

    out_msh = _Path(out_msh)
    out_msh.parent.mkdir(parents=True, exist_ok=True)

    nodes = result["nodes"]
    tets = result["tets"]
    surfaces = result["surfaces"]
    top = surfaces["top"]
    bot = surfaces["bottom"]
    sides = surfaces["sides"]
    tris = np.vstack([top, sides, bot]).astype(np.int64) if (
        len(top) + len(bot) + len(sides)
    ) else np.empty((0, 3), dtype=np.int64)

    cells = [("tetra", tets)]
    phys = [np.ones(len(tets), dtype=np.int32)]
    geom = [np.ones(len(tets), dtype=np.int32)]
    if len(tris):
        cells.append(("triangle", tris))
        # Physical Surface tag：1 = pad 体，2 = top，4 = sides，5 = bottom
        # （注意：gmsh 的 triangle tag 必须 ≥ 2；1 留给 volume 实体）
        tag = np.concatenate(
            [
                np.full(len(top), 2, dtype=np.int32),
                np.full(len(sides), 4, dtype=np.int32),
                np.full(len(bot), 5, dtype=np.int32),
            ]
        )
        phys.append(tag)
        geom.append(tag.copy())

    mesh = meshio.Mesh(
        points=nodes,
        cells=cells,
        cell_data={"gmsh:physical": phys, "gmsh:geometrical": geom},
        field_data={
            "pad_volume": np.array([3, 1]),
            "top": np.array([2, 2]),
            "sides": np.array([2, 4]),
            "bottom": np.array([2, 5]),
        },
    )
    meshio.write(str(out_msh), mesh, file_format="gmsh22", binary=False)


# ---------------------------------------------------------------------------
# 公开入口
# ---------------------------------------------------------------------------
def pad_block_tet(
    nx: int,
    ny: int,
    nz: int,
    size_mm: tuple[float, float, float] = (300.0, 300.0, 200.0),
    *,
    out_msh: str | None = None,
) -> dict:
    """矩形 crash-pad 的结构化四面体网格（mm）。

    Parameters
    ----------
    nx, ny, nz:
        沿 X/Y/Z 方向的 cell 数（≥ 1，整数）。
    size_mm:
        ``(sx, sy, sz)``，单位 mm，盒子从原点 ``(0,0,0)`` 张到 ``(sx, sy, sz)``。
        默认 ``(300, 300, 200)`` —— S5 §2.2 验收靶的双足接触尺寸 + 200 mm 厚。
    out_msh:
        若非空，写一份 ``.msh 2.2``（具 ``Physical Surface``：top/sides/bottom）。

    Returns
    -------
    dict
        与 :func:`climbing.coupling.meshing.vtp_to_tet` 同 schema。至少含：

        * ``nodes`` (N, 3) float64, mm
        * ``tets`` (M, 4) int64
        * ``surfaces`` ``{"top", "bottom", "sides"}`` → (K, 3) int
          （三角形面节点索引，外法向 CCW）
        * ``node_sets`` 同上 3 键 → (N,) int
        * ``volume_mm3``、``surface_area_mm2``
        * ``min_jacobian``、``jacobian_frac_lt_0_3``（裸方块恒为 0.0）
        * ``method``、``n_nodes``、``n_tets``、``n_boundary_faces``
        * ``size_mm``、``cell_size_mm``、``nx``/``ny``/``nz``

        额外字段 ``reference_volume_mm3``、``volume_error_pct`` 提供解析对照。
    """
    nx, ny, sx, sy, sz = _validate_inputs(nx, ny, nz, size_mm)
    dx, dy, dz = sx / nx, sy / ny, sz / nz

    # --- 1. 节点 + tet（确定性，无外部依赖）------------------------------
    nodes, tets = _structured_pad_mesh(nx, ny, nz, sx, sy, sz)

    # --- 2. 定向 + 退化剔除（复用 meshing.py）-------------------------------------
    nodes, tets = _orient_and_clean(nodes, tets)
    n_nodes = int(nodes.shape[0])
    n_tets = int(tets.shape[0])
    # 节点数 / tet 数在裸方块上恒等于公式值；任何残留偏差都意味着
    # _orient_and_clean 内部的退化解丢弃逻辑被触发了（结构化网格不应该退化）。
    expected_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    expected_tets = 6 * nx * ny * nz
    if n_nodes != expected_nodes or n_tets != expected_tets:
        raise RuntimeError(
            f"结构化网格退化：nodes={n_nodes}/{expected_nodes}, "
            f"tets={n_tets}/{expected_tets}；可能 dx/dy/dz 为 0 或数值溢出。"
        )

    # --- 3. 边界三角面 + 外法向 -----------------------------------------------
    # 直接按盒子拓扑枚举 6 面，避免 _boundary_faces 的 sorted-key 丢方向。
    surfaces, normals = _box_boundary_faces(nx, ny, nz)
    node_sets = {name: np.unique(surfaces[name], axis=None).astype(np.int64)
                 for name in ("top", "bottom", "sides")}

    # --- 4. 体积 / 面积 / Jacobian -------------------------------------------------
    tet_vols = _tet_volumes(nodes, tets)
    volume_mm3 = float(tet_vols.sum())
    expected_volume = sx * sy * sz

    # 面积 = ∑ |(v1-v0) × (v2-v0)| / 2 over (top, sides, bottom)
    ordered_tris = np.concatenate(
        [surfaces["top"], surfaces["sides"], surfaces["bottom"]], axis=0
    )
    pts = nodes[ordered_tris]
    face_area = 0.5 * np.linalg.norm(
        np.cross(pts[:, 1] - pts[:, 0], pts[:, 2] - pts[:, 0]),
        axis=1,
    )
    surface_area_mm2 = float(face_area.sum())

    jac = _scaled_jacobian(nodes, tets)
    min_jac = float(jac.min())
    frac_lt_03 = float((jac < 0.3).mean())
    if min_jac <= 0:
        raise RuntimeError(
            f"pad 网格出现非正 Jacobian（min={min_jac:.3e}）"
            "——结构化拆分理论上不可能退化，请检查输入。"
        )

    # --- 5. 拼装结果 dict（与 meshing.vtp_to_tet 同 schema）-------------------
    result: dict = {
        "nodes": nodes.astype(np.float64),
        "tets": tets.astype(np.int64),
        "surfaces": surfaces,
        "node_sets": node_sets,
        "boundary_faces": np.concatenate(
            [surfaces["top"], surfaces["sides"], surfaces["bottom"]], axis=0
        ).astype(np.int64),
        "boundary_normals": normals.astype(np.float64),
        "volume_mm3": volume_mm3,
        "reference_volume_mm3": float(expected_volume),
        "volume_error_pct": abs(volume_mm3 - expected_volume) / expected_volume * 100.0,
        "surface_area_mm2": surface_area_mm2,
        "min_jacobian": min_jac,
        "jacobian_frac_lt_0_3": frac_lt_03,
        "n_nodes": n_nodes,
        "n_tets": n_tets,
        "n_boundary_faces": int(sum(len(surfaces[k]) for k in ("top", "bottom", "sides"))),
        "n_top_faces": int(len(surfaces["top"])),
        "n_bottom_faces": int(len(surfaces["bottom"])),
        "n_sides_faces": int(len(surfaces["sides"])),
        "method": "structured-hex-split-6",
        "size_mm": (sx, sy, sz),
        "cell_size_mm": (dx, dy, dz),
        "nx": nx,
        "ny": ny,
        "nz": nz,
    }

    # --- 6. 可选写盘 ---------------------------------------------------------
    if out_msh is not None:
        _write_msh(result, _Path(out_msh))
        result["msh_path"] = str(_Path(out_msh))

    return result


# ---------------------------------------------------------------------------
# 六面体路径（hex8）
# ---------------------------------------------------------------------------
def _structured_pad_nodes(
    nx: int, ny: int, nz: int, sx: float, sy: float, sz: float
) -> np.ndarray:
    """结构化节点坐标（与 tet 路径**完全同一套**编号与坐标）。"""
    dx, dy, dz = sx / nx, sy / ny, sz / nz
    is_ = np.arange(nx + 1, dtype=np.float64)
    js = np.arange(ny + 1, dtype=np.float64)
    ks = np.arange(nz + 1, dtype=np.float64)
    K, J, I = np.meshgrid(ks, js, is_, indexing="ij")
    return np.stack([I.ravel() * dx, J.ravel() * dy, K.ravel() * dz], axis=1)


def _hex_node_ids(i: int, j: int, k: int, nx: int, ny: int) -> np.ndarray:
    """一个 cell 的 8 个顶点，顺序 = FEBio ``hex8``。

    底面（z=k）四顶点**从 +Z 看逆时针**：``v0(i,j) v1(i+1,j) v2(i+1,j+1)
    v3(i,j+1)``；顶面（z=k+1）同序为 ``v4..v7``。与 :func:`_hex_to_six_tets`
    用的是同一套 v0..v7 编号（hex8 与 6-tet 拆分共享顶点定义）。
    """
    return np.array(
        [
            _grid_node_id(i,     j,     k,     nx, ny),
            _grid_node_id(i + 1, j,     k,     nx, ny),
            _grid_node_id(i + 1, j + 1, k,     nx, ny),
            _grid_node_id(i,     j + 1, k,     nx, ny),
            _grid_node_id(i,     j,     k + 1, nx, ny),
            _grid_node_id(i + 1, j,     k + 1, nx, ny),
            _grid_node_id(i + 1, j + 1, k + 1, nx, ny),
            _grid_node_id(i,     j + 1, k + 1, nx, ny),
        ],
        dtype=np.int64,
    )


def _box_boundary_quads(
    nx: int, ny: int, nz: int
) -> dict[str, np.ndarray]:
    """盒子外壳的 **quad4** 面（外法向 CCW）：``top``/``bottom``/``sides``。

    每个面的顶点顺序使 ``(p1-p0) × (p2-p0)`` = 外法向（与 tet 路径同约定）。
    """
    top: list[list[int]] = []
    bot: list[list[int]] = []
    sides: list[list[int]] = []

    for j in range(ny):
        for i in range(nx):
            top.append([
                _grid_node_id(i,     j,     nz, nx, ny),
                _grid_node_id(i + 1, j,     nz, nx, ny),
                _grid_node_id(i + 1, j + 1, nz, nx, ny),
                _grid_node_id(i,     j + 1, nz, nx, ny),
            ])                                            # +Z 外法向
    for j in range(ny):
        for i in range(nx):
            bot.append([
                _grid_node_id(i,     j,     0, nx, ny),
                _grid_node_id(i,     j + 1, 0, nx, ny),
                _grid_node_id(i + 1, j + 1, 0, nx, ny),
                _grid_node_id(i + 1, j,     0, nx, ny),
            ])                                            # -Z 外法向
    for k in range(nz):
        for j in range(ny):                               # +X
            sides.append([
                _grid_node_id(nx, j,     k,     nx, ny),
                _grid_node_id(nx, j + 1, k,     nx, ny),
                _grid_node_id(nx, j + 1, k + 1, nx, ny),
                _grid_node_id(nx, j,     k + 1, nx, ny),
            ])
        for j in range(ny):                               # -X
            sides.append([
                _grid_node_id(0, j,     k,     nx, ny),
                _grid_node_id(0, j,     k + 1, nx, ny),
                _grid_node_id(0, j + 1, k + 1, nx, ny),
                _grid_node_id(0, j + 1, k,     nx, ny),
            ])
        for i in range(nx):                               # +Y
            sides.append([
                _grid_node_id(i,     ny, k,     nx, ny),
                _grid_node_id(i,     ny, k + 1, nx, ny),
                _grid_node_id(i + 1, ny, k + 1, nx, ny),
                _grid_node_id(i + 1, ny, k,     nx, ny),
            ])
        for i in range(nx):                               # -Y
            sides.append([
                _grid_node_id(i,     0, k,     nx, ny),
                _grid_node_id(i + 1, 0, k,     nx, ny),
                _grid_node_id(i + 1, 0, k + 1, nx, ny),
                _grid_node_id(i,     0, k + 1, nx, ny),
            ])

    return {
        "top": np.asarray(top, dtype=np.int64),
        "bottom": np.asarray(bot, dtype=np.int64),
        "sides": np.asarray(sides, dtype=np.int64),
    }


def pad_block_hex(
    nx: int,
    ny: int,
    nz: int,
    size_mm: tuple[float, float, float] = (300.0, 300.0, 200.0),
    *,
    out_msh: str | None = None,
) -> dict:
    """矩形 crash-pad 的结构化 **六面体（hex8）** 网格（mm）。

    与 :func:`pad_block_tet` **同 schema**，差别只有：

    * ``hexes`` (M, 8) 取代 ``tets``（``M = nx·ny·nz``，不再 ×6）；
    * ``surfaces`` / ``boundary_faces`` 是 **quad4** (K, 4) 而非 tri3；
    * ``n_tets`` 为 0，另给 ``n_hexes``、``method="structured-hex8"``。

    ``min_jacobian`` / ``jacobian_frac_lt_0_3`` 用**同一 cell 的 6-tet 分解**
    计算（目的是与 tet 路径的 Jacobian 定义逐位可比，而不是评价 hex8 本身
    的形状质量）—— 对规则长方体该值恒定。

    动机：``tet4`` 在名义应变 ``ξ ≳ 0.25`` 时反演（
    ``docs/S5_contact_convergence.md``），``hex8`` 在大应变下鲁棒得多。
    """
    nx, ny, sx, sy, sz = _validate_inputs(nx, ny, nz, size_mm)
    dx, dy, dz = sx / nx, sy / ny, sz / nz

    nodes = _structured_pad_nodes(nx, ny, nz, sx, sy, sz)
    n_nodes = int(nodes.shape[0])
    expected_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    if n_nodes != expected_nodes:
        raise RuntimeError(f"结构化 hex 网格节点数异常：{n_nodes}/{expected_nodes}")

    hexes = np.empty((nx * ny * nz, 8), dtype=np.int64)
    tets = np.empty((6 * nx * ny * nz, 4), dtype=np.int64)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = i + nx * (j + ny * k)
                hexes[c] = _hex_node_ids(i, j, k, nx, ny)
                tets[6 * c:6 * c + 6] = _hex_to_six_tets(i, j, k, nx, ny, nz)

    surfaces = _box_boundary_quads(nx, ny, nz)
    node_sets = {name: np.unique(surfaces[name], axis=None).astype(np.int64)
                 for name in ("top", "bottom", "sides")}

    tet_vols = _tet_volumes(nodes, tets)
    volume_mm3 = float(tet_vols.sum())
    expected_volume = sx * sy * sz

    pts = nodes[np.concatenate([surfaces["top"], surfaces["sides"], surfaces["bottom"]], axis=0)]
    face_area = 0.5 * (
        np.linalg.norm(np.cross(pts[:, 1] - pts[:, 0], pts[:, 2] - pts[:, 0]), axis=1)
        + np.linalg.norm(np.cross(pts[:, 2] - pts[:, 0], pts[:, 3] - pts[:, 0]), axis=1)
    )
    surface_area_mm2 = float(face_area.sum())

    jac = _scaled_jacobian(nodes, tets)
    min_jac = float(jac.min())
    if min_jac <= 0:
        raise RuntimeError(f"pad hex 网格 6-tet 分解出现非正 Jacobian（{min_jac:.3e}）")

    result: dict = {
        "nodes": nodes.astype(np.float64),
        "hexes": hexes.astype(np.int64),
        "tets": np.empty((0, 4), dtype=np.int64),
        "surfaces": surfaces,
        "node_sets": node_sets,
        "boundary_faces": np.concatenate(
            [surfaces["top"], surfaces["sides"], surfaces["bottom"]], axis=0).astype(np.int64),
        "volume_mm3": volume_mm3,
        "reference_volume_mm3": float(expected_volume),
        "volume_error_pct": abs(volume_mm3 - expected_volume) / expected_volume * 100.0,
        "surface_area_mm2": surface_area_mm2,
        "min_jacobian": min_jac,
        "jacobian_frac_lt_0_3": float((jac < 0.3).mean()),
        "n_nodes": n_nodes,
        "n_hexes": int(hexes.shape[0]),
        "n_tets": 0,
        "n_boundary_faces": int(sum(len(surfaces[k]) for k in ("top", "bottom", "sides"))),
        "n_top_faces": int(len(surfaces["top"])),
        "n_bottom_faces": int(len(surfaces["bottom"])),
        "n_sides_faces": int(len(surfaces["sides"])),
        "method": "structured-hex8",
        "size_mm": (sx, sy, sz),
        "cell_size_mm": (dx, dy, dz),
        "nx": nx,
        "ny": ny,
        "nz": nz,
    }
    if out_msh is not None:
        import meshio

        out_msh = _Path(out_msh)
        out_msh.parent.mkdir(parents=True, exist_ok=True)
        top, bot, sides = surfaces["top"], surfaces["bottom"], surfaces["sides"]
        quads = np.vstack([top, sides, bot]).astype(np.int64)
        meshio.write(
            str(out_msh),
            meshio.Mesh(
                points=nodes,
                cells=[("hexahedron", hexes), ("quad", quads)],
                cell_data={"gmsh:physical": [
                    np.ones(len(hexes), dtype=np.int32),
                    np.concatenate([np.full(len(top), 2, np.int32),
                                    np.full(len(sides), 4, np.int32),
                                    np.full(len(bot), 5, np.int32)]),
                ]},
            ),
            file_format="gmsh22", binary=False)
        result["msh_path"] = str(out_msh)
    return result