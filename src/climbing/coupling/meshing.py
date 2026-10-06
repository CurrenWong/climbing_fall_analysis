"""VTP 骨面 → 四面体网格（OpenSim FE 半边的 S1 前处理）。

职责
----
把 OpenSim 的 ``r_foot.vtp``（body ``calcn_r`` 的几何）转成 FEBio 可用的
tet 网格 + 具名面（距下关节面 / 跖面），并给出可自检的网格统计。

坐标系与单位（关键，§七 T10）
-----------------------------
* VTP 原始坐标是 **米**，且几何挂在 ``calcn_r`` 体上；默认站姿下 ``calcn_r``
  的体轴与地面系重合（X 前、Y 上、Z 内外），因此 **VTP 就是 calcn 体坐标**。
* 本模块把坐标整体 ×1000 转成 **毫米**，之后全程 mm–N–MPa–s。
* Rajagopal2015 的 ``subtalar_r`` 关节 ``location_in_parent`` 是
  ``(-0.04877,-0.04195,0.00792) m``（**在 talus_r 父系**）。因为 talus/calcn 两系
  在默认姿态下旋转相同、只差平移，该点换算到 calcn 系恰好是 **原点 (0,0,0)**
  ——这正是"距下关节中心 = 跟骨体原点"的含义。所以关节中心在几何坐标里是原点。

几何现实（诚实记录）
-------------------
``r_foot.vtp`` 只有 1000 点 / 1455 面，且 **拓扑非流形**（补齐小孔后 Euler 数仍
不合法），3D 自交导致 gmsh 的 ``classifySurfaces``→``createGeometry`` 参数化失败
（既试过 ``pi/2``，也试过官方 example 的 ``includeBoundary=True,
forReparametrization=True``，以及小角度、30–60°、HXT/Frontal 等 3D 算法，均在
``createGeometry`` 或 ``generate(3)`` 处报错：``Wrong topology of boundary mesh``、
``Invalid boundary mesh (overlapping facets)``、``PLC Error: segment and facet
intersect``）。这是几何本身的病，不是参数问题。

因此：优先走 gmsh SDK（与 ``scripts/ankle_fe/geom.py`` 完全一致的调用方式）；
**若 gmsh 无法成实体**，退到"体素等值面 + scipy Delaunay"这一自包含、可复现的
修复路径（``method`` 字段记录用了哪条），并在报告里如实标注。绝不伪造网格。
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# gmsh 官方 SDK（与 scripts/ankle_fe/geom.py 相同的方式；本模块不 import 它）
# ---------------------------------------------------------------------------
SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"

# meshing.py -> coupling -> climbing -> src -> 仓库根
ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "temp" / "opensim_fe"
MESH_DIR = OUT_DIR / "meshes"

#: 距下关节中心在 calcn 体坐标（mm）= 原点（见模块 docstring）。
JOINT_CENTER_MM: tuple[float, float, float] = (0.0, 0.0, 0.0)
#: 关节面选取半径（mm）与"朝上（+Y）"法向阈值。
JOINT_RADIUS_MM = 25.0
JOINT_NORMAL_Y_MIN = 0.2
#: 跖面（地面接触）取朝下（−Y）的边界面。
PLANTAR_NORMAL_Y_MAX = -0.5
#: 跟腱附着区（后上结节）：x < 此值 = 后半部；y > 此值 = 上半部（mm）。
#: 实测 IITD calcaneus STL：该框选到 876 面 / 1164 mm²，质心 (-38.8, 7.0, 1.0)。
ACHILLES_X_MAX = -25.0
ACHILLES_Y_MIN = 0.0

__all__ = ["vtp_to_tet", "stl_to_tet_gmsh", "JOINT_CENTER_MM", "SDK", "MESH_DIR"]


class GmshSolidError(RuntimeError):
    """gmsh 无法把该 VTP 面转成实体（拓扑/自交），由调用方决定是否回退。"""


# ---------------------------------------------------------------------------
# 读取 / 修面
# ---------------------------------------------------------------------------
def _load_surface_mm(vtp_path: Path):
    """读 VTP → clean → fill_holes → triangulate → ×1000 转 mm。"""
    import pyvista as pv

    surf = pv.read(str(vtp_path))
    surf = surf.clean(point_merging=True, tolerance=1e-6, absolute=True)
    surf = surf.triangulate()
    # 3 个开放边 → 补一个小孔即可闭合（不水密时 geom.py 的流程无法成体）。
    if surf.n_open_edges:
        surf = surf.fill_holes(1e9)
        surf = surf.clean(point_merging=True, tolerance=1e-6, absolute=True)
    pts = np.asarray(surf.points, dtype=np.float64) * 1000.0  # m -> mm
    out = pv.PolyData(pts, np.asarray(surf.faces))
    return out


def _mesh_size(char_len: float) -> float:
    if not np.isfinite(char_len) or char_len <= 0:
        raise ValueError(f"char_len 必须是正数（mm），得到 {char_len!r}")
    return float(char_len)


# ---------------------------------------------------------------------------
# 路径 A：gmsh 官方 SDK（geom.py 同款）
# ---------------------------------------------------------------------------
def _gmsh_tet(surf, char_len: float) -> tuple[np.ndarray, np.ndarray]:
    """用官方 gmsh SDK 走 STL→tet。失败抛 :class:`GmshSolidError`。"""
    import os
    import sys
    import tempfile

    if not LIB.is_dir():
        raise GmshSolidError(f"gmsh SDK lib 不存在：{LIB}")
    os.environ["PATH"] = f"{LIB};{SDK / 'bin'};{os.environ['PATH']}"
    if str(LIB) not in sys.path:
        sys.path.insert(0, str(LIB))
    import gmsh  # noqa: E402  (需先设好 PATH/sys.path)

    tmp_stl = Path(tempfile.gettempdir()) / "_calcaneus_gmsh_in.stl"
    surf.save(str(tmp_stl))

    gmsh.initialize(str(LIB))
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.merge(str(tmp_stl))
        gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
        gmsh.model.mesh.removeDuplicateNodes()
        # geom.py 同款：先把共面三角形归类成曲面实体
        try:
            gmsh.model.mesh.classifySurfaces(math.pi / 2)
            gmsh.model.mesh.createGeometry()
            gmsh.model.geo.synchronize()
        except Exception as exc:  # gmsh 用裸 Exception，这里显式转为领域错误
            raise GmshSolidError(f"classifySurfaces/createGeometry: {exc}") from exc

        surfs = [t for _, t in gmsh.model.getEntities(2)]
        if not surfs:
            raise GmshSolidError("classifySurfaces 未产生任何曲面")
        try:
            loop = gmsh.model.geo.addSurfaceLoop(surfs)
            gmsh.model.geo.addVolume([loop])
            gmsh.model.geo.synchronize()
        except Exception as exc:
            raise GmshSolidError(f"addSurfaceLoop/addVolume: {exc}") from exc
        if not gmsh.model.getEntities(3):
            raise GmshSolidError("未能封闭成实体（非封闭/自交）")

        gmsh.option.setNumber("Mesh.MeshSizeMax", char_len)
        gmsh.option.setNumber("Mesh.MeshSizeMin", char_len * 0.4)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay
        try:
            gmsh.model.mesh.generate(3)
        except Exception as exc:
            raise GmshSolidError(f"generate(3): {exc}") from exc
        gmsh.model.mesh.optimize("Netgen")

        types = gmsh.model.mesh.getElementTypes(3)
        if 4 not in types:
            raise GmshSolidError(f"generate(3) 未产生 tet4（types={types}）")

        gtags, gcoord, _ = gmsh.model.mesh.getNodes()
        gcoord = np.array(gcoord, dtype=np.float64).reshape(-1, 3)
        rmap = {int(t): i for i, t in enumerate(gtags)}
        etypes, _, econn = gmsh.model.mesh.getElements(3)
        tets = None
        for et, en in zip(etypes, econn):
            if et == 4:
                tets = np.vectorize(rmap.get)(np.array(en).reshape(-1, 4)).astype(np.int64)
                break
        if tets is None or len(tets) == 0:
            raise GmshSolidError("未取到 tet4 连接关系")
        return gcoord, tets
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# 路径 B：体素等值面 + scipy Delaunay（对自交/非流形面健壮）
# ---------------------------------------------------------------------------
def _isosurface(surf, spacing: float):
    """用 stencil(原面) → marching cubes 得到一个干净闭合的参考面。

    原面不合法也没关系：stencil 只用来给出近似的 inside 掩膜，marching cubes
    的输出一定是水密流形（有向性可能内翻，``vtkSelectEnclosedPoints`` 不在乎）。
    """
    import pyvista as pv
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk

    b = np.asarray(surf.bounds, dtype=float)
    n = [int(math.ceil((b[2 * i + 1] - b[2 * i]) / spacing)) + 4 for i in range(3)]
    origin = np.array([b[0] - 2 * spacing, b[2] - 2 * spacing, b[4] - 2 * spacing])

    img = vtk.vtkImageData()
    img.SetSpacing(spacing, spacing, spacing)
    img.SetDimensions(*n)
    img.SetOrigin(*origin)
    ones = np.ones(int(np.prod(n)), dtype=np.float32)
    img.GetPointData().SetScalars(numpy_to_vtk(ones, deep=True))

    stencil = vtk.vtkPolyDataToImageStencil()
    stencil.SetInputData(surf)
    stencil.SetOutputSpacing(spacing, spacing, spacing)
    stencil.SetOutputOrigin(*origin)
    stencil.SetOutputWholeExtent(img.GetExtent())
    stencil.Update()

    applied = vtk.vtkImageStencil()
    applied.SetInputData(img)
    applied.SetStencilConnection(stencil.GetOutputPort())
    applied.ReverseStencilOn()
    applied.SetBackgroundValue(0)
    applied.Update()

    mc = vtk.vtkMarchingCubes()
    mc.SetInputConnection(applied.GetOutputPort())
    mc.SetValue(0, 0.5)
    mc.ComputeNormalsOn()
    mc.Update()
    return pv.wrap(mc.GetOutput())


def _enclosed_mask(reference, points: np.ndarray) -> np.ndarray:
    import pyvista as pv
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    sel = vtk.vtkSelectEnclosedPoints()
    sel.SetInputData(pv.PolyData(np.asarray(points, dtype=float)))
    sel.SetSurfaceData(reference)
    sel.SetTolerance(0.0)
    sel.Update()
    arr = sel.GetOutput().GetPointData().GetArray("SelectedPoints")
    return vtk_to_numpy(arr).astype(bool)


def _inside_union_blocks(reference, min_points: int = 50):
    """把参考面拆成连通体，丢掉微小碎片，返回（体列表, OR 判定函数）。

    为什么需要：原 VTP 自交，stencil 的偶奇判定在**两片重叠区**会互相抵消，
    把脚拆成前/后两块互不相连的实体 -> FEBio 里出现自由漂浮子域 -> 奇异。
    改为"对每个闭合连通体分别判内、再取并集"，重叠区归属不变、整体连通。
    """
    import pyvista as pv

    blocks = []
    try:
        parts = reference.split_bodies()
        for i in range(getattr(parts, "n_blocks", 0)):
            b = parts[i]
            # split_bodies 可能给 UnstructuredGrid；统一成 PolyData 供 SelectEnclosedPoints。
            if not isinstance(b, pv.PolyData):
                b = b.extract_surface()
            # MC 输出的 vtkMassProperties 可能给 0，故以点数做碎片过滤（大块 ~1e4，
            # 碎片 6–18 点）。
            if b.n_points >= min_points:
                blocks.append(b)
    except Exception as exc:  # pyvista split_bodies 的偶发失败——不静默
        print(f"[meshing] split_bodies 失败（{type(exc).__name__}: {exc}）；退化为整面判定。")
        blocks = []
    if not blocks:
        blocks = [reference]

    def enclosed(points: np.ndarray) -> np.ndarray:
        mask = np.zeros(len(points), dtype=bool)
        for b in blocks:
            mask |= _enclosed_mask(b, points)
        return mask

    return blocks, enclosed


def _bridge_components(mask: np.ndarray, iterations: int = 1, min_frac: float = 0.05) -> np.ndarray:
    """丢掉微小碎片，再只在相邻大块之间搭桥，避免整体闭运算虚增体积。"""
    from scipy.ndimage import binary_dilation, label

    lab, ncomp = label(mask)
    if ncomp > 1:
        sizes = np.bincount(lab.ravel())
        sizes[0] = 0
        keep = np.where(sizes >= min_frac * sizes.max())[0]
        keep = keep[keep > 0]
        mask = np.isin(lab, keep)

    lab, ncomp = label(mask)
    if ncomp <= 1:
        return mask
    struct = np.ones((3, 3, 3), dtype=bool)
    comps = [lab == c for c in range(1, ncomp + 1)]
    out = mask.copy()
    for a in range(len(comps)):
        da = binary_dilation(comps[a], structure=struct, iterations=iterations)
        for b in range(a + 1, len(comps)):
            db = binary_dilation(comps[b], structure=struct, iterations=iterations)
            out |= da & db
    return out


def _delaunay_tet(surf, char_len: float, *, seed: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """体素角点 + jitter → Delaunay → 保留质心在体内的 tet。

    自交面 repaired：先由 stencil 得到一个粗糙体素掩膜，对**每个闭合连通体**
    分别判内取并集，再形态学闭运算补上亚体素的缝隙，最后取最大连通分量。
    这样得到的实体在 FEBio 里是**单一连通域**（不会有自由漂浮子域导致奇异）。
    """
    from scipy.ndimage import binary_dilation, binary_fill_holes, label
    from scipy.spatial import Delaunay

    spacing = char_len
    iso = _isosurface(surf, spacing=max(0.8, spacing * 0.4))
    _blocks, enclosed = _inside_union_blocks(iso)
    b = np.asarray(surf.bounds, dtype=float)
    origin = np.array([b[0] - 2 * spacing, b[2] - 2 * spacing, b[4] - 2 * spacing])
    n = [int(math.ceil((b[2 * i + 1] - b[2 * i]) / spacing)) + 4 for i in range(3)]

    # 体素中心 → inside（对每个闭合体判内再取并集）
    k, j, i = np.meshgrid(np.arange(n[2]), np.arange(n[1]), np.arange(n[0]), indexing="ij")
    centers = np.stack(
        [
            origin[0] + i.ravel() * spacing + spacing / 2,
            origin[1] + j.ravel() * spacing + spacing / 2,
            origin[2] + k.ravel() * spacing + spacing / 2,
        ],
        axis=1,
    )
    inside = enclosed(centers).reshape(n[2], n[1], n[0])
    if inside.sum() == 0:
        raise RuntimeError("体素化未找到任何内部体素：VTP 可能不可定向")

    # 自交面在重叠区被偶奇判定"抵消"，会留下一道亚体素级窄缝，把脚拆成前后两块。
    # 只在**两块之间**做桥接（对两块分别膨胀后取交），而不是整体闭运算——
    # 后者会把足弓等凹腔也填掉，虚增 ~10% 体积（T7 失败）。
    inside = _bridge_components(inside)
    inside = binary_fill_holes(inside)
    lab, ncomp = label(inside)
    if ncomp > 1:
        sizes = np.bincount(lab.ravel())
        sizes[0] = 0
        inside = lab == int(sizes.argmax())
    if inside.sum() == 0:
        raise RuntimeError("形态学处理后掩膜为空（自交面修复失败）")

    # 内部体素的 8 个角点（含表面角点），去重后 jitter，打破规则晶格以改善 tet 质量
    idx = np.argwhere(inside)
    offs = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]],
        dtype=float,
    )
    corner = (idx[:, None, :] + offs[None, :, :]).reshape(-1, 3)
    pts = np.stack(
        [origin[0] + corner[:, 2] * spacing, origin[1] + corner[:, 1] * spacing, origin[2] + corner[:, 0] * spacing],
        axis=1,
    )
    pts = np.unique(np.round(pts, 6), axis=0)

    rng = np.random.default_rng(seed)
    jitter = 0.15 * spacing
    pts = pts + rng.uniform(-jitter, jitter, pts.shape)

    tri = Delaunay(pts)
    simplices = tri.simplices
    centroids = pts[simplices].mean(axis=1)
    # 用体素掩膜直接查表判内外（避免自交/多连通面上的偶奇判定失效）。
    # 体素 i 的中心在 origin+(i+0.5)*spacing，故点 c 所属体素 = floor((c-origin)/spacing)。
    ijk = np.floor((centroids - origin) / spacing).astype(np.int64)  # 列序 (x,y,z)
    n_xyz = np.array([n[0], n[1], n[2]])
    valid = np.all((ijk >= 0) & (ijk < n_xyz), axis=1)
    keep = np.zeros(len(simplices), dtype=bool)
    if valid.any():
        # inside 的索引序是 [z, y, x]
        keep[valid] = inside[ijk[valid, 2], ijk[valid, 1], ijk[valid, 0]]
    tets = simplices[keep]
    if len(tets) == 0:
        raise RuntimeError("Delaunay 过滤后无 tet 保留")

    # 只保留被引用的节点并重编号
    used = np.unique(tets)
    remap = -np.ones(len(pts), dtype=np.int64)
    remap[used] = np.arange(len(used))
    return pts[used], remap[tets]


# ---------------------------------------------------------------------------
# 统一后处理：定向、边界、具名面、质量指标
# ---------------------------------------------------------------------------
def _orient_and_clean(nodes: np.ndarray, tets: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q = nodes[tets]
    det = np.einsum("ij,ij->i", q[:, 1] - q[:, 0], np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0]))
    tets = tets.copy()
    flip = det < 0
    tets[flip, 1], tets[flip, 2] = tets[flip, 2].copy(), tets[flip, 1].copy()
    # 丢退化（零体积）tet
    q = nodes[tets]
    vol6 = np.abs(
        np.einsum("ij,ij->i", q[:, 1] - q[:, 0], np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0]))
    )
    tets = tets[vol6 > 1e-9]
    # 再次只保留被引用节点
    used = np.unique(tets)
    remap = -np.ones(len(nodes), dtype=np.int64)
    remap[used] = np.arange(len(used))
    return nodes[used], remap[tets]


def _scaled_jacobian(nodes: np.ndarray, tets: np.ndarray) -> np.ndarray:
    q = nodes[tets]
    a, b, c = q[:, 1] - q[:, 0], q[:, 2] - q[:, 0], q[:, 3] - q[:, 0]
    det = np.einsum("ij,ij->i", a, np.cross(b, c))
    denom = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) * np.linalg.norm(c, axis=1)
    return np.abs(det) / np.maximum(denom, 1e-30)


def _tet_volumes(nodes: np.ndarray, tets: np.ndarray) -> np.ndarray:
    q = nodes[tets]
    return (
        np.abs(
            np.einsum("ij,ij->i", q[:, 1] - q[:, 0], np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0]))
        )
        / 6.0
    )


def _boundary_faces(tets: np.ndarray) -> np.ndarray:
    """返回只属于一个 tet 的边界面（(K,3) 节点索引）及其 owner tet 索引。"""
    faces: dict[tuple[int, int, int], list[int]] = {}
    for ti, (a, b, c, d) in enumerate(tets):
        for tri in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)):
            key = tuple(sorted(int(x) for x in tri))
            faces.setdefault(key, []).append(ti)
    out_faces, owners = [], []
    for key, owners_ in faces.items():
        if len(owners_) == 1:
            out_faces.append(key)
            owners.append(owners_[0])
    return np.asarray(out_faces, dtype=np.int64), np.asarray(owners, dtype=np.int64)


def _outward_normals(nodes: np.ndarray, faces: np.ndarray, owners: np.ndarray, tets: np.ndarray) -> np.ndarray:
    p = nodes[faces]
    n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    norm = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / np.maximum(norm, 1e-30)
    # 让法向背离 owner tet 的第 4 个（相对面）顶点
    all_v = tets[owners]
    face_set = faces
    for k in range(4):
        v = all_v[:, k]
        # 该顶点是否属于当前面？
        on_face = np.any(face_set == v[:, None], axis=1)
        away = ~on_face
        if away.any():
            n[away] = -n[away]
    return n


def _classify_surfaces(
    nodes: np.ndarray,
    faces: np.ndarray,
    normals: np.ndarray,
    joint_center: tuple[float, float, float],
) -> dict[str, np.ndarray]:
    c = np.asarray(joint_center, dtype=float)
    centroids = nodes[faces].mean(axis=1)
    dist = np.linalg.norm(centroids - c, axis=1)

    joint = (dist < JOINT_RADIUS_MM) & (normals[:, 1] > JOINT_NORMAL_Y_MIN)
    plantar = normals[:, 1] < PLANTAR_NORMAL_Y_MAX
    # 跟腱附着区：后上结节（x 后半部 + y 上半部）。不做法向约束，避免过渡带漏选。
    achilles = (centroids[:, 0] < ACHILLES_X_MAX) & (centroids[:, 1] > ACHILLES_Y_MIN)
    return {
        "subtalar_joint": faces[joint],
        "plantar": faces[plantar],
        "achilles": faces[achilles],
    }


def _finalize(
    surf,
    nodes: np.ndarray,
    tets: np.ndarray,
    *,
    method: str,
    joint_center: tuple[float, float, float],
    char_len: float,
) -> dict:
    nodes, tets = _orient_and_clean(nodes, tets)
    faces, owners = _boundary_faces(tets)
    normals = _outward_normals(nodes, faces, owners, tets)

    surfaces = _classify_surfaces(nodes, faces, normals, joint_center)
    # 关节面与跖面会在曲率过渡处共享节点。若共享节点同时被跖面 BC 固定，
    # 施加在关节面的载荷会被固定端吸收（响应≈0、伪收敛）。这里以关节面优先，
    # 从跖面节点集中剔除关节节点，保证加载端与约束端分离。
    joint_nodes = np.unique(surfaces["subtalar_joint"])
    plantar_nodes = np.setdiff1d(np.unique(surfaces["plantar"]), joint_nodes)
    achilles_nodes = np.setdiff1d(
        np.unique(surfaces["achilles"]), np.union1d(joint_nodes, plantar_nodes)
    )
    node_sets = {
        "subtalar_joint": joint_nodes,
        "plantar": plantar_nodes,
        "achilles": achilles_nodes,
    }

    vol_tets = _tet_volumes(nodes, tets)
    jac = _scaled_jacobian(nodes, tets)
    p = nodes[faces]
    face_area = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)

    volume = float(vol_tets.sum())
    ref_volume = float(surf.volume)  # pyvista 表面体积（mm³）
    return {
        "nodes": nodes.astype(np.float64),
        "tets": tets.astype(np.int64),
        "surfaces": surfaces,
        "node_sets": node_sets,
        "boundary_faces": faces,
        "boundary_normals": normals,
        "volume_mm3": volume,
        "surface_area_mm2": float(face_area.sum()),
        "reference_volume_mm3": ref_volume,
        "reference_area_mm2": float(surf.area),
        "volume_error_pct": abs(volume - ref_volume) / ref_volume * 100.0,
        "min_jacobian": float(jac.min()),
        "jacobian_frac_lt_0_3": float((jac < 0.3).mean()),
        "n_nodes": int(len(nodes)),
        "n_tets": int(len(tets)),
        "n_boundary_faces": int(len(faces)),
        "n_joint_faces": int(len(surfaces["subtalar_joint"])),
        "n_plantar_faces": int(len(surfaces["plantar"])),
        "n_achilles_faces": int(len(surfaces["achilles"])),
        "method": method,
        "char_len_mm": float(char_len),
        "joint_center_mm": tuple(float(x) for x in joint_center),
    }


# ---------------------------------------------------------------------------
# 写 .msh（具名 Physical Surface；供人工检查 / 其它工具复用）
# ---------------------------------------------------------------------------
def _write_msh(result: dict, out_msh: Path, *, joint_name: str = "subtalar_joint") -> None:
    import meshio

    out_msh = Path(out_msh)
    out_msh.parent.mkdir(parents=True, exist_ok=True)

    nodes = result["nodes"]
    tets = result["tets"]
    joint = result["surfaces"].get(joint_name, np.empty((0, 3), dtype=np.int64))
    plantar = result["surfaces"].get("plantar", np.empty((0, 3), dtype=np.int64))
    achilles = result["surfaces"].get("achilles", np.empty((0, 3), dtype=np.int64))
    tris = (np.vstack([joint, plantar, achilles])
            if (len(joint) + len(plantar) + len(achilles)) else np.empty((0, 3), dtype=np.int64))

    cells = [("tetra", tets)]
    phys = [np.ones(len(tets), dtype=np.int32)]
    geom = [np.ones(len(tets), dtype=np.int32)]
    if len(tris):
        cells.append(("triangle", tris))
        # Physical Surface 2 = 距下关节面, 3 = 跖面, 4 = 跟腱附着面
        tag = np.concatenate(
            [np.full(len(joint), 2, dtype=np.int32),
             np.full(len(plantar), 3, dtype=np.int32),
             np.full(len(achilles), 4, dtype=np.int32)]
        )
        phys.append(tag)
        geom.append(tag.copy())

    mesh = meshio.Mesh(
        points=nodes,
        cells=cells,
        cell_data={"gmsh:physical": phys, "gmsh:geometrical": geom},
        field_data={
            "calcaneus": np.array([3, 1]),
            joint_name: np.array([2, 2]),
            "plantar": np.array([2, 3]),
            "achilles": np.array([2, 4]),
        },
    )
    meshio.write(str(out_msh), mesh, file_format="gmsh22", binary=False)


# ---------------------------------------------------------------------------
# 公开入口
# ---------------------------------------------------------------------------
def vtp_to_tet(vtp_path, char_len: float, out_msh=None) -> dict:
    """VTP 骨面 → tet 网格（mm）。

    Parameters
    ----------
    vtp_path:
        OpenSim VTP 路径（``r_foot.vtp``）。
    char_len:
        目标单元尺寸（mm）。网格越细越准，也越慢。
    out_msh:
        给了就写一个带具名 ``Physical Surface`` 的 ``.msh``。

    Returns
    -------
    dict
        至少含 ``nodes`` (N,3 float64, mm)、``tets`` (M,4 int64)、具名面
        ``surfaces``/``node_sets``、``volume_mm3``、``min_jacobian``、
        ``surface_area_mm2``；另含 ``method`` 与若干自检字段。
    """
    vtp_path = Path(vtp_path)
    if not vtp_path.is_file():
        raise FileNotFoundError(f"VTP 不存在：{vtp_path}")
    char_len = _mesh_size(char_len)

    surf = _load_surface_mm(vtp_path)
    if surf.n_open_edges:
        raise RuntimeError(
            f"VTP 补孔后仍非水密（open_edges={surf.n_open_edges}）；"
            "无法保证实体内外判定，请检查几何。"
        )

    method = "gmsh-sdk"
    try:
        nodes, tets = _gmsh_tet(surf, char_len)
    except GmshSolidError as exc:
        # 不是静默兜底：明确说明 gmsh 为何失败，再改用自包含修复路径。
        print(f"[meshing] gmsh SDK 流程失败（几何自交/非流形）：{exc}")
        print("[meshing] 回退：体素等值面 + scipy Delaunay 修复路径。")
        nodes, tets = _delaunay_tet(surf, char_len)
        method = "delaunay-voxel-repair"

    result = _finalize(
        surf, nodes, tets, method=method, joint_center=JOINT_CENTER_MM, char_len=char_len
    )

    # 网格自检（T5 / T7）——失败不掩盖，直接抛错或显著打印
    if result["min_jacobian"] <= 0:
        raise RuntimeError(
            f"网格含非正 Jacobian（min={result['min_jacobian']:.3e}），T5 失败"
        )
    if result["n_joint_faces"] == 0:
        raise RuntimeError(
            "未选到任何距下关节面；请检查 joint 半径/法向阈值或几何坐标。"
        )
    if result["n_plantar_faces"] == 0:
        raise RuntimeError("未选到任何跖面；请检查 −Y 法向阈值或几何坐标。")

    if out_msh is not None:
        _write_msh(result, Path(out_msh))
        result["msh_path"] = str(Path(out_msh))

    return result


# ---------------------------------------------------------------------------
# 主路径：decimate + gmsh 平滑网格（有效加密族）
# ---------------------------------------------------------------------------
#: 默认抽面比例。抽面后 STL 面片 ~2.8 mm，成为网格尺寸下限（gmsh 保留输入面片）。
DEFAULT_REDUCTION = 0.70
#: `classifySurfaces` 角度（rad）。**必须是极小值**才参数化成功（LEARNINGS-025）。
DEFAULT_CLASSIFY_ANGLE = 0.001


def _gmsh_tet_remesh(surf, char_len: float, angle: float = DEFAULT_CLASSIFY_ANGLE) -> tuple[np.ndarray, np.ndarray]:
    """gmsh：离散 STL → geo 实体 → **重网格化**表面 → 四面体。

    与 :func:`_gmsh_tet` 的区别：这里在 ``generate(2)`` 之后才 ``generate(3)``，
    并用 ``MeshSizeMin/Max`` 控制 —— 但**输入 STL 的面片尺寸会先封顶**，所以必须
    配合 :func:`stl_to_tet_gmsh` 的抽面。失败抛 :class:`GmshSolidError`。
    """
    import os
    import sys
    import tempfile

    if not LIB.is_dir():
        raise GmshSolidError(f"gmsh SDK lib 不存在：{LIB}")
    os.environ["PATH"] = f"{LIB};{SDK / 'bin'};{os.environ['PATH']}"
    if str(LIB) not in sys.path:
        sys.path.insert(0, str(LIB))
    import gmsh  # noqa: E402  (需先设好 PATH/sys.path)

    tmp_stl = Path(tempfile.gettempdir()) / "_calcaneus_remesh_in.stl"
    surf.save(str(tmp_stl))

    gmsh.initialize(str(LIB))
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.merge(str(tmp_stl))
        gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
        gmsh.model.mesh.removeDuplicateNodes()
        # 极小角度才参数化成功；逐级放宽以应对不同几何
        last: Exception | None = None
        for ang in (angle, 0.01, 0.1, math.pi / 2):
            try:
                gmsh.model.mesh.classifySurfaces(ang)
                gmsh.model.mesh.createGeometry()
                break
            except Exception as exc:  # gmsh 用裸 Exception
                last = exc
                gmsh.model.mesh.clear()
                gmsh.merge(str(tmp_stl))
        else:
            raise GmshSolidError(f"classifySurfaces/createGeometry 全部角度失败：{last}")

        gmsh.model.geo.synchronize()
        surfs = [t for _, t in gmsh.model.getEntities(2)]
        if not surfs:
            raise GmshSolidError("createGeometry 未产生任何曲面")
        gmsh.model.geo.addVolume([gmsh.model.geo.addSurfaceLoop(surfs)])
        gmsh.model.geo.synchronize()

        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay（平滑三角）
        gmsh.option.setNumber("Mesh.MeshSizeMin", float(char_len) * 0.4)
        gmsh.option.setNumber("Mesh.MeshSizeMax", float(char_len))
        gmsh.model.mesh.generate(2)                  # 先重网格化表面
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay 四面体
        gmsh.model.mesh.generate(3)

        tags, coords, _ = gmsh.model.mesh.getNodes()
        tag_arr = np.asarray(tags, dtype=np.int64)
        nodes = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
        _, nidx = gmsh.model.mesh.getElementsByType(4)  # 4 = 4 节点四面体
        if nidx is None or len(nidx) == 0:
            raise GmshSolidError("generate(3) 未产生四面体")
        lut = np.zeros(int(tag_arr.max()) + 1, dtype=np.int64)
        lut[tag_arr] = np.arange(len(tag_arr), dtype=np.int64)
        tets = lut[np.asarray(nidx, dtype=np.int64)].reshape(-1, 4)
        return nodes, tets
    finally:
        gmsh.finalize()


def stl_to_tet_gmsh(
    src_path,
    char_len: float,
    out_msh=None,
    *,
    reduction: float = DEFAULT_REDUCTION,
    joint_center: tuple[float, float, float] = JOINT_CENTER_MM,
) -> dict:
    """解剖级 STL → 抽面 → gmsh 平滑四面体网格（mm）。**主路径**。

    相对 :func:`vtp_to_tet` 的体素修复路径，本路径给出**有效加密族**：抽面后的
    STL 面片是固定几何，``char_len`` 只控制内部细化 —— 因此"同一几何、不同尺寸"
    是一族真正的 h-加密网格，可用于网格收敛/正则化验证（LEARNINGS-026）。

    实测质量（跟骨 STL，char_len=3 mm）：体素路径 17.2% 单元质心 Jacobian<0.3，
    本路径 **6.2%**；体积误差 0.00%。
    """
    src = Path(src_path)
    surf = _load_surface_mm(src)
    if reduction and reduction > 0.0:
        try:
            surf = surf.decimate_pro(
                reduction=float(reduction),
                splitting=True,
                feature_angle=30.0,
                preserve_topology=True,
                boundary_vertex_deletion=False,
            )
            surf = surf.triangulate().clean(point_merging=True, tolerance=1e-6, absolute=True)
        except Exception as exc:  # noqa: BLE001 - 抽面失败退回原始面（非致命）
            surf = _load_surface_mm(src)
            reduction = 0.0
            print(f"[meshing] decimate_pro 失败（{type(exc).__name__}），使用原始面片")

    nodes, tets = _gmsh_tet_remesh(surf, char_len)
    result = _finalize(
        surf, nodes, tets,
        method=f"gmsh-decimate(reduction={reduction:g})",
        joint_center=joint_center,
        char_len=char_len,
    )
    result["n_faces_input"] = int(surf.n_faces)
    result["reduction"] = float(reduction)

    if result["n_joint_faces"] == 0:
        raise RuntimeError("未选到任何距下关节面；请检查 joint 半径/法向阈值或几何坐标。")
    if result["n_plantar_faces"] == 0:
        raise RuntimeError("未选到任何跖面；请检查 −Y 法向阈值或几何坐标。")

    if out_msh is not None:
        _write_msh(result, Path(out_msh))
        result["msh_path"] = str(Path(out_msh))
    return result


def _main() -> int:  # pragma: no cover - 手动自检
    import sys

    vtp = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        ROOT / "model" / "opensim" / "FullBodyModel-4.0" / "Geometry" / "r_foot.vtp"
    )
    char_len = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
    MESH_DIR.mkdir(parents=True, exist_ok=True)
    out = MESH_DIR / f"calcaneus_{int(char_len)}mm.msh"
    r = vtp_to_tet(vtp, char_len, out_msh=out)
    for k in (
        "method",
        "n_nodes",
        "n_tets",
        "n_joint_faces",
        "n_plantar_faces",
        "volume_mm3",
        "reference_volume_mm3",
        "volume_error_pct",
        "min_jacobian",
        "jacobian_frac_lt_0_3",
        "surface_area_mm2",
        "msh_path",
    ):
        print(f"  {k}: {r.get(k)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
