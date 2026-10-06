"""用 pyfebio 0.3.0 生成跟骨 FEBio 模型（.feb）。

单位制
------
**mm–N–MPa–s**（MSK 骨 FE 惯例，见方案 §七 T10 / 附 A.7）：
长度 mm、力 N、应力 MPa、密度 tonne/mm³（1e-9）。全程一致。

材料
----
``IsotropicElastic``，E=15000 MPa、ν=0.3 —— *皮质骨量级*的各向同性近似
（``bone.py`` 的跟骨强度参考是 100/150 MPa，对应皮质骨；E 取皮质骨量级）。
这是 S1 的第一刀，不做皮质/松质分区。

边界条件与载荷（§3.2）
----------------------
* **载荷面**：由 ``meshing.vtp_to_tet`` 给出的具名面 ``subtalar_joint``
  （距下关节面，≈原点附近朝上的面）。
* **约束**：跖面（``plantar``，−Y 朝下的地面接触面）三向固定，消除刚体位移。
* **加载**：优先"幽灵刚体"模式
  ``add_simple_rigid_body`` + ``BCRigid`` + ``RigidForceLoad``；
  若该模式在 FEBio 4.13 上不收敛，调用方可改用 ``use_rigid=False`` 的
  ``PressureLoad`` 面压兜底（p = F/A，A=关节面面积）。
* **分析类型**：``Control(analysis="STATIC")``。S1 只求"给定关节力下峰值应力"，
  冲击时程留给后续 DYNAMIC 切片；static 是这一刀的正确最小模型。

pyfebio 0.3.0 / FEBio 4.13 兼容性
--------------------------------
pyfebio 0.3.0 按 FEBio 4.0 语法序列化；本机 ``febio4.exe`` 是 **4.13.0**，
它对 ``<solver type="solid">``、``<time_stepper type=...>``、``<qn_method type=...>``
等 ``type`` 属性报 "invalid attribute"。因此本模块：
1. ``time_stepper=None``（不写）；
2. 落盘后把 ``<solver>…</solver>`` 替换成合法的空 ``<solver/>``（用 FEBio 默认
   线性/非线性求解设置）。
这是显式的格式适配，不是猜测——见 :func:`_fix_febio413`。
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .plantar_bc import apply_plantar_bc

__all__ = [
    "build_calcaneus_feb",
    "joint_area_mm2",
    "joint_face_frame",
    "moment_to_pressure_gradient",
    "linear_pressure_bands",
    "DEFAULT_JOINT_FACE",
]

DEFAULT_JOINT_FACE = "subtalar_joint"

#: pyfebio/FEBio 版本适配：把 pyfebio 写出的 solver 块换成合法的空 <solver/>
_SOLVER_RE = re.compile(r"<solver[^>]*>.*?</solver>", re.S)


def joint_area_mm2(mesh: dict, joint_face: str = DEFAULT_JOINT_FACE) -> float:
    """具名关节面的面积（mm²），用于 PressureLoad 的 p=F/A。"""
    nodes = np.asarray(mesh["nodes"], dtype=np.float64)
    tris = np.asarray(mesh["surfaces"][joint_face], dtype=np.int64)
    if len(tris) == 0:
        raise ValueError(f"关节面 {joint_face!r} 为空")
    p = nodes[tris]
    return float(0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1).sum())


# ---------------------------------------------------------------------------
# Phase-S4：3D wrench（力矢量 + 力矩）→ FE 边界的共享几何/载荷原语
# ---------------------------------------------------------------------------
def joint_face_frame(mesh: dict, joint_face: str = DEFAULT_JOINT_FACE) -> dict:
    """关节载荷面的**面积加权几何量**（mm-N-MPa 单位）。

    返回 ``area_mm2``、``centroid_mm``、``normal_unit``（外法向，取 n·(+Y)>0）、
    ``second_moment_mm4``（相对形心的二阶矩 ``I_c = ∫(x−c)(x−c)ᵀ dA``）、
    逐面片 ``facet_centroids_mm`` / ``facet_areas_mm2`` / ``facet_normals_unit``。

    面片法向**逐片朝外**校正（与曲面曲率无关，只是把不一致的三角形绕向统一）。
    这些量用于把 3D 力矩表达成关节面线性压力梯度，以及把力矢量分解为
    法向压力 + 切向牵引。
    """
    nodes = np.asarray(mesh["nodes"], dtype=np.float64)
    tris = np.asarray(mesh["surfaces"][joint_face], dtype=np.int64).reshape(-1, 3)
    if len(tris) == 0:
        raise ValueError(f"关节面 {joint_face!r} 为空，无法建立 wrench 几何")
    p = nodes[tris]
    cen = p.mean(axis=1)
    cross = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    area = 0.5 * np.linalg.norm(cross, axis=1)
    if float(area.sum()) <= 0.0:
        raise ValueError(f"关节面 {joint_face!r} 面积非正")
    nf = cross / np.linalg.norm(cross, axis=1, keepdims=True)

    nhat = (area[:, None] * nf).sum(axis=0)
    nhat = nhat / np.linalg.norm(nhat)
    sign = 1.0 if nhat[1] >= 0.0 else -1.0
    nf = nf * sign
    nhat = nhat * sign
    flip = (nf @ nhat) < 0.0
    nf[flip] *= -1.0

    total_area = float(area.sum())
    c = (area[:, None] * cen).sum(axis=0) / total_area
    # 二阶矩：逐三角形**精确积分** ∫x xᵀ dA = (A/12)[Σ p_k p_kᵀ + (Σp_k)(Σp_k)ᵀ]，
    # 再平行轴到形心 I_c = ∫x xᵀ dA − A c cᵀ。（用质心点质量近似对粗面片误差大。）
    M2 = np.zeros((3, 3), dtype=np.float64)
    for tri, af in zip(tris, area):
        q = nodes[tri]                    # (3,3)
        s = q.sum(axis=0)
        M2 += (af / 12.0) * (q.T @ q + np.outer(s, s))
    I = M2 - total_area * np.outer(c, c)
    return {
        "joint_face": joint_face,
        "area_mm2": total_area,
        "centroid_mm": c,
        "normal_unit": nhat,
        "second_moment_mm4": I,
        "facet_centroids_mm": cen,
        "facet_areas_mm2": area,
        "facet_normals_unit": nf,
        "n_facets": int(len(tris)),
    }


def _uniform_pressure_moment_nmm(frame: dict, pressure_mpa: float) -> np.ndarray:
    """均匀法向压力 ``p0`` 相对形心的力矩（N·mm）。

    曲面上逐片法向变化 → 均匀压力也带一个力矩（平面面上恒为 0）。这个量在
    :func:`moment_to_pressure_gradient` 里先扣除，使**总**载荷（均匀压力 + 梯度）
    的力矩等于目标（见 ``base_uniform_pressure_mpa``）。
    """
    cen = np.asarray(frame["facet_centroids_mm"], dtype=np.float64)
    nf = np.asarray(frame["facet_normals_unit"], dtype=np.float64)
    area = np.asarray(frame["facet_areas_mm2"], dtype=np.float64)
    c = np.asarray(frame["centroid_mm"], dtype=np.float64)
    d = cen - c
    return float(pressure_mpa) * np.einsum("f,fi->i", area, np.cross(d, -nf))


def moment_to_pressure_gradient(
    moment_nm, frame: dict, *, base_uniform_pressure_mpa: float = 0.0
) -> tuple[float, np.ndarray, dict]:
    """把力矩 (N·m,3) 解成关节面**线性压力梯度** ``p(x)=p0+g·x`` 的坡度与方向。

    线性压力场 ``p=g·(x−c)`` 对面片求矩关于形心是 ``g`` 的**精确线性函数**：::

        M(g) = J g,   J_ij = Σ_f A_f (d_f × (−n̂_f))_i d_f,j     (d_f = x_f − c)

    其中 ``n̂_f`` 是**逐面片真实外法向**（已含曲面曲率）。``base_uniform_pressure_mpa``
    给出随力矢量一起施加的**均匀法向压力 p0**；曲面下 p0 自带力矩
    ``M(p0)``（见 :func:`_uniform_pressure_moment_nmm`），故梯度目标取
    ``M_needed = M_target − M(p0)``，解最小二乘 ``J g = M_needed``，返回
    ``(slope, direction_unit, meta)``，``slope=|g|``、``direction_unit=g/|g|``
    （零力矩 → slope=0、方向零向量）。

    **边界（诚实）**：法向压力场对**平面**面的力矩恒 ⊥ 法向（扭转分量无法表达）；
    曲面上逐片法向变化使该限制放松但引入近似。完整/精确力矩请走 ``use_rigid=True``
    的 ``RigidMomentLoad`` 路径。
    """
    nhat = np.asarray(frame["normal_unit"], dtype=np.float64)
    cen = np.asarray(frame["facet_centroids_mm"], dtype=np.float64)
    nf = np.asarray(frame["facet_normals_unit"], dtype=np.float64)
    area = np.asarray(frame["facet_areas_mm2"], dtype=np.float64)
    c = np.asarray(frame["centroid_mm"], dtype=np.float64)
    M_nmm = np.asarray(moment_nm, dtype=np.float64).ravel() * 1000.0
    if M_nmm.shape != (3,):
        raise ValueError(f"moment_nm 必须是 3 分量，得到 {moment_nm!r}")
    if not np.all(np.isfinite(M_nmm)):
        raise ValueError(f"moment_nm 含 NaN/Inf：{moment_nm!r}")
    d = cen - c
    cross_term = np.cross(d, -nf)                       # (n,3)
    J = np.einsum("f,fi,fj->ij", area, cross_term, d)   # 3×3，M(g)=J g
    M_base = _uniform_pressure_moment_nmm(frame, base_uniform_pressure_mpa)
    M_needed = M_nmm - M_base
    g, *_ = np.linalg.lstsq(J, M_needed, rcond=None)
    slope = float(np.linalg.norm(g))
    direction = g / slope if slope > 0.0 else np.zeros(3, dtype=np.float64)
    achieved = M_base + J @ g
    residual = M_nmm - achieved
    M_ip = M_nmm - float(M_nmm @ nhat) * nhat
    meta = {
        "moment_nm": [float(x) for x in np.asarray(moment_nm, dtype=float).ravel()],
        "moment_inplane_nm": [float(x) for x in (M_ip / 1000.0)],
        "base_uniform_pressure_mpa": float(base_uniform_pressure_mpa),
        "base_uniform_moment_nm": [float(x) for x in (M_base / 1000.0)],
        "achievable_moment_nm": [float(x) for x in (achieved / 1000.0)],
        "dropped_normal_axis_nm": float(residual @ nhat) / 1000.0,
        "moment_residual_nm": float(np.linalg.norm(residual)) / 1000.0,
        "slope_mpa_per_mm": slope,
        "direction_unit": [float(x) for x in direction],
        "facet_exact": True,
    }
    return slope, direction, meta


def linear_pressure_bands(
    model,
    fmesh,
    mesh: dict,
    joint_face: str,
    frame: dict,
    *,
    p0_mpa: float = 0.0,
    slope: float = 0.0,
    direction=None,
    n_bands: int = 6,
    name_prefix: str | None = None,
) -> tuple[list[tuple[str, float]], dict]:
    """把压力场 ``p(x)=p0+slope·((x−c)·ĝ)`` 表达成关节面分带常压。

    沿面内单位方向 ``ĝ`` 把关节面三角片分到 ``n_bands`` 个带，每带取**面积加权
    平均压力**。梯度项按面积加权零点平移（``t_ref``），因此 ``Σ_b P_b A_b = p0·A``
    中梯度项贡献严格为零 → **合力精确守恒**（只由 p0 决定）。返回
    ``([(surface_name, pressure_MPa)], meta)``，并把每个带注册到 ``model.mesh_``。

    ``meta["realized_moment_nm"]`` 用**真实逐片法向**数值积分带化压力场相对形心的
    力矩（含曲面效应），供报告核对；与目标力矩之差即平面近似 + 分带离散误差。
    """
    tris = np.asarray(mesh["surfaces"][joint_face], dtype=np.int64).reshape(-1, 3)
    cen = np.asarray(frame["facet_centroids_mm"], dtype=np.float64)
    area = np.asarray(frame["facet_areas_mm2"], dtype=np.float64)
    nf = np.asarray(frame["facet_normals_unit"], dtype=np.float64)
    c = np.asarray(frame["centroid_mm"], dtype=np.float64)
    total_area = float(area.sum())

    n_bands = int(n_bands)
    if n_bands < 1:
        raise ValueError(f"n_bands 至少 1，得到 {n_bands}")

    direction = np.zeros(3) if direction is None else np.asarray(direction, dtype=float).ravel()
    slope = float(slope)
    if direction.shape != (3,) or not np.all(np.isfinite(direction)):
        raise ValueError(f"direction 必须是有限 3 向量，得到 {direction!r}")
    dnorm = float(np.linalg.norm(direction))
    if slope != 0.0 and dnorm > 0.0:
        de = direction / dnorm
        t = cen @ de
        t_ref = float((area * t).sum() / total_area)
        p_f = float(p0_mpa) + slope * (t - t_ref)
    else:
        de = np.zeros(3)
        t = np.zeros(len(cen))
        p_f = np.full(len(cen), float(p0_mpa))

    if len(np.unique(t)) > 1 and n_bands > 1:
        edges = np.linspace(float(t.min()), float(t.max()), n_bands + 1)
        bidx = np.clip(np.digitize(t, edges[1:-1]), 0, n_bands - 1)
    else:
        bidx = np.zeros(len(t), dtype=np.int64)

    prefix = name_prefix or f"{joint_face}_g"
    bands: list[tuple[str, float]] = []
    M_real = np.zeros(3)
    for b in range(int(bidx.max()) + 1):
        sel = np.flatnonzero(bidx == b)
        if len(sel) == 0:
            continue
        ab = float(area[sel].sum())
        pb = float((p_f[sel] * area[sel]).sum() / ab)
        name = f"{prefix}{b}"
        surf = fmesh.Surface(name=name)
        for i, fi in enumerate(sel):
            tri = tris[int(fi)]
            surf.add_tri3(fmesh.Tri3Element(id=i + 1, text=",".join(map(str, tri + 1))))
        model.mesh_.add_surface(surf)
        bands.append((name, pb))
        for fi in sel:
            M_real += np.cross(cen[int(fi)] - c, -pb * nf[int(fi)]) * area[int(fi)]

    net_force = float((p_f * area).sum())
    meta = {
        "joint_face": joint_face,
        "p0_mpa": float(p0_mpa),
        "slope_mpa_per_mm": slope,
        "direction_unit": [float(x) for x in de],
        "n_bands": len(bands),
        "bands": [{"surface": nm, "pressure_mpa": pb} for nm, pb in bands],
        "facet_pressure_min_mpa": float(p_f.min()),
        "facet_pressure_max_mpa": float(p_f.max()),
        "net_normal_force_n": float(net_force),
        "realized_moment_nm": [float(x) for x in (M_real / 1000.0)],
    }
    return bands, meta


def _fix_febio413(path: Path) -> list[str]:
    """把 FEBio 4.0 语法微调到本机 4.13.0 能读的形式，返回改动说明。"""
    text = path.read_text(encoding="ISO-8859-1")
    notes: list[str] = []
    if _SOLVER_RE.search(text):
        text = _SOLVER_RE.sub("<solver/>", text)
        notes.append("solver 块 → <solver/>（FEBio 4.13 拒绝其 type 属性，改用默认设置）")
    path.write_text(text, encoding="ISO-8859-1")
    return notes


def build_calcaneus_feb(
    mesh: dict,
    out_feb,
    *,
    E_mpa: float = 15000.0,
    nu: float = 0.3,
    load_n: float = 4000.0,
    load_dir=(0.0, -1.0, 0.0),
    joint_face=None,
    use_rigid: bool = True,
    time_steps: int = 1,
    step_size: float | None = None,
    achilles_n: float = 0.0,
    achilles_dir=(0.15, 0.985, 0.0),
    achilles_surface: str = "achilles",
    plantar_bc: str = "fixed",
    spring_k: float = 1000.0,
    subtalar_force=None,
    subtalar_moment=None,
    moment_bands: int = 6,
) -> Path:
    """由 :func:`climbing.coupling.meshing.vtp_to_tet` 的返回值写出 ``.feb``。

    Parameters
    ----------
    mesh:
        ``vtp_to_tet`` 的返回字典（mm）。
    out_feb:
        输出 ``.feb`` 路径。
    E_mpa, nu:
        各向同性线弹性参数（MPa，无量纲）。
    load_n:
        关节合力大小（N）。**这是一个参数，不是物理结论**；默认 4000 N 竖直。
    load_dir:
        载荷方向（单位向量，FE 坐标系）。
    joint_face:
        施加关节载荷的具名面；默认 ``"subtalar_joint"``。
    use_rigid:
        True=幽灵刚体 + RigidForceLoad；False=PressureLoad 兜底（p=F/A）。
    time_steps, step_size:
        STATIC 载荷步数与每步时长（默认 10×0.1，载荷线性加载到 1.0）。
    plantar_bc:
        跖面 BC 类型：``"fixed"``（默认，三向零位移，原行为）/ ``"roller"``
        （法向零位移 + 切向销）/ ``"roller_free"``（诊断变体）/ ``"spring"``
        （三向 Winkler 弹簧，刚度 ``spring_k`` N/mm）/ ``"contact"``（**波3 最小版**：
        跖面 ⇄ 对偶固定薄板（数值刚性壁）的真实 FEBio 接触 + Coulomb 摩擦，recipe
        取 ``pad_contact`` 的 udg NEW PENALTY 配方；**opt-in，默认不变**）。
        见 :mod:`climbing.coupling.plantar_bc`。
    spring_k:
        Winkler 弹簧刚度 k（N/mm），仅 ``plantar_bc="spring"`` 时使用。
    subtalar_force:
        **Phase-S4（opt-in，默认 ``None``）** 完整 3D 距下力矢量 (N,3)。默认 ``None``
        时用 ``load_n × load_dir``（旧标量口径，逐位不变）。给了就替换旧口径：
        ``use_rigid=True`` 时按分量走 ``RigidForceLoad``；``use_rigid=False`` 时按
        关节面法向分解为**法向压力**(``PressureLoad``) + **切向牵引**(``TractionLoad``)，
        完整传递含剪切的 3D 力。
    subtalar_moment:
        **Phase-S4（opt-in，默认 ``None``）** 完整 3D 距下力矩 (N·m,3)。默认 ``None``
        = 零力偶（旧行为）。``use_rigid=True`` 时用 ``RigidMomentLoad``（值按 N·mm，
        即 ×1000）；``use_rigid=False`` 时用关节面**线性压力梯度** ``p(x)=p0+g·x``
        （分带常压表达，见 :func:`linear_pressure_bands`）。压力场无法表达沿面法向的
        扭转分量，该分量会被丢弃并记录（见 :func:`moment_to_pressure_gradient`）。
    moment_bands:
        压力梯度分带数（默认 6），仅 ``use_rigid=False`` 且给了 ``subtalar_moment`` 时用。

    Returns
    -------
    pathlib.Path
        写好的 ``.feb``。附带同名 ``.feb.notes.txt``（若发生格式适配）——不，直接
        返回路径；适配说明同时打印。
    """
    # 延迟 import：避免在没装 pyfebio 的环境 import 本模块就炸。
    from pyfebio.model import Model
    from pyfebio import (
        boundary as fbc,
        control,
        material,
        mesh as fmesh,
        meshdomains,
        output as fout,
        rigid,
        step as fstep,
        loads as floads,
    )

    joint_face = joint_face or DEFAULT_JOINT_FACE
    out_feb = Path(out_feb)
    out_feb.parent.mkdir(parents=True, exist_ok=True)

    nodes = np.asarray(mesh["nodes"], dtype=np.float64)
    tets = np.asarray(mesh["tets"], dtype=np.int64)
    if joint_face not in mesh.get("surfaces", {}):
        raise KeyError(f"mesh 里没有具名面 {joint_face!r}；有：{list(mesh.get('surfaces', {}))}")
    if "plantar" not in mesh.get("node_sets", {}):
        raise KeyError("mesh 里没有 'plantar' 节点集（跖面约束所需）")

    d = np.asarray(load_dir, dtype=float)
    nrm = float(np.linalg.norm(d))
    if nrm == 0.0:
        raise ValueError("load_dir 不能是零向量")
    d = d / nrm

    # --- Phase-S4：可选 3D wrench（默认 None → 旧标量口径逐位不变）----------
    force_vec: np.ndarray | None = None
    if subtalar_force is not None:
        force_vec = np.asarray(subtalar_force, dtype=float).ravel()
        if force_vec.shape != (3,):
            raise ValueError(f"subtalar_force 必须是 3 分量 (N)，得到 {subtalar_force!r}")
        if not np.all(np.isfinite(force_vec)):
            raise ValueError(f"subtalar_force 含 NaN/Inf：{subtalar_force!r}")
    moment_vec: np.ndarray | None = None
    if subtalar_moment is not None:
        moment_vec = np.asarray(subtalar_moment, dtype=float).ravel()
        if moment_vec.shape != (3,):
            raise ValueError(f"subtalar_moment 必须是 3 分量 (N·m)，得到 {subtalar_moment!r}")
        if not np.all(np.isfinite(moment_vec)):
            raise ValueError(f"subtalar_moment 含 NaN/Inf：{subtalar_moment!r}")
    wrench = force_vec is not None or moment_vec is not None

    model = Model()
    # --- Mesh -------------------------------------------------------------
    model.mesh_.add_node_domain(fmesh.numpy_to_nodes(nodes, name="calcaneus"))
    model.mesh_.add_element_domain(fmesh.numpy_to_elements(tets + 1, "tet4", name="calcaneus"))

    joint_nodes = np.asarray(mesh["node_sets"][joint_face], dtype=np.int64)
    plantar_nodes = np.asarray(mesh["node_sets"]["plantar"], dtype=np.int64)
    if len(joint_nodes) == 0:
        raise ValueError("关节节点集为空")
    if len(plantar_nodes) == 0:
        raise ValueError("跖面节点集为空")
    model.mesh_.add_node_set(fmesh.NodeSet(name=joint_face, text=",".join(map(str, joint_nodes + 1))))
    model.mesh_.add_node_set(fmesh.NodeSet(name="plantar", text=",".join(map(str, plantar_nodes + 1))))

    # 具名关节面（PressureLoad 兜底用；也便于 .feb 里追溯）
    surf = fmesh.Surface(name=joint_face)
    for i, tri in enumerate(np.asarray(mesh["surfaces"][joint_face], dtype=np.int64)):
        surf.add_tri3(fmesh.Tri3Element(id=i + 1, text=",".join(map(str, tri + 1))))
    model.mesh_.add_surface(surf)

    # 跟腱附着面（可选）：后上结节。提供具名 node_set + surface，供 TractionLoad 用。
    has_achilles = (
        float(achilles_n) > 0.0
        and achilles_surface in mesh.get("surfaces", {})
        and len(np.asarray(mesh["surfaces"][achilles_surface])) > 0
    )
    if achilles_surface in mesh.get("surfaces", {}):
        ach_nodes = np.asarray(mesh["node_sets"][achilles_surface], dtype=np.int64)
        if len(ach_nodes) > 0:
            model.mesh_.add_node_set(
                fmesh.NodeSet(name=achilles_surface, text=",".join(map(str, ach_nodes + 1)))
            )
        ach_surf = fmesh.Surface(name=achilles_surface)
        for i, tri in enumerate(np.asarray(mesh["surfaces"][achilles_surface], dtype=np.int64)):
            ach_surf.add_tri3(fmesh.Tri3Element(id=i + 1, text=",".join(map(str, tri + 1))))
        model.mesh_.add_surface(ach_surf)

    # --- Material ---------------------------------------------------------
    model.material_.add_material(
        material.IsotropicElastic(
            name="bone",
            E=material.MaterialParameter(text=float(E_mpa)),
            v=material.MaterialParameter(text=float(nu)),
        )
    )
    model.meshdomains_.add_solid_domain(meshdomains.SolidDomain(name="calcaneus", mat="bone"))

    # --- Step -------------------------------------------------------------
    # 为什么默认只走 1 个载荷步：FEBio 的 STATIC 收敛判据用**相对**残差/位移；
    # 多步时后续步的增量会小到 ~1e-21，required 也等比缩到 ~1e-27，求解器在
    # 机器精度上反复迭代却"永不满足"，最终 Max nr of reformations。
    # 本切片只求"给定载荷下的峰值应力"，一次加载到位（总时间 1.0）是正确最小模型。
    if step_size is None:
        step_size = 1.0 / max(int(time_steps), 1)
    st = fstep.StepEntry(id=1, name="Step")
    st.control = control.Control(
        analysis="STATIC",
        time_steps=int(time_steps),
        step_size=float(step_size),
        time_stepper=None,  # FEBio 4.13 拒绝 <time_stepper type=...>
    )

    bcs: list = []
    plantar_meta = apply_plantar_bc(
        model,
        nodes=nodes,
        plantar_local=plantar_nodes,
        plantar_tris=np.asarray(
            mesh.get("surfaces", {}).get("plantar", []), dtype=np.int64
        ).reshape(-1, 3),
        bcs=bcs,
        kind=plantar_bc,
        spring_k=spring_k,
    )
    # ``plantar_bc="contact"`` (opt-in): the contact interface lives in the step's
    # ``<Contact>`` block.  ``apply_plantar_bc`` only sees ``model`` (not the step),
    # so it returns the pyfebio ``Contact`` object here and we attach it.  Default
    # ``"fixed"`` returns no ``"contact"`` key ⇒ bit-for-bit unchanged deck.
    if plantar_meta.get("contact") is not None:
        st.contact = plantar_meta["contact"]

    surface_loads = []
    if use_rigid:
        joint_center = tuple(float(x) for x in mesh.get("joint_center_mm", (0.0, 0.0, 0.0)))
        model.add_simple_rigid_body(joint_center, "ghost")
        bcs.append(fbc.BCRigid(name="joint_to_ghost", node_set=joint_face, rb="ghost"))
        st.boundary = fbc.Boundary(all_bcs=bcs)
        # 纯力时锁住刚体转动、放开平动（力由骨网格约束）。给了力矩时，必须**放开
        # 承载该力矩的转动自由度**，否则 RigidMomentLoad 会被固定的转动 DOF 吸收而无效；
        # 无力矩的转动 DOF 仍锁住，避免无力矩的自由转动模态。
        if moment_vec is None:
            rigid_bcs = [rigid.RigidFixed(rb="ghost", Ru_dof=1, Rv_dof=1, Rw_dof=1)]
        else:
            lock = {dof: 1 for i, dof in enumerate(("Ru", "Rv", "Rw")) if moment_vec[i] == 0.0}
            rigid_bcs = [rigid.RigidFixed(rb="ghost", **lock)] if lock else []
        rigid_loads = []
        fv = force_vec if force_vec is not None else load_n * d
        for comp, dof in enumerate(("Rx", "Ry", "Rz")):
            if fv[comp] != 0.0:
                rigid_loads.append(
                    rigid.RigidForceLoad(
                        rb="ghost", dof=dof, value=rigid.Value(text=float(fv[comp]))
                    )
                )
        if moment_vec is not None:
            # RigidMomentLoad.value 用 N·mm（力×长度）；输入是 N·m。
            for comp, dof in enumerate(("Ru", "Rv", "Rw")):
                if moment_vec[comp] != 0.0:
                    rigid_loads.append(
                        rigid.RigidMomentLoad(
                            rb="ghost", dof=dof,
                            value=rigid.Value(text=float(moment_vec[comp] * 1000.0)),
                        )
                    )
        st.rigid = rigid.Rigid(all_rigid_bcs=rigid_bcs, all_rigid_loads=rigid_loads)
        fmag = float(np.linalg.norm(fv))
        load_desc = f"幽灵刚体 RigidForceLoad |F|={fmag:.1f} N × ({fv[0]:.3f},{fv[1]:.3f},{fv[2]:.3f})"
        if moment_vec is not None:
            load_desc += (f" + RigidMomentLoad |M|={np.linalg.norm(moment_vec):.2f} N·m "
                          f"({moment_vec[0]:.3f},{moment_vec[1]:.3f},{moment_vec[2]:.3f})")
    elif wrench:
        # Phase-S4 压力路径：法向压力（合力 p0·A）+ 切向牵引 + 力矩线性压力梯度。
        frame = joint_face_frame(mesh, joint_face)
        fv = force_vec if force_vec is not None else load_n * d
        f_n = float(fv @ frame["normal_unit"])
        f_t = fv - f_n * frame["normal_unit"]
        p0 = -f_n / frame["area_mm2"]  # 压入面为正（FEBio: F=−p·n̂_out）
        st.boundary = fbc.Boundary(all_bcs=bcs)
        if moment_vec is None:
            surface_loads.append(
                floads.PressureLoad(surface=joint_face, pressure=floads.Scale(text=float(p0)))
            )
            band_meta = None
        else:
            slope, direction, gmeta = moment_to_pressure_gradient(
                moment_vec, frame, base_uniform_pressure_mpa=p0
            )
            bands, band_meta = linear_pressure_bands(
                model, fmesh, mesh, joint_face, frame,
                p0_mpa=p0, slope=slope, direction=direction,
                n_bands=moment_bands, name_prefix=f"{joint_face}_g",
            )
            for name, pb in bands:
                surface_loads.append(
                    floads.PressureLoad(surface=name, pressure=floads.Scale(text=float(pb)))
                )
            band_meta["dropped_normal_axis_nm"] = gmeta["dropped_normal_axis_nm"]
            band_meta["moment_inplane_nm"] = gmeta["moment_inplane_nm"]
        ft_mag = float(np.linalg.norm(f_t))
        if ft_mag > 1e-12:
            ta = f_t / frame["area_mm2"]
            surface_loads.append(
                floads.TractionLoad(
                    surface=joint_face,
                    scale=floads.Scale(text="1.0"),
                    traction=f"{ta[0]:.6f},{ta[1]:.6f},{ta[2]:.6f}",
                )
            )
        load_desc = (
            f"3D wrench: F=({fv[0]:.1f},{fv[1]:.1f},{fv[2]:.1f}) N，法向压力 p0={p0:.4f} MPa "
            f"+ 切向牵引 |F_t|={ft_mag:.1f} N"
            + (f" + 力矩梯度 {band_meta['n_bands']} 带" if band_meta else "")
        )
    else:
        area = joint_area_mm2(mesh, joint_face)
        pressure = float(load_n) / area  # N / mm² = MPa
        st.boundary = fbc.Boundary(all_bcs=bcs)
        surface_loads.append(
            floads.PressureLoad(surface=joint_face, pressure=floads.Scale(text=pressure))
        )
        load_desc = f"PressureLoad p={pressure:.4f} MPa（F={load_n:.1f} N / A={area:.2f} mm²）"

    # 跟腱牵引（后上结节）：traction = F/A × 单位方向（MPa）。把大小直接写进向量，
    # 不用 scale 相乘，避免 FEBio 的 scale 语义歧义。
    if has_achilles:
        ad = np.asarray(achilles_dir, dtype=float)
        nn = float(np.linalg.norm(ad))
        if nn == 0.0:
            raise ValueError("achilles_dir 不能是零向量")
        ad = ad / nn
        a_area = joint_area_mm2(mesh, achilles_surface)
        trac = float(achilles_n) / a_area  # MPa
        vec = trac * ad
        surface_loads.append(
            floads.TractionLoad(
                surface=achilles_surface,
                scale=floads.Scale(text="1.0"),  # lc 缺省=常量；默认 Scale(lc=1) 会引用不存在的曲线
                traction=f"{vec[0]:.6f},{vec[1]:.6f},{vec[2]:.6f}",
            )
        )
        load_desc += (f" + 跟腱牵引 F={achilles_n:.1f} N / A={a_area:.1f} mm² = {trac:.4f} MPa "
                      f"× ({ad[0]:.2f},{ad[1]:.2f},{ad[2]:.2f})")

    if surface_loads:
        st.loads = floads.Loads(all_surface_loads=surface_loads)

    # --- Output -----------------------------------------------------------
    st.control.plot_level = "PLOT_MAJOR_ITRS"
    model.step_.add_step(st)
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
        print(f"[febio_model] FEBio 4.13 适配：{n}")
    print(f"[febio_model] 载荷：{load_desc}")
    return out_feb
