"""FEBio 后处理：``.xplt`` → 峰值 von Mises 应力。

链路（方案 §3.5 / 附 A.4）
--------------------------
``.feb`` --(run_febio)--> ``.xplt`` --(pyfebio.xplt.to_hdf5)--> ``.hdf5``
--> 按单元读 Voigt ``[sxx,syy,szz,sxy,syz,sxz]`` --> von Mises --> 峰值。

von Mises（Voigt 约定 sxy=σ_xy, syz=σ_yz, sxz=σ_xz）::

    σ_vm = sqrt( 0.5[(σxx−σyy)² + (σyy−σzz)² + (σzz−σxx)²]
                 + 3(σxy² + σyz² + σxz²) )

单位：``.feb`` 是 mm–N–MPa–s，所以读出的应力直接是 **MPa**。

``pyfebio.xplt.to_hdf5`` 把每个 state 存成独立 dataset（``ERRORS.md
[ERR-...-007]`` 的布局）：``states/{id}/element_data/stress/{domain}``。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = [
    "peak_von_mises",
    "run_and_peak",
    "von_mises_from_voigt",
    "von_mises_stats",
    "von_mises_gauge",
    "gauge_von_mises",
    "element_centroids_and_volumes",
    # --- S5 §3.3（T3）：位移 / 支反力 / 接触合力 ---
    "build_hdf5",
    "read_displacement",
    "read_reaction_forces",
    "contact_force_sum",
]


def von_mises_from_voigt(s: np.ndarray) -> np.ndarray:
    """(...,6) Voigt 应力 → von Mises（MPa）。列序 [xx,yy,zz,xy,yz,xz]。"""
    s = np.asarray(s, dtype=np.float64)
    if s.shape[-1] != 6:
        raise ValueError(f"需要 6 分量 Voigt 应力，得到 shape={s.shape}")
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(
        0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2)
        + 3.0 * (xy**2 + yz**2 + xz**2)
    )


def _domain_datasets(state_group, domain: str) -> list[str]:
    """返回该 state 下 element_data/stress 里可用的 domain 名。"""
    if "element_data" not in state_group or "stress" not in state_group["element_data"]:
        return []
    keys = list(state_group["element_data"]["stress"].keys())
    if domain in keys:
        return [domain]
    # 兜底：取第一个非刚体域
    rest = [k for k in keys if k.lower() != "ghost"]
    return rest[:1]


def _read_last_stress(h5: "object", domain: str) -> tuple[np.ndarray, str]:
    """取"最后一个 state"的单元应力（STATIC 末步=满载荷）。"""
    state_ids = sorted(h5["states"].keys(), key=int)
    if not state_ids:
        raise RuntimeError("HDF5 里没有 states：FEBio 可能未写出结果")
    for sid in reversed(state_ids):
        grp = h5["states"][sid]
        names = _domain_datasets(grp, domain)
        if names:
            return np.asarray(grp["element_data"]["stress"][names[0]][:], dtype=np.float64), names[0]
    raise RuntimeError(f"任何 state 都找不到 {domain!r} 的单元应力")


def peak_von_mises(
    xplt_path,
    *,
    hdf5_path=None,
    domain: str = "calcaneus",
    mesh: dict | None = None,
    write_vtu=None,
) -> float:
    """读 ``.xplt``（或直接给 ``.feb``）→ 返回峰值 von Mises（MPa）。

    Parameters
    ----------
    xplt_path:
        ``.xplt`` 结果文件；若给 ``.feb``，先经 :func:`run_and_peak` 的
        ``run_febio`` 求解再读。
    hdf5_path:
        中间 HDF5 路径；默认与 xplt 同名的 ``.hdf5``。
    domain:
        element block 名（本切片是 ``"calcaneus"``）。
    mesh:
        给了且 ``write_vtu`` 非空时，写一个**节点平均**的 ``.vtu``。
    write_vtu:
        ``.vtu`` 输出路径。

    Returns
    -------
    float
        峰值 von Mises（MPa）。所有 state、所有单元取最大。
    """
    from pyfebio import xplt as xplt_mod
    import h5py

    p = Path(xplt_path)
    if p.suffix.lower() == ".feb":
        from climbing.coupling.febio_run import run_febio

        run_febio(p)  # 失败会抛错（绝不静默 rc!=0）
        p = p.with_suffix(".xplt")
    elif not p.is_file():
        # 便利：只给了 .xplt 名但还没求解时，找同名 .feb 先跑（S1 入口就是这么调的）。
        sibling_feb = p.with_suffix(".feb")
        if sibling_feb.is_file():
            from climbing.coupling.febio_run import run_febio

            run_febio(sibling_feb)
        else:
            raise FileNotFoundError(f"结果文件不存在：{p}（也没有同名 .feb 可求解）")
    if not p.is_file():
        raise FileNotFoundError(f"求解后仍未生成结果文件：{p}")

    h5p = Path(hdf5_path) if hdf5_path is not None else p.with_suffix(".hdf5")
    xplt_mod.to_hdf5(str(p), str(h5p))

    with h5py.File(h5p, "r") as h5:
        stress, used_domain = _read_last_stress(h5, domain)
    vm = von_mises_from_voigt(stress)

    if write_vtu is not None and mesh is not None:
        _write_node_averaged_vtu(Path(write_vtu), mesh, stress, vm)
    return float(vm.max())


def von_mises_stats(
    xplt_path,
    *,
    hdf5_path=None,
    domain: str = "calcaneus",
    mesh: dict | None = None,
    write_vtu=None,
    percents: tuple[float, ...] = (95.0, 99.0),
) -> dict:
    """读结果 → von Mises 分布统计（max / p95 / p99 / mean / median）。

    **为什么需要它**：峰值 σ_vm 受**边界奇异**支配——面压载荷只加到关节面边缘、
    跖面又全约束，棱边处应力随网格加密无界增长。实测同一模型 4→3→2 mm 网格峰值
    126→180→223 MPa（未收敛）；因此**峰值不能直接当强度判据**。稳健做法取
    **p95/p99 或体平均**（避开奇异带）。
    """
    from pyfebio import xplt as xplt_mod
    import h5py

    p = Path(xplt_path)
    if not p.is_file():
        raise FileNotFoundError(f"结果文件不存在：{p}")
    h5p = Path(hdf5_path) if hdf5_path is not None else p.with_suffix(".hdf5")
    xplt_mod.to_hdf5(str(p), str(h5p))
    with h5py.File(h5p, "r") as h5:
        stress, _ = _read_last_stress(h5, domain)
    vm = von_mises_from_voigt(stress)

    if write_vtu is not None and mesh is not None:
        _write_node_averaged_vtu(Path(write_vtu), mesh, stress, vm)

    out = {
        "max": float(vm.max()),
        "mean": float(vm.mean()),
        "median": float(np.median(vm)),
        "n_elem": int(vm.size),
    }
    for pc in percents:
        out[f"p{pc:g}"] = float(np.percentile(vm, pc))
    return out


# ---------------------------------------------------------------------------
# 过程区正则化（把边界奇异峰值 → 网格无关的稳健应力）
# ---------------------------------------------------------------------------
def element_centroids_and_volumes(mesh: dict) -> tuple[np.ndarray, np.ndarray]:
    """tet 网格 → 单元质心 ``(n,3)`` 与体积 ``(n,)``（mm）。"""
    nodes = np.asarray(mesh["nodes"], dtype=np.float64)
    tets = np.asarray(mesh["tets"], dtype=np.int64)
    c = nodes[tets].mean(axis=1)
    a, b, cc, d = (nodes[tets[:, k]] for k in range(4))
    vol = np.abs(np.einsum("ij,ij->i", b - a, np.cross(cc - a, d - a))) / 6.0
    return c, vol


def gauge_von_mises(vm: np.ndarray, mesh: dict, radius_mm: float) -> np.ndarray:
    """对每个单元，取半径 ``radius_mm`` 内单元的**体积加权平均** σ_vm。

    这是**正则化**：把边界奇异处随网格发散的点峰值，替换成"过程区尺度"上
    的局部平均应力 —— 当半径 > 过程区长度时，该量**与网格无关**，可作强度判据
    （见 LEARNINGS-20261003-023）。返回逐单元的 gauge 值。
    """
    from scipy.spatial import cKDTree

    vm = np.asarray(vm, dtype=np.float64)
    c, vol = element_centroids_and_volumes(mesh)
    if len(c) != len(vm):
        raise ValueError(f"单元数不匹配：vm={len(vm)} vs tets={len(c)}")
    tree = cKDTree(c)
    nb = tree.query_ball_point(c, float(radius_mm))
    out = np.empty(len(vm), dtype=np.float64)
    for i, idx in enumerate(nb):
        w = vol[idx]
        out[i] = float((vm[idx] * w).sum() / w.sum())
    return out


def von_mises_gauge(
    xplt_path,
    mesh: dict,
    *,
    radius_mm: float,
    domain: str = "calcaneus",
    hdf5_path=None,
) -> dict:
    """读结果 → **过程区正则化** von Mises（体积加权、半径 ``radius_mm`` mm）。

    Returns
    -------
    dict
        ``gauge_max`` —— 正则化峰值（**作强度判据**）；
        ``gauge_p99`` / ``gauge_p95``；``raw_max`` —— 未正则化的裸峰值（诊断）；
        ``radius_mm``、``n_elem``。
    """
    from pyfebio import xplt as xplt_mod
    import h5py

    p = Path(xplt_path)
    if not p.is_file():
        raise FileNotFoundError(f"结果文件不存在：{p}")
    h5p = Path(hdf5_path) if hdf5_path is not None else p.with_suffix(".hdf5")
    xplt_mod.to_hdf5(str(p), str(h5p))
    with h5py.File(h5p, "r") as h5:
        stress, _ = _read_last_stress(h5, domain)
    vm = von_mises_from_voigt(stress)
    g = gauge_von_mises(vm, mesh, radius_mm)
    return {
        "gauge_max": float(g.max()),
        "gauge_p99": float(np.percentile(g, 99)),
        "gauge_p95": float(np.percentile(g, 95)),
        "raw_max": float(vm.max()),
        "radius_mm": float(radius_mm),
        "n_elem": int(g.size),
    }


def _write_node_averaged_vtu(path: Path, mesh: dict, stress: np.ndarray, vm: np.ndarray) -> None:
    """把 per-element σ_vm 平均到节点，写 ParaView 可读的 ``.vtu``。"""
    import meshio

    nodes = np.asarray(mesh["nodes"], dtype=np.float64)
    tets = np.asarray(mesh["tets"], dtype=np.int64)
    if len(stress) != len(tets):
        raise ValueError(
            f"单元数不匹配：stress={len(stress)} vs tets={len(tets)}；"
            "请确认 domain 与网格对得上。"
        )
    acc = np.zeros(len(nodes))
    cnt = np.zeros(len(nodes))
    np.add.at(acc, tets.ravel(), np.repeat(vm, 4))
    np.add.at(cnt, tets.ravel(), 1)
    node_vm = acc / np.maximum(cnt, 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(
        str(path),
        meshio.Mesh(
            points=nodes,
            cells=[("tetra", tets)],
            point_data={"von_mises_MPa": node_vm},
            cell_data={"von_mises_MPa": [vm]},
        ),
        file_format="vtu",
    )


def run_and_peak(feb_path, *, workdir=None, domain: str = "calcaneus", **peak_kwargs) -> float:
    """先 :func:`febio_run.run_febio` 求解，再 :func:`peak_von_mises`。"""
    from climbing.coupling.febio_run import run_febio

    feb_path = Path(feb_path)
    run_febio(feb_path, workdir=workdir)
    return peak_von_mises(feb_path, domain=domain, **peak_kwargs)


# ---------------------------------------------------------------------------
# S5 §3.3（T3）位移 / 支反力 / 接触合力
# ---------------------------------------------------------------------------
# 单位：``.feb`` 是 mm–N–MPa–s，故读出的位移 = **mm**，反力 = **N**。
#
# HDF5 布局（``pyfebio.xplt.parse_state`` 写出来的）：
#     states/{state_id}/node_data/<var>/<set_name>
# 其中 ``<var>`` 是 ``PlotDataVariables`` 字符串（"displacement" /
# "reaction forces"）；``<set_name>`` 是节点集名 —— ``pyfebio.xplt`` 把全局
# 节点集固定写成 ``"nodes"``（``XpltTag 0x01041000 → pyname="nodes"``），
# 故 dataset path 在本仓库里实际是 ``displacement/nodes`` 与
# ``reaction forces/nodes``，但 spec 也提到 ``.../all``，所以下面做**完全
# 容忍**（dataset / group 都吃，inner 名是 ``nodes`` / ``all`` / 别的都接受）。
#
# 反应力的 key 还要做大小写/空格容忍（找子串 "reaction"），以防将来 FEBio
# 改名（比如 "Reaction Forces" / "reaction_forces" 之类）。
#
# ``build_hdf5`` 复刻了 ``scripts/opensim_fe/bone_feb._xplt_to_hdf5_patched``
# 的 FEBio 4.13 兼容性修补：FEBio 4.13 在 shell 域上**省略 PLT_DOM_NAME**，
# ``pyfebio.parse_mesh`` 又按 ``domain['name']`` 给 dataset 命名 —— 多个
# 无名 shell 域会撞 ``h5py.create_dataset``。修法：给每个无名域按解析顺序
# 补一个 ``__unnamed_<i+1>`` 名字（同时把 ``id`` 改成位置号，保证
# ``parse_mesh`` / ``parse_state`` 共享同一份 ``mesh_dict['domains']``）。
# 这段逻辑只在 shell 域出错时有用，对纯体域（如本测试）是 no-op。
# ---------------------------------------------------------------------------


def build_hdf5(xplt_path, *, hdf5_path=None) -> Path:
    """``xplt`` → HDF5（FEBio 4.13 兼容版）。

    内部补丁（**只换 ``_parse_domain_section``，不修改 pyfebio 源码**，finally
    里还原）：无名 shell 域 → ``__unnamed_<pos>``。该修补是 S5 在 FEBio 4.13 上
    跑通的必要前置（``l3_shell_feb`` / ``parietal_r`` 已实测）；纯体域模型下
    是 no-op。

    Returns
    -------
    Path
        写出的 ``.hdf5`` 绝对路径。
    """
    import pyfebio.xplt as xp_mod

    p = Path(xplt_path)
    if not p.is_file():
        raise FileNotFoundError(f"xplt 文件不存在：{p}")
    h5p = Path(hdf5_path) if hdf5_path is not None else p.with_suffix(".hdf5")
    h5p.parent.mkdir(parents=True, exist_ok=True)

    original_section = xp_mod._parse_domain_section

    def _patched_domain_section(buf):
        domains = original_section(buf)
        for i, d in enumerate(domains):
            d["id"] = np.array([i + 1], dtype=np.int32)
            nm = d.get("name")
            empty = (not nm) or (
                isinstance(nm, bytes) and not nm.strip(b"\x00")
            ) or (isinstance(nm, str) and nm.strip() == "")
            if empty:
                d["name"] = f"__unnamed_{i + 1}"
        return domains

    xp_mod._parse_domain_section = _patched_domain_section
    try:
        xp_mod.to_hdf5(str(p), str(h5p))
    finally:
        xp_mod._parse_domain_section = original_section
    return h5p


def _resolve_xplt_hdf5(xplt_path, *, hdf5_path=None) -> tuple[Path, Path]:
    """解析 xplt + hdf5 路径，并确保 hdf5 存在（必要时调用 :func:`build_hdf5`）。

    若给 ``.feb`` 或只给了 ``.xplt`` 名但还没求解，会先经
    :func:`febio_run.run_febio` 求解。
    """
    p = Path(xplt_path)
    if p.suffix.lower() == ".feb":
        from climbing.coupling.febio_run import run_febio

        run_febio(p)
        p = p.with_suffix(".xplt")
    elif not p.is_file():
        sibling_feb = p.with_suffix(".feb")
        if sibling_feb.is_file():
            from climbing.coupling.febio_run import run_febio

            run_febio(sibling_feb)
        else:
            raise FileNotFoundError(f"结果文件不存在：{p}（也没有同名 .feb 可求解）")
    if not p.is_file():
        raise FileNotFoundError(f"求解后仍未生成结果文件：{p}")

    h5p = Path(hdf5_path) if hdf5_path is not None else p.with_suffix(".hdf5")
    if not h5p.is_file() or h5p.stat().st_mtime < p.stat().st_mtime:
        build_hdf5(p, hdf5_path=h5p)
    return p, h5p


def _find_node_dataset(node_data_group, var_substring: str) -> tuple[str, "np.ndarray"]:
    """在 ``states/{sid}/node_data`` 下找指定子串的变量（容忍大小写/空格）。

    pyfebio 把 ``<var>`` 写成单独的 dataset（``.../displacement/...``）或一个
    group（``.../displacement/nodes``）—— 都吃；返回第一个匹配。
    """
    import h5py

    target = var_substring.lower().replace(" ", "")
    for name in node_data_group:
        if target in name.lower().replace(" ", ""):
            sub = node_data_group[name]
            if isinstance(sub, h5py.Dataset):
                return name, np.asarray(sub)
            if isinstance(sub, h5py.Group):
                # pyfebio 的实际布局：``<var>/<node_set>``；group 下通常只有一个 dataset
                for inner in sub:
                    return f"{name}/{inner}", np.asarray(sub[inner])
    raise KeyError(
        f"在 node_data 下找不到变量（substring={var_substring!r}）。"
        f"现有变量：{list(node_data_group.keys())}"
    )


def _read_node_var(
    xplt_path,
    var_substring: str,
    *,
    hdf5_path=None,
    last: bool = True,
) -> np.ndarray:
    """读 ``node_data/<var_substring>``。

    Parameters
    ----------
    last:
        True  → 仅末 state，shape ``(N, 3)``；
        False → 所有 state 堆叠，shape ``(S, N, 3)``。
    """
    import h5py

    _, h5p = _resolve_xplt_hdf5(xplt_path, hdf5_path=hdf5_path)
    with h5py.File(h5p, "r") as h5:
        if "states" not in h5:
            raise RuntimeError("HDF5 缺 states group；FEBio 可能未写结果")
        sids = sorted(h5["states"].keys(), key=int)
        if not sids:
            raise RuntimeError("HDF5 states 为空")
        selected = [sids[-1]] if last else sids
        out: list[np.ndarray] = []
        for sid in selected:
            nd = h5[f"states/{sid}/node_data"]
            _, arr = _find_node_dataset(nd, var_substring)
            if arr.ndim != 2 or arr.shape[-1] != 3:
                raise ValueError(
                    f"期望 (N, 3) 形状，得到 {arr.shape}（{var_substring}）"
                )
            out.append(arr.astype(np.float64, copy=False))
    if last:
        return out[0]
    return np.stack(out, axis=0)


def read_displacement(xplt_path, *, hdf5_path=None, last: bool = True) -> np.ndarray:
    """读节点位移 ``(N, 3)`` mm（``last=False`` 时 → ``(S, N, 3)``）。

    接受 ``.feb`` / ``.xplt`` / 已有的 ``.hdf5``；缺 HDF5 时通过
    :func:`build_hdf5` 现场补建。
    """
    return _read_node_var(xplt_path, "displacement", hdf5_path=hdf5_path, last=last)


def read_reaction_forces(xplt_path, *, hdf5_path=None, last: bool = True) -> np.ndarray:
    """读节点反力 ``(N, 3)`` N（``last=False`` 时 → ``(S, N, 3)``）。

    Key 匹配容忍大小写 / 空格：找 ``node_data`` 下任何名字含子串
    ``"reaction"`` 的变量。这样无论 FEBio 写 ``"reaction forces"`` /
    ``"ReactionForces"`` / ``"reaction_forces"`` 都能命中。
    """
    return _read_node_var(xplt_path, "reaction", hdf5_path=hdf5_path, last=last)


def contact_force_sum(
    xplt_path, *, hdf5_path=None, last: bool = True
) -> np.ndarray:
    """反应力向量在所有节点上的求和 → ``(3,)`` N。

    即 **刚性支承 / 接触面对体的合力**（S5 G3 守恒检查要的就是这个）。
    FEBio 4.13 的反力符号约定：

    * 约束节点 → 反应力 = **约束力本身**（支承推 / 拉体的力）。
    * 被加载节点 → 反应力 = **−施加荷载**（FEBio 在内部为平衡荷载所施加的力）。

    因此对自平衡系统（含 ``BCZeroDisplacement`` 满约束）的求和会整体抵消，
    G3 应比较 **``|contact_force_sum|`` 与 ``|applied load|``**（或
    ``contact_force_sum`` 与 ``-applied load`` 的分量 —— 取决于看哪个方向）。

    Parameters
    ----------
    last:
        True  → 末 state，shape ``(3,)``；
        False → 所有 state，shape ``(S, 3)``。
    """
    rx = read_reaction_forces(xplt_path, hdf5_path=hdf5_path, last=last)
    return rx.sum(axis=-2)
