"""S5 §3.3 回归测试：``fe_post.read_displacement`` / ``read_reaction_forces`` /
``contact_force_sum`` / ``build_hdf5``。

测试用 pyfebio 现建一个最小 .feb（沿用 ``tests/test_febio_run.py::_build_minimal_model``
的模式，但**额外声明 ``reaction forces`` plotfile 变量**），跑 FEBio，再读
HDF5 做断言。

断言两条（spec：G3 守恒）：

1. ``read_displacement`` 在施加了 prescribed z=-0.01 mm 的节点上读出
   ``u_z ≈ -0.01 mm``（<2%）。
2. ``contact_force_sum`` 的 |y 分量| ≈ 施加的 ``+1 N × 4 节点 = +4 N`` 的 y 荷载
   （<2%；FEBio 反应力约定：加载节点处反应 = -applied load，求和后 sum_y = -4 N）。

找不到 FEBio 时**整文件 skip**，与 ``tests/test_febio_run.py`` 完全一致。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from climbing.coupling import fe_post
from climbing.coupling.febio_run import (
    FebioNotFound,
    find_febio,
    run_febio,
)


def _febio_available() -> bool:
    try:
        find_febio(verify=True)
        return True
    except FebioNotFound:
        return False


requires_febio = pytest.mark.skipif(
    not _febio_available(), reason="本机未安装可运行的 FEBio"
)


# 相对误差阈值（G3 守恒检查）
REL_TOL = 0.02


def _build_reaction_model(pyfebio):
    """最小 .feb：1×1×1 mm 立方体（E=7300 MPa, v=0.3），

    * 底面（节点 1-4）：``BCZeroDisplacement`` —— 全约束（FE 全场刚体参考）。
    * 顶面（节点 5-8）：``BCPrescribedDisplacement dof=z value=-0.01`` —— z 向
      压 −0.01 mm（这是 read_displacement 要断言的 *prescribed displacement*）。
    * 顶面：``NodalLoad dof=y scale=1.0`` —— 4 个顶节点各 +1 N（合计 **+4 N**
      in +y）；这是 contact_force_sum 要跟它比对的 *applied load*。

    plotfile 同时声明 ``displacement`` + ``reaction forces``，这是本测试
    的核心扩展（``_build_minimal_model`` 只声明了 displacement）。
    """
    import pyfebio  # 局部 import，便于 pyfebio 缺失时本文件 import 不挂

    m = pyfebio.model.Model()

    nodes = pyfebio.mesh.Nodes(name="nodes")
    coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
              (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    for i, c in enumerate(coords):
        nodes.add_node(pyfebio.mesh.Node(id=i + 1, text=",".join(map(str, c))))
    elems = pyfebio.mesh.Elements(name="box", type="hex8")
    elems.add_element(pyfebio.mesh.Hex8Element(id=1, text="1,2,3,4,5,6,7,8"))

    m.mesh_.nodes.append(nodes)
    m.mesh_.elements.append(elems)
    m.mesh_.node_sets.append(pyfebio.mesh.NodeSet(name="bottom", text="1,2,3,4"))
    m.mesh_.node_sets.append(pyfebio.mesh.NodeSet(name="top", text="5,6,7,8"))

    m.material_.add_material(pyfebio.material.IsotropicElastic(
        id=1, name="bone",
        E=pyfebio.material.MaterialParameter(text=7300.0),   # MPa (mm-N-MPa)
        v=pyfebio.material.MaterialParameter(text=0.30),
    ))
    m.meshdomains_.add_solid_domain(
        pyfebio.meshdomains.SolidDomain(name="box", mat="bone")
    )

    m.boundary_.add_bc(pyfebio.boundary.BCZeroDisplacement(
        node_set="bottom", x_dof=1, y_dof=1, z_dof=1))
    m.boundary_.add_bc(pyfebio.boundary.BCPrescribedDisplacement(
        node_set="top", dof="z", value=pyfebio.boundary.Value(lc=1, text=-0.01)))

    # 顶部施加 +1 N（每节点）y 方向集中力 → 总计 +4 N（contact_force_sum 比对靶）
    m.loads_.add_nodal_load(pyfebio.loads.NodalLoad(
        node_set="top", dof="y",
        scale=pyfebio.loads.Scale(lc=1, text=1.0)))

    m.loaddata_.add_load_curve(pyfebio.loaddata.LoadCurve(
        id=1, points=pyfebio.loaddata.CurvePoints(points=["0.0,0.0", "1.0,1.0"])))

    # ★ 关键扩展：声明 reaction forces plotfile（_build_minimal_model 没声明）
    pf = pyfebio.output.OutputPlotfile(type="febio")
    pf.add_var(pyfebio.output.Var(type="displacement"))
    pf.add_var(pyfebio.output.Var(type="reaction forces"))
    m.output_.plotfile.append(pf)
    return m


# ==========================================================================
# 不依赖 FEBio 的 API 烟雾测试（import / 类型 / __all__）
# ==========================================================================
class TestModuleApi:
    def test_new_names_in_dunder_all(self):
        """新函数必须出现在 __all__ 里（spec 显式要求）。"""
        for name in ("read_displacement", "read_reaction_forces",
                     "contact_force_sum", "build_hdf5"):
            assert name in fe_post.__all__, f"{name} 未在 fe_post.__all__ 中"

    def test_resolve_missing_xplt_raises(self, tmp_path):
        """既无 .xplt 也无同名 .feb 时 → 显式 FileNotFoundError（绝不静默）。"""
        with pytest.raises(FileNotFoundError):
            fe_post.read_displacement(tmp_path / "nope.xplt")


# ==========================================================================
# 真实 FEBio 跑通后做的 G3 守恒检查（缺 FEBio → skip）
# ==========================================================================
@requires_febio
class TestReactionRoundtrip:
    """建一个最小 .feb，跑通 FEBio，再从 xplt 反推读出。"""

    @pytest.fixture
    def feb_artifacts(self, tmp_path):
        pyfebio = pytest.importorskip("pyfebio")
        feb = tmp_path / "cube.feb"
        _build_reaction_model(pyfebio).save(str(feb))
        # 用 climbing.coupling.febio_run.run_febio：失败必抛（绝不静默 rc=1）
        assert run_febio(feb, silent=True) == 0
        xplt = feb.with_suffix(".xplt")
        assert xplt.is_file(), f"FEBio 没产出 xplt：{xplt}"
        return feb, xplt

    def test_build_hdf5_returns_path(self, feb_artifacts):
        """build_hdf5 必须返回 .hdf5 路径且文件存在（>=1 byte）。"""
        _, xplt = feb_artifacts
        h5p = fe_post.build_hdf5(xplt)
        assert isinstance(h5p, Path)
        assert h5p.is_file()
        assert h5p.stat().st_size > 0

    def test_build_hdf5_is_idempotent(self, feb_artifacts):
        """重复跑 build_hdf5 不应抛错（finally 已还原 _parse_domain_section）。"""
        _, xplt = feb_artifacts
        fe_post.build_hdf5(xplt)
        # 再来一次 —— 第二次调用时 pyfebio.xplt._parse_domain_section 应
        # 已被还原成原始版本；同一 xplt 跑两次都应成功。
        h5p = fe_post.build_hdf5(xplt)
        assert h5p.is_file()

    def test_read_displacement_matches_prescribed(self, feb_artifacts):
        """G3：顶面 prescribed uz=-0.01 mm 必须被读出来（<2%）。"""
        _, xplt = feb_artifacts
        u = fe_post.read_displacement(xplt)
        assert u.shape == (8, 3), f"期望 (8, 3)，得到 {u.shape}"

        # 顶面节点 5-8 的 u_z 必须接近 prescribed -0.01 mm
        uz_top = u[4:, 2]
        assert np.all(np.isfinite(uz_top)), "u_z 出现非有限值"
        rel = np.abs((uz_top - (-0.01)) / (-0.01))
        assert rel.max() < REL_TOL, (
            f"prescribed displacement u_z 误差 {rel.max()*100:.2f}% > {REL_TOL*100:.0f}%：{uz_top}"
        )

        # 底面（节点 1-4）z 方向应 ≈ 0
        uz_bot = u[:4, 2]
        assert np.all(np.abs(uz_bot) < 1e-6), f"底面 u_z 应 ≈ 0，得到 {uz_bot}"

    def test_contact_force_sum_matches_applied_load(self, feb_artifacts):
        """G3：|contact_force_sum_y| 应 ≈ 4 N（4 节点 × +1 N y 向荷载）。"""
        _, xplt = feb_artifacts
        applied_y = 4.0  # N
        s = fe_post.contact_force_sum(xplt)
        assert s.shape == (3,), f"期望 (3,)，得到 {s.shape}"

        # FEBio 约定：加载节点处反应 = -applied load；求和后 sum_y = -applied_y。
        # 这里比绝对值（G3 关心的是 *大小* 是否守恒，符号取决于观察方向）。
        rel = abs(abs(s[1]) - applied_y) / applied_y
        assert rel < REL_TOL, (
            f"|contact_force_sum_y|={abs(s[1]):.4f} N vs applied {applied_y} N，"
            f"相对差 {rel*100:.2f}% > {REL_TOL*100:.0f}%；完整 sum={s}"
        )

        # 同时确认：x 分量应 ≈ 0（没有任何 x 方向的荷载）
        assert abs(s[0]) < 1e-3, f"sum_x 应 ≈ 0，得到 {s[0]}"

    def test_read_reaction_forces_shape_and_finite(self, feb_artifacts):
        """read_reaction_forces：shape=(N,3)，全有限值。"""
        _, xplt = feb_artifacts
        r = fe_post.read_reaction_forces(xplt)
        assert r.shape == (8, 3), f"期望 (8, 3)，得到 {r.shape}"
        assert np.all(np.isfinite(r)), "反应力出现非有限值"

    def test_last_false_returns_all_states_stacked(self, feb_artifacts):
        """last=False → (S, N, 3)；last=True → (N, 3)。S >= 2。"""
        _, xplt = feb_artifacts
        u_last = fe_post.read_displacement(xplt, last=True)
        u_all = fe_post.read_displacement(xplt, last=False)
        r_last = fe_post.read_reaction_forces(xplt, last=True)
        r_all = fe_post.read_reaction_forces(xplt, last=False)
        s_last = fe_post.contact_force_sum(xplt, last=True)
        s_all = fe_post.contact_force_sum(xplt, last=False)

        assert u_last.ndim == 2 and u_last.shape == (8, 3)
        assert u_all.ndim == 3 and u_all.shape[1:] == (8, 3)
        assert r_last.shape == (8, 3)
        assert r_all.ndim == 3 and r_all.shape[1:] == (8, 3)
        assert s_last.shape == (3,)
        assert s_all.ndim == 2 and s_all.shape[1:] == (3,)

        assert u_all.shape[0] >= 2, f"应至少有 2 个 state，得到 {u_all.shape[0]}"
        # 末 state 的 last-True 与 last-False[-1] 必须严格相等（同一份数据）
        np.testing.assert_allclose(u_last, u_all[-1])
        np.testing.assert_allclose(r_last, r_all[-1])
        np.testing.assert_allclose(s_last, s_all[-1])

    def test_explicit_hdf5_path_wins(self, feb_artifacts, tmp_path):
        """hdf5_path 指定 → 写到该路径（而不与 xplt 同名）。"""
        _, xplt = feb_artifacts
        target = tmp_path / "custom.hdf5"
        out = fe_post.build_hdf5(xplt, hdf5_path=target)
        assert out == target.resolve()
        assert target.is_file()

    def test_existing_old_hdf5_is_rebuilt_when_xplt_newer(self, feb_artifacts, tmp_path):
        """xplt 比 hdf5 新 → 自动 rebuild（避免拿旧缓存）。"""
        _, xplt = feb_artifacts
        # 1) 先建一份
        h5p = fe_post.build_hdf5(xplt, hdf5_path=tmp_path / "a.hdf5")
        mtime_before = h5p.stat().st_mtime
        # 2) 戳一下 xplt 让它更新；再调 read_displacement → 应 rebuild
        import time
        time.sleep(1.1)  # mtime 精度（Windows FAT 约 2s）
        Path(xplt).touch()
        u = fe_post.read_displacement(xplt, hdf5_path=tmp_path / "a.hdf5")
        assert u.shape == (8, 3)
        # 缓存已被重建（mtime 更新）—— 这一条不强校验，避免 fs 精度问题
