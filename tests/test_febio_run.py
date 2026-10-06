"""S0.2 回归测试：FEBio 定位与安全调用。

锁死的核心失效：FEBio 不在 PATH 时，``pyfebio.model.run_model()``
静默返回 rc=1 而不报错 —— 见 ``.learnings/ERRORS.md [ERR-...-007]``。
本模块的用例保证"找不到 / 跑挂"一定会**抛错**。

测试模型用 pyfebio 现建（不手写 XML）—— 避免 FEBio 4.x 的 ``<Control>``
schema 漂移（手写版会因缺 ``<solver>`` 直接解析失败）。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from climbing.coupling.febio_run import (
    DEFAULT_CANDIDATES,
    ENV_VAR,
    FebioNotFound,
    FebioRunError,
    find_febio,
    probe_febio,
    run_febio,
)


def _febio_available() -> bool:
    try:
        find_febio(verify=True)
        return True
    except FebioNotFound:
        return False


#: 需要真实 FEBio 的用例统一跳过条件（本机已装则应全部执行）
requires_febio = pytest.mark.skipif(
    not _febio_available(), reason="本机未安装可运行的 FEBio"
)


def _build_minimal_model(pyfebio):
    """最小可解模型：单个 hex8 立方体，底面固定、顶面压 1%（弹性骨）。"""
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

    m.loaddata_.add_load_curve(pyfebio.loaddata.LoadCurve(
        id=1, points=pyfebio.loaddata.CurvePoints(points=["0.0,0.0", "1.0,1.0"])))
    return m


# ==========================================================================
# find_febio
# ==========================================================================
class TestFindFebio:
    def test_explicit_bogus_path_raises(self):
        """显式给了不可用路径 -> 立即抛错，不静默回退到别的候选。"""
        with pytest.raises(FebioNotFound):
            find_febio(r"Z:\definitely\not\here\febio4.exe")

    def test_env_var_override_wins(self, monkeypatch, tmp_path):
        """FEBIO_EXE 应优先于默认候选（哪怕文件只是存在、未校验）。"""
        fake = tmp_path / "febio4.exe"
        fake.write_bytes(b"")
        monkeypatch.setenv(ENV_VAR, str(fake))
        assert find_febio() == fake.resolve()

    @requires_febio
    def test_finds_known_default_path(self):
        exe = find_febio(verify=True)
        assert exe.is_file()
        assert any(
            exe == Path(c).resolve()
            for c in DEFAULT_CANDIDATES
            if Path(c).is_file()
        )

    @requires_febio
    def test_probe_returns_version(self):
        ver = probe_febio(find_febio())
        assert ver and ver != "unknown"
        assert ver.split(".")[0].isdigit()

    def test_not_found_message_mentions_env_var(self, monkeypatch, tmp_path):
        """错误信息要能指路（提示 FEBIO_EXE / PATH）。"""
        monkeypatch.delenv(ENV_VAR, raising=False)
        with pytest.raises(FebioNotFound) as ei:
            find_febio(tmp_path / "nope.exe")
        assert ENV_VAR in str(ei.value)


# ==========================================================================
# run_febio
# ==========================================================================
class TestRunFebio:
    def test_missing_input_raises(self, tmp_path):
        with pytest.raises((FileNotFoundError, FebioNotFound)):
            run_febio(tmp_path / "does_not_exist.feb")

    @requires_febio
    def test_invalid_feb_raises_not_silent(self, tmp_path):
        """回归核心：坏输入必须抛 FebioRunError，而不是静默返回 rc=1。

        这正是 pyfebio.run_model() 的反例 —— 它会把 rc=1 当普通返回值吞掉。
        """
        bad = tmp_path / "bad.feb"
        bad.write_text("this is not a valid febio model\n")
        with pytest.raises(FebioRunError):
            run_febio(bad)

    @requires_febio
    def test_valid_model_returns_zero(self, tmp_path):
        """一个最小合法模型应成功求解并返回 0，且产出 .xplt。"""
        pyfebio = pytest.importorskip("pyfebio")
        feb = tmp_path / "unit.feb"
        _build_minimal_model(pyfebio).save(str(feb))

        assert run_febio(feb, silent=True) == 0
        assert (tmp_path / "unit.xplt").is_file()
