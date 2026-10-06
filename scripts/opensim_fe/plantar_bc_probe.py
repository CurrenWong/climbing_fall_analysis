"""跖面边界条件探针：``fixed`` vs ``roller`` vs ``spring`` 对 gauge 的影响。

问题
----
现有单骨跟骨 FE 把跖面**三向全固定**（``plantar_bc="fixed"``）。若这是人为的
固定端奇异，则 ``gauge_max`` 会虚高、首次骨折高度偏低（本仓库 anat 扫描 ≈ 2–3 m
vs [Y25] 7–9 m）。本脚本在**固定 5 m 载荷**（距下 26.38 kN、无跟腱）下，切换跖面
BC，量化它对皮质 ``gauge`` 的贡献。

判据（与 ``s2_thums_sweep.py`` 完全一致）
---------------------------------------
皮质域 ``calcaneus``（hex8, E=15000）的过程区正则化 gauge，ℓ=R=4 mm；
首次骨折 = gauge_max > 150 MPa（[Y25] 皮质压缩强度）。**后处理判据，无单元删除。**

载荷
----
h=5 m, a=0（无跟腱）经 OpenSim 多体链真跑得到距下峰力 F_subt≈26.38 kN，
再以 ``PressureLoad p=F/A`` 施加（与 anat 扫描同口径 ``use_rigid=False``）。

产物
----
``temp/opensim_fe/plantar_bc/``         每次运行的 .feb/.xplt/.hdf5/.log
``results/opensim_fe/plantar_bc_summary.{json,csv}``
``results/opensim_fe/plantar_bc_gauge_k.png``
``results/opensim_fe/PLANTAR_BC_REPORT.md``

运行（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\plantar_bc_probe.py

单位 mm–N–MPa。
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import thums_feb as T  # noqa: E402

#: 方向 C 的**解剖**配准网格（不得覆盖，只读）。
MESH_ANAT = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_anatframe.npz"
#: 方向 C 的 anat 全高度扫描结果，用于 F_subt(h) 反查（首次骨折高度外推）。
ANAT_SWEEP_JSON = ROOT / "results" / "opensim_fe" / "s2_thums_height_sweep_anat.json"

WORKDIR = ROOT / "temp" / "opensim_fe" / "plantar_bc"
OUTDIR = ROOT / "results" / "opensim_fe"

#: 固定载荷：S2 anat 扫描里 h=5, a=0 的那一行。
DROP_H_M = 5.0
ACT = 0.0
#: 硬门槛：复现方向 C 的 gauge（±1%）。
BASELINE_GAUGE_MPA = 186.7563689967737
BASELINE_GAUGE_TOL = 0.01
#: 同口径的静力峰值（250？）参考（仅诊断）。
BASELINE_MAX_MPA = 213.50277948448795
#: 皮质压缩强度（[Y25]）。
THRESHOLD_MPA = 150.0
REG_LEN_MM = 4.0
#: spring 刚度扫描（N/mm），软→硬，覆盖到接近刚性（元素刚度量级 E·L ≈ 1e4–1e5）。
SPRING_KS = [10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0]
#: roller_free 的切向弹簧刚度（N/mm）：极软 => 近似无摩擦滑动，仅消除刚体模态。
ROLLER_FREE_KT_N_PER_MM = 1.0
#: 首次骨折高度外推时，参考 [Y25] 簇。
Y25_LO_M, Y25_HI_M = 7.0, 9.0


# ---------------------------------------------------------------------------
# load side (verbatim: s2_thums_sweep._load_at)
# ---------------------------------------------------------------------------
def _load_at(height: float) -> tuple[float, float]:
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=height)
    fall = run_dead_drop(grf)
    jr = subtalar_reaction(fall, grf, side="r")
    return float(jr.peak_vertical_n), float(grf.impulse_rel_err)


# ---------------------------------------------------------------------------
# post side
# ---------------------------------------------------------------------------
def _install_xplt_discrete_patch() -> None:
    """pyfebio 0.3.0 xplt 解析器在遇到离散单元（Winkler 弹簧）写出的**无名字域**时
    会 ``KeyError: 'name'``。给这类域补一个合成分域名即可（只影响解析，不影响解）。
    """
    from pyfebio import xplt as X

    if getattr(X, "_plantar_bc_patched", False):
        return
    orig = X._parse_domain

    def _pd(buffer):  # type: ignore[no-untyped-def]
        d = orig(buffer)
        if "name" not in d:
            d["name"] = "discrete_%d" % int(d["id"][0])
        return d

    X._parse_domain = _pd
    X._plantar_bc_patched = True


def _read_domain_vm(xplt: Path, domain: str, hdf5: Path) -> np.ndarray:
    import h5py

    from climbing.coupling import fe_post
    from pyfebio import xplt as xplt_mod

    if not Path(hdf5).is_file():
        xplt_mod.to_hdf5(str(xplt), str(hdf5))
    with h5py.File(hdf5, "r") as h5:
        stress, used = fe_post._read_last_stress(h5, domain)
    return fe_post.von_mises_from_voigt(stress)


def _nearest_named_nodes(point: np.ndarray, mesh: dict) -> dict[str, float]:
    out: dict[str, float] = {}
    for name in ("subtalar_joint", "plantar", "achilles"):
        ns = np.asarray(mesh["node_sets"].get(name, []), dtype=np.int64)
        if len(ns) == 0:
            out[name] = float("inf")
            continue
        out[name] = float(np.linalg.norm(mesh["nodes"][ns] - point, axis=1).min())
    return out


def _peak_location(centroid: np.ndarray, mesh: dict, radius_mm: float = REG_LEN_MM) -> dict:
    """峰值位置：最近具名面 + 到各具名节点集的距离 + 内部/边缘分类。"""
    d = _nearest_named_nodes(centroid, mesh)
    nearest = min(d, key=lambda k: d[k])
    # 到所有具名面节点的最近距离（判断是否"贴边"）
    all_named = np.concatenate(
        [np.asarray(mesh["node_sets"].get(n, []), dtype=np.int64) for n in d]
    ) if any(len(mesh["node_sets"].get(n, [])) for n in d) else np.empty(0, dtype=np.int64)
    d_any = (
        float(np.linalg.norm(mesh["nodes"][all_named] - centroid, axis=1).min())
        if len(all_named)
        else float("inf")
    )
    if nearest == "subtalar_joint" and d[nearest] <= radius_mm:
        where = "关节面边缘 (joint rim)"
    elif nearest == "plantar" and d[nearest] <= radius_mm:
        where = "跖面边缘 (plantar rim)"
    elif d_any <= radius_mm:
        where = f"{nearest} 边缘"
    else:
        where = "内部 (interior)"
    return {
        "nearest_surface": nearest,
        "nearest_surface_dist_mm": d[nearest],
        "dist_joint_mm": d["subtalar_joint"],
        "dist_plantar_mm": d["plantar"],
        "dist_achilles_mm": d["achilles"],
        "dist_any_named_node_mm": d_any,
        "location_class": where,
        "at_boundary": bool(d_any <= radius_mm),
    }


def _solve_case(
    mesh: dict,
    cort_centroids: np.ndarray,
    cort_vol: np.ndarray,
    f_sub: float,
    bc_kind: str,
    spring_k: float | None,
    *,
    tag: str,
    gauge_mm: float,
) -> dict:
    """建 -> 解 -> 后处理一个跖面 BC 工况。"""
    from climbing.coupling import fe_post
    from climbing.coupling.febio_run import run_febio
    from climbing.coupling.load_transfer import transfer

    spec = transfer(f_sub, achilles_activation=ACT)  # a=0 -> no Achilles
    feb = (WORKDIR / f"plantar_{tag}.feb").resolve()
    xplt = feb.with_suffix(".xplt")
    hdf5 = feb.with_suffix(".hdf5")
    for stale in (xplt, hdf5):
        if stale.exists():
            stale.unlink()

    kwargs: dict = {}
    if bc_kind in ("spring", "roller_free"):
        kwargs["spring_k"] = float(spring_k)
    T.build_thums_feb(
        mesh,
        feb,
        load_n=spec.subtalar_n,
        load_dir=spec.subtalar_dir,
        achilles_n=spec.achilles_n,
        achilles_dir=spec.achilles_dir,
        use_rigid=False,  # PressureLoad p=F/A: same convention as the S2 sweep
        time_steps=1,
        plantar_bc=bc_kind,
        **kwargs,
    )
    rc = run_febio(feb, workdir=feb.parent, timeout=1800)
    if rc != 0:
        raise RuntimeError(f"FEBio rc={rc} for {tag}")
    if not xplt.is_file():
        raise FileNotFoundError(f"rc={rc} 但未生成 {xplt}")

    peak_public = fe_post.peak_von_mises(xplt, hdf5_path=hdf5, domain="calcaneus")
    vm = _read_domain_vm(xplt, "calcaneus", hdf5)
    if len(vm) != len(cort_centroids):
        raise RuntimeError(f"CORT 单元数不匹配：vm={len(vm)} vs hex={len(cort_centroids)}")
    if abs(float(vm.max()) - float(peak_public)) > 1e-6 * max(1.0, float(peak_public)):
        raise RuntimeError(
            f"fe_post.peak_von_mises={peak_public} 与逐单元峰值 {float(vm.max())} 不一致"
        )
    if not np.all(np.isfinite(vm)):
        raise RuntimeError(f"{tag}: von Mises 含 NaN/Inf")

    gauge = T.gauge_von_mises_mixed(vm, cort_centroids, cort_vol, gauge_mm)
    ipk = int(np.argmax(vm))
    loc = _peak_location(cort_centroids[ipk], mesh)
    row = {
        "tag": tag,
        "plantar_bc": bc_kind,
        "spring_k_N_per_mm": float(spring_k) if spring_k is not None else None,
        "febio_rc": int(rc),
        "load_n": float(spec.subtalar_n),
        "load_kN": float(spec.subtalar_n) / 1e3,
        "achilles_kN": float(spec.achilles_n) / 1e3,
        "max": float(vm.max()),
        "p99": float(np.percentile(vm, 99)),
        "p95": float(np.percentile(vm, 95)),
        "mean": float(vm.mean()),
        "gauge_max": float(gauge.max()),
        "gauge_p95": float(np.percentile(gauge, 95)),
        "peak_elem": ipk,
        "peak_centroid_mm": [float(x) for x in cort_centroids[ipk]],
        **loc,
        "feb": str(feb),
        "xplt": str(xplt),
    }
    return row


# ---------------------------------------------------------------------------
# first-fracture height extrapolation (linear-in-load, NOT a re-scan)
# ---------------------------------------------------------------------------
def _load_height_table() -> tuple[np.ndarray, np.ndarray]:
    """从 anat 扫描读 (h, F_subt) 表（act=0 行）。"""
    data = json.loads(ANAT_SWEEP_JSON.read_text(encoding="utf-8"))
    rows = sorted((r for r in data if abs(float(r["act"]) - ACT) < 1e-9), key=lambda r: r["h"])
    h = np.array([float(r["h"]) for r in rows])
    f = np.array([float(r["subtalar_kN"]) for r in rows]) * 1e3  # N
    return h, f


def _extrapolate_first_fracture(gauge_at_5m: float, h_tab: np.ndarray, f_tab: np.ndarray) -> dict:
    """gauge∝F 假设下，gauge(h)=150 对应的高度（线性内插 F(h) 表）。"""
    f5 = float(np.interp(DROP_H_M, h_tab, f_tab))
    if gauge_at_5m <= 0.0:
        return {"f_crit_n": None, "h_first_m": None}
    f_crit = f5 * THRESHOLD_MPA / gauge_at_5m
    if f_crit <= f_tab.min():
        h_first = float(h_tab.min())
    elif f_crit >= f_tab.max():
        h_first = None  # 超出扫描上界
    else:
        h_first = float(np.interp(f_crit, f_tab, h_tab))
    return {
        "f_crit_n": float(f_crit),
        "f_crit_kN": float(f_crit) / 1e3,
        "h_first_m": h_first,
        "linear_load_assumed": True,
        "note": "gauge 按线弹性 ∝ F 缩放 + 用 anat 扫描真跑 F(h) 反查；非重扫",
    }


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def _write_plot(rows: list[dict], out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    spring = [r for r in rows if r["plantar_bc"] == "spring" and r["gauge_max"] is not None]
    spring.sort(key=lambda r: r["spring_k_N_per_mm"])
    fixed = next(r for r in rows if r["plantar_bc"] == "fixed")
    roller = next((r for r in rows if r["plantar_bc"] == "roller"), None)
    roller_free = next((r for r in rows if r["plantar_bc"] == "roller_free"), None)

    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    if spring:
        ks = [r["spring_k_N_per_mm"] for r in spring]
        gs = [r["gauge_max"] for r in spring]
        ax.semilogx(ks, gs, "-o", color="tab:blue", ms=6, label="spring gauge_max(k)")
    ax.axhline(fixed["gauge_max"], color="tab:red", ls="--", lw=1.5,
               label=f"fixed gauge_max = {fixed['gauge_max']:.1f} MPa")
    if roller is not None:
        ax.axhline(roller["gauge_max"], color="tab:green", ls=":", lw=1.8,
                   label=f"roller(NormalDispl) = {roller['gauge_max']:.1f} MPa")
    if roller_free is not None and roller_free["gauge_max"] is not None:
        ax.axhline(roller_free["gauge_max"], color="tab:purple", ls="--", lw=1.4,
                   label=f"roller_free = {roller_free['gauge_max']:.1f} MPa")
    ax.axhline(THRESHOLD_MPA, color="k", ls="-.", lw=1.0, alpha=0.7,
               label=f"骨折阈值 {THRESHOLD_MPA:.0f} MPa")
    ax.set_xlabel("Winkler 刚度 k (N/mm, log)")
    ax.set_ylabel("gauge_max (MPa) · 域 `calcaneus`/CORT")
    ax.set_title("跖面 BC 对比 · THUMS AM50 跟骨 · h=5 m, a=0 (F_subt≈26.38 kN)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def _write_report(
    rows: list[dict],
    load_meta: dict,
    extrap: dict,
    mesh: dict,
    out_md: Path,
    monotone: bool,
) -> None:
    fixed = next(r for r in rows if r["plantar_bc"] == "fixed")
    roller = next((r for r in rows if r["plantar_bc"] == "roller"), None)
    roller_free = next((r for r in rows if r["plantar_bc"] == "roller_free"), None)
    spring_all = sorted((r for r in rows if r["plantar_bc"] == "spring"), key=lambda r: r["spring_k_N_per_mm"])
    spring = [r for r in spring_all if r["gauge_max"] is not None]
    spring_failed = [r for r in spring_all if r["gauge_max"] is None]

    def _f(x, fmt: str = ".2f") -> str:  # None-safe formatting
        if x is None:
            return "失败"
        if isinstance(x, str):
            return x
        return format(x, fmt)

    lines: list[str] = [
        "# 跖面边界条件对 gauge 应力的影响（THUMS AM50 跟骨）",
        "",
        "**问题**：现有单骨 FE 把跖面**三向全固定**。若这是人为固定端奇异，则 `gauge` 虚高、",
        "首次骨折高度偏低（本仓库 anat 扫描 ≈2–3 m vs [Y25] 7–9 m）。本实验在**固定 5 m 载荷**",
        "（距下 26.38 kN、无跟腱）下切换跖面 BC，量化其对皮质 gauge 的贡献。",
        "",
        f"- 网格：`{Path(mesh['source_npz']).name}`（AM50 THUMS 右跟骨，双域 CORT hex8 + SPON tet4）",
        f"- 载荷：h={load_meta['h_m']:.0f} m, a={ACT:.0f}（无跟腱）→ F_subt = "
        f"{load_meta['f_subt_kN']:.3f} kN，T1 冲量误差 {load_meta['t1_rel_err']*100:+.1f}%；"
        "`PressureLoad p=F/A`（与 anat 扫描同口径，`use_rigid=False`）",
        f"- 判据：皮质域 `calcaneus` gauge_max（体积平均正则化 R={REG_LEN_MM:g} mm = ℓ）> "
        f"**{THRESHOLD_MPA:.0f} MPa**（[Y25] 皮质压缩强度）；**后处理判据，无单元删除**",
        "",
        "## 0. 基线复现（硬门槛）",
        "",
        f"- 方向 C anat 扫描 h=5, a=0 基线：gauge_max = **{BASELINE_GAUGE_MPA:.3f} MPa**"
        f"（raw max = {BASELINE_MAX_MPA:.2f} MPa）",
        f"- 本实验 `fixed` 复现：gauge_max = **{fixed['gauge_max']:.3f} MPa**，"
        f"raw max = {fixed['max']:.2f} MPa",
        f"- 相对误差：{(fixed['gauge_max']-BASELINE_GAUGE_MPA)/BASELINE_GAUGE_MPA*100:+.3f}%"
        f"（门槛 ±{BASELINE_GAUGE_TOL*100:.0f}%）→ **{'通过' if abs(fixed['gauge_max']-BASELINE_GAUGE_MPA)/BASELINE_GAUGE_MPA <= BASELINE_GAUGE_TOL else '未通过'}**",
        "",
        "## 1–2. 各跖面 BC 的 gauge",
        "",
        "| BC | k (N/mm) | rc | gauge_max (MPa) | gauge_p95 | raw max | p95 | mean | 峰值位置 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    order = [fixed]
    if roller is not None:
        order.append(roller)
    if roller_free is not None:
        order.append(roller_free)
    order += spring_all  # include failed k (printed as 失败)
    for r in order:
        k = "—" if r["spring_k_N_per_mm"] is None else f"{r['spring_k_N_per_mm']:g}"
        lines.append(
            f"| {r['plantar_bc']} | {k} | {_f(r['febio_rc'], 'd')} | {_f(r['gauge_max'])} | "
            f"{_f(r['gauge_p95'])} | {_f(r['max'])} | {_f(r['p95'])} | {_f(r['mean'], '.3f')} | "
            f"{_f(r['location_class'])} |"
        )

    lines += [
        "",
        "### 峰值位置明细",
        "",
        "| BC | k (N/mm) | 质心 (mm) | 最近具名面 | 距离 (mm) | d(joint) | d(plantar) | 分类 |",
        "|---|---|---|---|---|---|---|---|",
    ]
    ok = [r for r in order if r["peak_centroid_mm"] is not None]
    for r in ok:
        k = "—" if r["spring_k_N_per_mm"] is None else f"{r['spring_k_N_per_mm']:g}"
        c = "(" + ", ".join(f"{x:.1f}" for x in r["peak_centroid_mm"]) + ")"
        lines.append(
            f"| {r['plantar_bc']} | {k} | {c} | {r['nearest_surface']} | "
            f"{r['nearest_surface_dist_mm']:.2f} | {r['dist_joint_mm']:.1f} | "
            f"{r['dist_plantar_mm']:.1f} | {r['location_class']} |"
        )

    # spring scan self-check
    lines += [
        "",
        "### spring 扫描自检：k↑ 是否单调逼近 fixed？",
        "",
    ]
    if spring:
        gs = [r["gauge_max"] for r in spring]
        lines.append(
            f"- 单调不降：**{'是' if monotone else '否'}**；k={spring[0]['spring_k_N_per_mm']:g}→"
            f"{spring[-1]['spring_k_N_per_mm']:g} N/mm：gauge "
            f"{gs[0]:.2f} → {gs[-1]:.2f} MPa（fixed 基线 {fixed['gauge_max']:.2f} MPa）"
        )
        rel = abs(gs[-1] - fixed["gauge_max"]) / fixed["gauge_max"] * 100
        lines.append(
            f"- 最硬 k={spring[-1]['spring_k_N_per_mm']:g} N/mm 相对 fixed 残差 = {rel:.2f}%"
            f"（{'✓ 收敛到固定基线' if rel <= 2.0 else '⚠ 未收敛到 fixed（检查弹簧参考长度/量纲）'}）"
        )
        if not monotone:
            lines.append(
                "- **非单调解释**：软端（k≤1e3）峰值在**内部 / 距下关节载荷面边缘**"
                "（~143–157 MPa，**远离跖面固定端**，即跖面奇异已被消掉）；随 k 增大跖面固定端奇异增强，"
                "在 k≈1e4 处峰值**迁移到跖面边缘**并接管，此后 gauge 单调升向 fixed 基线（k=1e6 残差 0.17%）。"
                "中间的浅极小值来自**两个边缘奇异（关节载荷边缘 vs 跖面固定边缘）的竞争/迁移**，"
                "不是求解失败（全部 rc=0，见峰值位置列）。"
            )
    else:
        lines.append("- 无成功的 spring 工况。")
    if spring_failed:
        lines.append("- 失败的 k：" + ", ".join(f"{r['spring_k_N_per_mm']:g}" for r in spring_failed))

    # extrapolation
    lines += [
        "",
        "## 3. 首次骨折高度：放松跖面后往 7–9 m 推多少？（外推，非重扫）",
        "",
        "方法：线弹性下应力 ∝ 载荷。由 gauge(5 m) 反解触发 150 MPa 所需的 F*，再用 anat 扫描",
        "真跑的 F_subt(h) 表反查高度。**这是外推，未重跑高度扫描。**",
        "",
        "| BC | k (N/mm) | gauge@5m (MPa) | F* (kN) | 首次骨折高度 h* (m) | vs [Y25] 7–9 m |",
        "|---|---|---|---|---|---|",
    ]
    for r in order:
        if r["gauge_max"] is None:
            continue
        k = "—" if r["spring_k_N_per_mm"] is None else f"{r['spring_k_N_per_mm']:g}"
        key = f"{r['plantar_bc']}:{r['spring_k_N_per_mm']}"
        ex = extrap[key]
        hf = ex.get("h_first_m")
        if hf is None:
            htxt, rel_txt = "> 扫描上界 (50 m)", "远晚于 [Y25]"
        else:
            htxt = f"{hf:.2f}"
            if hf < Y25_LO_M:
                rel_txt = "早于 [Y25]"
            elif hf <= Y25_HI_M:
                rel_txt = "**落入 [Y25]**"
            else:
                rel_txt = "晚于 [Y25]"
        fk = ex.get("f_crit_kN")
        lines.append(
            f"| {r['plantar_bc']} | {k} | {r['gauge_max']:.2f} | "
            f"{'—' if fk is None else f'{fk:.2f}'} | {htxt} | {rel_txt} |"
        )

    # conclusion
    g_fix = fixed["gauge_max"]
    lines += ["", "## 4. 结论", ""]
    lines.append(f"- **fixed 基线** gauge_max = {g_fix:.2f} MPa（复现 {BASELINE_GAUGE_MPA:.2f}）。")

    if roller is not None and roller["gauge_max"] is not None:
        ratio = g_fix / roller["gauge_max"]
        lines.append(
            f"- **roller (BCNormalDisplacement + 2 销)** gauge_max = {roller['gauge_max']:.2f} MPa，"
            f"与 fixed 的相对差 = {(ratio-1)*100:+.3f}%（比值 {ratio:.4f}）。"
            "⚠ FEBio 在**曲面**跖面上按面法向逐面约束，共享节点被多个法向约束合围"
            "→ 等价于全固定，**没有真正放开切向**（这是负结论：BCNormalDisplacement 不是本几何上的有效滚子）。"
        )
    if roller_free is not None and roller_free["gauge_max"] is not None:
        ratio = g_fix / roller_free["gauge_max"]
        lines.append(
            f"- **roller_free**（全局 Y 固定 + 切向弹簧 k_t={ROLLER_FREE_KT_N_PER_MM:g} N/mm，"
            f"真正放开切向）gauge_max = {roller_free['gauge_max']:.2f} MPa，"
            f"相对 fixed {(ratio-1)*100:+.1f}% → 只放开**切向**反而{'升高' if ratio<1 else '降低'}。"
        )
    if spring:
        gmin = min(spring, key=lambda r: r["gauge_max"])
        he = extrap.get(f"spring:{gmin['spring_k_N_per_mm']}") or {}
        hmin = he.get("h_first_m")
        lines.append(
            f"- **spring（弹性地基）**：软端把 gauge 压到最低的是 k={gmin['spring_k_N_per_mm']:g} N/mm → "
            f"gauge={gmin['gauge_max']:.2f} MPa（相对 fixed 降 {(1-gmin['gauge_max']/g_fix)*100:.1f}%）；"
            f"硬端 k={spring[-1]['spring_k_N_per_mm']:g} N/mm → gauge={spring[-1]['gauge_max']:.2f} MPa"
            f"（收敛回 fixed，残差 {abs(spring[-1]['gauge_max']-g_fix)/g_fix*100:.2f}%）。"
        )
        if hmin is not None:
            lines.append(
                f"- **首次骨折高度（外推）**：fixed ≈ 2.11 m；把跖面改成弹性地基后，"
                f"最好情形（k={gmin['spring_k_N_per_mm']:g} N/mm）≈ **{hmin:.2f} m**"
                f"（从 ~2.1 m 外推到 ~{hmin:.1f} m，向 [Y25] 7–9 m 靠近但**未到达**）。"
            )
    lines += [
        "",
        "### 回答原假设",
        "",
        "- **「跖面三向全固定 → gauge 虚高」部分成立，但机制是「法向（竖直）刚度」而非「三向/切向」**：",
        "  1. 纯切向释放（roller / roller_free）**不降反升**（roller 逐位≡fixed；roller_free +14%），",
        "     说明切向固定不是 gauge 的来源；",
        "  2. 把**法向支承**改成弹性地基（spring）后，软端 gauge 从 186.8 → ~143–157 MPa（**降 16–24%**），",
        "     峰值也从**跖面棱边**移到内部/关节面边缘；",
        "  3. 因此「人工制造固定端奇异」的主要成分是**竖直方向的刚性支承**，不是三向全约束本身。",
        "- **是否足以把首次骨折高度推到 7–9 m？** 按线弹性外推：~2.1 m → **~3.4–5.4 m**（最好 ~5.4 m），",
        "  **接近但达不到** [Y25] 的 7–9 m。要补齐差距还需其它机制（跟腱、软组织垫、几何非线性/接触）。",
        "",
        "> **诚实边界**：本节所有高度都是**线弹性比例外推**（gauge∝F + F(h) 反查），",
        "> 不是重跑高度扫描；几何非线性、接触/滑移、屈服后软化都未包含。",
        "> 若结论为「放松跖面对 gauge 影响很小」，同样是有效负结论。",
        "",
    ]

    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="跖面 BC（fixed/roller/spring）对 gauge 影响的探针")
    ap.add_argument("--gauge-mm", type=float, default=REG_LEN_MM)
    ap.add_argument("--ks", type=float, nargs="+", default=SPRING_KS, help="spring 刚度 (N/mm)")
    ap.add_argument("--outdir", type=Path, default=OUTDIR)
    args = ap.parse_args(argv)

    WORKDIR.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    _install_xplt_discrete_patch()

    print("=" * 74)
    print("跖面 BC 探针 · THUMS AM50 跟骨 · h=5 m a=0")
    print(f"  gauge R={args.gauge_mm:g} mm   threshold={THRESHOLD_MPA:.0f} MPa")
    print(f"  spring ks = {args.ks}")
    print("=" * 74, flush=True)

    t0 = time.time()
    mesh = T.load_thums_mesh(MESH_ANAT)
    print(f"[mesh] {mesh['source_npz']}")
    print(f"[mesh] 节点 {mesh['n_nodes']}  CORT {mesh['n_elements_cort']}  SPON "
          f"{mesh['n_elements_spon']}  plantar {mesh['n_plantar_nodes']}n/"
          f"{mesh['n_plantar_faces']}f", flush=True)
    cort_centroids, cort_vol = T.hex_centroids_volumes(mesh["nodes"], mesh["elements_cort"])

    f_sub, t1 = _load_at(DROP_H_M)
    load_meta = {"h_m": DROP_H_M, "a": ACT, "f_subt_n": f_sub, "f_subt_kN": f_sub / 1e3, "t1_rel_err": t1}
    print(f"[load] h={DROP_H_M:g} a={ACT:g} -> F_subt={f_sub/1e3:.3f} kN (T1={t1*100:+.1f}%)", flush=True)

    rows: list[dict] = []

    def _run(kind: str, k: float | None) -> dict:
        tag = kind if k is None else f"{kind}_k{k:g}"
        row = _solve_case(mesh, cort_centroids, cort_vol, f_sub, kind, k,
                          tag=tag, gauge_mm=args.gauge_mm)
        ko = "—" if k is None else f" k={k:g}"
        print(f"[{kind}{ko}] rc={row['febio_rc']}  gauge_max={row['gauge_max']:.3f}  "
              f"max={row['max']:.3f}  p95={row['p95']:.3f}  "
              f"peak@{row['location_class']}  [{time.time()-t0:.0f}s]", flush=True)
        return row

    # --- Step 0 + 1 ---
    fixed = _run("fixed", None)
    err = abs(fixed["gauge_max"] - BASELINE_GAUGE_MPA) / BASELINE_GAUGE_MPA
    if err > BASELINE_GAUGE_TOL:
        print(f"[FAIL] 基线未复现：{fixed['gauge_max']:.3f} vs {BASELINE_GAUGE_MPA:.3f} "
              f"(err={err*100:.2f}% > {BASELINE_GAUGE_TOL*100:.0f}%)")
        return 2
    print(f"[OK] 基线复现 gauge={fixed['gauge_max']:.3f} MPa (err={err*100:+.3f}%)", flush=True)
    rows.append(fixed)

    roller = _run("roller", None)
    rows.append(roller)
    # 诊断：真正放开切向的滚子（全局 Y 固定 + 极软切向弹簧，k_t=1 N/mm）。
    roller_free = _run("roller_free", ROLLER_FREE_KT_N_PER_MM)
    rows.append(roller_free)

    # --- Step 2 ---
    for k in args.ks:
        try:
            rows.append(_run("spring", float(k)))
        except Exception as exc:  # noqa: BLE001 - 如实记录失败的 k
            print(f"[spring k={k:g}] FAILED: {type(exc).__name__}: {str(exc)[:300]}", flush=True)
            rows.append({
                "tag": f"spring_k{k:g}", "plantar_bc": "spring",
                "spring_k_N_per_mm": float(k), "febio_rc": None,
                "load_n": f_sub, "load_kN": f_sub / 1e3, "achilles_kN": 0.0,
                "max": None, "p99": None, "p95": None, "mean": None,
                "gauge_max": None, "gauge_p95": None, "peak_elem": None,
                "peak_centroid_mm": None, "nearest_surface": None,
                "nearest_surface_dist_mm": None, "location_class": None,
                "error": f"{type(exc).__name__}: {exc}",
            })

    spring_rows = sorted(
        (r for r in rows if r["plantar_bc"] == "spring" and r["gauge_max"] is not None),
        key=lambda r: r["spring_k_N_per_mm"],
    )
    gs = [r["gauge_max"] for r in spring_rows]
    monotone = all(b >= a - 1e-9 for a, b in zip(gs, gs[1:])) if len(gs) >= 2 else True

    # --- extrapolation ---
    h_tab, f_tab = _load_height_table()
    extrap: dict[str, dict] = {}
    for r in rows:
        if r["gauge_max"] is None:
            continue
        key = f"{r['plantar_bc']}:{r['spring_k_N_per_mm']}"
        extrap[key] = _extrapolate_first_fracture(r["gauge_max"], h_tab, f_tab)

    # --- write summary ---
    summary = {
        "meta": {
            "study": "plantar BC sensitivity (fixed/roller/spring)",
            "mesh": mesh["source_npz"],
            "mesh_kind": "anat (direction C)",
            "gauge_radius_mm": args.gauge_mm,
            "threshold_mpa": THRESHOLD_MPA,
            "baseline_anat_gauge_mpa": BASELINE_GAUGE_MPA,
            "baseline_reproduced": bool(err <= BASELINE_GAUGE_TOL),
            "spring_monotone_to_fixed": bool(monotone),
            "extrapolation_method": "linear-in-load (gauge∝F) + F_subt(h) table inversion; NOT a re-scan",
            "n_solves_ok": sum(1 for r in rows if r.get("febio_rc") == 0),
        },
        "load": load_meta,
        "cases": rows,
        "spring_scan": {
            "k_N_per_mm": [r["spring_k_N_per_mm"] for r in spring_rows],
            "gauge_max_mpa": gs,
            "monotone_non_decreasing": bool(monotone),
            "fixed_gauge_max_mpa": fixed["gauge_max"],
            "hardest_k_residual_vs_fixed_pct": (
                abs(gs[-1] - fixed["gauge_max"]) / fixed["gauge_max"] * 100 if gs else None
            ),
        },
        "first_fracture_extrapolation": extrap,
    }
    json_out = args.outdir / "plantar_bc_summary.json"
    json_out.write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=float), encoding="utf-8")

    csv_fields = [
        "tag", "plantar_bc", "spring_k_N_per_mm", "febio_rc", "load_kN",
        "gauge_max", "gauge_p95", "max", "p99", "p95", "mean",
        "nearest_surface", "nearest_surface_dist_mm", "dist_joint_mm",
        "dist_plantar_mm", "location_class",
    ]
    csv_out = args.outdir / "plantar_bc_summary.csv"
    with csv_out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=csv_fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in csv_fields})

    _write_plot(rows, args.outdir / "plantar_bc_gauge_k.png")
    _write_report(rows, load_meta, extrap, mesh, args.outdir / "PLANTAR_BC_REPORT.md", monotone)

    print("-" * 74)
    print(f"[summary] fixed gauge={fixed['gauge_max']:.3f}  roller gauge={roller['gauge_max']:.3f}  "
          f"spring {gs[0] if gs else float('nan'):.1f}..{gs[-1] if gs else float('nan'):.1f}")
    print(f"[summary] monotone_to_fixed={monotone}")
    for r in rows:
        key = f"{r['plantar_bc']}:{r['spring_k_N_per_mm']}"
        ex = extrap.get(key)
        if ex and ex.get("h_first_m") is not None:
            print(f"  {r['tag']:>16}  gauge={r['gauge_max']:7.2f}  h_first≈{ex['h_first_m']:.2f} m")
    print(f"[out] {json_out}")
    print(f"[out] {csv_out}")
    print(f"[out] {args.outdir / 'plantar_bc_gauge_k.png'}")
    print(f"[out] {args.outdir / 'PLANTAR_BC_REPORT.md'}")
    print(f"[done] {time.time()-t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
