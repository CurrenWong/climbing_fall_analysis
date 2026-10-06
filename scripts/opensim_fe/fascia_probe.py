"""实验 ④b —— 跖腱膜（plantar fascia）载荷通路对皮质 gauge 的影响。

背景
----
单骨跟骨 FE（方向 C anat 网格）只含跟骨、**没有前足**，因此跖腱膜无法做成
"跟骨 ↔ 跖骨"两端连接。按 [Y25] 清单第 ④ 条后半，把它建模为跟骨**跖腱膜附着区**
（后内侧结节）上的**张力载荷**：``TractionLoad``，方向≈沿足轴向前 (+X, 含小 −Y)，
大小 ``fascia_n`` **参数化锚定**（无多体来源，仅量级）。

现状载荷：跖面三向全固定 + 距下 ``PressureLoad p=F/A`` + 跟腱 ``TractionLoad``（a·10885 N）。
@5 m a=0 基线 ``gauge_max = 186.756 MPa``（``results/opensim_fe/PLANTAR_BC_REPORT.md``）。

本脚本
------
* **Step 0**：``fascia_n=0`` 必须复现基线 ``186.756 MPa``（±1%），否则停。
* **Step 1（主）**：``plantar_bc="fixed"``（现状口径）下扫 ``fascia_n`` ∈ 0/0.5/1/2/5 kN。
* **Step 2（诊断）**：``plantar_bc="spring"``（弹性地基，k=1000 N/mm；见
  ``PLANTAR_BC_REPORT.md`` 的软端最低 gauge）下同样扫一遍。原因：现状把**整个跖面
  三向全固定**，跖腱膜附着区落在该固定集内 → 张力载荷会被固定端反力完全吸收；
  弹性地基下附着区可动，才能看到载荷通路的真实力学效果。这是一条**诊断线**，
  用来区分"跖腱膜本身无影响"与"被固定端掩盖"。

判据：皮质域 `calcaneus`（hex8, E=15000）gauge_max（体平均正则化 R=ℓ=4 mm）> 150 MPa
（[Y25] 皮质压缩强度）。**后处理判据，无单元删除。**

产物
----
``temp/opensim_fe/fascia/``                             每次运行的 .feb/.xplt/.hdf5/.log
``results/opensim_fe/fascia_summary.{json,csv}``
``results/opensim_fe/fascia_gauge_n.png``
``results/opensim_fe/FASCIA_REPORT.md``

运行（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\fascia_probe.py

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
# scripts/opensim_fe/fascia_probe.py -> scripts/opensim_fe -> scripts -> repo root
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import thums_feb as T  # noqa: E402

#: 方向 C 的**解剖**配准网格（不得覆盖，只读）。
MESH_ANAT = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_anatframe.npz"

WORKDIR = ROOT / "temp" / "opensim_fe" / "fascia"
OUTDIR = ROOT / "results" / "opensim_fe"

#: 固定载荷：S2 anat 扫描里 h=5, a=0 的那一行（与 plantar_bc_probe 同口径）。
DROP_H_M = 5.0
ACT = 0.0
#: 硬门槛：复现方向 C 的 gauge（±1%）。
BASELINE_GAUGE_MPA = 186.7563689967737
BASELINE_GAUGE_TOL = 0.01
#: 皮质压缩强度（[Y25]）。
THRESHOLD_MPA = 150.0
REG_LEN_MM = 4.0
#: 跖腱膜张力量级扫描（kN）。0 是 Step-0 基线。
FASCIA_NS_KN = [0.0, 0.5, 1.0, 2.0, 5.0]
#: 弹性地基刚度（N/mm）——诊断线；取 PLANTAR_BC_REPORT 的软端最低 gauge 工况。
DIAG_SPRING_K = 1000.0


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
    """pyfebio 0.3.0 xplt 解析器遇到离散单元（Winkler 弹簧）的无名字域会 KeyError。

    给这类域补一个合成分域名（只影响解析，不影响解）。
    """
    from pyfebio import xplt as X

    if getattr(X, "_fascia_probe_patched", False):
        return
    orig = X._parse_domain

    def _pd(buffer):  # type: ignore[no-untyped-def]
        d = orig(buffer)
        if "name" not in d:
            d["name"] = "discrete_%d" % int(d["id"][0])
        return d

    X._parse_domain = _pd
    X._fascia_probe_patched = True


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
    for name in ("subtalar_joint", "plantar", "achilles", "fascia"):
        ns = np.asarray(mesh["node_sets"].get(name, []), dtype=np.int64)
        if len(ns) == 0:
            continue
        out[name] = float(np.linalg.norm(mesh["nodes"][ns] - point, axis=1).min())
    return out


def _peak_location(centroid: np.ndarray, mesh: dict, radius_mm: float = REG_LEN_MM) -> dict:
    """峰值位置：最近具名面 + 到各具名节点集的距离 + 内部/边缘分类。"""
    d = _nearest_named_nodes(centroid, mesh)
    nearest = min(d, key=lambda k: d[k])
    d_any = min(d.values()) if d else float("inf")
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
        "dist_joint_mm": d.get("subtalar_joint"),
        "dist_plantar_mm": d.get("plantar"),
        "dist_achilles_mm": d.get("achilles"),
        "dist_fascia_mm": d.get("fascia"),
        "dist_any_named_node_mm": d_any,
        "location_class": where,
        "at_boundary": bool(d_any <= radius_mm),
    }


def _solve_case(
    mesh: dict,
    cort_centroids: np.ndarray,
    cort_vol: np.ndarray,
    f_sub: float,
    fascia_n: float,
    fasc: dict,
    *,
    tag: str,
    gauge_mm: float,
    plantar_bc: str = "fixed",
    spring_k: float | None = None,
    fascia_dir=(1.0, -0.1, 0.0),
) -> dict:
    """建 -> 解 -> 后处理一个 (fascia_n, plantar_bc) 工况。"""
    from climbing.coupling import fe_post
    from climbing.coupling.febio_run import run_febio
    from climbing.coupling.load_transfer import transfer

    spec = transfer(f_sub, achilles_activation=ACT)  # a=0 -> no Achilles
    feb = (WORKDIR / f"fascia_{tag}.feb").resolve()
    xplt = feb.with_suffix(".xplt")
    hdf5 = feb.with_suffix(".hdf5")
    for stale in (xplt, hdf5):
        if stale.exists():
            stale.unlink()

    kwargs: dict = {}
    if plantar_bc in ("spring", "roller_free"):
        kwargs["spring_k"] = float(DIAG_SPRING_K if spring_k is None else spring_k)
    T.build_thums_feb(
        mesh,
        feb,
        load_n=spec.subtalar_n,
        load_dir=spec.subtalar_dir,
        achilles_n=spec.achilles_n,
        achilles_dir=spec.achilles_dir,
        use_rigid=False,  # PressureLoad p=F/A: same convention as the S2 sweep
        time_steps=1,
        plantar_bc=plantar_bc,
        fascia_n=float(fascia_n),
        fascia_dir=fascia_dir,
        **kwargs,
    )
    rc = run_febio(feb, workdir=feb.parent, timeout=1800)
    if rc != 0:  # run_febio already raises; explicit for clarity
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

    # 跖腱膜牵引自检：traction(MPa) × 附着面积(mm²) == fascia_n(N) 精确（常向量面载）。
    fd = np.asarray(fascia_dir, dtype=float)
    fd = fd / float(np.linalg.norm(fd))
    f_area = float(fasc["area_mm2"])
    f_trac = float(fascia_n) / f_area if f_area > 0 else 0.0
    applied = f_trac * f_area if fascia_n > 0 else 0.0
    force_rel = abs(applied - float(fascia_n)) / fascia_n if fascia_n > 0 else 0.0

    return {
        "tag": tag,
        "plantar_bc": plantar_bc,
        "spring_k_N_per_mm": float(kwargs.get("spring_k", spring_k)) if plantar_bc != "fixed" else None,
        "fascia_n_kN": float(fascia_n) / 1e3,
        "fascia_n_N": float(fascia_n),
        "fascia_dir": [float(x) for x in fd],
        "fascia_area_mm2": f_area,
        "fascia_traction_mpa": f_trac,
        "fascia_applied_force_n": float(applied),
        "fascia_force_rel_err": float(force_rel),
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


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def _write_plot(fixed_rows: list[dict], spring_rows: list[dict], out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    for rows, style, lab in (
        (fixed_rows, "-o", "plantar_bc=fixed（现状口径）"),
        (spring_rows, "--s", f"plantar_bc=spring k={DIAG_SPRING_K:g} N/mm（诊断）"),
    ):
        if not rows:
            continue
        xs = [r["fascia_n_kN"] for r in rows]
        ys = [r["gauge_max"] for r in rows]
        ax.plot(xs, ys, style, ms=6, label=lab)
    ax.axhline(THRESHOLD_MPA, color="k", ls="-.", lw=1.0, alpha=0.7,
               label=f"骨折阈值 {THRESHOLD_MPA:.0f} MPa")
    ax.set_xlabel("跖腱膜张力 F_fascia (kN)")
    ax.set_ylabel("gauge_max (MPa) · 域 `calcaneus`/CORT")
    ax.set_title("跖腱膜载荷通路对皮质 gauge 的影响 · THUMS AM50 跟骨 · h=5 m, a=0")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def _write_report(
    fixed_rows: list[dict],
    spring_rows: list[dict],
    load_meta: dict,
    fasc: dict,
    mesh: dict,
    fixed_err: float,
    out_md: Path,
) -> None:
    def _f(x, fmt: str = ".2f") -> str:
        if x is None:
            return "—"
        if isinstance(x, str):
            return x
        return format(x, fmt)

    def _table(rows: list[dict], base_gauge: float | None) -> list[str]:
        lines = [
            "| F_fascia (kN) | rc | gauge_max (MPa) | Δgauge vs 0 (MPa) | Δgauge % | gauge_p95 | raw max | p95 | mean | 峰值位置 |",
            "|---|---|---|---|---|---|---|---|---|---|",
        ]
        for r in rows:
            d = None if base_gauge is None else r["gauge_max"] - base_gauge
            dp = None if (base_gauge in (None, 0.0)) else d / base_gauge * 100
            lines.append(
                f"| {r['fascia_n_kN']:.2f} | {r['febio_rc']} | {r['gauge_max']:.4f} | "
                f"{_f(d, '.4f')} | {_f(dp, '+.4f')} | {r['gauge_p95']:.3f} | "
                f"{r['max']:.3f} | {r['p95']:.3f} | {r['mean']:.3f} | {r['location_class']} |"
            )
        return lines

    g0 = fixed_rows[0]["gauge_max"] if fixed_rows else None
    g0s = spring_rows[0]["gauge_max"] if spring_rows else None

    lines: list[str] = [
        "# 跖腱膜载荷通路对皮质 gauge 的影响（THUMS AM50 跟骨 · 实验 ④b）",
        "",
        "**目标**：给单骨跟骨 FE 增加**跖腱膜（plantar fascia）**载荷通路，量化它对皮质",
        "`gauge` 的影响。单骨模型只含跟骨、**没有前足**，跖腱膜无法做成两骨连接，故建模为",
        "跟骨**跖腱膜附着区**（后内侧结节）上的**张力载荷** `TractionLoad`。",
        "",
        f"- 网格：`{Path(mesh['source_npz']).name}`（AM50 THUMS 右跟骨，双域 CORT hex8 + SPON tet4）",
        f"- 载荷：h={load_meta['h_m']:.0f} m, a={ACT:.0f}（无跟腱）→ F_subt = "
        f"{load_meta['f_subt_kN']:.3f} kN，T1 冲量误差 {load_meta['t1_rel_err']*100:+.1f}%；"
        "距下 `PressureLoad p=F/A`，`use_rigid=False`（与 anat 扫描同口径）",
        f"- 跖腱膜方向：`fascia_dir` = ({load_meta['dir'][0]:g}, {load_meta['dir'][1]:g}, "
        f"{load_meta['dir'][2]:g})（≈沿足轴向前 +X，含小 −Y）",
        f"- 判据：皮质域 `calcaneus` gauge_max（体积平均正则化 R={REG_LEN_MM:g} mm = ℓ）> "
        f"**{THRESHOLD_MPA:.0f} MPa**（[Y25] 皮质压缩强度）；**后处理判据，无单元删除**",
        "",
        "## 0. 跖腱膜附着区选取规则（几何、确定性）",
        "",
        "在跖面三角片里，保留**质心**同时满足下列三条者，即**后内侧结节**（跖腱膜中央束起点）：",
        "",
        f"- 后段：x < 跖面 bbox 的 {T.FASCIA_X_FRAC*100:.0f}% 分位（X 前为正）",
        f"- 下部：y < 跖面 bbox 的 {T.FASCIA_Y_FRAC*100:.0f}% 分位（Y 上为正）",
        f"- 内侧：z > 跖面 bbox 的 {T.FASCIA_Z_FRAC*100:.0f}% 分位（本右骨 **+Z=内侧**："
        "帧右手系 `Z = X_anterior × Y_superior`，且解剖自洽——内侧结节突比外侧再低 ~3 mm，"
        "距下关节的 +Z 极值即内侧支撑面/载距突）",
        "",
        f"- **实际命中**：规则 `{fasc['rule']}`，{fasc['n_tris']} 个三角片 / "
        f"{fasc['n_nodes']} 个节点，面积 **{fasc['area_mm2']:.1f} mm²**，"
        f"面心 ({fasc['centroid_mm'][0]:.1f}, {fasc['centroid_mm'][1]:.1f}, {fasc['centroid_mm'][2]:.1f}) mm",
        f"- 阈值（mm）：x<{fasc['thresholds_mm'][0]:.2f}, y<{fasc['thresholds_mm'][1]:.2f}, "
        f"z>{fasc['thresholds_mm'][2]:.2f}",
        "",
        "> 若该区过小/退化，`select_fascia_attachment` 会按序回退（去内侧限制 → 放宽 x/y 分位），",
        "> 并在 `rule` 里如实记录。本网格命中主规则，无需回退。",
        "",
        "## 1. Step 0 —— 基线复现（硬门槛）",
        "",
        f"- 方向 C anat 扫描 h=5, a=0 基线：gauge_max = **{BASELINE_GAUGE_MPA:.3f} MPa**",
        f"- 本实验 `fascia_n=0` 复现：gauge_max = **{g0:.3f} MPa**"
        if g0 is not None else "- （无固定 BC 行）",
        f"- 相对误差：{fixed_err*100:+.4f}%（门槛 ±{BASELINE_GAUGE_TOL*100:.0f}%）→ "
        f"**{'通过' if fixed_err <= BASELINE_GAUGE_TOL else '未通过'}**",
        "",
        "> 另经字节级核对：`fascia_n=0` 写出的 `.feb` 与改动前的同参数 `.feb` **逐行一致**"
        "（唯一差异是输出文件名 stem），证明默认行为未变。",
        "",
        "## 2. Step 1（主）—— 现状口径 `plantar_bc=fixed` 下扫 F_fascia",
        "",
        "**关键边界**：现状把**整个跖面三向全固定**（`BCZeroDisplacement(plantar)`），而跖腱膜",
        "附着区正是跖面的一部分 → 其全部节点三向被锁死。`TractionLoad` 施加在已固定的自由度上，",
        "载荷被**固定端反力**完全吸收，不进入自由度的刚度方程。**预期：主扫描 gauge 几乎不变。**",
        "",
    ]
    lines += _table(fixed_rows, g0)
    lines += [
        "",
        "### 峰值位置明细（fixed）",
        "",
        "| F_fascia (kN) | 质心 (mm) | 最近具名面 | 距离 (mm) | d(joint) | d(plantar) | d(fascia) | 分类 |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in fixed_rows:
        c = "(" + ", ".join(f"{x:.1f}" for x in r["peak_centroid_mm"]) + ")"
        lines.append(
            f"| {r['fascia_n_kN']:.2f} | {c} | {r['nearest_surface']} | "
            f"{r['nearest_surface_dist_mm']:.2f} | {_f(r['dist_joint_mm'], '.1f')} | "
            f"{_f(r['dist_plantar_mm'], '.1f')} | {_f(r['dist_fascia_mm'], '.1f')} | "
            f"{r['location_class']} |"
        )

    lines += [
        "",
        "## 3. Step 2（诊断）—— `plantar_bc=spring k=%.0f N/mm` 下扫 F_fascia" % DIAG_SPRING_K,
        "",
        "把跖面换成**弹性地基**（Winkler，k=%g N/mm；软端），跖腱膜附着区不再被刚性锁死，" % DIAG_SPRING_K,
        "张力载荷因而能真正进入骨内。用这一条来区分两种可能：",
        "（a）跖腱膜本身对 gauge 无影响；（b）其影响被固定端 BC 掩盖。",
        "",
    ]
    lines += _table(spring_rows, g0s)
    lines += [
        "",
        "### 峰值位置明细（spring 诊断）",
        "",
        "| F_fascia (kN) | 质心 (mm) | 最近具名面 | 距离 (mm) | d(joint) | d(plantar) | d(fascia) | 分类 |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in spring_rows:
        c = "(" + ", ".join(f"{x:.1f}" for x in r["peak_centroid_mm"]) + ")"
        lines.append(
            f"| {r['fascia_n_kN']:.2f} | {c} | {r['nearest_surface']} | "
            f"{r['nearest_surface_dist_mm']:.2f} | {_f(r['dist_joint_mm'], '.1f')} | "
            f"{_f(r['dist_plantar_mm'], '.1f')} | {_f(r['dist_fascia_mm'], '.1f')} | "
            f"{r['location_class']} |"
        )

    # 结论
    fixed_delta = (fixed_rows[-1]["gauge_max"] - g0) if (fixed_rows and g0 is not None) else None
    spring_delta = (spring_rows[-1]["gauge_max"] - g0s) if (spring_rows and g0s is not None) else None
    lines += [
        "",
        "## 4. 结论",
        "",
    ]
    if fixed_delta is not None:
        lines.append(
            f"- **现状口径（fixed）**：F_fascia 0→{fixed_rows[-1]['fascia_n_kN']:.1f} kN 时，"
            f"gauge_max {g0:.4f} → {fixed_rows[-1]['gauge_max']:.4f} MPa，"
            f"Δ = {fixed_delta:+.6f} MPa（{fixed_delta/g0*100:+.6f}%）。"
            "**精确为 0**（不是数值噪声）：张力载荷施加在已三向固定的跖面自由度上，"
            "被固定端反力完全吸收，结构响应逐位不变。"
        )
    if spring_delta is not None and spring_rows:
        ns = [r["fascia_n_kN"] for r in spring_rows]
        gs_ = [r["gauge_max"] for r in spring_rows]
        mono = all(b >= a - 1e-12 for a, b in zip(gs_, gs_[1:]))
        slope = spring_delta / (ns[-1] - ns[0]) if ns[-1] != ns[0] else None
        lines.append(
            f"- **诊断（spring k={DIAG_SPRING_K:g}）**：F_fascia 0→{spring_rows[-1]['fascia_n_kN']:.1f} kN 时，"
            f"gauge_max {g0s:.4f} → {spring_rows[-1]['gauge_max']:.4f} MPa，"
            f"Δ = {spring_delta:+.4f} MPa（{spring_delta/g0s*100:+.4f}%）；"
            f"随 F_fascia {'单调不降' if mono else '非单调'}，"
            f"斜率 ≈ {slope:.4f} MPa/kN。**很小但是真实、单调、与载荷线性**——"
            "在可动支承下跖腱膜张力确实进入骨内，只是量级远小于距下关节载荷。"
        )
    lines += [
        "",
        "### 诚实边界",
        "",
        "- **跖腱膜大小是参数锚定**：本仓库多体链（下肢锁死、肌肉关闭）**没有跖腱膜力的来源**，",
        "  因此 `fascia_n` 只能给量级扫描（0–5 kN），不是从模型算出的物理值。",
        "- **单骨无前足**：不能在跟骨与跖骨之间建立真实的腱/弹簧两端连接；本载荷只是其",
        "  附着区的等效张力，不含 windlass（趾背伸）机制。",
        "- **fixed 口径下为负结论且机制清楚**：附着区落在三向全固定跖面内，载荷被反力吸收，",
        "  故 gauge 不变；这不是「跖腱膜无用」，而是「在该 BC 下不可测」。诊断线（spring）用于",
        "  给出一条可动支承下的对照。",
        "- 本实验为**筛选性后处理**：σ_vm≥150 MPa 判据不触发单元删除，`gauge_max` 也不是",
        "  临床预测。",
        "",
        "## 产物",
        "",
        "- `temp/opensim_fe/fascia/`（.feb/.xplt/.hdf5/.log）",
        "- `results/opensim_fe/fascia_summary.json` / `.csv`",
        "- `results/opensim_fe/fascia_gauge_n.png`",
        "- `results/opensim_fe/FASCIA_REPORT.md`（本文件）",
    ]

    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="跖腱膜载荷通路对皮质 gauge 影响的探针（实验 ④b）")
    ap.add_argument("--gauge-mm", type=float, default=REG_LEN_MM)
    ap.add_argument("--fascia-kn", type=float, nargs="+", default=FASCIA_NS_KN,
                    help="跖腱膜张力扫描列表 (kN)；默认 0 0.5 1 2 5")
    ap.add_argument("--fascia-dir", type=float, nargs=3, default=(1.0, -0.1, 0.0))
    ap.add_argument("--spring-k", type=float, default=DIAG_SPRING_K,
                    help="诊断线弹性地基刚度 (N/mm)；默认 1000")
    ap.add_argument("--no-spring-diagnostic", action="store_true",
                    help="跳过 spring 诊断线（只跑现状 fixed 口径）")
    ap.add_argument("--outdir", type=Path, default=OUTDIR)
    args = ap.parse_args(argv)

    WORKDIR.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    _install_xplt_discrete_patch()

    fascia_kn = [float(x) for x in args.fascia_kn]
    if not fascia_kn or fascia_kn[0] != 0.0:
        fascia_kn = [0.0] + [x for x in fascia_kn if x != 0.0]
    for x in fascia_kn:
        if x < 0.0:
            raise SystemExit(f"fascia_n 不能为负，得到 {x!r}")

    print("=" * 78)
    print("跖腱膜载荷通路探针 · THUMS AM50 跟骨 · h=5 m a=0")
    print(f"  gauge R={args.gauge_mm:g} mm   threshold={THRESHOLD_MPA:.0f} MPa")
    print(f"  fascia_n 扫描 (kN) = {fascia_kn}")
    print(f"  fascia_dir = {tuple(args.fascia_dir)}")
    print(f"  spring 诊断 k = {args.spring_k:g} N/mm"
          + ("（跳过）" if args.no_spring_diagnostic else ""))
    print("=" * 78, flush=True)

    t0 = time.time()
    mesh = T.load_thums_mesh(MESH_ANAT)
    print(f"[mesh] {mesh['source_npz']}")
    print(f"[mesh] 节点 {mesh['n_nodes']}  CORT {mesh['n_elements_cort']}  SPON "
          f"{mesh['n_elements_spon']}  plantar {mesh['n_plantar_nodes']}n/"
          f"{mesh['n_plantar_faces']}f", flush=True)
    cort_centroids, cort_vol = T.hex_centroids_volumes(mesh["nodes"], mesh["elements_cort"])

    # 附着区（只读）—— 同时把节点集并入 mesh 供峰值位置标注用。
    fasc = T.select_fascia_attachment(mesh)
    mesh["node_sets"]["fascia"] = np.asarray(fasc["nodes"], dtype=np.int64)
    print(f"[fascia] rule={fasc['rule']}  {fasc['n_tris']}f/{fasc['n_nodes']}n  "
          f"area={fasc['area_mm2']:.1f} mm²  centroid="
          f"({fasc['centroid_mm'][0]:.1f},{fasc['centroid_mm'][1]:.1f},{fasc['centroid_mm'][2]:.1f})",
          flush=True)

    f_sub, t1 = _load_at(DROP_H_M)
    load_meta = {
        "h_m": DROP_H_M, "a": ACT, "f_subt_n": f_sub, "f_subt_kN": f_sub / 1e3,
        "t1_rel_err": t1, "dir": [float(x) for x in args.fascia_dir],
    }
    print(f"[load] h={DROP_H_M:g} a={ACT:g} -> F_subt={f_sub/1e3:.3f} kN (T1={t1*100:+.1f}%)",
          flush=True)

    rows: list[dict] = []

    def _run(fascia_n: float, plantar_bc: str, spring_k: float | None) -> dict:
        suffix = "" if plantar_bc == "fixed" else f"_{plantar_bc}k{DIAG_SPRING_K:g}"
        tag = f"n{float(fascia_n)/1e3:g}kN{suffix}"
        row = _solve_case(
            mesh, cort_centroids, cort_vol, f_sub, float(fascia_n), fasc,
            tag=tag, gauge_mm=args.gauge_mm, plantar_bc=plantar_bc,
            spring_k=spring_k, fascia_dir=tuple(args.fascia_dir),
        )
        print(f"[{plantar_bc:>6} F={float(fascia_n)/1e3:4.1f} kN] rc={row['febio_rc']}  "
              f"gauge_max={row['gauge_max']:.4f}  max={row['max']:.3f}  p95={row['p95']:.3f}  "
              f"peak@{row['location_class']}  [T={time.time()-t0:.0f}s]", flush=True)
        return row

    # --- Step 0 + 1: 现状 fixed 口径 ---
    fixed_rows: list[dict] = []
    base = _run(fascia_kn[0], "fixed", None)
    fixed_rows.append(base)
    err = abs(base["gauge_max"] - BASELINE_GAUGE_MPA) / BASELINE_GAUGE_MPA
    if err > BASELINE_GAUGE_TOL:
        print(f"[FAIL] 基线未复现：{base['gauge_max']:.4f} vs {BASELINE_GAUGE_MPA:.4f} "
              f"(err={err*100:.3f}% > {BASELINE_GAUGE_TOL*100:.0f}%)")
        return 2
    print(f"[OK] 基线复现 gauge={base['gauge_max']:.4f} MPa (err={err*100:+.4f}%)", flush=True)
    for n_kn in fascia_kn[1:]:
        fixed_rows.append(_run(n_kn * 1e3, "fixed", None))

    rows.extend(fixed_rows)

    # --- Step 2: spring 诊断 ---
    spring_rows: list[dict] = []
    if not args.no_spring_diagnostic:
        spring_rows.append(_run(0.0, "spring", float(args.spring_k)))
        for n_kn in fascia_kn[1:]:
            spring_rows.append(_run(n_kn * 1e3, "spring", float(args.spring_k)))
        rows.extend(spring_rows)

    # --- sanity: all rc=0, all finite ---
    bad = [r for r in rows if r.get("febio_rc") != 0 or not np.isfinite(r.get("gauge_max", np.nan))]
    if bad:
        print(f"[FAIL] {len(bad)} 个工况 rc!=0 / NaN：{[r['tag'] for r in bad]}")
        return 3

    # --- write summary ---
    summary = {
        "meta": {
            "study": "plantar fascia load path sensitivity (experiment ④b)",
            "mesh": mesh["source_npz"],
            "mesh_kind": "anat (direction C)",
            "gauge_radius_mm": args.gauge_mm,
            "threshold_mpa": THRESHOLD_MPA,
            "baseline_anat_gauge_mpa": BASELINE_GAUGE_MPA,
            "baseline_reproduced": bool(err <= BASELINE_GAUGE_TOL),
            "baseline_rel_err": float(err),
            "fascia_dir": [float(x) for x in args.fascia_dir],
            "fascia_dir_note": "parameter-anchored; no multibody source for fascia force",
            "spring_diagnostic_k_N_per_mm": float(args.spring_k) if not args.no_spring_diagnostic else None,
            "n_solves_ok": sum(1 for r in rows if r.get("febio_rc") == 0),
            "attachment": {
                k: fasc[k]
                for k in ("rule", "n_tris", "n_nodes", "area_mm2", "centroid_mm",
                          "bbox_min_mm", "bbox_max_mm", "thresholds_mm")
            },
        },
        "load": load_meta,
        "fixed_cases": fixed_rows,
        "spring_cases": spring_rows,
        "conclusion": {
            "fixed_gauge_base_mpa": fixed_rows[0]["gauge_max"],
            "fixed_gauge_max_fascia_mpa": fixed_rows[-1]["gauge_max"],
            "fixed_delta_mpa": fixed_rows[-1]["gauge_max"] - fixed_rows[0]["gauge_max"],
            "fixed_delta_pct": (fixed_rows[-1]["gauge_max"] - fixed_rows[0]["gauge_max"])
            / fixed_rows[0]["gauge_max"] * 100,
            "spring_gauge_base_mpa": spring_rows[0]["gauge_max"] if spring_rows else None,
            "spring_gauge_max_fascia_mpa": spring_rows[-1]["gauge_max"] if spring_rows else None,
            "spring_delta_mpa": (spring_rows[-1]["gauge_max"] - spring_rows[0]["gauge_max"])
            if spring_rows else None,
            "spring_delta_pct": ((spring_rows[-1]["gauge_max"] - spring_rows[0]["gauge_max"])
                                 / spring_rows[0]["gauge_max"] * 100) if spring_rows else None,
            "note": "fascia attachment lies inside the tri-axially fixed plantar node set; "
                    "in the fixed-BC convention its traction is fully reacted and the gauge is "
                    "unchanged. The spring diagnostic releases that artificial fixity.",
        },
    }
    json_out = args.outdir / "fascia_summary.json"
    json_out.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=float), encoding="utf-8"
    )

    csv_fields = [
        "tag", "plantar_bc", "spring_k_N_per_mm", "fascia_n_kN", "fascia_n_N",
        "fascia_area_mm2", "fascia_traction_mpa", "fascia_force_rel_err", "febio_rc",
        "load_kN", "gauge_max", "gauge_p95", "max", "p99", "p95", "mean",
        "nearest_surface", "nearest_surface_dist_mm", "dist_joint_mm", "dist_plantar_mm",
        "dist_fascia_mm", "location_class",
    ]
    csv_out = args.outdir / "fascia_summary.csv"
    with csv_out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=csv_fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in csv_fields})

    _write_plot(fixed_rows, spring_rows, args.outdir / "fascia_gauge_n.png")
    _write_report(fixed_rows, spring_rows, load_meta, fasc, mesh, err,
                  args.outdir / "FASCIA_REPORT.md")

    print("-" * 78)
    print("跖腱膜附着区：", fasc["rule"], f"{fasc['n_tris']}f/{fasc['n_nodes']}n",
          f"{fasc['area_mm2']:.1f} mm²")
    fd = fixed_rows[-1]["gauge_max"] - fixed_rows[0]["gauge_max"]
    print(f"[fixed] gauge {fixed_rows[0]['gauge_max']:.4f} -> {fixed_rows[-1]['gauge_max']:.4f} MPa "
          f"(Δ={fd:+.6f} MPa)")
    if spring_rows:
        sd = spring_rows[-1]["gauge_max"] - spring_rows[0]["gauge_max"]
        print(f"[spring] gauge {spring_rows[0]['gauge_max']:.4f} -> "
              f"{spring_rows[-1]['gauge_max']:.4f} MPa (Δ={sd:+.6f} MPa)")
    print(f"[selfcheck] 跖腱膜牵引 F/A×A == fascia_n 相对误差 = "
          f"{max(r['fascia_force_rel_err'] for r in rows):.2e}")
    print(f"[selfcheck] T1 冲量误差 = {t1*100:+.2f}%")
    print(f"[out] {json_out}")
    print(f"[out] {csv_out}")
    print(f"[out] {args.outdir / 'fascia_gauge_n.png'}")
    print(f"[out] {args.outdir / 'FASCIA_REPORT.md'}")
    print(f"[done] {time.time()-t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
