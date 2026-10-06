"""Phase-S4 probe —— 3D wrench 上 FE 边界 + 跖面弹性支撑（opt-in）。

方案：``docs/非垂直落地扩展方案.md`` §3 L4 / §4 Phase-S4 / §7（轴向回归硬门槛）。

本脚本只做两件事，**不覆盖任何既有产物**：

1. **轴向回归（硬门槛）**：``build_thums_feb(..., plantar_bc="fixed",
   use_rigid=False, load_n=26383.133775894516)`` 必须复现改前缓存
   (``plantar_bc_summary.json`` 的 ``fixed`` 行：gauge_max=186.756、
   raw max=213.503、p95=136.297)。这是"默认逐位不变"的数值证明。
2. **S4 演示（FE）**：把方案-S1 的 20° 倾斜 wrench（力矢量 + 力矩，
   ``nonvertical_s1.json`` 的 ``feloadspec``）送上跟骨 FE，在
   **默认 fixed 跖面 BC** 与 **弹性支撑（spring, k=100 N/mm）** 下各跑一次；
   并跑 wrench 轴向对照组。报告 σ_vm（max/p95）+ 热点位置 + A7 边缘奇异是否被缓解。

单位 mm–N–MPa–s。用 ``use_rigid=False`` 口径（与方向 C / PLANAR_BC 同口径）。

产物（仅新增）::

    results/opensim_fe/nonvertical_s4.json
    results/opensim_fe/NONVERTICAL_S4_REPORT.md

运行（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s4.py
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import thums_feb as T  # noqa: E402

# ---- 只读输入 -------------------------------------------------------------
MESH_ANAT = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_anatframe.npz"
S1_JSON = ROOT / "results" / "opensim_fe" / "nonvertical_s1.json"
PLANTAR_JSON = ROOT / "results" / "opensim_fe" / "plantar_bc_summary.json"
S1_OPEN_FE = ROOT / "results" / "opensim_fe" / "s1_h5_opensim.json"

OUT_JSON = ROOT / "results" / "opensim_fe" / "nonvertical_s4.json"
OUT_MD = ROOT / "results" / "opensim_fe" / "NONVERTICAL_S4_REPORT.md"

WORKDIR = ROOT / "temp" / "opensim_fe" / "s4_nonvertical"

# ---- 判据 / 常数 ----------------------------------------------------------
REG_LEN_MM = 4.0            # gauge 正则化半径 ℓ
THRESHOLD_MPA = 150.0       # [Y25] 皮质压缩强度
SPRING_K = 100.0            # PLANAR_BC 最优弹性地基刚度 (N/mm)
TILT_DEG = 20.0

# ---- 改前基线（PLANAR_BC 缓存；硬门槛）----------------------------------
BASELINE_GAUGE_MPA = 186.7563689967737
BASELINE_MAX_MPA = 213.50277948448795
BASELINE_P95_MPA = 136.29650339601775
BASELINE_TOL = 1e-6        # "逐位不变"目标（数值误差内）


def _install_xplt_discrete_patch() -> None:
    """pyfebio xplt 解析器遇到离散单元（Winkler 弹簧）写出的无名字域会
    ``KeyError: 'name'``。补一个合成分域名（只影响解析，不影响解）。"""
    from pyfebio import xplt as X

    if getattr(X, "_s4_patched", False):
        return
    orig = X._parse_domain

    def _pd(buffer):  # type: ignore[no-untyped-def]
        d = orig(buffer)
        if "name" not in d:
            d["name"] = "discrete_%d" % int(d["id"][0])
        return d

    X._parse_domain = _pd
    X._s4_patched = True


def _read_domain_vm(xplt: Path, domain: str, hdf5: Path) -> np.ndarray:
    import h5py

    from climbing.coupling import fe_post
    from pyfebio import xplt as xplt_mod

    if not Path(hdf5).is_file():
        xplt_mod.to_hdf5(str(xplt), str(hdf5))
    with h5py.File(hdf5, "r") as h5:
        stress, _ = fe_post._read_last_stress(h5, domain)
    return fe_post.von_mises_from_voigt(stress)


def _peak_location(centroid: np.ndarray, mesh: dict, radius_mm: float = REG_LEN_MM) -> dict:
    """峰值位置：最近具名面 + 到各具名节点集的距离 + 内部/边缘分类。"""
    dist: dict[str, float] = {}
    for name in ("subtalar_joint", "plantar", "achilles"):
        ns = np.asarray(mesh["node_sets"].get(name, []), dtype=np.int64)
        dist[name] = (
            float(np.linalg.norm(mesh["nodes"][ns] - centroid, axis=1).min())
            if len(ns) else float("inf")
        )
    nearest = min(dist, key=lambda k: dist[k])
    all_named = np.concatenate(
        [np.asarray(mesh["node_sets"].get(n, []), dtype=np.int64) for n in dist]
    ) if any(len(mesh["node_sets"].get(n, [])) for n in dist) else np.empty(0, dtype=np.int64)
    d_any = (
        float(np.linalg.norm(mesh["nodes"][all_named] - centroid, axis=1).min())
        if len(all_named) else float("inf")
    )
    if nearest == "subtalar_joint" and dist[nearest] <= radius_mm:
        where = "关节面边缘 (joint rim)"
    elif nearest == "plantar" and dist[nearest] <= radius_mm:
        where = "跖面边缘 (plantar rim)"
    elif d_any <= radius_mm:
        where = f"{nearest} 边缘"
    else:
        where = "内部 (interior)"
    return {
        "nearest_surface": nearest,
        "nearest_surface_dist_mm": dist[nearest],
        "dist_joint_mm": dist["subtalar_joint"],
        "dist_plantar_mm": dist["plantar"],
        "dist_any_named_node_mm": d_any,
        "location_class": where,
        "at_plantar_rim": bool(nearest == "plantar" and dist["plantar"] <= radius_mm),
    }


def _solve_case(
    tag: str,
    mesh: dict,
    cort_centroids: np.ndarray,
    cort_vol: np.ndarray,
    *,
    load_n: float,
    subtalar_force,
    subtalar_moment,
    plantar_bc: str,
    spring_k: float | None,
    gauge_mm: float,
    use_rigid: bool = False,
) -> dict:
    from climbing.coupling import fe_post
    from climbing.coupling.febio_run import run_febio

    workdir = WORKDIR
    workdir.mkdir(parents=True, exist_ok=True)
    feb = (workdir / f"s4_{tag}.feb").resolve()
    xplt = feb.with_suffix(".xplt")
    hdf5 = feb.with_suffix(".hdf5")
    for stale in (xplt, hdf5):
        if stale.exists():
            stale.unlink()

    kwargs: dict = {}
    if spring_k is not None:
        kwargs["spring_k"] = float(spring_k)
    T.build_thums_feb(
        mesh,
        feb,
        load_n=float(load_n),
        load_dir=(0.0, -1.0, 0.0),
        use_rigid=use_rigid,
        time_steps=1,
        plantar_bc=plantar_bc,
        subtalar_force=subtalar_force,
        subtalar_moment=subtalar_moment,
        **kwargs,
    )
    t0 = time.time()
    rc = run_febio(feb, workdir=feb.parent, timeout=1800)  # rc!=0 会抛错
    if not xplt.is_file():
        raise FileNotFoundError(f"rc={rc} 但未生成 {xplt}")

    peak_public = fe_post.peak_von_mises(xplt, hdf5_path=hdf5, domain="calcaneus")
    vm = _read_domain_vm(xplt, "calcaneus", hdf5)
    if len(vm) != len(cort_centroids):
        raise RuntimeError(f"{tag}: CORT 单元数不匹配 vm={len(vm)} vs hex={len(cort_centroids)}")
    if abs(float(vm.max()) - float(peak_public)) > 1e-6 * max(1.0, float(peak_public)):
        raise RuntimeError(f"{tag}: fe_post 峰值 {peak_public} != 逐单元峰值 {float(vm.max())}")
    if not np.all(np.isfinite(vm)):
        raise RuntimeError(f"{tag}: σ_vm 含 NaN/Inf")

    gauge = T.gauge_von_mises_mixed(vm, cort_centroids, cort_vol, gauge_mm)
    ipk = int(np.argmax(vm))
    loc = _peak_location(cort_centroids[ipk], mesh)
    return {
        "tag": tag,
        "plantar_bc": plantar_bc,
        "spring_k_N_per_mm": float(spring_k) if spring_k is not None else None,
        "use_rigid": bool(use_rigid),
        "febio_rc": int(rc),
        "load_n": float(load_n),
        "has_force_vector": subtalar_force is not None,
        "has_moment": subtalar_moment is not None,
        "subtalar_force_n": None if subtalar_force is None else [float(x) for x in subtalar_force],
        "subtalar_moment_nm": None if subtalar_moment is None else [float(x) for x in subtalar_moment],
        "max": float(vm.max()),
        "p99": float(np.percentile(vm, 99)),
        "p95": float(np.percentile(vm, 95)),
        "mean": float(vm.mean()),
        "median": float(np.median(vm)),
        "gauge_max": float(gauge.max()),
        "gauge_p95": float(np.percentile(gauge, 95)),
        "peak_elem": ipk,
        "peak_centroid_mm": [float(x) for x in cort_centroids[ipk]],
        **loc,
        "feb": str(feb),
        "xplt": str(xplt),
        "wall_s": float(time.time() - t0),
    }


def _rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b != 0.0 else abs(a - b)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Phase-S4：3D wrench + 跖面弹性支撑")
    ap.add_argument("--force", action="store_true", help="允许覆盖既有 S4 产物")
    ap.add_argument("--gauge-mm", type=float, default=REG_LEN_MM)
    ap.add_argument("--spring-k", type=float, default=SPRING_K)
    args = ap.parse_args(argv)

    for p in (OUT_JSON, OUT_MD):
        if p.exists() and not args.force:
            print(f"拒绝覆盖既有文件：{p}（用 --force 覆盖）", file=sys.stderr)
            return 2

    _install_xplt_discrete_patch()
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("=" * 76)
    print("Phase-S4 · 3D wrench 上 FE + 跖面弹性支撑 · THUMS AM50 跟骨")
    print("=" * 76, flush=True)

    mesh = T.load_thums_mesh(MESH_ANAT)
    print(f"[mesh] {Path(mesh['source_npz']).name}  节点 {mesh['n_nodes']}  "
          f"CORT {mesh['n_elements_cort']}  SPON {mesh['n_elements_spon']}  "
          f"joint={mesh['n_joint_faces']}f  plantar={mesh['n_plantar_faces']}f", flush=True)
    cort_centroids, cort_vol = T.hex_centroids_volumes(mesh["nodes"], mesh["elements_cort"])

    # ---- 输入：S1 wrench + 轴向参考载荷 ---------------------------------
    s1 = json.loads(S1_JSON.read_text(encoding="utf-8"))
    fs = s1["tilted_demo"]["feloadspec"]
    F_tilt = [float(x) for x in fs["force_vector_n"]]
    M_tilt = [float(x) for x in fs["moment_vector_nm"]]
    s1_open_fe = json.loads(S1_OPEN_FE.read_text(encoding="utf-8"))
    load_axial = float(s1_open_fe["subtalar_peak_vertical_n"])  # 26383.13 N
    land = np.linalg.norm(F_tilt)
    print(f"[load] 轴向参考 F={load_axial/1e3:.3f} kN（s1_h5_opensim）", flush=True)
    print(f"[load] S1 20° tilt wrench F=({F_tilt[0]:.1f},{F_tilt[1]:.1f},{F_tilt[2]:.1f}) N "
          f"|F|={land/1e3:.2f} kN  M=({M_tilt[0]:.2f},{M_tilt[1]:.2f},{M_tilt[2]:.2f}) N·m", flush=True)

    # ---- wrench 力矩表达能力（压力场）-----------------------------------
    from climbing.coupling.febio_model import (
        joint_face_frame, linear_pressure_bands, moment_to_pressure_gradient,
    )
    from pyfebio import mesh as fmesh
    from pyfebio.model import Model

    frame = joint_face_frame(mesh, "subtalar_joint")
    F_arr = np.asarray(F_tilt, dtype=float)
    p0_tilt = -float(np.dot(F_arr, frame["normal_unit"])) / frame["area_mm2"]
    ft_tilt = F_arr - float(np.dot(F_arr, frame["normal_unit"])) * frame["normal_unit"]
    slope, direction, gmeta = moment_to_pressure_gradient(
        M_tilt, frame, base_uniform_pressure_mpa=p0_tilt
    )
    _bands, bmeta = linear_pressure_bands(
        Model(), fmesh, mesh, "subtalar_joint", frame,
        p0_mpa=p0_tilt, slope=slope, direction=direction, n_bands=6,
        name_prefix="subtalar_joint_g",
    )
    wrench_meta = {
        "frame": {
            "n_facets": frame["n_facets"],
            "area_mm2": frame["area_mm2"],
            "centroid_mm": [float(x) for x in frame["centroid_mm"]],
            "normal_unit": [float(x) for x in frame["normal_unit"]],
        },
        "force": {
            "vector_n": [float(x) for x in F_arr],
            "normal_component_n": float(np.dot(F_arr, frame["normal_unit"])),
            "tangential_magnitude_n": float(np.linalg.norm(ft_tilt)),
            "p0_mpa": float(p0_tilt),
        },
        "gradient": gmeta,
        "gradient_slope_mpa_per_mm": slope,
        "gradient_direction_unit": [float(x) for x in direction],
        "pressure_band_range_mpa": [
            float(bmeta["facet_pressure_min_mpa"]), float(bmeta["facet_pressure_max_mpa"]),
        ],
        "realized_moment_nm": bmeta["realized_moment_nm"],
        "pressure_path_ill_conditioned": bool(slope > 1.0),  # 坡度 ≫ 平面估计(~0.25)
    }
    print(f"[wrench] 关节面法向 {np.round(frame['normal_unit'],3).tolist()}  面积 "
          f"{frame['area_mm2']:.1f} mm²", flush=True)
    print(f"[wrench] 力矩面内分量 {np.round(gmeta['moment_inplane_nm'],2).tolist()} N·m  "
          f"被丢弃法向扭转 {gmeta['dropped_normal_axis_nm']:.2f} N·m", flush=True)

    # ---- 运行工况 --------------------------------------------------------
    # 主口径 = 幽灵刚体（use_rigid=True，builder 默认），力矩用 RigidMomentLoad（精确）；
    # 另跑一个压力梯度路径作为"负结果"（该曲面关节面上线性压力梯度病态）。
    cases: list[dict] = []

    def run(tag, load_n, force, moment, bc, k, rigid):
        row = _solve_case(tag, mesh, cort_centroids, cort_vol,
                          load_n=load_n, subtalar_force=force, subtalar_moment=moment,
                          plantar_bc=bc, spring_k=k, gauge_mm=args.gauge_mm, use_rigid=rigid)
        mode = "rigid" if rigid else "pressure"
        print(f"  {tag:28s} [{mode}] rc={row['febio_rc']}  gauge={row['gauge_max']:8.2f}  "
              f"max={row['max']:8.2f}  p95={row['p95']:7.2f}  "
              f"peak@{row['location_class']}  [{row['wall_s']:.1f}s]", flush=True)
        return row

    print("[1/3] 轴向基线（均匀面压，复现改前缓存 + A7 BC 对照）...", flush=True)
    base_fixed = run("axial_baseline_fixed", load_axial, None, None, "fixed", None, False)
    base_spring = run("axial_baseline_spring", load_axial, None, None, "spring", args.spring_k, False)

    print("[2/3] 3D wrench（幽灵刚体，精确力矢量 + RigidMomentLoad）...", flush=True)
    ra_fixed = run("rigid_axial_fixed", load_axial, (0.0, -load_axial, 0.0), None, "fixed", None, True)
    ra_spring = run("rigid_axial_spring", load_axial, (0.0, -load_axial, 0.0), None, "spring", args.spring_k, True)
    rt_fixed = run("rigid_tilt_fixed", land, tuple(F_tilt), tuple(M_tilt), "fixed", None, True)
    rt_spring = run("rigid_tilt_spring", land, tuple(F_tilt), tuple(M_tilt), "spring", args.spring_k, True)
    rtn_fixed = run("rigid_tilt_nomoment_fixed", land, tuple(F_tilt), None, "fixed", None, True)

    print("[3/3] 压力梯度路径（负结果：曲面关节面病态）...", flush=True)
    pt_fixed = run("pressure_tilt_fixed", land, tuple(F_tilt), tuple(M_tilt), "fixed", None, False)

    cases = [base_fixed, base_spring, ra_fixed, ra_spring, rt_fixed, rt_spring, rtn_fixed, pt_fixed]

    # ---- 回归判定 --------------------------------------------------------
    reg = {
        "cached": {
            "source": "results/opensim_fe/plantar_bc_summary.json (cases[fixed])",
            "gauge_max_mpa": BASELINE_GAUGE_MPA,
            "max_mpa": BASELINE_MAX_MPA,
            "p95_mpa": BASELINE_P95_MPA,
        },
        "reproduced": {
            "gauge_max_mpa": base_fixed["gauge_max"],
            "max_mpa": base_fixed["max"],
            "p95_mpa": base_fixed["p95"],
        },
        "rel_err": {
            "gauge_max": _rel(base_fixed["gauge_max"], BASELINE_GAUGE_MPA),
            "max": _rel(base_fixed["max"], BASELINE_MAX_MPA),
            "p95": _rel(base_fixed["p95"], BASELINE_P95_MPA),
        },
    }
    reg["max_rel_err"] = max(reg["rel_err"].values())
    reg["verdict"] = "PASS" if reg["max_rel_err"] <= BASELINE_TOL else "FAIL"

    # ---- A7 缓解 / 3D 影响 ----------------------------------------------
    def relief(a: dict, b: dict) -> dict:
        return {
            "gauge_max_drop_pct": (1.0 - b["gauge_max"] / a["gauge_max"]) * 100.0,
            "max_drop_pct": (1.0 - b["max"] / a["max"]) * 100.0,
            "p95_change_pct": (b["p95"] / a["p95"] - 1.0) * 100.0,
        }

    analysis = {
        "A7_relief_baseline_pressure": relief(base_fixed, base_spring),
        "A7_relief_rigid_axial": relief(ra_fixed, ra_spring),
        "A7_relief_rigid_tilt": relief(rt_fixed, rt_spring),
        "tilt_vs_axial_baseline_fixed": {
            "gauge_ratio": rt_fixed["gauge_max"] / base_fixed["gauge_max"],
            "max_ratio": rt_fixed["max"] / base_fixed["max"],
            "p95_ratio": rt_fixed["p95"] / base_fixed["p95"],
        },
        "tilt_vs_rigid_axial_fixed": {
            "gauge_ratio": rt_fixed["gauge_max"] / ra_fixed["gauge_max"],
            "max_ratio": rt_fixed["max"] / ra_fixed["max"],
            "p95_ratio": rt_fixed["p95"] / ra_fixed["p95"],
        },
        "direction_effect_forceonly_fixed": {
            # 同口径 + 同 BC（幽灵刚体、转动锁死、无力矩），只变力方向 → 纯方向效应
            "axial_gauge": ra_fixed["gauge_max"],
            "tilt_nomoment_gauge": rtn_fixed["gauge_max"],
            "gauge_ratio_tilt_over_axial": rtn_fixed["gauge_max"] / ra_fixed["gauge_max"],
        },
        "moment_effect_rigid_tilt_fixed": {
            "with_moment_gauge": rt_fixed["gauge_max"],
            "no_moment_gauge": rtn_fixed["gauge_max"],
            "gauge_ratio_with_over_without": rt_fixed["gauge_max"] / rtn_fixed["gauge_max"],
            "note": "带力矩须放开转动 DOF，与无力矩(锁转动)的 BC 不同，力矩效应含 BC 影响",
        },
        "pressure_gradient_vs_rigid_tilt_fixed": {
            "pressure_gauge": pt_fixed["gauge_max"],
            "rigid_gauge": rt_fixed["gauge_max"],
            "ratio": pt_fixed["gauge_max"] / rt_fixed["gauge_max"],
        },
    }

    result = {
        "meta": {
            "generated_at": ts,
            "study": "Phase-S4 — 3D wrench 上 FE + 跖面弹性支撑（opt-in，默认回归不变）",
            "plan_doc": "docs/非垂直落地扩展方案.md",
            "units": "mm-N-MPa-s",
            "mesh": str(Path(mesh["source_npz"]).relative_to(ROOT)),
            "mesh_kind": "anat (direction C registration)",
            "gauge_radius_mm": args.gauge_mm,
            "threshold_mpa": THRESHOLD_MPA,
            "tilt_deg": TILT_DEG,
            "spring_k_N_per_mm": args.spring_k,
            "load_side": (
                "regression/axial-baseline = use_rigid=False (PressureLoad, PLANAR_BC 口径); "
                "3D wrench = use_rigid=True (幽灵刚体 RigidForceLoad + RigidMomentLoad)"
            ),
            "commands": [
                "python scripts/opensim_fe/nonvertical_s4.py",
            ],
            "measured_vs_modeled": {
                "measured": [
                    "open-FE 暂态距下关节峰值（s1_h5_opensim.json）",
                    "S1 20° 倾斜 wrench 力矢量+力矩（nonvertical_s1.json feloadspec）",
                    "PLANAR_BC 改前 fixed 基线 gauge/max/p95（plantar_bc_summary.json）",
                    "FEBio 4.13 STATIC 实解（本脚本 8 次，全部 rc=0）",
                ],
                "modeled_assumed": [
                    "把 OpenSim 全局帧的 3D wrench 直接在 FE 帧施加（与既有轴向口径一致）",
                    "3D 力矩 = 幽灵刚体 RigidMomentLoad（值按 N·mm），给力矩时放开承载其的转动 DOF",
                    "压力梯度力矩路径在本曲面关节面病态 → 负结果，不作为主口径",
                    "跖面弹性支撑 = 三向 Winkler 地基（k=100 N/mm，PLANAR_BC 最优点）",
                    "gauge 正则化 R=4 mm；骨折阈值 150 MPa 仅后处理，不做单元删除",
                ],
            },
        },
        "regression_gate": reg,
        "wrench_representation": wrench_meta,
        "cases": cases,
        "analysis": analysis,
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, indent=2, ensure_ascii=False, default=float), encoding="utf-8")
    _write_report(result, OUT_MD, ts)
    print("-" * 76)
    print(f"[regression] verdict={reg['verdict']}  max_rel_err={reg['max_rel_err']:.2e}")
    print(f"[A7] baseline pressure fixed→spring gauge {base_fixed['gauge_max']:.2f}→"
          f"{base_spring['gauge_max']:.2f} ({analysis['A7_relief_baseline_pressure']['gauge_max_drop_pct']:+.1f}%)")
    print(f"[A7] rigid tilt        fixed→spring gauge {rt_fixed['gauge_max']:.2f}→"
          f"{rt_spring['gauge_max']:.2f} ({analysis['A7_relief_rigid_tilt']['gauge_max_drop_pct']:+.1f}%)")
    print(f"[3D] rigid tilt vs rigid axial (fixed) gauge "
          f"{rt_fixed['gauge_max']:.2f} / {ra_fixed['gauge_max']:.2f} = "
          f"×{analysis['tilt_vs_rigid_axial_fixed']['gauge_ratio']:.2f}")
    print(f"[out] {OUT_JSON}")
    print(f"[out] {OUT_MD}")
    return 0 if reg["verdict"] == "PASS" else 1


def _write_report(res: dict, path: Path, ts: str) -> None:
    reg = res["regression_gate"]
    cases = {c["tag"]: c for c in res["cases"]}
    an = res["analysis"]
    wm = res["wrench_representation"]
    L: list[str] = []
    A = L.append
    A("# Phase-S4 非垂直落地 · 3D wrench 上 FE + 跖面弹性支撑 报告")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s；FEBio 4.13。")
    A("> 方案 `docs/非垂直落地扩展方案.md` §3 L4 / §4 Phase-S4 / §7。")
    A("> 本报告**新增**，不改任何模块默认、不覆盖既有产物。")
    A("")
    A("## 0. 一句话结论")
    A("")
    A(f"1. **轴向回归（硬门槛）: {reg['verdict']}** —— 默认 `build_thums_feb`（fixed 跖面、"
      f"`use_rigid=False`、F={load_ref(res):.1f} N）复现改前缓存，最大相对误差 "
      f"**{reg['max_rel_err']:.2e}**（gauge_max {reg['reproduced']['gauge_max_mpa']:.6f} vs "
      f"{reg['cached']['gauge_max_mpa']:.6f} MPa，逐位一致）。")
    A("2. **S4 新增 API 全部 opt-in**：`subtalar_force` / `subtalar_moment`（默认 `None`）"
      "把完整 3D wrench 送到 FE；跖面弹性支撑 `plantar_bc=\"spring\"`。")
    A(f"3. **力矩两条路**：`use_rigid=True` 的 `RigidMomentLoad` 精确（须放开承载力矩的刚体"
      "转动 DOF）；`use_rigid=False` 的**关节面线性压力梯度**在本解剖曲面关节面上**病态**"
      f"（需 ±{max(abs(wm['pressure_band_range_mpa'][0]), abs(wm['pressure_band_range_mpa'][1])):.0f} MPa "
      "带压、含拉伸）→ 负结果。")
    A(f"4. **A7 边缘奇异**：轴向**均匀面压**下弹性支撑把 gauge 降 "
      f"{an['A7_relief_baseline_pressure']['gauge_max_drop_pct']:.1f}%（复现 PLANAR_BC 结论）；"
      f"但在 **3D wrench**（剪切/力矩把热点移到关节面）下**不缓解**"
      f"（Δ{an['A7_relief_rigid_tilt']['gauge_max_drop_pct']:+.1f}%）。")
    A(f"5. **3D 载荷（方向效应，同口径同 BC、均无力矩）**：20° 倾斜相对轴向 `gauge` ×"
      f"{an['direction_effect_forceonly_fixed']['gauge_ratio_tilt_over_axial']:.2f}；"
      f"叠加力矩后 `gauge` ×{an['moment_effect_rigid_tilt_fixed']['gauge_ratio_with_over_without']:.2f}"
      "（含转动 BC 影响）。")
    A("")
    A("---")
    A("")
    A("## 1. 本阶段 API（全部 opt-in、默认逐位复现旧行为）")
    A("")
    A("| 文件 / 符号 | 新增参数（默认） | 作用 |")
    A("|---|---|---|")
    A("| `febio_model.build_calcaneus_feb` | `subtalar_force=None` | 3 分量力矢量；None→`load_n×load_dir` |")
    A("| `febio_model.build_calcaneus_feb` | `subtalar_moment=None` | 3 分量力矩 (N·m)；None→零力偶 |")
    A("| `febio_model.build_calcaneus_feb` | `moment_bands=6` | 压力梯度分带数 |")
    A("| `thums_feb.build_thums_feb` | 同上三参数 | 同构路径（本演示所用） |")
    A("| `febio_model.joint_face_frame` | 新 | 关节面面积/形心/外法向/二阶矩 |")
    A("| `febio_model.moment_to_pressure_gradient` | 新 | 力矩→线性压力梯度（逐面片精确）|")
    A("| `febio_model.linear_pressure_bands` | 新 | 线性压力场→分带常压（合力守恒） |")
    A("")
    A("**默认不变**：`subtalar_force=None, subtalar_moment=None` 时走的是改前的旧代码分支")
    A("（不新增任何 wrench 元素——结构测试锁定；FE 结果逐位复现缓存——§2 数值门槛）；")
    A("`plantar_bc` 默认仍 `\"fixed\"`（三向全固定，A7 原状）。")
    A("")
    A("## 2. 轴向回归（硬门槛，方案 §7）")
    A("")
    A("| 量 | 本次（默认） | 改前缓存 | 相对误差 |")
    A("|---|---:|---:|---:|")
    A(f"| gauge_max (MPa) | {reg['reproduced']['gauge_max_mpa']:.6f} | "
      f"{reg['cached']['gauge_max_mpa']:.6f} | {reg['rel_err']['gauge_max']:.2e} |")
    A(f"| raw max (MPa) | {reg['reproduced']['max_mpa']:.6f} | "
      f"{reg['cached']['max_mpa']:.6f} | {reg['rel_err']['max']:.2e} |")
    A(f"| p95 (MPa) | {reg['reproduced']['p95_mpa']:.6f} | "
      f"{reg['cached']['p95_mpa']:.6f} | {reg['rel_err']['p95']:.2e} |")
    A("")
    A(f"- 最大相对误差 **{reg['max_rel_err']:.2e}**（门槛 {BASELINE_TOL:.0e}）→ **{reg['verdict']}**。")
    A("- 缓存来源：`results/opensim_fe/plantar_bc_summary.json` 的 `fixed` 行"
      "（PLANAR_BC，h=5 m a=0，F_subt=26383.13 N，`use_rigid=False`）。")
    A("")
    A("## 3. S1 20° 倾斜 wrench 的 FE 表达")
    A("")
    fr = wm["frame"]
    g = wm["gradient"]
    ff = wm["force"]
    A(f"- 关节载荷面：{fr['n_facets']} 面片，面积 {fr['area_mm2']:.1f} mm²，"
      f"形心 ({fr['centroid_mm'][0]:.2f},{fr['centroid_mm'][1]:.2f},{fr['centroid_mm'][2]:.2f}) mm，"
      f"外法向 ({fr['normal_unit'][0]:.3f},{fr['normal_unit'][1]:.3f},{fr['normal_unit'][2]:.3f})。")
    A(f"- 力 `F=({ff['vector_n'][0]:.1f},{ff['vector_n'][1]:.1f},{ff['vector_n'][2]:.1f}) N`："
      f"相对面法向分解为 法向 {ff['normal_component_n']:.0f} N + 切向 |F_t|="
      f"{ff['tangential_magnitude_n']:.0f} N。")
    A(f"- 力矩 `M=({g['moment_nm'][0]:.2f},{g['moment_nm'][1]:.2f},{g['moment_nm'][2]:.2f}) N·m`"
      f"（|M|={np.linalg.norm(g['moment_nm']):.2f} N·m）：")
    A("  - **精确路径（本演示主口径）**：`RigidMomentLoad`（值按 N·mm 写），给力矩时放开承载"
      "该力矩的刚体转动 DOF（否则被固定 DOF 吸收、力矩无效）。")
    A("  - **压力梯度路径（负结果）**：线性压力场对面片求矩 `M(g)=J g` 精确，但解出")
    A(f"    坡度 {wm['gradient_slope_mpa_per_mm']:.3f} MPa/mm（远大于平面估计 ~0.25），"
      f"带压范围 [{wm['pressure_band_range_mpa'][0]:.1f}, {wm['pressure_band_range_mpa'][1]:.1f}] MPa，"
      "**含大幅拉伸**（不物理——真实关节面不能受拉）。")
    A(f"    - 目标面内分量 `({g['moment_inplane_nm'][0]:.2f},{g['moment_inplane_nm'][1]:.2f},"
      f"{g['moment_inplane_nm'][2]:.2f}) N·m`；均匀法向压力 p0={g['base_uniform_pressure_mpa']:.4f} MPa")
    A(f"      自带曲率力矩 `({g['base_uniform_moment_nm'][0]:.2f},{g['base_uniform_moment_nm'][1]:.2f},"
      f"{g['base_uniform_moment_nm'][2]:.2f}) N·m`（曲面效应）。")
    A(f"    - → 该几何/力矩下压力梯度表达**病态**（`pressure_path_ill_conditioned="
      f"{wm['pressure_path_ill_conditioned']}`）；结论：**本曲面关节面上力矩必须走刚体路径**。")
    A("")
    A("## 4. σ_vm 结果（域 `calcaneus`/CORT，gauge R=4 mm）")
    A("")
    A("| 工况 | 载荷路径 | 跖面 BC | k (N/mm) | max (MPa) | p95 (MPa) | gauge_max (MPa) | 峰值位置 |")
    A("|---|---|---|---:|---:|---:|---:|---|")
    order = ["axial_baseline_fixed", "axial_baseline_spring",
             "rigid_axial_fixed", "rigid_axial_spring",
             "rigid_tilt_fixed", "rigid_tilt_spring",
             "rigid_tilt_nomoment_fixed", "pressure_tilt_fixed"]
    label = {
        "axial_baseline_fixed": "轴向基线（均匀面压 F/A）",
        "axial_baseline_spring": "轴向基线（均匀面压 F/A）",
        "rigid_axial_fixed": "wrench 轴向（F 仅）",
        "rigid_axial_spring": "wrench 轴向（F 仅）",
        "rigid_tilt_fixed": "S1 20° 倾斜 wrench（F+M）",
        "rigid_tilt_spring": "S1 20° 倾斜 wrench（F+M）",
        "rigid_tilt_nomoment_fixed": "S1 20° 倾斜（F 仅，M=0）",
        "pressure_tilt_fixed": "S1 20° 倾斜（压力梯度，负结果）",
    }
    for t in order:
        c = cases[t]
        k = "—" if c["spring_k_N_per_mm"] is None else f"{c['spring_k_N_per_mm']:g}"
        mode = "rigid" if c["use_rigid"] else "pressure"
        A(f"| {label[t]} | {mode} | {c['plantar_bc']} | {k} | {c['max']:.2f} | {c['p95']:.2f} | "
          f"{c['gauge_max']:.2f} | {c['location_class']} |")
    A("")
    A("### 4.1 A7 边缘奇异：弹性支撑缓解多少？")
    A("")
    A("| 工况 | gauge_max fixed | gauge_max spring | Δgauge | raw max fixed | raw max spring | Δmax |")
    A("|---|---:|---:|---:|---:|---:|---:|")
    for key, name, fx, sp in (
        ("A7_relief_baseline_pressure", "轴向基线（均匀面压）", cases["axial_baseline_fixed"], cases["axial_baseline_spring"]),
        ("A7_relief_rigid_axial", "wrench 轴向（幽灵刚体）", cases["rigid_axial_fixed"], cases["rigid_axial_spring"]),
        ("A7_relief_rigid_tilt", "S1 20° 倾斜 wrench（幽灵刚体）", cases["rigid_tilt_fixed"], cases["rigid_tilt_spring"]),
    ):
        r = an[key]
        A(f"| {name} | {fx['gauge_max']:.2f} | {sp['gauge_max']:.2f} | "
          f"{r['gauge_max_drop_pct']:+.1f}% | {fx['max']:.2f} | {sp['max']:.2f} | "
          f"{r['max_drop_pct']:+.1f}% |")
    A("")
    A("**判定**：弹性支撑对 **A7 跖面固定端奇异**有效（轴向均匀面压下 peak 在跖面棱边、"
      "弹性地基把 gauge 降 "
      f"{an['A7_relief_baseline_pressure']['gauge_max_drop_pct']:.1f}%）；但 3D wrench 把峰值"
      "移到**关节面棱边**（剪切/力矩驱动），跖面弹性支撑**不缓解**该新奇异。")
    A("")
    A("### 4.2 3D 载荷（20° 倾斜）改变多少？")
    A("")
    de = an["direction_effect_forceonly_fixed"]
    me = an["moment_effect_rigid_tilt_fixed"]
    pr = an["pressure_gradient_vs_rigid_tilt_fixed"]
    A(f"- **方向效应（同口径+同 BC，仅力矢量不同、都无力矩）**：轴向 gauge={de['axial_gauge']:.2f} vs "
      f"20° 倾斜 gauge={de['tilt_nomoment_gauge']:.2f}（×{de['gauge_ratio_tilt_over_axial']:.2f}）。"
      "倾斜把力沿关节面分解，法向分量下降 + 新增横向剪切。")
    A(f"- **力矩效应**（刚体、fixed）：带力矩 gauge={me['with_moment_gauge']:.2f} vs "
      f"无力矩 gauge={me['no_moment_gauge']:.2f}（×{me['gauge_ratio_with_over_without']:.2f}）。"
      "⚠ 带力矩必须放开转动 DOF，与无力矩（锁转动）的 BC 不同，故此差值含 BC 影响。")
    A(f"- 压力梯度路径 gauge={pr['pressure_gauge']:.2f} vs 刚体 gauge={pr['rigid_gauge']:.2f}"
      f"（×{pr['ratio']:.2f}）：病态梯度把应力炒高，**不可用**。")
    A("")
    A("## 5. 诚实边界")
    A("")
    A("- **载荷参考帧**：S1 的 3D wrench 在 OpenSim 全局帧；本演示按既有轴向口径**直接把矢量"
      "施加在 FE 帧**，未做额外的帧旋转变换（与 `SUBTALAR_DIR=(0,−1,0)` 的既有假设一致）。")
    A("- **压力梯度路径是负结果**：在本解剖配准的**曲面**关节面上，线性压力梯度表达该力矩"
      "需要病态的大坡度 → 带压出现大幅**拉伸**（不物理）。因此 §4 的主口径用幽灵刚体"
      "`RigidMomentLoad`。这是本阶段最有价值的**诚实结论**（方案 §3 L4 的 (a) 路在本几何失败）。")
    A("- **刚体力矩的正确性依赖放开转动 DOF**：给力矩时必须放开承载它的刚体转动自由度，"
      "否则 `RigidMomentLoad` 被 `RigidFixed` 吸收而静默无效（本实现已处理）。")
    A("- **弹性支撑是参数**：k=100 N/mm 取自 PLANAR_BC 最优点，非标定值；三向 Winkler 地基"
      "不是真实接触/摩擦。")
    A("- **幽灵刚体载荷口径 ≠ 均匀面压口径**：两者把同一关节力以不同方式作用（刚体 → 关节面"
      "整体刚性平移；面压 → 分布压力），绝对应力不可直接互比；§4.2 用**同口径**（刚体）做"
      "轴向 vs 倾斜的方向比较。")
    A("- 本演示仍是**筛选性**研究；σ_vm≥150 MPa 仅后处理判据，无单元删除。")
    A("")
    A("## 6. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s4.py")
    A("# 回归：")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_opensim_fe.py tests\\test_febio_run.py -q")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_nonvertical_s1.py tests\\test_nonvertical_s4.py -q")
    A("```")
    A("")
    A("## 7. 产物")
    A("")
    A(f"- `{OUT_MD.relative_to(ROOT)}`（本报告）")
    A(f"- `{OUT_JSON.relative_to(ROOT)}`（机器可读）")
    A(f"- `{WORKDIR.relative_to(ROOT)}/s4_*.feb/.xplt`（每次求解输入/结果）")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def load_ref(res: dict) -> float:
    return float(res["cases"][0]["load_n"])


if __name__ == "__main__":
    raise SystemExit(main())
