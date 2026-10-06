"""THUMS AM50 calcaneus two-material-domain FEBio driver (Option C).

Builds the registered THUMS calcaneus (mm, OpenSim ``calcn_r`` frame) into a
two-domain FEBio model -- cortical (CORT, hex8, ``"calcaneus"``) + trabecular
(SPON, tet4, ``"trabecular"``) -- solves it, and reports the peak von Mises in
the fracture-relevant **cortical** domain.

Run from the repo root with the project venv::

    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\thums_feb_run.py

Outputs (under ``temp/opensim_fe/thums_calcaneus/``)::

    thums_calcaneus.feb          two-domain FEBio model
    thums_calcaneus.xplt         FEBio results
    thums_feb_result.json        peak sigma_vm, load, domain, counts, stats
    THUMS_FEB_REPORT.md          human-readable report

The load path mirrors ``febio_model.build_calcaneus_feb``: ghost rigid body +
``RigidForceLoad`` on ``subtalar_joint`` (PressureLoad fallback), optional
Achilles ``TractionLoad``, ``plantar`` fixed, single STATIC step.

This is a *screening* study: the absolute stress is not a clinical prediction.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import thums_feb as T  # noqa: E402

OUT_DIR = ROOT / "temp" / "opensim_fe" / "thums_calcaneus"
FEB = OUT_DIR / "thums_calcaneus.feb"
XPLT = OUT_DIR / "thums_calcaneus.xplt"
HDF5 = OUT_DIR / "thums_calcaneus.hdf5"
RESULT_JSON = OUT_DIR / "thums_feb_result.json"
REPORT_MD = OUT_DIR / "THUMS_FEB_REPORT.md"

#: cortical yield (compression) reference from the [Y25] strength table.
CORT_YIELD_MPA = 150.0


def _read_domain_vm(xplt: Path, domain: str, hdf5: Path | None = None) -> np.ndarray:
    """Read the last-state element von Mises stress for ``domain`` (MPa).

    Converts ``.xplt`` -> HDF5 first if ``hdf5`` is missing.  Reusing the HDF5
    already written by :func:`fe_post.peak_von_mises` avoids a second parse.
    """
    import h5py

    from climbing.coupling import fe_post
    from pyfebio import xplt as xplt_mod

    h5p = Path(hdf5) if hdf5 is not None else xplt.with_suffix(".hdf5")
    if not h5p.is_file():
        xplt_mod.to_hdf5(str(xplt), str(h5p))
    with h5py.File(h5p, "r") as h5:
        stress, used = fe_post._read_last_stress(h5, domain)
    print(f"[post] domain={domain!r} -> xplt block {used!r} "
          f"({stress.shape[0]} elements)")
    return fe_post.von_mises_from_voigt(stress)


def _nearest_named_surface(
    point: np.ndarray, mesh: dict, names: tuple[str, ...] = ("subtalar_joint", "plantar")
) -> tuple[str, float]:
    best_name, best_d = "other", float("inf")
    for name in names:
        ns = np.asarray(mesh["node_sets"].get(name, []), dtype=np.int64)
        if len(ns) == 0:
            continue
        d = float(np.linalg.norm(mesh["nodes"][ns] - point, axis=1).min())
        if d < best_d:
            best_name, best_d = name, d
    return best_name, best_d


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="THUMS 双材料域跟骨 FEBio（CORT+SPON）")
    ap.add_argument("--load-n", type=float, default=4000.0, help="距下关节合力 (N)")
    ap.add_argument("--load-dir", type=float, nargs=3, default=(0.0, -1.0, 0.0),
                    metavar=("X", "Y", "Z"), help="载荷方向（FE 坐标）")
    ap.add_argument("--no-rigid", action="store_true", help="直接用 PressureLoad (p=F/A)")
    ap.add_argument("--achilles-n", type=float, default=0.0, help="跟腱牵引力 (N)")
    ap.add_argument(
        "--plantar-bc",
        choices=("fixed", "roller", "roller_free", "spring", "contact"),
        default="fixed",
        help="跖面 BC（默认 fixed，逐位保留历史行为）；contact = 波3 最小版刚性壁接触",
    )
    ap.add_argument("--time-steps", type=int, default=1, help="STATIC 载荷步数")
    ap.add_argument("--gauge-radius", type=float, default=5.0, help="正则化半径 (mm)")
    args = ap.parse_args(argv)

    from climbing.coupling import fe_post
    from climbing.coupling.febio_run import FebioRunError, run_febio

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 68)
    print("THUMS AM50 跟骨 · 双材料域 FEBio（CORT hex8 + SPON tet4）")
    print("=" * 68)

    mesh = T.load_thums_mesh()
    print(f"[mesh] 节点 {mesh['n_nodes']}  "
          f"CORT hex8={mesh['n_elements_cort']} (E=15000)  "
          f"SPON tet4={mesh['n_elements_spon']} (E=73.4)")
    print(f"[mesh] 体积 CORT {mesh['volume_cort_mm3']:.1f} mm³ + "
          f"SPON {mesh['volume_spon_mm3']:.1f} mm³ = "
          f"{mesh['volume_cort_mm3']+mesh['volume_spon_mm3']:.1f} mm³")
    print(f"[mesh] 边界面 {mesh['n_boundary_faces']} 片（三角化 {mesh['n_boundary_tris']}）  "
          f"关节 {mesh['n_joint_faces']} / 跖面 {mesh['n_plantar_faces']} / "
          f"跟腱 {mesh['n_achilles_faces']}")
    print(f"[mesh] node_sets: subtalar_joint={mesh['n_joint_nodes']} "
          f"plantar={mesh['n_plantar_nodes']} achilles={mesh['n_achilles_nodes']}")
    if mesh["n_achilles_nodes"] == 0:
        print("[mesh] 注意：THUMS-calcnframe 下 meshing.py 的跟腱框选 "
              "(x<-25 & y>0) 为空（后结节仅到 x≈-9 mm）-> 不施加跟腱牵引。")

    use_rigid = not args.no_rigid
    for _stale in (XPLT, HDF5):
        if _stale.exists():
            _stale.unlink()

    feb = T.build_thums_feb(
        mesh, FEB,
        load_n=args.load_n, load_dir=tuple(args.load_dir),
        use_rigid=use_rigid, time_steps=args.time_steps,
        achilles_n=args.achilles_n, plantar_bc=args.plantar_bc,
    )
    load_mode = "rigid" if use_rigid else "pressure"
    print(f"[bc] plantar_bc={args.plantar_bc}")

    try:
        rc = run_febio(feb, workdir=feb.parent, timeout=900)
    except FebioRunError as exc:
        if args.no_rigid:
            print(f"[feb] 面压模式仍失败：\n{str(exc)[-1500:]}")
            raise
        print("[feb] ⚠ 幽灵刚体模式未收敛，回退到 PressureLoad (p=F/A) 重试。")
        print(str(exc)[-800:])
        load_mode = "pressure"
        if XPLT.exists():
            XPLT.unlink()
        feb = T.build_thums_feb(
            mesh, FEB,
            load_n=args.load_n, load_dir=tuple(args.load_dir),
            use_rigid=False, time_steps=args.time_steps,
            achilles_n=args.achilles_n, plantar_bc=args.plantar_bc,
        )
        rc = run_febio(feb, workdir=feb.parent, timeout=900)

    print(f"[feb] FEBio 退出码 {rc}  ({load_mode})")
    if not XPLT.is_file():
        raise FileNotFoundError(f"FEBio 退出码 {rc} 但未生成 {XPLT}")

    # --- cortical (fracture-relevant) post-processing ---------------------
    # Public headline number via the shared fe_post entry point (domain
    # "calcaneus"); its HDF5 is then reused to read the full per-element field.
    peak_public = fe_post.peak_von_mises(XPLT, hdf5_path=HDF5, domain="calcaneus")
    print(f"[post] fe_post.peak_von_mises(domain='calcaneus') = {peak_public:.3f} MPa")

    vm = _read_domain_vm(XPLT, "calcaneus", hdf5=HDF5)
    vm_spon = _read_domain_vm(XPLT, "trabecular", hdf5=HDF5)
    cort_centroids, cort_vol = T.hex_centroids_volumes(
        mesh["nodes"], mesh["elements_cort"]
    )
    if len(vm) != len(cort_centroids):
        raise RuntimeError(
            f"CORT 单元数不匹配：vm={len(vm)} vs hex={len(cort_centroids)}"
        )

    peak = float(vm.max())
    if abs(peak - float(peak_public)) > 1e-6 * max(1.0, peak):
        raise RuntimeError(
            f"fe_post.peak_von_mises={peak_public} 与逐单元峰值 {peak} 不一致"
        )
    ipk = int(np.argmax(vm))
    peak_centroid = cort_centroids[ipk]
    peak_region, peak_region_dist = _nearest_named_surface(peak_centroid, mesh)
    gauge = T.gauge_von_mises_mixed(vm, cort_centroids, cort_vol, args.gauge_radius)

    stats = {
        "max": peak,
        "p99": float(np.percentile(vm, 99)),
        "p95": float(np.percentile(vm, 95)),
        "mean": float(vm.mean()),
        "median": float(np.median(vm)),
        "n_elem": int(vm.size),
    }
    stats_spon = {
        "max": float(vm_spon.max()),
        "p95": float(np.percentile(vm_spon, 95)),
        "mean": float(vm_spon.mean()),
        "n_elem": int(vm_spon.size),
    }

    # Matched-load comparison against the old homogeneous whole-foot pipeline.
    # S1 @ h=5 m: load 26383 N -> max 179.8 / p95 14.97 (E=15000 whole bone).
    OLD = {"load_n": 26383.13, "max": 179.807, "p95": 14.972, "p99": 30.269}
    scale = args.load_n / OLD["load_n"]
    comparison = {
        "old_pipeline": {
            "note": "S1 h=5 m, homogeneous E=15000 whole-foot mesh, load via PressureLoad",
            **OLD,
        },
        "matched_load_scaling": {
            "load_n": args.load_n,
            "old_max_at_this_load_mpa": OLD["max"] * scale,
            "old_p95_at_this_load_mpa": OLD["p95"] * scale,
            "new_max_mpa": peak,
            "new_p95_mpa": stats["p95"],
            "max_ratio_new_over_old": peak / (OLD["max"] * scale),
            "p95_ratio_new_over_old": stats["p95"] / (OLD["p95"] * scale),
        },
    }

    result = {
        "domain": "calcaneus",
        "domain_role": "CORT (cortical, fracture-relevant)",
        "load_n": float(args.load_n),
        "load_dir": [float(x) for x in args.load_dir],
        "load_mode": load_mode,
        "achilles_n": float(args.achilles_n),
        "febio_rc": int(rc),
        "materials": {
            "bone_cort": {"E_mpa": 15000.0, "nu": 0.3},
            "bone_spon": {"E_mpa": 73.4, "nu": 0.45},
        },
        "counts": {
            "n_nodes": mesh["n_nodes"],
            "n_elements_cort_hex8": mesh["n_elements_cort"],
            "n_elements_spon_tet4": mesh["n_elements_spon"],
            "n_nodes_cort": mesh["n_nodes_cort"],
            "n_nodes_spon": mesh["n_nodes_spon"],
            "volume_cort_mm3": mesh["volume_cort_mm3"],
            "volume_spon_mm3": mesh["volume_spon_mm3"],
            "n_joint_faces": mesh["n_joint_faces"],
            "n_plantar_faces": mesh["n_plantar_faces"],
            "n_achilles_faces": mesh["n_achilles_faces"],
            "n_joint_nodes": mesh["n_joint_nodes"],
            "n_plantar_nodes": mesh["n_plantar_nodes"],
            "n_achilles_nodes": mesh["n_achilles_nodes"],
        },
        "peak_von_mises_mpa": peak,
        "peak_occurrence": {
            "element_index": ipk,
            "centroid_mm": [float(x) for x in peak_centroid],
            "nearest_surface": peak_region,
            "nearest_surface_dist_mm": float(peak_region_dist),
        },
        "cort_stats_mpa": stats,
        "spon_stats_mpa": stats_spon,
        "gauge": {
            "radius_mm": float(args.gauge_radius),
            "gauge_max_mpa": float(gauge.max()),
            "gauge_p99_mpa": float(np.percentile(gauge, 99)),
            "raw_max_mpa": peak,
        },
        "cort_yield_mpa": CORT_YIELD_MPA,
        "peak_over_yield": peak / CORT_YIELD_MPA,
        "gauge_over_yield": float(gauge.max()) / CORT_YIELD_MPA,
        "comparison": comparison,
        "artifacts": {
            "feb": str(FEB),
            "xplt": str(XPLT),
            "hdf5": str(HDF5),
            "source_npz": mesh["source_npz"],
        },
        "note": (
            "screening study; peak is affected by the joint-load/plantar-fix edge "
            "singularity -- use gauge_max as the mesh-robust strength indicator"
        ),
    }
    RESULT_JSON.write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # --- report -----------------------------------------------------------
    _write_report(mesh, result)

    print("-" * 68)
    print(f"峰值 σ_vm (域 'calcaneus'/CORT) : {peak:.3f} MPa")
    print(f"  位置: 单元 #{ipk} 质心 {np.round(peak_centroid,2).tolist()} mm, "
          f"最近面 {peak_region} ({peak_region_dist:.1f} mm)")
    print(f"  p99={stats['p99']:.2f}  p95={stats['p95']:.2f}  "
          f"mean={stats['mean']:.2f} MPa")
    print(f"  正则化 gauge_max @ r={args.gauge_radius:g} mm = "
          f"{float(gauge.max()):.2f} MPa")
    print(f"  皮质屈服参考 {CORT_YIELD_MPA:.0f} MPa "
          f"-> peak/yield={peak/CORT_YIELD_MPA:.2f}, "
          f"gauge/yield={float(gauge.max())/CORT_YIELD_MPA:.2f}")
    print(f"[out] {FEB}")
    print(f"[out] {XPLT}")
    print(f"[out] {RESULT_JSON}")
    print(f"[out] {REPORT_MD}")
    return 0


def _write_report(mesh: dict, r: dict) -> None:
    c = r["counts"]
    st = r["cort_stats_mpa"]
    cmp_ = r["comparison"]["matched_load_scaling"]
    cmp_ratio = cmp_["p95_ratio_new_over_old"]
    if 0.5 <= cmp_ratio <= 2.0:
        verdict = "同一量级（比值 0.5–2）"
    elif cmp_ratio > 2.0:
        verdict = "偏高"
    else:
        verdict = "偏低"
    lines = [
        "# THUMS AM50 跟骨双材料域 FEBio 报告（Option C）",
        "",
        f"- 输入网格：`{Path(mesh['source_npz']).name}`（mm，OpenSim `calcn_r` 帧，"
        f"原点 = 距下关节中心）",
        f"- 求解：FEBio 4.13 STATIC，退出码 **{r['febio_rc']}**，加载方式 **{r['load_mode']}**",
        f"- 载荷：{r['load_n']:.1f} N × ({r['load_dir'][0]:.3f},"
        f"{r['load_dir'][1]:.3f},{r['load_dir'][2]:.3f})"
        + (f" + 跟腱 {r['achilles_n']:.1f} N" if r["achilles_n"] else ""),
        "",
        "## 材料域",
        "",
        "| 域 | 名称 | 类型 | 单元 | E (MPa) | ν |",
        "|---|---|---|---|---|---|",
        f"| CORT（皮质，断裂相关） | `calcaneus` | hex8 | {c['n_elements_cort_hex8']} | "
        f"15000 | 0.30 |",
        f"| SPON（松质） | `trabecular` | tet4 | {c['n_elements_spon_tet4']} | 73.4 | 0.45 |",
        "",
        f"节点 {c['n_nodes']}（CORT 用 {c['n_nodes_cort']}，SPON 用 {c['n_nodes_spon']}，"
        "两域共用节点/共形界面）。"
        f"体积：CORT {c['volume_cort_mm3']:.1f} + SPON {c['volume_spon_mm3']:.1f} = "
        f"{c['volume_cort_mm3']+c['volume_spon_mm3']:.1f} mm³。",
        "",
        "## 具名面 / 节点集（复用 `meshing.py` 阈值）",
        "",
        f"- `subtalar_joint`：{c['n_joint_faces']} 面 / {c['n_joint_nodes']} 节点"
        "（dist<25 mm 且 n·(+Y)>0.2）",
        f"- `plantar`：{c['n_plantar_faces']} 面 / {c['n_plantar_nodes']} 节点"
        "（n·Y<−0.5，已剔除关节节点）",
        f"- `achilles`：{c['n_achilles_faces']} 面 / {c['n_achilles_nodes']} 节点"
        "（x<−25 且 y>0）",
        "",
        "> **跟腱面为空**：配准后的 THUMS 跟骨后结节最小 x ≈ −9 mm（见 "
        "`REGISTRATION_REPORT.md`），因此 `meshing.py` 的 x<−25 框选不到任何面，"
        "跟腱牵引量为 0（与参考 builder 默认 `achilles_n=0` 一致）。这是坐标现实，"
        "不是缺省省略。",
        "",
        "## 结果：皮质域（`domain=\"calcaneus\"`）von Mises",
        "",
        f"- **峰值 σ_vm = {r['peak_von_mises_mpa']:.3f} MPa**",
        f"  - 出现于单元 #{r['peak_occurrence']['element_index']}，质心 "
        f"{[round(x,2) for x in r['peak_occurrence']['centroid_mm']]} mm，"
        f"最近具名面 `{r['peak_occurrence']['nearest_surface']}`"
        f"（{r['peak_occurrence']['nearest_surface_dist_mm']:.1f} mm）",
        f"- 分布：p99={st['p99']:.2f}，p95={st['p95']:.2f}，mean={st['mean']:.2f}，"
        f"median={st['median']:.2f} MPa（{st['n_elem']} 单元）",
        f"- 过程区正则化 gauge_max @ r={r['gauge']['radius_mm']:g} mm = "
        f"**{r['gauge']['gauge_max_mpa']:.2f} MPa**（gauge_p99="
        f"{r['gauge']['gauge_p99_mpa']:.2f}）",
        "",
        f"屈服参考（[Y25] 皮质压缩）{r['cort_yield_mpa']:.0f} MPa → "
        f"peak/yield={r['peak_over_yield']:.2f}，gauge/yield={r['gauge_over_yield']:.2f}。",
        "",
        "## 与旧管线（整体足网格，均匀 E=15000）对比",
        "",
        f"- 旧 S1 @ h=5 m：载荷 {r['comparison']['old_pipeline']['load_n']:.0f} N → "
        f"max={r['comparison']['old_pipeline']['max']:.1f}、"
        f"p95={r['comparison']['old_pipeline']['p95']:.1f} MPa（均匀 E=15000 整足网格）",
        f"- 同一载荷 {cmp_['load_n']:.0f} N 按线弹性线性折算：旧 max≈"
        f"{cmp_['old_max_at_this_load_mpa']:.2f}、p95≈{cmp_['old_p95_at_this_load_mpa']:.2f} MPa。",
        f"- 新模型：max={cmp_['new_max_mpa']:.2f}、p95={cmp_['new_p95_mpa']:.2f} MPa；"
        f"比值 max×{cmp_['max_ratio_new_over_old']:.2f}、p95×{cmp_['p95_ratio_new_over_old']:.2f}。",
        f"- 判定：p95 比值 {cmp_ratio:.2f} → **{verdict}**。几何（真实 AM50 跟骨 vs "
        "粗化整足）与材料分区（CORT/SPON）都改变绝对值，故只作量级核对，"
        "不作为一致性命中。",
        "",
        "## 诚实说明",
        "",
        "- 峰值受**加载面/固定面棱边奇异**支配，网格加密会发散；本网格为 THUMS 原生",
        "  分辨率（3957 单元），未做加密，故 `max` 仅作诊断，`gauge_max`（过程区正则化）",
        "  才是有物理意义的强度指标。",
        "- 本报告不做单元删除 / 断裂判据，只求解并报告应力；σ_vm≥150 MPa 的判据在",
        "  后处理中应用。",
        "- 载荷大小是**参数**（默认 4000 N），不是物理结论。",
        "",
        "## 产物",
        "",
        f"- `{Path(r['artifacts']['feb']).name}`",
        f"- `{Path(r['artifacts']['xplt']).name}`",
        f"- `{RESULT_JSON.name}`",
        f"- `{REPORT_MD.name}`",
    ]
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
