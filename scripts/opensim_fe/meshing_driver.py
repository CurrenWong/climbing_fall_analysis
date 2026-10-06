"""S1 · 跟骨 FE 半边可复现入口：VTP → tet → .feb → FEBio → 峰值 σ_vm。

用法（从仓库根，用项目 venv）::

    python scripts/opensim_fe/meshing_driver.py --load-n 4000
    python scripts/opensim_fe/meshing_driver.py --load-n 4000 --char-len 3.0
    python scripts/opensim_fe/meshing_driver.py --load-n 4000 --use-pressure  # 面压兜底

流程
----
1. ``climbing.coupling.meshing.vtp_to_tet``：读 ``r_foot.vtp``（calcn_r 几何）→
   mm tet 网格 + 具名面（``subtalar_joint`` / ``plantar``），写
   ``temp/opensim_fe/meshes/calcaneus_<char>mm.msh``。
2. 网格自检：T7（体积误差 < 2%）、T5（min Jacobian > 0）。失败即报，不带病求解。
3. ``build_calcaneus_feb``：写 ``temp/opensim_fe/calcaneus.feb``（mm–N–MPa–s）。
   默认"幽灵刚体"加载；若 FEBio 不收敛则自动改用面压 p=F/A。
4. ``run_febio`` 求解（非 0 退出码会抛错，绝不静默）。
5. ``peak_von_mises`` 读 ``.xplt`` → 峰值 σ_vm（MPa），并写节点平均 ``.vtu``。

⚠️ ``--load-n`` 是**参数**（默认 4000 N 竖直），不是物理结论。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

DEFAULT_VTP = ROOT / "model" / "opensim" / "FullBodyModel-4.0" / "Geometry" / "r_foot.vtp"
OUT_DIR = ROOT / "temp" / "opensim_fe"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="跟骨 S1 FE 半边（VTP→tet→feb→σ_vm）")
    ap.add_argument("--load-n", type=float, default=4000.0, help="距下关节合力 (N)，默认 4000")
    ap.add_argument("--char-len", type=float, default=3.0, help="目标单元尺寸 (mm)")
    ap.add_argument("--vtp", type=Path, default=DEFAULT_VTP, help="输入 VTP 路径")
    ap.add_argument("--E-mpa", type=float, default=15000.0, help="杨氏模量 (MPa)，皮质骨量级")
    ap.add_argument("--nu", type=float, default=0.3, help="泊松比")
    ap.add_argument("--time-steps", type=int, default=1,
                    help="STATIC 载荷步数（默认 1：多步会触发 FEBio 相对收敛判据的机器精度问题）")
    ap.add_argument("--use-pressure", action="store_true",
                    help="强制用面压 p=F/A（跳过幽灵刚体）")
    ap.add_argument("--no-vtu", action="store_true", help="不写节点平均 .vtu")
    args = ap.parse_args(argv)

    from climbing.coupling.febio_run import FebioRunError, run_febio
    from climbing.coupling.febio_model import build_calcaneus_feb, joint_area_mm2
    from climbing.coupling.fe_post import von_mises_stats
    from climbing.coupling.meshing import MESH_DIR, vtp_to_tet

    MESH_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"=== 跟骨 S1 FE 半边 · 载荷 {args.load_n:.0f} N (竖直) ===")

    # --- 1) 网格 ---------------------------------------------------------
    msh_path = MESH_DIR / f"calcaneus_{args.char_len:g}mm.msh"
    mesh = vtp_to_tet(args.vtp, args.char_len, out_msh=msh_path)

    # --- 2) 自检 ---------------------------------------------------------
    t7_ok = mesh["volume_error_pct"] < 2.0
    t5_ok = mesh["min_jacobian"] > 0.0
    print(f"[mesh] 方法={mesh['method']}  节点={mesh['n_nodes']:,}  tet={mesh['n_tets']:,}")
    print(f"[mesh] 关节面 {mesh['n_joint_faces']} 片 / 跖面 {mesh['n_plantar_faces']} 片"
          f"  （面积 {joint_area_mm2(mesh):.1f} mm²）")
    print(f"[T7 ] 体积 {mesh['volume_mm3']:.1f} vs 参考 {mesh['reference_volume_mm3']:.1f} mm³  "
          f"误差 {mesh['volume_error_pct']:.2f}%  {'✅' if t7_ok else '❌'}")
    print(f"[T5 ] min scaled Jacobian {mesh['min_jacobian']:.2e}  "
          f"（质心<0.3 占比 {mesh['jacobian_frac_lt_0_3']*100:.1f}%）"
          f"  {'✅>0' if t5_ok else '❌<=0'}")
    if not (t7_ok and t5_ok):
        print("[STOP] 网格自检未过；按方案 §七 先修网格，不进入求解。")
        return 1

    # --- 3) .feb ---------------------------------------------------------
    feb_path = OUT_DIR / "calcaneus.feb"
    xplt_path = feb_path.with_suffix(".xplt")
    if xplt_path.exists():
        xplt_path.unlink()  # 防呆：避免读到上一次的结果

    used = "pressure" if args.use_pressure else "rigid"
    feb = build_calcaneus_feb(
        mesh,
        feb_path,
        E_mpa=args.E_mpa,
        nu=args.nu,
        load_n=args.load_n,
        use_rigid=not args.use_pressure,
        time_steps=args.time_steps,
    )

    # --- 4) 求解（刚体不收敛则回退面压）---------------------------------
    try:
        rc = run_febio(feb, workdir=feb.parent)
    except FebioRunError as exc:
        if args.use_pressure:
            print(f"[feb] ❌ 面压模式仍失败：{str(exc)[:400]}")
            raise
        print("[feb] ⚠ 幽灵刚体模式未收敛，回退到 PressureLoad (p=F/A) 重试。")
        used = "pressure"
        feb = build_calcaneus_feb(
            mesh, feb_path, E_mpa=args.E_mpa, nu=args.nu, load_n=args.load_n,
            use_rigid=False, time_steps=args.time_steps,
        )
        rc = run_febio(feb, workdir=feb.parent)

    # --- 5) 后处理 -------------------------------------------------------
    if not xplt_path.is_file():
        raise FileNotFoundError(f"FEBio 退出码 {rc} 但未生成 {xplt_path}")
    vtu = None if args.no_vtu else OUT_DIR / "calcaneus_vm.vtu"
    stats = von_mises_stats(xplt_path, mesh=mesh, write_vtu=vtu)
    sigma = stats["max"]  # 峰值仅作诊断；**不是**强度判据（见下）

    # 阈值参考（[Y25] 材料强度，皮质/压缩 150 MPa）
    try:
        from climbing.bone import MATERIAL_STRENGTH_MPA

        thr_c = MATERIAL_STRENGTH_MPA["calcaneus"][1]
    except Exception:  # noqa: BLE001 - 阈值只是参考，不阻断报告
        thr_c = float("nan")

    print("=== 报告 ===")
    print(f"VTP            : {args.vtp}")
    print(f"方法           : mesh={mesh['method']}  load={used}")
    print(f"命令           : python scripts/opensim_fe/meshing_driver.py --load-n {args.load_n:.0f}")
    print(f"文件           : {msh_path}")
    print(f"                 {feb}")
    print(f"                 {xplt_path}")
    if vtu:
        print(f"                 {vtu}")
    print(f"网格           : {mesh['n_nodes']:,} 节点 / {mesh['n_tets']:,} tet"
          f"  vol_err {mesh['volume_error_pct']:.2f}%  minJ {mesh['min_jacobian']:.2e}")
    print(f"FEBio 退出码   : {rc}")
    print(f"σ_vm 分布      : 峰值(诊断)={stats['max']:.2f}  p99={stats['p99']:.2f}  "
          f"p95={stats['p95']:.2f}  mean={stats['mean']:.2f}  median={stats['median']:.2f} MPa")
    print(f"稳健判据(p95)  : {stats['p95']:.2f} MPa  vs 阈值 {thr_c:.0f} MPa (压缩, [Y25])  ->  "
          f"{'低于阈值' if stats['p95'] < thr_c else '超过阈值'}")
    print(f"峰值 σ_vm      : {stats['max']:.3f} MPa")  # 峰值：仅供诊断/向后兼容解析
    print(f"[stats] p95={stats['p95']:.3f} p99={stats['p99']:.3f} "
          f"mean={stats['mean']:.3f} max={stats['max']:.3f}")
    failed = [n for n, ok in (("T5", t5_ok), ("T7", t7_ok)) if not ok]
    print(f"失败不变量     : {', '.join(failed) if failed else '无'}")
    print("说明           : 关节面取 calcn 原点(距下关节中心; talus 系 "
          "(-48.77,-41.95,7.92)mm 换算) 附近朝上的上关节面；")
    print("                 载荷 --load-n 是**参数**（默认 4000 N 竖直），不是物理结论；")
    print("                 关节面为 voxel 修复网格，T5 仅保证 Jacobian>0，非高质量网格。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
