"""S1 垂直切片入口：GRF → OpenSim 正动力学 → 距下关节反力 → FEBio von Mises。

方案 ``docs/OpenSim_FE复现方案.md`` §五 S1。用法::

    python scripts/opensim_fe/s1_vertical_slice.py --height 5
    python scripts/opensim_fe/s1_vertical_slice.py --height 5 --skip-fe   # 只跑多体半段

产物写到 ``results/opensim_fe/s1_h<H>.json``（多体段）与
``results/opensim_fe/s1_h<H>_fe.json``（有限元段）。

⚠️ 这是**研究管线**：绝对应力值不是临床预测（见方案 §六 D5、README 已知局限）。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))


def _dump(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=float), encoding="utf-8")


def run_opensim_half(height: float, *, accuracy: float = 1e-4, max_step_s: float = 2e-5):
    """S1.1–S1.3：GRF → FD → 距下关节反力。"""
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=height)
    if not grf.impulse_ok:
        raise RuntimeError(
            f"T1 冲量不变量失败：∫Fdt/(m·v0)-1 = {grf.impulse_rel_err*100:+.1f}% "
            f"(阈值 10%)。停，先查 GRF 标定/左右分配。"
        )
    fall = run_dead_drop(grf, accuracy=accuracy, max_step_s=max_step_s)
    jr = subtalar_reaction(fall, grf, side="r")
    return grf, fall, jr


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="S1 垂直切片（跟骨 @ 刚性地面）")
    ap.add_argument("--height", type=float, default=5.0, help="释放高度 (m)")
    ap.add_argument("--skip-fe", action="store_true", help="只跑多体半段")
    ap.add_argument("--load-n", type=float, default=None,
                    help="覆盖传给 FE 的载荷 (N)；默认取距下关节峰值纵向力")
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "opensim_fe")
    args = ap.parse_args(argv)

    print(f"=== S1 垂直切片：跟骨 @ {args.height:g} m 刚性地面 ===")
    grf, fall, jr = run_opensim_half(args.height)

    print(f"[S1.1] GRF  峰值 {grf.peak_total_n/1e3:6.1f} kN  "
          f"窗口 {grf.window_s[-1]*1e3:5.1f} ms  v0 {grf.impact_speed_ms:5.2f} m/s")
    print(f"[T1 ] 冲量  ∫Fdt={grf.impulse_ns:7.1f} N·s  m·v0={grf.momentum_ns:7.1f} N·s  "
          f"err {grf.impulse_rel_err*100:+5.1f}%  {'✅' if grf.impulse_ok else '❌'}")
    print(f"[S1.2] FD   {fall.n_states} 帧 / {fall.t_s[-1]*1e3:.1f} ms")
    print(f"[S1.3] 距下关节 peak|F|={jr.peak_force_n/1e3:.2f} kN  "
          f"纵向={jr.peak_vertical_n/1e3:.2f} kN @ {jr.peak_vertical_time_s*1e3:.2f} ms  "
          f"peak|M|={jr.peak_moment_nm:.2f} N·m")

    multibody = {
        "height_m": args.height,
        "grf_peak_n": grf.peak_total_n,
        "grf_window_s": list(grf.window_s),
        "impact_speed_ms": grf.impact_speed_ms,
        "impulse_ns": grf.impulse_ns,
        "momentum_ns": grf.momentum_ns,
        "impulse_rel_err": grf.impulse_rel_err,
        "t1_ok": bool(grf.impulse_ok),
        "fd_frames": fall.n_states,
        "subtalar_peak_force_n": jr.peak_force_n,
        "subtalar_peak_vertical_n": jr.peak_vertical_n,
        "subtalar_peak_vertical_time_s": jr.peak_vertical_time_s,
        "subtalar_peak_moment_nm": jr.peak_moment_nm,
    }
    mb_path = args.outdir / f"s1_h{args.height:g}_opensim.json"
    _dump(mb_path, multibody)
    print(f"      -> {mb_path}")

    load_n = args.load_n if args.load_n is not None else jr.peak_vertical_n
    if args.skip_fe:
        print(f"[S1.4/S1.5] 跳过 FE。载荷 = {load_n/1e3:.2f} kN（传给 meshing_driver --load-n）")
        return 0

    # ---- S1.4/S1.5：有限元半段 --------------------------------------------
    # 直接复用经过验证的 meshing_driver（VTP→tet→.feb→FEBio→σ_vm），
    # 避免在本文件重复/漂移 FE 接口。
    import re
    import subprocess

    driver = ROOT / "scripts" / "opensim_fe" / "meshing_driver.py"
    geom = ROOT / "model" / "Model" / "Geometry" / "calcaneus_r.stl"
    cmd = [sys.executable, str(driver), "--vtp", str(geom),
           "--use-pressure", "--load-n", f"{load_n:.0f}"]
    print(f"[S1.4/S1.5] {' '.join(str(c) for c in cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT))
    tail = "\n".join((proc.stdout or "").splitlines()[-14:])
    print(tail)
    if proc.returncode != 0:
        print(f"[S1.5] ❌ FE 半边失败 (rc={proc.returncode})")
        if proc.stderr:
            print(proc.stderr[-800:])
        return proc.returncode

    m = re.search(r"\[stats\] p95=([0-9.]+) p99=([0-9.]+) mean=([0-9.]+) max=([0-9.]+)",
                  proc.stdout or "")
    p95 = p99 = mean = mx = None
    if m:
        p95, p99, mean, mx = (float(x) for x in m.groups())
    print(f"[S1.5] load={load_n/1e3:.2f} kN -> σ_vm p95={p95} MPa"
          f"  (p99={p99}, mean={mean}, 峰值={mx})")
    _dump(args.outdir / f"s1_h{args.height:g}_fe.json",
          {"height_m": args.height, "load_n": load_n,
           "von_mises_p95_mpa": p95, "von_mises_p99_mpa": p99,
           "von_mises_mean_mpa": mean, "von_mises_max_mpa": mx,
           "note": "p95 为稳健判据；峰值受面压/固定边奇异支配，未收敛"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
