"""S2 · 跟骨全高度扫描（刚性地面）：h = 1..50 m，取 σ_vm(h)。

方案 ``docs/OpenSim_FE复现方案.md`` §五 S2。用法::

    python scripts/opensim_fe/s2_height_sweep.py                    # 默认高度网格, cl=3
    python scripts/opensim_fe/s2_height_sweep.py --heights 5        # 冒烟测试
    python scripts/opensim_fe/s2_height_sweep.py --char-len 4.0     # 粗网格（峰值为正则化峰值）
    python scripts/opensim_fe/s2_height_sweep.py --heights 1 2 3 5 7 10 15 20 30 50

设计
----
* **网格只算一次**（几何与 E 不随 h 变），每个高度只改载荷 `p=F/A` 重解 ——
  省掉每高度重网格化。
* 每个高度仍**真跑** OpenSim 正动力学（不靠 √h 外推），以核对 T9 单调性。
* 判据同时给 **max / p95 / p99**：max 是边界奇异（见 LEARNINGS-023），
  p95/p99 是稳健量；首次骨折高度对三者分别算，并如实标注对网格的依赖。
* 材料是**力控线弹性** → σ_vm 与 E 无关，只随载荷（∝ 冲击速度 ∝ √h）走。

产物：``results/opensim_fe/s2_height_sweep[_cl<L>].json`` / ``.csv``。
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

#: 默认高度网格（1–50 m，靠近 [Y25] 的 7–9 m 边界处加密）。
DEFAULT_GRID = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50]
#: 跟骨压缩强度阈值（MPa，[Y25] 量级）。
THRESHOLD_MPA = 150.0
#: **正则化长度 ℓ（mm）** —— 骨内应力的物理下限尺度，不是任意网格参数：
#: 松质骨 RVE ≈ 5 mm（Harrigan 1988）；皮质骨断裂过程区 ≈ 1–5 mm。
#: 网格单元尺寸取 ℓ，使裸峰值成为"ℓ 尺度正则化峰值"（见 LEARNINGS-023）。
REG_LEN_MM = 4.0


def _mesh_once(char_len: float):
    from climbing.coupling import meshing

    stl = ROOT / "model" / "Model" / "Geometry" / "calcaneus_r.stl"
    mesh = meshing.vtp_to_tet(stl, char_len, out_msh=meshing.MESH_DIR / f"calcaneus_{char_len:g}mm.msh")
    print(f"[mesh] method={mesh['method']}  nodes={mesh['n_nodes']:,}  tet={mesh['n_tets']:,}  "
          f"vol_err={mesh['volume_error_pct']:.2f}%  minJ={mesh['min_jacobian']:.2e}", flush=True)
    return mesh


def _load_at(height: float) -> tuple[float, float]:
    """跑多体半段 → 返回 (距下关节纵向峰值力 N, T1 相对误差)。"""
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=height)
    fall = run_dead_drop(grf)
    jr = subtalar_reaction(fall, grf, side="r")
    return float(jr.peak_vertical_n), float(grf.impulse_rel_err)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="S2 跟骨全高度扫描（刚性地面）")
    ap.add_argument("--heights", type=float, nargs="+", default=None, help="高度列表 (m)")
    ap.add_argument("--char-len", type=float, default=REG_LEN_MM,
                    help=f"FE 单元尺寸 (mm)；默认 = 正则化长度 ℓ={REG_LEN_MM} mm")
    ap.add_argument("--gauge-mm", type=float, default=REG_LEN_MM,
                    help=f"gauge 体积平均半径 (mm)；默认 = ℓ={REG_LEN_MM}")
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "opensim_fe")
    args = ap.parse_args(argv)
    heights = args.heights if args.heights else DEFAULT_GRID

    from climbing.coupling import febio_model
    from climbing.coupling.fe_post import von_mises_gauge, von_mises_stats
    from climbing.coupling.febio_run import run_febio

    args.outdir.mkdir(parents=True, exist_ok=True)
    feb = (ROOT / "temp" / "opensim_fe" / "s2.feb").resolve()
    feb.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    mesh = _mesh_once(args.char_len)
    print(f"[S2.1] 扫 {len(heights)} 个高度：{heights}", flush=True)

    rows: list[dict] = []
    for h in heights:
        f_n, t1 = _load_at(h)
        febio_model.build_calcaneus_feb(mesh, feb, load_n=f_n, use_rigid=False, time_steps=1)
        run_febio(feb, workdir=feb.parent)
        xplt = feb.with_suffix(".xplt")
        st = von_mises_stats(xplt)
        gg = von_mises_gauge(xplt, mesh, radius_mm=args.gauge_mm)
        row = {"h": float(h), "load_kn": f_n / 1e3, "t1_rel_err": t1,
               "max": st["max"], "p99": st["p99"], "p95": st["p95"], "mean": st["mean"],
               "gauge_max": gg["gauge_max"], "gauge_p99": gg["gauge_p99"]}
        rows.append(row)
        print(f"h={h:5.1f} m  F={f_n/1e3:6.2f} kN  max={st['max']:7.2f}  "
              f"p99={st['p99']:6.2f}  p95={st['p95']:6.2f}  "
              f"gauge(R{args.gauge_mm:g})={gg['gauge_max']:6.2f} MPa  T1={t1*100:+.1f}%  "
              f"({time.time()-t0:.0f}s)", flush=True)

    tag = f"_cl{args.char_len:g}"
    (args.outdir / f"s2_height_sweep{tag}.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False, default=float), encoding="utf-8")
    with (args.outdir / f"s2_height_sweep{tag}.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("=== S2 报告 ===", flush=True)
    for key in ("max", "p99", "p95", "gauge_max"):
        frac = [r["h"] for r in rows if r[key] > THRESHOLD_MPA]
        # 单调性（容差 2%）：σ_vm(h) 应该随 h 非降
        drops = sum(1 for a, b in zip(rows, rows[1:]) if b[key] < a[key] * 0.98)
        print(f"[{key:>9}] 首次 >{THRESHOLD_MPA:.0f} MPa 高度 = "
              f"{frac[0] if frac else '未触发(≤50 m)'} m  |  单调性 {len(rows)-1-drops}/{len(rows)-1} 段非降")
    print(f"[meta] mesh cl={args.char_len:g}mm (=ℓ={REG_LEN_MM}mm 正则化长度)  tet={mesh['n_tets']:,}  "
          f"minJ={mesh['min_jacobian']:.2e}  vol_err={mesh['volume_error_pct']:.2f}%")
    print(f"[meta] gauge 半径 R={args.gauge_mm:g}mm（体积加权）。max = ℓ 尺度正则化峰值（操作判据）；"
          f"gauge(R) = 更强平滑的下界；裸峰值随网格发散（LEARNINGS-023）。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
