"""S3.2b · 跟腱力对跟骨首次骨折高度的影响（刚性地面，参数化激活 a∈[0,1]）。

方案 ``docs/OpenSim_FE复现方案.md`` §五 S3.2b。用法::

    python scripts/opensim_fe/s3_achilles_sweep.py --heights 5 10 20 --acts 0 1
    python scripts/opensim_fe/s3_achilles_sweep.py                 # 默认网格

设计
----
* **网格只算一次**；每个 (高度, 激活) 只改载荷重解。
* 高度 → 跑 OpenSim 正动力学 → 取距下关节峰值纵向力 ``F_subt``（真跑，不外推）。
* 跟腱力 ``F_ach = a · F_max(triceps)``（见 ``load_transfer``），a=0 即无跟腱基线。
* 判据：``gauge_max``（过程区 ℓ=4 mm 正则化峰值）越过 150 MPa 的最低高度 = 首次骨折高度。
* 报告**每个激活 a 的首次骨折高度**，直接检验「补跟腱能否把高度拉回 [Y25] 7–9 m」。

⚠️ 网格走 ``vtp_to_tet``（体素修复）—— gmsh 平滑路径单骨 >5 万 tet，会触发
``pyfebio.to_hdf5`` 的 HDF5 attribute 上限（见 LEARNINGS-027 / S3.0）。
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

DEFAULT_GRID = [1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 40, 50]
THRESHOLD_MPA = 150.0
REG_LEN_MM = 4.0


def _mesh_once(char_len: float):
    from climbing.coupling import meshing

    stl = ROOT / "model" / "Model" / "Geometry" / "calcaneus_r.stl"
    mesh = meshing.vtp_to_tet(stl, char_len, out_msh=meshing.MESH_DIR / f"s3_{char_len:g}mm.msh")
    print(f"[mesh] method={mesh['method']} nodes={mesh['n_nodes']:,} tet={mesh['n_tets']:,} "
          f"achilles_faces={mesh['n_achilles_faces']} vol_err={mesh['volume_error_pct']:.2f}% "
          f"minJ={mesh['min_jacobian']:.2e}", flush=True)
    return mesh


def _load_at(height: float) -> tuple[float, float]:
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=height)
    fall = run_dead_drop(grf)
    jr = subtalar_reaction(fall, grf, side="r")
    return float(jr.peak_vertical_n), float(grf.impulse_rel_err)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="S3.2b 跟腱力 × 高度扫描")
    ap.add_argument("--heights", type=float, nargs="+", default=None)
    ap.add_argument("--acts", type=float, nargs="+", default=[0.0, 1.0],
                    help="三头肌激活水平列表 (0..1)")
    ap.add_argument("--char-len", type=float, default=REG_LEN_MM)
    ap.add_argument("--gauge-mm", type=float, default=REG_LEN_MM)
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "opensim_fe")
    args = ap.parse_args(argv)
    heights = args.heights if args.heights else DEFAULT_GRID

    from climbing.coupling.fe_post import von_mises_stats, von_mises_gauge
    from climbing.coupling.febio_run import run_febio
    from climbing.coupling.load_transfer import transfer

    args.outdir.mkdir(parents=True, exist_ok=True)
    workdir = (ROOT / "temp" / "opensim_fe")
    workdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    mesh = _mesh_once(args.char_len)
    print(f"[S3.2b] 高度×激活 = {len(heights)}×{len(args.acts)}  a∈{args.acts}", flush=True)

    rows: list[dict] = []
    for h in heights:
        f_sub, t1 = _load_at(h)
        for a in args.acts:
            spec = transfer(f_sub, achilles_activation=a)
            feb = (workdir / f"s3_h{h:g}_a{a:g}.feb").resolve()
            spec.build(mesh, feb)
            run_febio(feb)
            xplt = feb.with_suffix(".xplt")
            st = von_mises_stats(xplt)
            gg = von_mises_gauge(xplt, mesh, radius_mm=args.gauge_mm)
            row = {"h": float(h), "act": float(a), "subtalar_kn": f_sub / 1e3,
                   "achilles_kn": spec.achilles_n / 1e3, "t1_rel_err": t1,
                   "max": st["max"], "p99": st["p99"], "p95": st["p95"], "mean": st["mean"],
                   "gauge_max": gg["gauge_max"]}
            rows.append(row)
            print(f"h={h:5.1f} a={a:3.1f}  F_sub={f_sub/1e3:6.2f}kN F_ach={spec.achilles_n/1e3:6.2f}kN  "
                  f"max={st['max']:7.1f} p95={st['p95']:6.1f} mean={st['mean']:5.2f}  "
                  f"gauge={gg['gauge_max']:6.2f} MPa  ({time.time()-t0:.0f}s)", flush=True)

    (args.outdir / "s3_achilles_sweep.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False, default=float), encoding="utf-8")
    with (args.outdir / "s3_achilles_sweep.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("\n=== S3.2b 报告：各激活 a 的首次骨折高度（gauge_max > 150 MPa）===", flush=True)
    for a in args.acts:
        hs = [r["h"] for r in rows if r["act"] == a and r["gauge_max"] > THRESHOLD_MPA]
        print(f"  a={a:3.1f} (F_ach_max≈{a*10885/1e3:.1f}kN): 首次骨折高度 = "
              f"{hs[0] if hs else '未触发(≤%.0fm)' % max(heights)} m", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
