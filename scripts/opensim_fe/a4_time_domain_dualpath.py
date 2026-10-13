# -*- coding: utf-8 -*-
"""A4 时间域双路径度量（预注册实现）—— 见 docs/methods/A4_时间域双路径_预注册.md

判据：Δt(h) = t_arr(lumbar, h) − t_arr(ankle_r, h) 是否随高度【变号】。
  Δt > 0 ⇒ 脊柱晚于踝 ⇒ 自下而上传导
  Δt < 0 ⇒ 脊柱早于踝 ⇒ 头惯量（自上而下）路径

链路/预设与 g4_load_chain_compliant.py 逐字相同（T6 冻结），只改「取首达时刻」。
纯新增；不改 coupling/、不改任何冻结输入。
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

OUT_JSON = ROOT / "results" / "opensim_fe" / "a4_time_domain_dualpath.json"

MODEL_MASS_KG = 75.337          # 与 g4_load_chain_compliant.py 一致
HEIGHTS = list(range(1, 51))
JOINTS = ["ankle_r", "knee_r", "hip_r", "lumbar"]   # T7 冻结
ALPHA = 0.05                    # T2 冻结
ALPHA_SENS = [0.02, 0.10]       # 只报，不判


def _t_arr(t: np.ndarray, f: np.ndarray, alpha: float) -> float | None:
    """首个 |f(t)| >= alpha * max|f| 的时刻（T2）。"""
    a = np.abs(np.asarray(f, float))
    m = float(a.max())
    if not np.isfinite(m) or m <= 0:
        return None
    idx = int(np.argmax(a >= alpha * m))
    return float(t[idx]) if a[idx] >= alpha * m else None


def _one_height(h: int) -> dict:
    from climbing.coupling.opensim_fall import run_dead_drop, PassiveStiffness
    from climbing.coupling.opensim_grf import ground_reaction
    from climbing.coupling.joint_reactions import joint_reaction

    grf = ground_reaction(height_m=float(h), mass_kg=MODEL_MASS_KG)
    fall = run_dead_drop(grf, passive=PassiveStiffness.nominal_landing())

    out: dict = {"h_m": int(h), "grf_peak_n": float(grf.peak_total_n)}
    for j in JOINTS:
        jr = joint_reaction(fall, grf, joint=j)
        t = np.asarray(jr.t_s, float)
        fz = np.asarray(jr.force_n, float)[:, 1]         # 竖直分量
        out[f"peak_{j}_n"] = float(np.max(np.abs(fz)))
        out[f"t_arr_{j}_s"] = _t_arr(t, fz, ALPHA)
        for a in ALPHA_SENS:
            out[f"t_arr_{j}_a{a:g}"] = _t_arr(t, fz, a)
    ta, tl = out.get("t_arr_ankle_r_s"), out.get("t_arr_lumbar_s")
    out["dt_s"] = (tl - ta) if (ta is not None and tl is not None) else None
    for a in ALPHA_SENS:
        taa, tla = out.get(f"t_arr_ankle_r_a{a:g}"), out.get(f"t_arr_lumbar_a{a:g}")
        out[f"dt_a{a:g}_s"] = (tla - taa) if (taa is not None and tla is not None) else None
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--heights", default="", help="逗号分隔，默认 1..50")
    ap.add_argument("--out-suffix", default="")
    ap.add_argument("--reuse", action="store_true")
    args = ap.parse_args(argv)

    hs = ([int(x) for x in args.heights.split(",") if x.strip()]
          if args.heights else HEIGHTS)
    out_path = OUT_JSON if not args.out_suffix else \
        OUT_JSON.with_name(OUT_JSON.stem + args.out_suffix + OUT_JSON.suffix)

    done: dict[int, dict] = {}
    if args.reuse and out_path.exists():
        prev = json.loads(out_path.read_text(encoding="utf-8"))
        done = {int(r["h_m"]): r for r in prev.get("rows", [])}
        print(f"reuse: {len(done)} 行")

    for k, h in enumerate(hs, 1):
        if h in done:
            continue
        t0 = time.time()
        r = _one_height(h)
        done[h] = r
        dt = r["dt_s"]
        print(f"[{k}/{len(hs)}] h={h:2d} m  GRF={r['grf_peak_n']/1e3:7.1f} kN  "
              f"t_arr ankle={(r['t_arr_ankle_r_s'] or float('nan'))*1e3:7.3f} ms  "
              f"lumbar={(r['t_arr_lumbar_s'] or float('nan'))*1e3:7.3f} ms  "
              f"dt={(dt*1e3 if dt is not None else float('nan')):+8.3f} ms  "
              f"({time.time()-t0:4.1f}s)", flush=True)
        rows = [done[k2] for k2 in sorted(done)]
        out_path.write_text(json.dumps(
            {"alpha": ALPHA, "alpha_sens": ALPHA_SENS, "joints": JOINTS,
             "chain": "ground_reaction -> run_dead_drop(passive=nominal_landing) -> joint_reaction",
             "rows": rows}, indent=2, ensure_ascii=False), encoding="utf-8")

    rows = [done[k2] for k2 in sorted(done)]
    dts = np.array([r["dt_s"] for r in rows if r["dt_s"] is not None], float)
    if dts.size:
        print(f"\n=== Δt 汇总（n={dts.size}）===")
        print(f"  范围 {dts.min()*1e3:+.3f} … {dts.max()*1e3:+.3f} ms")
        print(f"  变号: {'有' if (dts.min() < 0 < dts.max()) else '无'}"
              f"   （Δt>0: {(dts>0).sum()}  Δt<0: {(dts<0).sum()}  =0: {(dts==0).sum()}）")
    print(f"\n产物: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
