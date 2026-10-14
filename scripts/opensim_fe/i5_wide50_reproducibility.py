# -*- coding: utf-8 -*-
"""I4b · `wide50` 非唯一性红旗的重复解检查。

同一组合（`share=0.2241` + `wide50`，以及作为对照的 `share=0.10` + `wide50`）
在 **3 个独立进程**里各解一次，逐次比对 `sigma_h_median_mpa` 与 G3 OR。

判据（事先写明）：
* 3 次**逐值一致** ⇒ 红旗清除 ⇒ L2 可去掉
* 出现分叉（相对差 ≫ 数值噪声）⇒ `wide50` 作废 ⇒ 退回 (a) 结论、L2 保留
"""
from __future__ import annotations

import json
import pathlib

import numpy as np

RES = pathlib.Path("results/opensim_fe")
RUNS = [
    ("bandscan", "g9_fibula_fe_support_bandscan.json"),
    ("rep2", "g9_fibula_fe_support_rep2.json"),
    ("rep3", "g9_fibula_fe_support_rep3.json"),
]
WIN = [1.1774, 2.1866]
# 数值噪声基准：G9 自己在 long95 上的两次同条件重跑差 0.3%（I″ 实测 10.53 vs 10.5）
NOISE = 0.005  # 0.5%


def main() -> int:
    docs = {}
    for tag, fn in RUNS:
        p = RES / fn
        if not p.is_file():
            print(f"[{tag}] 缺失：{fn}")
            continue
        docs[tag] = json.loads(p.read_text(encoding="utf-8"))
    if len(docs) < 2:
        print("可用样本不足，无法比对")
        return 1

    keys = sorted({k for d in docs.values() for k in d["combos"]})
    print(f"{'combo':<26} {'n':>2} {'σ 首':>9} {'σ 末':>9} "
          f"{'最大相对差':>11} {'OR 极差':>10}  判据")
    verdict = {}
    for k in keys:
        sigs, ors = [], []
        for tag, d in docs.items():
            c = d["combos"].get(k)
            if not c or not c.get("sigma_h_median_mpa"):
                continue
            sigs.append((tag, np.asarray(c["sigma_h_median_mpa"], float)))
            g3 = (d.get("g3") or {}).get(k, {}).get("fibula", {})
            if g3.get("OR") is not None:
                ors.append(g3["OR"])
        if len(sigs) < 2:
            print(f"{k:<26} {len(sigs):>2}   （样本不足）")
            continue
        ref = sigs[0][1]
        worst = 0.0
        for _, s in sigs[1:]:
            rel = np.abs(s - ref) / np.maximum(np.abs(ref), 1e-9)
            worst = max(worst, float(rel.max()))
        orr = (max(ors) - min(ors)) if len(ors) >= 2 else None
        ok = worst <= NOISE
        verdict[k] = {"n_samples": len(sigs), "max_rel_diff_pct": worst * 100,
                      "or_range": orr, "reproducible": bool(ok)}
        print(f"{k:<26} {len(sigs):>2} {ref[0]:>9.3f} {ref[-1]:>9.3f} "
              f"{worst * 100:>10.3f}% "
              f"{'--' if orr is None else f'{orr:.4f}':>10}  "
              f"{'✅ 可复现' if ok else '❌ 分叉'}")

    out = RES / "I5_WIDE50_REPRODUCIBILITY.json"
    out.write_text(json.dumps({
        "piece": "I4b · wide50 重复解可复现性检查",
        "runs": [t for t, _ in RUNS if t in docs],
        "noise_baseline_pct": NOISE * 100,
        "noise_basis": "G9 long95 同条件重跑差约 0.3%（I″：10.53 vs 10.5 MPa）",
        "combos": verdict,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n写出 {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())