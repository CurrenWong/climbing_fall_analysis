# -*- coding: utf-8 -*-
"""A1 聚类边界复核：Δ=2.0 是「真差异」还是「报边界口径」产物？

判据（冻结）：|Δ| < 2.0 m，目标边界 = [7.0, 16.0] m（[Y25]）。
本脚本**不改判据、不改数据**，只把「边界如何从离散相报成一个数」的三种约定并列算 Δ，
看 A1 的 FAIL 是否对约定敏感。纯只读（读 s4_alignment_bc2.json）。
"""
from __future__ import annotations

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC = ROOT / "results" / "opensim_fe" / "s4_alignment_bc2.json"

TARGETS = [7.0, 16.0]   # [Y25]
TOL = 2.0               # 判据：|Δ| < 2.0 m


def main() -> int:
    data = json.loads(SRC.read_text(encoding="utf-8"))
    a1 = data["A1"]
    phases = a1["phases"]

    # 网格分辨率 = 相邻高度差（应为 1 m）
    hs = sorted(phases[0]["heights"])
    step = min(b - a for a, b in zip(hs, hs[1:])) if len(hs) > 1 else None

    n = len(phases)
    print(f"判据: |Δ| < {TOL} m;  目标边界 = {TARGETS} m  ([Y25])")
    print(f"相数 = {n};  高度网格步长 = {step} m")
    for p in phases:
        r = p["range_m"]
        print(f"  P{p['phase']}: {r[0]}–{r[1]} m  ({p.get('label','')})")
    print()

    # 第 i 个相变发生在 phases[i] 与 phases[i+1] 之间
    rows = []
    for i in range(n - 1):
        lo = phases[i]["range_m"][1]        # 相内最后一格
        hi = phases[i + 1]["range_m"][0]    # 次相第一格
        mid = 0.5 * (lo + hi)               # 相变中点
        tgt = TARGETS[i]
        convs = {
            "末格 (冻结口径, 脚本现用)": lo,
            "相变中点": mid,
            "次相首格": hi,
        }
        print(f"=== 边界 #{i+1}  目标 {tgt} m   相变区间 [{lo}, {hi}] m ===")
        for name, val in convs.items():
            d = abs(val - tgt)
            verdict = "PASS" if d < TOL else "FAIL"
            print(f"  {name:<26} 报 {val:>5.1f} m   Δ = {d:>4.1f}   -> {verdict}")
        d_lo = abs(lo - tgt)
        d_mid = abs(mid - tgt)
        rows.append((i + 1, tgt, lo, hi, d_lo, d_mid, abs(hi - tgt)))
        print()

    print("=== 结论判定 ===")
    for idx, tgt, lo, hi, d_lo, d_mid, d_hi in rows:
        span = hi - lo
        if d_lo >= TOL and d_mid < TOL:
            tag = ("**口径产物**：Δ 从 %.1f 降到 %.1f ⇒ 相变位置在判据容差内，"
                   "FAIL 系「取末格」约定所致" % (d_lo, d_mid))
        elif d_lo >= TOL and d_mid >= TOL:
            excess = d_mid - TOL
            if excess <= span:
                tag = ("**真差异但仅超 %.1f m（≤ 1 格）**：三种约定下 Δ ∈ [%.1f, %.1f]，"
                       "均 ≥ %.1f" % (excess, d_hi, d_lo, TOL))
            else:
                tag = ("**真差异**：三种约定下 Δ ∈ [%.1f, %.1f] 全 ≥ %.1f，"
                       "口径挪不动" % (d_hi, d_lo, TOL))
        else:
            tag = "PASS（所有约定）"
        print(f"  边界 #{idx}（目标 {tgt:.0f} m）：{tag}")

    # 汇总：任一约定下能否双双 PASS
    print()
    for name, pick in (("末格", lambda lo, hi: lo),
                       ("中点", lambda lo, hi: (lo + hi) / 2),
                       ("次相首格", lambda lo, hi: hi)):
        ds = [abs(pick(phases[i]["range_m"][1], phases[i + 1]["range_m"][0]) - TARGETS[i])
              for i in range(n - 1)]
        ok = all(d < TOL for d in ds)
        print(f"  约定「{name}」⇒ Δ = {[round(d,1) for d in ds]} ⇒ A1 = "
              f"{'PASS' if ok else 'FAIL'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
