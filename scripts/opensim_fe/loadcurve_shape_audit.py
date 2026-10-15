# -*- coding: utf-8 -*-
"""载荷链变更的**形状**审计：归档 vs 当前，载荷-高度曲线是否只差一个常数？

**动机**：`STALE_LOADCHAIN_AUDIT.md` 发现 `now/存档` 的 gauge 比 = **1.0431 … 1.1202**，
且 **h=8 = 0.9746 (<1)** —— **若标定只是常数缩放，比值必恒定** ⇒ 形状变了。

用法::

    .venv\\Scripts\\python.exe scripts\\opensim_fe\\loadcurve_shape_audit.py
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(ROOT / "src"))

import s2_thums_sweep as S2                            # noqa: E402

R = ROOT / "results" / "opensim_fe"
ARCH = ['1', '2', '3', '4', '5', '6', '8']


def main() -> int:
    b = json.load(io.open(R / 's2_thums_height_sweep.json', encoding='utf-8'))
    arch = {float(x['h']): x for x in b if float(x.get('act', 0)) == 0.0}
    print('=== 归档（标定前）vs 当前 `_load_at` ===')
    print('%-6s %-14s %-14s %-10s %-14s %-14s' % (
        'h', '归档 subtalar_kN', '当前 subtalar_kN', '比', '归档 F_sub(N)?', '当前 F_sub(N)'))
    rows = []
    for h in sorted(arch):
        a = float(arch[h]['subtalar_kN'])
        f_sub, t1 = S2._load_at(h)
        cur = f_sub / 1000.0
        rows.append((h, a, cur, cur / a))
        extra = {k: arch[h][k] for k in ('achilles_kN', 't1_rel_err')
                 if k in arch[h]}
        print('%-6g %-14.4f %-14.4f %-10.4f %-14s %s' % (
            h, a, cur, (cur / a), '', extra))
    print()
    print('--- 比值是否恒定？ ---')
    rr = [r[3] for r in rows]
    print('  比值:', ['%.4f' % x for x in rr])
    print('  极差: %.4f   相对跨度: %.2f%%' % (
        max(rr) - min(rr), 100 * (max(rr) - min(rr)) / (sum(rr) / len(rr))))
    print('  ⇒ %s' % ('**非常数** ⇒ 载荷表形状变了（不是缩放）' if
                      (max(rr) - min(rr)) > 0.02 else '近似常数 ⇒ 疑似缩放'))
    print()
    print('--- 曲线形状（归一化到 h=2）---')
    print('%-6s %-14s %-14s' % ('h', '归档/归档@2', '当前/当前@2'))
    a2 = arch[2.0]['subtalar_kN']
    c2 = next(r[2] for r in rows if r[0] == 2.0)
    for h, a, c, _ in rows:
        print('%-6g %-14.4f %-14.4f' % (h, a / a2, c / c2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
