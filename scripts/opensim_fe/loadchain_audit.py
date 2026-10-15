# -*- coding: utf-8 -*-
"""载荷链审计 CLI（闸门 10 的日常工具）。

扫描 ``results/opensim_fe/`` 下所有 ``*sweep*.json``，把每份的**载荷-高度表**与
**当前链**（``s2_thums_sweep._load_at``）比指纹，报告：

* **✅ 同链**（指纹一致）
* **⚠️ 异链·常数缩放**（比值恒定）—— 说明只差一个整体系数
* **⛔ 异链·形状不同**（比值随高度变）—— 说明是**另一张表**（本仓库的实际情况）

用法::

    .venv\\Scripts\\python.exe scripts\\opensim_fe\\loadchain_audit.py
    .venv\\Scripts\\python.exe scripts\\opensim_fe\\loadchain_audit.py --json out.json
"""
from __future__ import annotations

import argparse
import glob
import io
import json
import os
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for p in (str(_HERE), str(ROOT / "src")):
    if p not in sys.path:
        sys.path.insert(0, p)

import s2_thums_sweep as S2                            # noqa: E402
from climbing.coupling.loadchain import (              # noqa: E402
    compare, fingerprint, table_from_json)

RESULTS = ROOT / "results" / "opensim_fe"
#: 当前链载荷表的缓存（``_load_at`` 每高度跑一次 OpenSim，≈ 数十秒 ⇒ 必须缓存）
CACHE = RESULTS / "LOADCHAIN_CURRENT.json"


def current_table(heights, *, refresh: bool = False) -> dict[float, float]:
    """当前链的 ``{h: F_sub}``（N）。**带磁盘缓存**（``LOADCHAIN_CURRENT.json``）。

    缓存是增量写入的：算完一个高度就落盘一次 ⇒ 中断后重跑可续。
    """
    cache: dict[float, float] = {}
    if CACHE.is_file() and not refresh:
        try:
            cache = {float(k): float(v)
                     for k, v in json.load(io.open(CACHE, encoding="utf-8")).items()}
        except Exception:                              # noqa: BLE001
            cache = {}
    todo = [float(h) for h in heights if float(h) not in cache]
    for h in todo:
        v = float(S2._load_at(h)[0])
        cache[h] = v
        io.open(CACHE, "w", encoding="utf-8").write(
            json.dumps({("%g" % k): cache[k] for k in sorted(cache)},
                       ensure_ascii=False, indent=2))
        print('  [cache] h=%g -> %.4f N  (%d/%d)' % (
            h, v, len(cache), len(cache) + len(todo) - todo.index(h) - 1))
    return cache


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", type=str, default="*sweep*.json")
    ap.add_argument("--act", type=float, default=0.0)
    ap.add_argument("--json", type=str, default=None)
    args = ap.parse_args(argv)

    files = sorted(glob.glob(str(RESULTS / args.glob)))
    print('扫描 %d 份（%s）' % (len(files), args.glob))

    # 先读出各文件的表，收集高度并集；再算当前链
    tabs = {}
    for f in files:
        try:
            tabs[os.path.basename(f)] = table_from_json(f, act=args.act)
        except Exception as exc:                       # noqa: BLE001
            print('  [skip] %s : %s' % (os.path.basename(f), exc))
    hs = sorted({h for t in tabs.values() for h in t})
    print('高度并集: %s' % ['%g' % h for h in hs])
    cur = current_table(hs)
    fp_cur = fingerprint(cur)
    print('当前链指纹 = %s\n' % fp_cur)

    rows = []
    print('%-52s %-18s %-14s %s' % ('文件', '载荷指纹', '比值跨度', '判定'))
    print('-' * 108)
    for name, tab in sorted(tabs.items()):
        c = compare(tab, cur)
        if c['same_fingerprint']:
            verd = '✅ 同链'
        elif c['is_constant_scale']:
            verd = '⚠️ 异链·常数缩放 ×%.4f' % (c['ratio_min'] or 0)
        else:
            verd = '⛔ 异链·**形状不同**（比值 %.4f→%.4f，跨度 %.1f%%）' % (
                c['ratio_min'] or 0, c['ratio_max'] or 0, c['ratio_span_pct'] or 0)
        print('%-52s %-18s %-14s %s' % (
            name[:52], c['fp_a'], ('%.4f' % (c['ratio_span_pct'] or 0)) + '%'
            if c['ratio_span_pct'] is not None else '—', verd))
        rows.append(dict(file=name, fp=c['fp_a'],
                         same=c['same_fingerprint'],
                         const_scale=c['is_constant_scale'],
                         ratio_min=c['ratio_min'], ratio_max=c['ratio_max'],
                         span_pct=c['ratio_span_pct'],
                         n_common=c['n_common']))

    n_same = sum(r['same'] for r in rows)
    print()
    print('合计: %d 份，同链 %d 份，异链 %d 份' % (len(rows), n_same,
                                                  len(rows) - n_same))
    if args.json:
        io.open(args.json, 'w', encoding='utf-8').write(json.dumps(
            {'current_fingerprint': fp_cur, 'heights': hs, 'files': rows},
            ensure_ascii=False, indent=2))
        print('[out]', args.json)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
