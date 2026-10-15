# -*- coding: utf-8 -*-
"""**V1 能否修**：`max`/`p95`（源侧 [Y25] 要的家族）在 `contact` BC 下是否还退化？

**背景**：V1 判 ❌【口径不可复现】的核心论据之一是「`max`/`p95` 是**伪影签名**
（其口径下**每个骨都 1 m**）」。而**探针已证**：该伪影的物理位置＝**跖面固定端的夹持边缘**
（热点恒为单元 #378、距跖面 2.99 mm），且**换成接触 BC 后热点会挪走**。

⇒ 本件直接查：**同一网格、同一链、只换跖面 BC** 时，`max`/`p95`/`mean`/`gauge_max`
各列的**首折高度**是否变化；若 `max` 不再是「1 m 全通过」，则 **V1 有可复现路径**。

数据（都已落盘、同链）：
* `s2_thums_height_sweep_anat_refixed_now.json`（`fixed`，同链）
* `s2_thums_height_sweep_anat_contact_hold100.json`（`contact + hold k=100`）
* `s2_thums_height_sweep_refixed_now.json` / `..._contact_hold100.json`（`calcnframe` 对照）

用法::

    .venv\\Scripts\\python.exe scripts\\opensim_fe\\v1_max_caliber_check.py
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
R = ROOT / "results" / "opensim_fe"
THR = 150.0          # 皮质骨强度阈值（跟骨，[Y25] Table 1）
COLS = ('max', 'p99', 'p95', 'mean', 'gauge_max')


def load(f: str) -> dict:
    b = json.load(io.open(R / f, encoding='utf-8'))
    return {float(x['h']): x for x in b if float(x.get('act', 0)) == 0.0}


def first_over(rows: dict, col: str, thr: float = THR):
    over = sorted(h for h, r in rows.items() if float(r.get(col, -1)) > thr)
    return over[0] if over else None


def main() -> int:
    cases = [
        ('anat', 's2_thums_height_sweep_anat_refixed_now.json',
         's2_thums_height_sweep_anat_contact_hold100.json'),
        ('calcnframe', 's2_thums_height_sweep_refixed_now.json',
         's2_thums_height_sweep_contact_hold100.json'),
    ]
    for lab, ff, fc in cases:
        a, c = load(ff), load(fc)
        print('=' * 100)
        print('%s  —— 首折（阈值 %.0f MPa）按口径列' % (lab, THR))
        print('%-12s %-16s %-16s %-10s %s' % (
            '口径列', 'fixed 首折', 'contact+hold 首折', '是否改变', '说明'))
        for col in COLS:
            fa, fc_ = first_over(a, col), first_over(c, col)
            chg = '—' if fa == fc_ else '★ 变'
            note = ''
            if fa == 1.0 and fc_ == 1.0:
                note = '**两者都退化到 h=1**'
            print('%-12s %-16s %-16s %-10s %s' % (
                col, ('%s m' % fa) if fa else '(无)', ('%s m' % fc_) if fc_ else '(无)',
                chg, note))
        print()
        print('--- 逐高度原值（%s）---' % lab)
        hs = sorted(set(a) & set(c))
        hdr = '%-6s' % 'h'
        for col in COLS:
            hdr += '| %-11s %-11s ' % ('fixed.' + col, 'cont.' + col)
        print(hdr)
        for h in hs:
            line = '%-6g' % h
            for col in COLS:
                line += '| %-11.2f %-11.2f ' % (float(a[h].get(col, -1)),
                                                 float(c[h].get(col, -1)))
            print(line)
        print()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
