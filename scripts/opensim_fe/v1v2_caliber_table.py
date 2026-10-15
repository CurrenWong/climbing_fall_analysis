# -*- coding: utf-8 -*-
"""V1/V2 口径审计：建「口径 × 部位 × 首折」全量表（**零新增解算**，只读冻结矩阵）。

用途：V1（跟骨首折）/ V2（<7 m 无骨折）的**口径与部位集合**未裁定 —— 本工具产出
全量口径表，供 `results/opensim_fe/V1V2_CALIBER_PREREG_DRAFT.md` 引用。

用法：.venv\\Scripts\\python.exe scripts\\opensim_fe\\v1v2_caliber_table.py

输出：results/opensim_fe/V1V2_CALIBER_TABLE.json

注意（闸门 9）：**报任何首折必须带口径 + 部位集合 + 模型/载荷链三标签**；
本项目实测 `max`/`p95` 口径给出「每个骨都 1 m」的伪影签名（网格太粗），
`median` 是为 A2 的阶跃对齐另立的（见 `A2_AGGREGATION_RULING.md` §2.1）。
"""
from __future__ import annotations

import io
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
fm = json.load(io.open(os.path.join(ROOT, 'results', 'opensim_fe',
                                    'fracture_matrix_bc2.json'), encoding='utf-8'))
H = np.asarray(fm['heights_m'], dtype=float)
parts = {p['bone']: float(p.get('sigma_c_mpa') or 150.0) for p in fm['parts']}
print('heights:', H.min(), '..', H.max(), 'n =', len(H))
print('stats 可用:', list(fm['sigma_vm_mpa'].keys()))
print('部位 (sigma_c):', parts)
print()

ORDER = ['calcaneus_r', 'fibula_r', 'tibia_r', 'femur_r', 'T6', 'L3', 'C5', 'parietal_r']
STATS = ['max', 'p95', 'mean', 'median']

print('=' * 118)
print('表 A · 各部位首折高度（σ_c 来自 [Y25] Table 1；► = <7 m 违反 V2）')
print('=' * 118)
print('%-13s %8s | %10s %10s %10s %10s' % ('bone', 'σ_c', *STATS))
rows = {}
for b in ORDER:
    if b not in parts:
        continue
    sc = parts[b]
    cells = []
    for st in STATS:
        sig = np.asarray(fm['sigma_vm_mpa'][st].get(b, []), dtype=float)
        i1 = next((i for i, v in enumerate(sig) if v > sc), None) if sig.size else None
        hv = H[i1] if i1 is not None else None
        cells.append(hv)
    rows[b] = cells
    print('%-13s %8.0f | %s' % (b, sc, ' '.join(
        ('%-9s' % ('%.0f m' % c)) if c is not None else ('%-9s' % '>50 m') for c in cells)))

print()
print('=' * 118)
print('表 B · V2（「<7 m 无骨折」）在各口径下的违反部位')
print('=' * 118)
for j, st in enumerate(STATS):
    bad = [b for b in ORDER if b in rows and rows[b][j] is not None and rows[b][j] < 7.0]
    ok = [b for b in ORDER if b in rows and (rows[b][j] is None or rows[b][j] >= 7.0)]
    print('  %-8s 违反 %d/%d：%-42s ｜ 满足：%s'
          % (st, len(bad), len([b for b in ORDER if b in rows]),
             '、'.join(bad) or '无', '、'.join(ok)))

print()
print('=' * 118)
print('表 C · V1（跟骨首折，目标 7–9 m）在各口径下的值')
print('=' * 118)
for st in STATS:
    c = rows.get('calcaneus_r', [None] * 4)[STATS.index(st)]
    print('  %-8s %s' % (st, ('%.0f m' % c) if c is not None else '>50 m'))
print('  %-8s %s' % ('1D 截面', '>40 m（validate_vs_fem；自述高估 17–34×）'))
print('  %-8s %s' % ('S2 gauge', '3 m 旧链 → 2 m 当前链（单骨 AM50 跟骨 + ℓ=4mm）'))
print('  %-8s %s' % ('1D 子区域', '足 8 m / 腓 7 m（纯轴向 + 子区域强度）'))
print('  %-8s %s' % ('旧整足 FE', '26.5 m（平滑网格 + ℓ=4mm）'))

# 落盘
out = {'heights_m': H.tolist(),
       'sigma_c': parts,
       'first_fracture_m': {b: rows[b] for b in rows},
       'stats': STATS}
dst = os.path.join(ROOT, 'results', 'opensim_fe', 'V1V2_CALIBER_TABLE.json')
io.open(dst, 'w', encoding='utf-8').write(
    json.dumps(out, ensure_ascii=False, indent=1))
print('\n产物:', os.path.relpath(dst, ROOT))
