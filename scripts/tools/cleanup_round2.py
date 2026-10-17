"""第二轮清理（**只 move**）：① .bak/_backup ② results 无消费者产物 ③ rerun 旧副本
④ docs/rigid_vtarget ⑤ 10 个同前缀脚本群（引用守卫筛选）。

守卫：被 `archive/`、`temp/`、LISTERS（仅"列出名字"的审计/清理产物）之外的文件引用 ⇒ **拒移**。
"""
import fnmatch
import io
import os
import re
import shutil

BASE = r'D:\Project\climbing_fall_analysis'
DEST = os.path.join(BASE, 'archive', '2026-10-17_cleanup')
SKIP = {'.git', '.venv', '__pycache__', '.pytest_cache', '.ruff_cache', 'node_modules', 'archive'}
TEXT_EXT = ('.py', '.md', '.json', '.txt', '.csv', '.toml', '.cfg', '.yml', '.yaml', '.ini')
LISTERS = {'results/AUDIT_2026-10-17.md', 'results/SWEEP2_2026-10-17.md',
           'results/STALE_REFS.md', 'CANONICAL.md',
           'scripts/tools/audit_2026_10_17.py', 'scripts/tools/sweep2_2026_10_17.py',
           'scripts/tools/scan_stale_refs.py', 'scripts/tools/cleanup_parallel.py',
           'scripts/tools/archive_cleanup.py', 'scripts/tools/archive_doc_clusters.py',
           'scripts/tools/cleanup_round2.py'}
GROUPS = ('g4', 'make', 'nonvertical', 'verify', 'fe', 'a1', 'm1', 'p2', 'pad', 'thums')
SCAN = ('src', 'scripts', 'tests', 'docs', 'results', 'paper')
round2 = {}


def rel(p):
    return os.path.relpath(p, BASE).replace('\\', '/')


# ---- 收集可扫描文本（作为"引用方"） -----------------------------------
files = []
for sub in list(SCAN) + ['']:
    root = os.path.join(BASE, sub) if sub else BASE
    for dp, dns, fs in os.walk(root):
        dns[:] = [d for d in dns if d not in SKIP]
        if sub == '':
            dns[:] = []          # 根目录只看文件
        for f in fs:
            if f.endswith(TEXT_EXT):
                files.append(os.path.join(dp, f))
files = list(dict.fromkeys(files))
texts = {}
for p in files:
    try:
        texts[rel(p)] = io.open(p, encoding='utf-8', errors='ignore').read()
    except OSError:
        pass
print(f"引用方文本 {len(texts)} 个")

# ---- 候选 --------------------------------------------------------------
cands = {}

# ① .bak / _backup —— ⚠️ **不能走 `texts`**：这类文件名以**日期结尾**
#   （如 `.bak_precal_20261009`），会被 `TEXT_EXT` 扩展名过滤漏掉（首版只收到 3/12）。
for dp_, dns_, fs_ in os.walk(os.path.join(BASE, 'results')):
    dns_[:] = [d for d in dns_ if d not in SKIP]
    for f_ in fs_:
        if re.search(r'(\.bak|_backup)', f_, re.I):
            rr = rel(os.path.join(dp_, f_))
            if rr not in cands.get('bak_backup', []):
                cands.setdefault('bak_backup', []).append(rr)
for r in texts:
    b = os.path.basename(r)
    if re.search(r'(\.bak|_backup)', b, re.I):
        cands.setdefault('bak_backup', []).append(r)

# ③ rerun 旧副本（同名，上级已有）
rd = 'results/opensim_fe/rerun_20261008/'
for r in texts:
    if r.startswith(rd):
        up = 'results/opensim_fe/' + os.path.basename(r)
        if up in texts:
            cands.setdefault('rerun_dupes', []).append(r)

# ② results 无消费者产物
for r in texts:
    if not r.startswith('results/') or not r.endswith(('.json', '.txt', '.csv')):
        continue
    n = os.path.basename(r)
    users = [q for q in texts if q != r and n in texts[q]]
    if not users and r not in LISTERS:
        cands.setdefault('dead_results', []).append(r)

# ④ docs 早期一次性
p = 'docs/rigid_vtarget_2026-10-08.md'
if os.path.isfile(os.path.join(BASE, p.replace('/', os.sep))):
    cands.setdefault('docs', []).append(p)

# ⑤ 脚本群
for r in texts:
    if not r.startswith('scripts/') or not r.endswith('.py'):
        continue
    if r.split('/')[-1].split('_')[0] in GROUPS:
        cands.setdefault('script_groups', []).append(r)

for k, v in sorted(cands.items()):
    print(f"  候选 {k}: {len(v)}")

# ---- 守卫 + 移动 -------------------------------------------------------
moved, blocked = [], []
allc = {c for v in cands.values() for c in v}
for sub, items in sorted(cands.items()):
    for c in items:
        src = os.path.join(BASE, c.replace('/', os.sep))
        if not os.path.isfile(src):
            continue
        n = os.path.basename(c)
        refs = [q for q in texts
                if q != c and q not in allc and q not in LISTERS
                and not q.startswith(('archive/', 'temp/'))
                and n in texts[q]]
        if refs:
            blocked.append((c, refs[:3]))
            continue
        dst = os.path.join(DEST, f'round2_{sub}', n)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.move(src, dst)
        moved.append((c, rel(dst)))

print(f"\n▲ 移动 {len(moved)}｜拒移 {len(blocked)}")
for c, r in blocked[:15]:
    print(f"   拒 {c} ← {', '.join(r[:2])}")
lines = ['# 第二轮清理 MANIFEST（2026-10-17）', '', '> **只 move，不删。** 守卫排除清理/审计产物。', '',
         f'## 已移动（{len(moved)}）', '']
lines += [f'- `{x}` → `{y}`' for x, y in moved]
lines += ['', f'## 拒移（{len(blocked)}）', '']
lines += [f'- `{x}` ← {", ".join(r)}' for x, r in blocked]
io.open(os.path.join(DEST, 'MANIFEST_round2.md'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('MANIFEST_round2.md 已写')
