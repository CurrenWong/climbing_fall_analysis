"""**多线程**整理入口 v2 —— 修正候选口径（目录级单元）与移动方式（整目录 rename）。

v1 的问题（已实测）：
  * 我把 `temp/pyfebio_demo/**` 的 **1,891 个文件**逐个当候选 ⇒ 候选 ~2,100
    ⇒ 比对 3191×2100 ≈ **670 万次** ⇒ 外推 **~26 s**；
  * 移动走 1,966 次逐文件 `shutil.move`。
  * 线程**无增益**（GIL；实测 16 线程读 vs 串行 = 0.052 vs 0.047 s）⇒ 不要再指望线程，
    **要靠减少工作量**。

v2 的三处改动：
  1. **目录级单元**：`temp/pyfebio_demo/` 是**一个单元**（不再拆成 1,891 个）⇒ 候选 ~150；
  2. **前置门剪枝**：目录单元先做一次 `目录名 in txt` 的廉价检查，命中才逐 basename 细查
     ⇒ 昂贵的比较只在极少数文件上发生；
  3. **整目录 `shutil.move`（同盘即一次 rename）**，不是 1,966 次。

用法：`python scripts/tools/cleanup_parallel.py [--jobs N] [--dry-run]`
"""
from __future__ import annotations

import argparse
import fnmatch
import io
import os
import shutil
import time
from concurrent.futures import ThreadPoolExecutor

BASE = r'D:\Project\climbing_fall_analysis'
STAMP = '2026-10-17_cleanup'
DEST = os.path.join(BASE, 'archive', STAMP)
SKIP_DIRS = {'.git', '.venv', '__pycache__', '.pytest_cache', '.ruff_cache',
             'node_modules', 'archive'}
TEXT_EXT = ('.py', '.md', '.json', '.txt', '.toml', '.cfg', '.yml', '.yaml', '.ini')
TODO_CUTOFF = time.mktime(time.strptime('2026-10-10', '%Y-%m-%d'))
DOCS_CUTOFF = time.mktime(time.strptime('2026-10-08', '%Y-%m-%d'))
DOCS_KEEP = ('*结论*', 'README.md', 'FE_BC_GOTCHAS.md', 'FE_PRESOLVE_GATES.md',
             'PAD_GOTCHAS.md', 'FEBIO_GOTCHAS.md', '阶段总结.md')
SELF = {'scripts/tools/cleanup_parallel.py', 'scripts/tools/archive_cleanup.py',
        'scripts/tools/archive_doc_clusters.py', 'scripts/tools/scan_stale_refs.py'}
#: ⛔ 绝不纳入：FE 产物目录（AGENTS.md：temp 子目录产物不动）
NEVER = {'temp/opensim_fe', 'temp/thums', 'temp/parts_gmsh'}
INVALID_FILES = {
    'results_invalid': [
        'results/UNLOAD_PROCESS_REPORT.md',
        'results/second_contact_probe.json',
        'results/b2_report.txt', 'results/smoke.txt',
        'results/calibration_ff_scan.csv', 'results/calibration_report.txt',
        'results/calibration_rope.csv', 'results/rope_report.txt',
    ],
    'scripts_invalid': ['scripts/verify_flex_observable.py'],
    'docs_stale': ['SECOND_CONTACT_TODO.md'],
}


def rel(p):
    return os.path.relpath(p, BASE).replace('\\', '/')


def build_units():
    """返回 unit 列表：{'disp','kind','path','sub','keys','gate'}"""
    units = []
    for sub, items in INVALID_FILES.items():
        for c in items:
            p = os.path.join(BASE, c.replace('/', os.sep))
            units.append({'disp': c, 'kind': 'file', 'path': p, 'sub': sub,
                          'keys': [os.path.basename(c)], 'gate': None})
    nd = os.path.join(BASE, 'docs', 'notes')
    if os.path.isdir(nd):
        for f in sorted(os.listdir(nd)):
            p = os.path.join(nd, f)
            if f.endswith('_TODO.md') and os.stat(p).st_mtime < TODO_CUTOFF:
                units.append({'disp': rel(p), 'kind': 'file', 'path': p,
                              'sub': 'docs_stale', 'keys': [f], 'gate': None})
    for dp, dns, fs in os.walk(os.path.join(BASE, 'docs')):
        dns[:] = [d for d in dns if d not in SKIP_DIRS and d != 'archive']
        for f in fs:
            if not f.endswith('.md') or any(fnmatch.fnmatch(f, k) for k in DOCS_KEEP):
                continue
            p = os.path.join(dp, f)
            if os.stat(p).st_mtime < DOCS_CUTOFF:
                units.append({'disp': rel(p), 'kind': 'file', 'path': p,
                              'sub': 'docs_old', 'keys': [f], 'gate': None})
    # ★ 目录级单元：temp/pyfebio_demo 整体
    pdd = os.path.join(BASE, 'temp', 'pyfebio_demo')
    if os.path.isdir(pdd) and 'temp/pyfebio_demo' not in NEVER:
        inner = [f for _, _, fs in os.walk(pdd) for f in fs]
        units.append({'disp': 'temp/pyfebio_demo', 'kind': 'dir', 'path': pdd,
                      'sub': 'temp_scratch',
                      'keys': sorted({os.path.basename(x) for x in inner}) + ['pyfebio_demo'],
                      'gate': 'pyfebio_demo'})
    td = os.path.join(BASE, 'temp')
    if os.path.isdir(td):
        for f in sorted(os.listdir(td)):
            p = os.path.join(td, f)
            if os.path.isfile(p) and f.endswith('.txt'):
                units.append({'disp': rel(p), 'kind': 'file', 'path': p,
                              'sub': 'temp_scratch', 'keys': [f], 'gate': None})
    return units


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=min(16, os.cpu_count() or 8))
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args(argv)
    T = time.perf_counter()

    units = build_units()
    unit_keys = [(u, u['keys'], u['gate']) for u in units]
    print(f"单元 {len(units)} 个（含目录单元 "
          f"{sum(1 for u in units if u['kind'] == 'dir')}）｜jobs={a.jobs}")

    files = []
    for dp, dns, fs in os.walk(BASE):
        dns[:] = [d for d in dns if d not in SKIP_DIRS]
        for f in fs:
            if f.endswith(TEXT_EXT):
                files.append(os.path.join(dp, f))
    print(f"文本文件 {len(files)} 个｜收集耗时 {time.perf_counter() - T:.2f}s")
    t1 = time.perf_counter()

    # ★ 成本模型（实测校正，2026-10-17）：O(候选数 × 总文本量)。
    #   v1（逐候选 `in`）：约 2.1e3 × 328 MB ≈ 690 GB ⇒ 估计 ~9 min；
    #   v2（单元化 121 候选）：121 × 328 MB ≈ 40 GB ⇒ **实测 32 s ✓**；
    #   v3（本版）：把**文件单元**的全部 basename 合成**一个编译正则** ⇒ 每份文本
    #   只扫 **一遍** ⇒ 328 MB ≈ 1–3 s。**目录单元**的 861 个内部 basename 不进正则
    #   （否则正则分支过多反而慢），改用**前置门**：只有文本里出现目录名时才细查（极少）。
    import re
    _file_names = sorted({k for u in units if u['kind'] == 'file' for k in u['keys']})
    _PAT = re.compile('|'.join(re.escape(n) for n in _file_names)) if _file_names else None

    def scan(path):
        r = rel(path)
        try:
            txt = io.open(path, encoding='utf-8', errors='ignore').read()
        except OSError:
            return r, []
        found = {m.group(0) for m in _PAT.finditer(txt)} if _PAT is not None else set()
        out = []
        for u, keys, gate in unit_keys:
            if gate is not None:                       # 目录单元：前置门
                if gate in txt and any(k in txt for k in keys):
                    out.append(u['disp'])
            elif found.intersection(keys):
                out.append(u['disp'])
        return r, out

    refs = {u['disp']: [] for u in units}
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for r, found in ex.map(scan, files, chunksize=64):
            if r.startswith('archive/') or r.startswith('temp/') or r in SELF:
                continue
            if any(not os.path.isfile(os.path.join(BASE, f.replace('/', os.sep)))
                   or rel(os.path.join(BASE, f.replace('/', os.sep))) == r
                   for f in found):
                continue
            for d in found:
                refs[d].append(r)
    print(f"扫描耗时 {time.perf_counter() - t1:.2f}s（共 {time.perf_counter() - T:.2f}s）")
    t2 = time.perf_counter()

    moved, blocked, missing = [], [], []
    for u in units:
        if not os.path.exists(u['path']):
            missing.append(u['disp'])
            continue
        who = refs.get(u['disp'], [])
        if who:
            blocked.append((u['disp'], who[:4]))
            continue
        dst = os.path.join(DEST, u['sub'], os.path.basename(u['path']))
        if not a.dry_run:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.move(u['path'], dst)     # 目录 ⇒ 同盘一次 rename
        moved.append((u['disp'], rel(dst)))
    print(f"移动耗时 {time.perf_counter() - t2:.2f}s（总 {time.perf_counter() - T:.2f}s）")

    lines = [f"# 并行整理清单 v2（{STAMP}）", "",
             f"> `scripts/tools/cleanup_parallel.py`（jobs={a.jobs}，单元级候选）；"
             f"**只 move，不删**。", "",
             f"## 已移动（{len(moved)}）", ""]
    lines += [f"- `{x}` → `{y}`" for x, y in moved]
    lines += ["", f"## 拒移（{len(blocked)}，仍被引用）", ""]
    lines += [f"- `{x}` ← {', '.join(w)}" for x, w in blocked]
    lines += ["", f"## 已不在原位（{len(missing)}）", ""]
    lines += [f"- `{x}`" for x in missing]
    if not a.dry_run:
        os.makedirs(DEST, exist_ok=True)
        io.open(os.path.join(DEST, 'MANIFEST_parallel.md'), 'w',
                encoding='utf-8').write('\n'.join(lines) + '\n')
    print(f"\n移动 {len(moved)}｜拒移 {len(blocked)}｜已不在原位 {len(missing)}"
          f"｜{'DRY-RUN' if a.dry_run else 'MANIFEST_parallel.md 已写'}")
    for x, y in moved[:12]:
        print(f"   移 {x} → {y}")
    for x, w in blocked[:8]:
        print(f"   拒 {x}  ← {', '.join(w[:2])}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
