"""归档**簇模式**：若某目录下的全部 md 都是「过期」且**只被该目录内部引用**，
则整个目录可作为一个**簇**一起归档（无需逐个人工定夺）。

规则（见 `CANONICAL.md` §3）：
  1. 目录内 md 全部 mtime < CUTOFF；
  2. 目录内 md 的引文**全部来自该目录内部**（外部零引用）；
  3. 白名单文件名（`*结论*` / `README.md` / 三本 gotcha / `阶段总结.md`）不入簇。
满足即整目录 move 到 `archive/<STAMP>/docs_clusters/`。**只 move，不删。**
"""
import io
import os
import shutil
import time

BASE = r'D:\Project\climbing_fall_analysis'
STAMP = '2026-10-17_cleanup'
DEST = os.path.join(BASE, 'archive', STAMP, 'docs_clusters')
SKIP_DIRS = {'.git', '.venv', '__pycache__', '.pytest_cache', '.ruff_cache',
             'node_modules', 'archive'}
TEXT_EXT = ('.py', '.md', '.json', '.txt', '.toml', '.cfg', '.yml', '.yaml', '.ini')
CUTOFF = time.mktime(time.strptime('2026-10-08', '%Y-%m-%d'))
KEEP = ('结论',)          # 文件名含「结论」⇒ 永不入簇


def rel(p):
    return os.path.relpath(p, BASE).replace('\\', '/')


def all_text_files():
    out = []
    for dp, dns, fs in os.walk(BASE):
        dns[:] = [d for d in dns if d not in SKIP_DIRS]
        for f in fs:
            if f.endswith(TEXT_EXT):
                out.append(os.path.join(dp, f))
    return out


def main():
    files = all_text_files()
    cache = {}

    def read(p):
        if p not in cache:
            try:
                cache[p] = io.open(p, encoding='utf-8', errors='ignore').read()
            except OSError:
                cache[p] = ''
        return cache[p]

    # 按「最近一层 docs/ 子目录」分组
    groups = {}
    dd = os.path.join(BASE, 'docs')
    for dp, dns, fs in os.walk(dd):
        dns[:] = [d for d in dns if d not in SKIP_DIRS and d != 'archive']
        mds = [os.path.join(dp, f) for f in fs if f.endswith('.md')]
        if mds:
            groups[rel(dp)] = mds

    moved, blocked, skipped = [], [], []
    for g, paths in sorted(groups.items()):
        if g == 'docs':
            continue
        names = [os.path.basename(p) for p in paths]
        if len(names) < 2:
            skipped.append((g, '文件数 < 2'))
            continue
        if any(any(k in n for k in KEEP) or n == 'README.md' for n in names):
            skipped.append((g, '含白名单（结论 / README）'))
            continue
        old = [n for n, p in zip(names, paths) if os.stat(p).st_mtime < CUTOFF]
        if len(old) != len(names):
            skipped.append((g, f'仅 {len(old)}/{len(names)} 过期'))
            continue
        pset = set(paths)
        ext = []
        for p in paths:
            n = os.path.basename(p)
            for q in files:
                if q in pset:
                    continue
                qr = rel(q)
                if qr.startswith('archive/') or qr.startswith('temp/'):
                    continue
                if qr == 'scripts/archive_cleanup.py' or qr == 'scripts/archive_doc_clusters.py':
                    continue
                if n in read(q):
                    ext.append((n, qr))
        if ext:
            blocked.append((g, len(names), ext[:5]))
            continue
        dst = os.path.join(DEST, g.replace('docs/', '', 1).replace('/', '__'))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.move(os.path.join(BASE, g.replace('/', os.sep)), dst)
        moved.append((g, len(names), rel(dst)))

    print("=" * 92)
    print("可整簇归档（外部零引用）")
    print("=" * 92)
    for g, n, d in moved:
        print(f"  ✅ {g}  （{n} 个 md）→ {d}")
    if not moved:
        print("  （无）")
    print()
    print("=" * 92)
    print("被外部引用挡住（需人工：先把引用改指向新件）")
    print("=" * 92)
    for g, n, ext in blocked:
        print(f"  ⛔ {g}（{n} 个）← 外部引用 {len(ext)}+："
              + "；".join(f"{a} ← {b}" for a, b in ext[:3]))
    print()
    print("=" * 92)
    print("不入簇（保留）")
    print("=" * 92)
    for g, why in skipped:
        print(f"  ➖ {g}：{why}")
    print(f"\n簇归档 {len(moved)} 个目录 / {sum(n for _, n, _ in moved)} 个文件"
          f"；受阻 {len(blocked)}；保留 {len(skipped)}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
