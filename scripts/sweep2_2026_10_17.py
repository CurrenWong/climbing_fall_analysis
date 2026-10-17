"""第二次扫描（只读）：从 **文件名 / 日期 / 功能** 三个信号找可清理的代码与过期文档。

范围：`src/` `scripts/` `tests/` `docs/` `results/` + 根目录文件。
排除：`archive/` `temp/`（草稿区）`.venv/.git` 等。
产物：`results/SWEEP2_2026-10-17.md`
"""
import collections
import io
import os
import re
import time

BASE = r'D:\Project\climbing_fall_analysis'
SKIP = {'.git', '.venv', '__pycache__', '.pytest_cache', '.ruff_cache',
        'node_modules', 'archive', 'temp', 'parts_gmsh', 'parts'}
SCAN = ('src', 'scripts', 'tests', 'docs', 'results')
TEXT_EXT = ('.py', '.md', '.json', '.txt')
TODAY = time.mktime(time.strptime('2026-10-17', '%Y-%m-%d'))

NAME_BAD = {
    '旧版本后缀': r'(_old|_OLD|_bak|\.bak|~$|_copy|_backup|_v\d+$|_old2|_backup\d*)',
    '临时/草稿名': r'(^tmp_|^tmp|_tmp$|_temp$|^scratch|^test_|^_test|^demo|^untitled|未命名)',
    '序号尾巴': r'(\d\.py$|-\d\.py$|_\d\.py$)',
    '名字含日期': r'20\d\d-?\d?\d-?\d?\d',
}


def rel(p):
    return os.path.relpath(p, BASE).replace('\\', '/')


def main():
    files = []
    for sub in SCAN:
        root = os.path.join(BASE, sub)
        for dp, dns, fs in os.walk(root):
            dns[:] = [d for d in dns if d not in SKIP]
            for f in fs:
                files.append(os.path.join(dp, f))
    for f in os.listdir(BASE):
        p = os.path.join(BASE, f)
        if os.path.isfile(p):
            files.append(p)

    out = [f"# 第二次扫描（文件名 / 日期 / 功能）{time.strftime('%Y-%m-%d %H:%M')}", "",
           f"- 范围 {len(files)} 个文件（`archive/`、`temp/` 除外）", ""]

    # ---------- A. 文件名 ----------
    print("=" * 92)
    print("A. 文件名可疑")
    print("=" * 92)
    out += ["## A. 文件名可疑", ""]
    for label, pat in NAME_BAD.items():
        hits = sorted(rel(p) for p in files
                      if re.search(pat, os.path.basename(p), re.I))
        if not hits:
            continue
        print(f"\n▼ {label}（{len(hits)}）")
        out += [f"### {label}（{len(hits)}）", ""]
        for h in hits[:25]:
            print(f"    {h}")
            out.append(f"- `{h}`")
        if len(hits) > 25:
            out.append(f"- …另 {len(hits) - 25} 个")

    # ---------- B. 日期 ----------
    print("\n" + "=" * 92)
    print("B. 日期分层（mtime）")
    print("=" * 92)
    buckets = collections.Counter()
    stale = collections.defaultdict(list)
    for p in files:
        try:
            age = (TODAY - os.path.getmtime(p)) / 86400.0
        except OSError:
            continue
        b = ('≤3 天（现行）' if age <= 3 else '4–7 天' if age <= 7
             else '8–12 天' if age <= 12 else '>12 天（早期）')
        buckets[b] += 1
        if age > 12:
            stale[os.path.dirname(rel(p)) or '.'].append(rel(p))
    out += ["", "## B. 日期分层", "", "| 桶 | 文件数 |", "|---|---:|"]
    for k in ('≤3 天（现行）', '4–7 天', '8–12 天', '>12 天（早期）'):
        print(f"  {k:<16} {buckets[k]:>4}")
        out.append(f"| {k} | {buckets[k]} |")
    print("\n  >12 天（早期）按目录聚集 top8：")
    out += ["", "### 早期文件（>12 天）按目录", ""]
    for d, lst in sorted(stale.items(), key=lambda kv: -len(kv[1]))[:8]:
        print(f"    {len(lst):>3} 个  {d}")
        out.append(f"- **{len(lst)} 个** `{d}`")

    # ---------- C. 功能 ----------
    print("\n" + "=" * 92)
    print("C. 功能重叠 / 无消费者的输出")
    print("=" * 92)
    out += ["", "## C. 功能", ""]
    # C1 脚本前缀聚集
    sc = [rel(p) for p in files if p.endswith('.py') and rel(p).startswith('scripts/')]
    prefix = collections.Counter(os.path.basename(s).split('_')[0] for s in sc)
    print("\n▼ 同一前缀的脚本数（≥5 ⇒ 试探性脚本群）")
    out += ["### 同前缀脚本群（≥5）", "", "| 前缀 | 个数 |", "|---|---:|"]
    for k, v in prefix.most_common():
        if v >= 5:
            print(f"    {k:<14} {v}")
            out.append(f"| `{k}` | {v} |")
    # C2 同名多份
    base = collections.defaultdict(list)
    for p in files:
        base[os.path.basename(p)].append(rel(p))
    dup = {k: v for k, v in base.items() if len(v) > 1}
    print(f"\n▼ 同名多份：{len(dup)} 组")
    out += ["", f"### 同名多份（{len(dup)} 组）", ""]
    for k, v in sorted(dup.items())[:20]:
        print(f"    {k}: {', '.join(v)}")
        out.append(f"- `{k}` → {', '.join(v)}")
    # C3 无消费者的 results 输出
    texts = {}
    for p in files:
        if p.endswith(TEXT_EXT):
            try:
                texts[rel(p)] = io.open(p, encoding='utf-8', errors='ignore').read()
            except OSError:
                pass
    dead = []
    for r, t in texts.items():
        if not r.startswith('results/') or not r.endswith(('.json', '.txt', '.csv')):
            continue
        n = os.path.basename(r)
        users = [q for q in texts if q != r and 'archive/' not in q and n in texts[q]]
        if not users:
            dead.append(r)
    print(f"\n▼ results/ 下**无任何文档/脚本引用**的产物：{len(dead)}")
    out += ["", f"### results/ 无消费者的产物（{len(dead)}）", ""]
    for d in sorted(dead):
        print(f"    {d}")
        out.append(f"- `{d}`")

    io.open(os.path.join(BASE, 'results', 'SWEEP2_2026-10-17.md'), 'w',
            encoding='utf-8').write('\n'.join(out) + '\n')
    print("\nOUT = results/SWEEP2_2026-10-17.md")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
