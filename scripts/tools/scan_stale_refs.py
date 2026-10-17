"""扫描**过期/悬空引用**（只读）。

口径：
  * 目标 = `archive/**` 下所有文件的 basename（已被归档 ⇒ 指向它的引用已**悬空**）
    ＋ 显式补充清单（如已作废结论件、已封档路线文档）。
  * 源 = 全仓文本文件（**排除** `archive/` 自身与 `.git` 等）。
产物：`results/STALE_REFS.md`（按目标分组，逐条给出 `源:行: 原文`）。
"""
import io
import os
import time

BASE = r'D:\Project\climbing_fall_analysis'
SKIP_DIRS = {'.git', '.venv', '__pycache__', '.pytest_cache', '.ruff_cache',
             'node_modules', 'archive'}
TEXT_EXT = ('.py', '.md', '.json', '.txt', '.toml', '.cfg', '.yml', '.yaml', '.ini')
#: 显式补充的「过期目标」（即便文件还在原位，其内容/口径已作废）
EXTRA_STALE = {
    'UNLOAD_PROCESS_REPORT.md': '数字全部作废（见 CANONICAL.md §2）',
    'second_contact_probe.json': '旧后处理口径产出',
    'verify_flex_observable.py': '表格式违反契约（事故件）',
}


def rel(p):
    return os.path.relpath(p, BASE).replace('\\', '/')


def main():
    archived = {}
    ar = os.path.join(BASE, 'archive')
    for dp, dns, fs in os.walk(ar):
        dns[:] = [d for d in dns if d not in SKIP_DIRS]
        for f in fs:
            archived.setdefault(f, []).append(rel(os.path.join(dp, f)))
    for k, v in EXTRA_STALE.items():
        archived.setdefault(k, []).append('（原位，内容' + v + '）')

    sources = []
    for dp, dns, fs in os.walk(BASE):
        dns[:] = [d for d in dns if d not in SKIP_DIRS]
        for f in fs:
            if f.endswith(TEXT_EXT):
                sources.append(os.path.join(dp, f))

    hits = {}          # target -> [(src, lineno, text)]
    for p in sources:
        sr = rel(p)
        try:
            lines = io.open(p, encoding='utf-8', errors='ignore').read().splitlines()
        except OSError:
            continue
        for i, line in enumerate(lines):
            for tgt in archived:
                if tgt in line:
                    hits.setdefault(tgt, []).append((sr, i + 1, line.strip()[:150]))

    out = [f"# 过期/悬空引用扫描（{time.strftime('%Y-%m-%d %H:%M')}）", "",
           f"- 归档文件 {len(archived)} 个（含显式补充 {len(EXTRA_STALE)}）",
           f"- 命中目标 {len(hits)} 个，引用条目 {sum(len(v) for v in hits.values())} 条", ""]
    print("=" * 96)
    print("过期/悬空引用（按目标分组）")
    print("=" * 96)
    for tgt in sorted(hits, key=lambda t: -len(hits[t])):
        rows = hits[tgt]
        where = archived[tgt][0]
        print(f"\n▼ {tgt}   （{len(rows)} 条；已归档 → {where}）")
        out += [f"## `{tgt}`  （{len(rows)} 条）", f"归档位置：`{where}`", ""]
        for sr, ln, tx in rows[:40]:
            print(f"    {sr}:{ln}  |  {tx[:110]}")
            out.append(f"- `{sr}:{ln}` — {tx}")
        if len(rows) > 40:
            out.append(f"- …另 {len(rows) - 40} 条")
        out.append("")
    dst = os.path.join(BASE, 'results', 'STALE_REFS.md')
    io.open(dst, 'w', encoding='utf-8').write('\n'.join(out) + '\n')
    print(f"\n命中目标 {len(hits)}｜条目 {sum(len(v) for v in hits.values())}")
    print(f"OUT = results/STALE_REFS.md")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
