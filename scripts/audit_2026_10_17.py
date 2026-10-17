"""审计（只读，多线程）：① 剩下的代码是否都在流程上 ② 文档是否仍按已作废口径陈述。

产物：`results/AUDIT_2026-10-17.md`
"""
import io
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor

BASE = r'D:\Project\climbing_fall_analysis'
SKIP = {'.git', '.venv', '__pycache__', '.pytest_cache', '.ruff_cache', 'node_modules', 'archive'}
TEXT_EXT = ('.py', '.md', '.json', '.txt', '.toml', '.cfg', '.yml', '.yaml', '.ini')

#: 已作废口径（CANONICAL.md §2）—— 出现这些应带作废标注
VOID = {
    'h*=1.34-3.66': r'1\.34\s*[–~-]\s*3\.66|3\.66\s*m',
    'h*=>9 / 7.17 / 2.09': r'7\.17\s*m|2\.09\s*m|>\s*9\s*m',
    'h*=<0.5 m': r'h\*\s*<\s*0\.5|<\s*0\.5\s*m',
    '旧“与流行病学 1.1–4 一致”': r'与流行病学|1\.1\s*[–~-]\s*4\s*m',
    '“k_flex 是死参数”': r'死参数',
    'O1 荒谬值 673.42': r'673\.42',
    '“废表”叙事（已推翻）': r'废表',
    'M(θ) 表格式误述 (x_rad,y)': r'\(TAB_TH,\s*TAB_M\)|x_rad_array',
    '力台锚定在门路径仍有效': r'力台锚定三档|k\s*=\s*6\.3\s*/\s*8\.2\s*/\s*18\.5',
}


def rel(p):
    return os.path.relpath(p, BASE).replace('\\', '/')


def collect():
    out = []
    for dp, dns, fs in os.walk(BASE):
        dns[:] = [d for d in dns if d not in SKIP]
        for f in fs:
            if f.endswith(TEXT_EXT):
                out.append(os.path.join(dp, f))
    return out


def main():
    T = time.perf_counter()
    files = collect()

    def read(path):
        try:
            return path, io.open(path, encoding='utf-8', errors='ignore').read()
        except OSError:
            return path, ''

    with ThreadPoolExecutor(max_workers=16) as ex:
        texts = dict(ex.map(read, files, chunksize=64))
    print(f"读入 {len(texts)} 个文件，{time.perf_counter() - T:.1f}s")

    # ---------- ① 代码：入向引用 ----------
    code = {p: t for p, t in texts.items()
            if p.endswith('.py') and not rel(p).startswith(('temp/', 'archive/'))}
    src = {rel(p): t for p, t in code.items() if rel(p).startswith('src/')}
    scr = {rel(p): t for p, t in code.items() if rel(p).startswith('scripts/')}
    tst = {rel(p): t for p, t in code.items() if rel(p).startswith('tests/')}
    others = [rel(p) for p in texts
              if not rel(p).startswith(('temp/', 'archive/')) and not p.endswith('.py')]
    live_src = list(src) + list(scr) + list(tst) + others

    # 归档/作废目标（引用它们是坏的）
    archived = set()
    for dp, _, fs in os.walk(os.path.join(BASE, 'archive')):
        for f in fs:
            archived.add(f)

    rows_src, rows_scr = [], []
    for r, t in src.items():
        mod = os.path.basename(r)[:-3]
        inbound = [q for q in live_src if q != r and (f'import {mod}' in texts[os.path.join(BASE, q.replace("/", os.sep))]
                                                      or mod in texts[os.path.join(BASE, q.replace("/", os.sep))])]
        rows_src.append((r, len(inbound), inbound[:3]))
    for r, t in scr.items():
        base = os.path.basename(r)
        inbound = [q for q in live_src if q != r and base in texts[os.path.join(BASE, q.replace("/", os.sep))]]
        dead = sorted({a for a in archived if a in t and a.endswith(('.json', '.md'))})
        rows_scr.append((r, len(inbound), inbound[:2], dead[:3]))

    # ---------- ② 文档：过时口径命中 ----------
    doc_hits = {}
    for p, t in texts.items():
        r = rel(p)
        if r.startswith(('archive/', 'temp/')):
            continue
        if not (r.endswith(('.md', '.txt')) or r.startswith('results/')):
            continue
        for name, pat in VOID.items():
            for i, line in enumerate(t.splitlines()):
                if re.search(pat, line):
                    marked = any(k in line for k in ('作废', '禁引', '已推翻', '不得引', 'void', '§0b', '§0c', '§0d'))
                    doc_hits.setdefault(name, []).append((r, i + 1, line.strip()[:130], marked))

    out = ["# 审计 2026-10-17 ① 代码是否都在流程上 ② 文档是否合新口径", "",
           f"- 读入 {len(texts)} 个文件；`src/` {len(src)}｜`scripts/` {len(scr)}｜`tests/` {len(tst)}", ""]

    print("\n" + "=" * 92)
    print("① src/climbing 模块入向引用")
    print("=" * 92)
    out += ["## ① src/ 模块入向引用", "", "| 模块 | 入向引用数 | 示例 |", "|---|---:|---|"]
    for r, n, ex in sorted(rows_src, key=lambda x: x[1]):
        flag = '  ← ⚠️ 零引用' if n == 0 else ''
        print(f"  {r:<48} {n:>4}{flag}")
        out.append(f"| `{r}` | {n} | {', '.join(ex)} |")

    print("\n" + "=" * 92)
    print("② scripts/：被引用数 + 是否引用已归档文件")
    print("=" * 92)
    out += ["", "## ② scripts/", "", "| 脚本 | 入向引用 | 引用已归档文件 |", "|---|---:|---|"]
    for r, n, ex, dead in sorted(rows_scr, key=lambda x: (len(x[3]) == 0, x[1])):
        flag = '  ← ⚠️ 零引用' if n == 0 else ''
        dl = ('⚠️ ' + ', '.join(dead)) if dead else ''
        print(f"  {r:<44} {n:>4}{flag}  {dl}")
        out.append(f"| `{r}` | {n} | {dl} |")

    print("\n" + "=" * 92)
    print("③ 文档过时口径命中（★=未带作废标注）")
    print("=" * 92)
    out += ["", "## ③ 文档过时口径", ""]
    for name, hits in sorted(doc_hits.items(), key=lambda kv: -len(kv[1])):
        unmarked = [h for h in hits if not h[3]]
        print(f"\n▼ {name}：{len(hits)} 条（未标注 {len(unmarked)}）")
        out += [f"### {name}（{len(hits)} 条，未标注 {len(unmarked)}）", ""]
        for r, ln, tx, mk in hits[:12]:
            print(f"    {'  ' if mk else '★ '}{r}:{ln}  | {tx[:100]}")
            out.append(f"- {'`[已标注]`' if mk else '**★ 未标注**'} `{r}:{ln}` — {tx}")
        if len(hits) > 12:
            out.append(f"- …另 {len(hits) - 12} 条")

    dst = os.path.join(BASE, 'results', 'AUDIT_2026-10-17.md')
    io.open(dst, 'w', encoding='utf-8').write('\n'.join(out) + '\n')
    print(f"\nOUT = results/AUDIT_2026-10-17.md  （{time.perf_counter() - T:.1f}s）")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
