"""归档整理（**只 move，不删除**）+ 引用守卫 + MANIFEST —— 约定见 `CANONICAL.md` §3。

守卫：若某候选文件的名字在 `archive/` 与 `temp/` 之外**仍被引用**（且引用方不是本批
候选），则**拒移**并打印引用方 —— 防止把别人读的输入挪走（如 `fp_anchor_probe.json`
被 7 个脚本当 J0 正对照的期望值读）。
"""
import io
import os
import shutil
import time

BASE = r'D:\Project\climbing_fall_analysis'
STAMP = '2026-10-17_cleanup'
DEST = os.path.join(BASE, 'archive', STAMP)
SKIP_DIRS = {'.git', '.venv', '__pycache__', '.pytest_cache', '.ruff_cache',
             'node_modules', 'archive'}
TEXT_EXT = ('.py', '.md', '.json', '.txt', '.toml', '.cfg', '.yml', '.yaml', '.ini')

#: 明确作废 / 过期的产物（**不删，只移**）
CAND = {
    'results_invalid': [
        'results/UNLOAD_PROCESS_REPORT.md',      # 数字经 §0b/§0c/§0d 全部作废
        'results/second_contact_probe.json',     # 旧后处理口径产出
        'results/b2_report.txt', 'results/smoke.txt',
        'results/calibration_ff_scan.csv', 'results/calibration_report.txt',
        'results/calibration_rope.csv', 'results/rope_report.txt',
    ],
    'scripts_invalid': [
        'scripts/verify_flex_observable.py',     # ★ 表格式违反契约（事故件）
    ],
    'docs_stale': [
        'SECOND_CONTACT_TODO.md',                # 根目录会话草稿（应 gitignore）
    ],
}

#: 自动候选：`docs/notes/*_TODO.md` 中早于该日期的旧会话草稿
TODO_CUTOFF = time.mktime(time.strptime('2026-10-10', '%Y-%m-%d'))
#: 自动候选：`docs/` 下早于该日期的阶段性产物（阶段报告 / 一次性说明）
DOCS_CUTOFF = time.mktime(time.strptime('2026-10-08', '%Y-%m-%d'))
#: ⛔ 白名单：**永不**纳入自动归档（项目规则保护 / 活跃踩坑手册 / 入口）
DOCS_KEEP = ('*结论*', 'README.md', 'FE_BC_GOTCHAS.md', 'FE_PRESOLVE_GATES.md',
             'PAD_GOTCHAS.md', 'FEBIO_GOTCHAS.md', '阶段总结.md')
#: 自动候选：`temp/` **顶层** .py 草稿区（子目录产物不动）
TEMP_TOP = True


def rel(p):
    return os.path.relpath(p, BASE).replace('\\', '/')


def walk_all():
    out = []
    for dp, dns, fs in os.walk(BASE):
        dns[:] = [d for d in dns if d not in SKIP_DIRS]
        for f in fs:
            if f.endswith(TEXT_EXT):
                out.append(os.path.join(dp, f))
    return out


def build_candidates():
    cand = {k: list(v) for k, v in CAND.items()}
    # docs/notes 旧 TODO
    nd = os.path.join(BASE, 'docs', 'notes')
    if os.path.isdir(nd):
        for f in sorted(os.listdir(nd)):
            if f.endswith('_TODO.md'):
                p = os.path.join(nd, f)
                if os.stat(p).st_mtime < TODO_CUTOFF:
                    cand.setdefault('docs_stale', []).append(rel(p))
    # docs/ 下早于 DOCS_CUTOFF 的阶段性产物（白名单 + 引用守卫双重保护）
    import fnmatch
    for dp, dns, fs in os.walk(os.path.join(BASE, 'docs')):
        dns[:] = [d for d in dns if d not in SKIP_DIRS and d != 'archive']
        for f in fs:
            if not f.endswith('.md'):
                continue
            relp = rel(os.path.join(dp, f))
            if any(fnmatch.fnmatch(f, pat) for pat in DOCS_KEEP):
                continue
            if os.stat(os.path.join(dp, f)).st_mtime < DOCS_CUTOFF:
                cand.setdefault('docs_old', []).append(relp)
    # temp/ 顶层草稿
    td = os.path.join(BASE, 'temp')
    if TEMP_TOP and os.path.isdir(td):
        for f in sorted(os.listdir(td)):
            if f.endswith('.py') and os.path.isfile(os.path.join(td, f)):
                cand.setdefault('temp_drafts', []).append(rel(os.path.join(td, f)))
    return cand


def main():
    cand = build_candidates()
    all_cand = {c for v in cand.values() for c in v}
    files = walk_all()
    print(f"候选 {len(all_cand)} 个；扫描 {len(files)} 个文本文件做引用守卫…\n")
    moved, skipped = [], []
    for sub, items in cand.items():
        for c in items:
            src = os.path.join(BASE, c.replace('/', os.sep))
            if not os.path.isfile(src):
                print(f"  [skip ] 不存在: {c}")
                continue
            name = os.path.basename(c)
            refs = []
            for p in files:
                pr = rel(p)
                if pr == c or pr.startswith('temp/') or pr.startswith('archive/'):
                    continue
                if pr == 'scripts/tools/archive_cleanup.py':
                    continue                     # 本脚本自身含有候选名单 ⇒ 不是真引用
                if pr in all_cand:          # 同批候选（会一起移动）
                    continue
                try:
                    if name in io.open(p, encoding='utf-8', errors='ignore').read():
                        refs.append(pr)
                except OSError:
                    continue
            if refs:
                skipped.append((c, refs))
                print(f"  [拒移 ] {c}  ← 仍被引用: {', '.join(refs[:4])}")
                continue
            dst = os.path.join(DEST, sub, name)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.move(src, dst)
            moved.append((c, rel(dst)))
            print(f"  [移动 ] {c}  →  {rel(dst)}")
    # MANIFEST
    os.makedirs(DEST, exist_ok=True)
    lines = [f"# 归档清单 {STAMP}", "",
             "> 仅 **move**（可恢复）。事由：清理 2026-10-16~17 返工事故后的过期产物与草稿；",
             "> 约定见 `CANONICAL.md` §3。", "", "## 已移动", ""]
    lines += [f"- `{a}` → `{b}`" for a, b in moved]
    lines += ["", "## 拒移（仍被引用，需人工处置）", ""]
    lines += [f"- `{a}` ← {', '.join(r)}" for a, r in skipped]
    lines += ["", "## 未列入本批（仍为现行，勿动）", "",
              "- `results/fp_anchor_probe.json`（7 个脚本的输入锚）",
              "- `results/*_REPORT.md` 中被 `pad.py` / `tests` 引用的 3 个",
              "- `docs/reports/最终结论_2026-10-07.md`（项目规则：旧版结论件保留不改）"]
    io.open(os.path.join(DEST, 'MANIFEST.md'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
    print(f"\n移动 {len(moved)} 个，拒移 {len(skipped)} 个")
    print(f"MANIFEST = {rel(os.path.join(DEST, 'MANIFEST.md'))}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
