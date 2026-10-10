"""行号内容核对：出处索引 表A 的 `path:start-end` 与其右侧代码片段是否真的对得上。

只查表 A（方程 → 代码），因为该表同时给出**行号**与**代码文本**，可互相验证。
判据：把片段里的第一个特征语句拿出来，在 [start, end]±8 行窗口内搜到 ⇒ OK。
"""
from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
IDX = ROOT / "docs/paper/论文v2_出处索引.md"

WINDOW = 12   # 允许上下漂移行数


def norm(s: str) -> str:
    """归一化：去 LaTeX、反引号、空白，便于字面匹配。"""
    s = re.sub(r"\$[^$]*\$", " ", s)
    s = s.replace("`", "").replace("\\", "")
    s = re.sub(r"\s+", "", s)
    return s


def main() -> int:
    txt = IDX.read_text(encoding="utf-8")
    m = re.search(r"## 2\. 表 A.*?(?=\n---|\n## 3\.)", txt, re.S)
    if not m:
        print("FAIL: cannot locate 表 A")
        return 1
    section = m.group(0)

    rows = re.findall(r"^\|(.+)\|\s*$", section, re.M)
    print("=" * 76)
    print(f"表 A 行号 <-> 代码片段 核对  ({len(rows)-2} data rows expected)")
    print("=" * 76)

    ok = bad = skipped = 0
    problems: list[str] = []
    for row in rows:
        cells = [c.strip() for c in row.split("|")]
        if len(cells) < 3:
            continue
        if cells[0].startswith("---") or cells[0] in ("描述", ""):
            continue
        desc, loc = cells[0], cells[1]
        code_cell = cells[2] if len(cells) > 2 else ""
        lm = re.search(r"`([^`]+?):(\d+)(?:-(\d+))?`", loc)
        if not lm:
            skipped += 1
            continue
        path = lm.group(1)
        a = int(lm.group(2))
        b = int(lm.group(3) or lm.group(2))
        f = ROOT / path
        if not f.exists():
            # 短路径惯例
            cand = None
            for d in ["", "src", "scripts", "docs", "docs/reports", "docs/notes",
                      "src/climbing", "scripts/opensim_fe"]:
                c = (ROOT / d / pathlib.PurePosixPath(path).name) if d else (
                    ROOT / pathlib.PurePosixPath(path).name)
                if c.exists():
                    cand = c
                    break
            if cand is None:
                bad += 1
                problems.append(f"FILE MISSING  {path}  ({desc[:30]})")
                continue
            f = cand
        lines = f.read_text(encoding="utf-8",
                            errors="replace").splitlines()

        # 取代码片段中第一个"像语句"的 token
        stmts = [s for s in re.split(r"[;`]", code_cell) if s.strip()]
        target = None
        for s in stmts:
            s2 = norm(s)
            if len(s2) >= 8:
                target = s2
                break
        if not target:
            skipped += 1
            continue

        lo, hi = max(1, a - WINDOW), min(len(lines), b + WINDOW)
        window = norm("\n".join(lines[lo - 1:hi]))
        if target in window:
            ok += 1
        else:
            bad += 1
            problems.append(
                f"NOT AT LINE  {path}:{a}-{b}  ({desc[:34]})  "
                f"missing -> {target[:52]!r}")
            # 全文件找一下真正位置
            whole = norm("\n".join(lines))
            if target in whole:
                idx = whole.find(target)
                # 反推近似行
                pos = 0
                real_line = -1
                for i, ln in enumerate(lines, 1):
                    pos += len(norm(ln))
                    if pos >= idx:
                        real_line = i
                        break
                problems[-1] += f"  [actual ~line {real_line}]"

    for p in problems:
        print("  " + p)
    print()
    print(f"  matched={ok}  mismatch={bad}  skipped={skipped}")
    print(f"  => {'PASS' if bad == 0 else 'FAIL'}")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
