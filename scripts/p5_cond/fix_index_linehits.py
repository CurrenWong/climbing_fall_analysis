"""重新定位出处索引 表A 每行的行号（空白不敏感精确匹配 + 标识符窗口回退）。

用法：
  python fix_index_linehits.py            # 只打印 old -> new
  python fix_index_linehits.py --apply    # 打印并写回索引
输出末尾列出**无法自动判定**的行（需人工看），这些行不会被改。
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
IDX = ROOT / "docs/paper/论文v2_出处索引.md"
WIN = 14
APPLY = "--apply" in sys.argv


def norm(s: str) -> str:
    return re.sub(r"\s+", "", s.replace("`", "").replace("\\", ""))


def tokens(s: str) -> set[str]:
    s = re.sub(r"\$[^$]*\$", " ", s).replace("`", "")
    return set(re.findall(r"[A-Za-z_][A-Za-z0-9_]{3,}", s)) - {
        "self", "float", "none", "true", "false", "return", "numpy"}


def locate(lines: list[str], code: str) -> tuple[int, int, str] | None:
    """返回 (start, end, how) 或 None。"""
    frags = [norm(x) for x in re.split(r";", code) if len(norm(x)) >= 10]
    # 策略 A：逐 fragment 精确搜（取命中最多的连续区）
    best = None
    for frag in frags:
        for i, ln in enumerate(lines):
            if frag in norm(ln):
                # 收集该 fragment 附近命中的其它 fragment
                lo, hi = i + 1, i + 1
                for f2 in frags:
                    for j in range(max(0, i - 6), min(len(lines), i + 7)):
                        if f2 in norm(lines[j]):
                            lo = min(lo, j + 1)
                            hi = max(hi, j + 1)
                            break
                if best is None or (hi - lo) < (best[1] - best[0]):
                    best = (lo, hi)
        if best:
            break
    if best:
        return best[0], best[1], "exact"
    # 策略 B：标识符窗口
    want = tokens(code)
    if not want:
        return None
    scored = []
    for i in range(len(lines)):
        w = "\n".join(lines[i:i + WIN])
        h = len(want & set(re.findall(r"[A-Za-z_][A-Za-z0-9_]+", w)))
        scored.append((h, i))
    h, i = max(scored)
    if h / len(want) < 0.6:
        return None
    return i + 1, min(len(lines), i + WIN), f"ident({h}/{len(want)})"


def main() -> int:
    txt = IDX.read_text(encoding="utf-8")
    sec_m = re.search(r"## 2\. 表 A.*?(?=\n---|\n## 3\.)", txt, re.S)
    sec = sec_m.group(0)
    rows = list(re.finditer(r"^\|(.+)\|\s*$", sec, re.M))

    print("=" * 96)
    print("表 A 行号重定位   old(indexed) -> new(measured)")
    print("=" * 96)
    print(f"{'old':<16} {'new':<16} {'how':<10} 描述")
    print("-" * 96)

    edits: list[tuple[str, str]] = []   # (old_loc_cell, new_loc_cell)
    unresolved: list[str] = []
    for m in rows:
        row = m.group(1)
        cells = [c.strip() for c in row.split("|")]
        if len(cells) < 3 or cells[0].startswith("---") or not cells[0]:
            continue
        desc, loc, code = cells[0], cells[1], cells[2]
        lm = re.search(r"`([^`]+?):(\d+)(?:-(\d+))?`", loc)
        if not lm:
            continue
        path = lm.group(1)
        a, b = int(lm.group(2)), int(lm.group(3) or lm.group(2))
        f = ROOT / path
        if not f.exists():
            unresolved.append(f"{path}  FILE MISSING  ({desc[:40]})")
            continue
        lines = f.read_text(encoding="utf-8", errors="replace").splitlines()
        got = locate(lines, code)
        if got is None:
            unresolved.append(f"{path}:{a}-{b}  NOT LOCATED  ({desc[:40]})")
            continue
        na, nb, how = got
        # 置信度：exact 直接信；ident 需覆盖率高、token 多、且不在文件头文档串区
        cov = 1.0
        if how.startswith("ident("):
            k, n = how[6:-1].split("/")
            cov = int(k) / int(n)
        confident = (how == "exact"
                     or (cov >= 0.8 and tokens(code).__len__() >= 6 and na > 80))
        # 已知假匹配：文件头部（docstring 引用同一段代码）
        if na < 80 and a > 80:
            confident = False
        tag = ("same" if (na <= a and nb >= b)
               else ("MOVE" if confident else "UNSURE"))
        print(f"{a}-{b:<11} {na}-{nb:<11} {how:<10} {tag:<6} {desc[:44]}")
        if tag == "UNSURE":
            unresolved.append(
                f"{path}:{a}-{b} -> {na}-{nb} ({how})  ({desc[:36]})")
            continue
        if tag == "MOVE":
            old_cell = f"`{path}:{a}`" if a == b else f"`{path}:{a}-{b}`"
            new_cell = f"`{path}:{na}`" if na == nb else f"`{path}:{na}-{nb}`"
            edits.append((old_cell, new_cell))

    print()
    if unresolved:
        print("-- 无法自动判定（不会改动，需人工核）--")
        for u in unresolved:
            print("   " + u)
    print(f"\n  total edits to apply: {len(edits)}  unresolved: {len(unresolved)}")

    if APPLY and edits:
        t = txt
        for old, new in edits:
            if t.count(old) == 1:
                t = t.replace(old, new)
            else:
                print(f"  !! ambiguous ({t.count(old)}x), skipped: {old}")
        IDX.write_text(t, encoding="utf-8")
        print(f"  APPLIED -> {IDX.relative_to(ROOT)}")
    elif APPLY:
        print("  nothing to apply")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
