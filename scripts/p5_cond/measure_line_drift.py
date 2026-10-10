"""表 A 行号漂移精确测量。

方法：对每个片段抽 identifier token（len>=4），在整文件里找"命中 token 数最多"的
连续 12 行窗口，报告 indexed 区间 vs best-match 区间。
判据：
  HIT    best 区间与 indexed 区间重叠（或偏移 <= WINDOW）=> 行号可用
  DRIFT  best 找得到但明显偏移  => 行号已漂移（需修）
  NONE   找不到任何合理位置     => 片段是转写/已改写（行号本身或文本需复核）
"""
from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
IDX = ROOT / "docs/paper/论文v2_出处索引.md"
WIN = 12


def tokens(s: str) -> set[str]:
    s = re.sub(r"\$[^$]*\$", " ", s).replace("`", "")
    toks = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]{3,}", s))
    for noise in ("self", "float", "none", "true", "false", "return",
                  "numpy", "float64", "shape", "dtype"):
        toks.discard(noise)
    return toks


def best_window(lines: list[str], want: set[str]) -> tuple[int, int, int]:
    """返回 (start, end, hits) —— 命中 token 最多的窗口。"""
    best = (0, 0, -1)
    for i in range(len(lines)):
        w = "\n".join(lines[i:i + WIN])
        h = len(want & set(re.findall(r"[A-Za-z_][A-Za-z0-9_]+", w)))
        if h > best[2]:
            best = (i + 1, min(len(lines), i + WIN), h)
    return best


def main() -> int:
    txt = IDX.read_text(encoding="utf-8")
    sec = re.search(r"## 2\. 表 A.*?(?=\n---|\n## 3\.)", txt, re.S).group(0)

    drift = hit = none = 0
    rows_out = []
    for row in re.findall(r"^\|(.+)\|\s*$", sec, re.M):
        cells = [c.strip() for c in row.split("|")]
        if len(cells) < 3 or cells[0].startswith("---") or cells[0] in ("",):
            continue
        desc, loc, code = cells[0], cells[1], cells[2]
        lm = re.search(r"`([^`]+?):(\d+)(?:-(\d+))?`", loc)
        if not lm:
            continue
        path, a, b = lm.group(1), int(lm.group(2)), int(lm.group(3) or lm.group(2))
        f = ROOT / path
        if not f.exists():
            for d in ("", "src", "scripts", "src/climbing"):
                c = (ROOT / d / pathlib.PurePosixPath(path).name) if d else \
                    (ROOT / pathlib.PurePosixPath(path).name)
                if c.exists():
                    f = c
                    break
            else:
                rows_out.append(("NOFILE", desc, path, a, b, 0, 0))
                none += 1
                continue
        lines = f.read_text(encoding="utf-8", errors="replace").splitlines()
        want = tokens(code)
        if not want:
            continue
        bs, be, bh = best_window(lines, want)
        cover = bh / max(1, len(want))
        overlap = not (be < a or bs > b)
        shift = min(abs(bs - a), abs(be - b))
        if overlap or shift <= WIN:
            hit += 1
            tag = "HIT"
        elif cover >= 0.5:
            drift += 1
            tag = "DRIFT"
        else:
            none += 1
            tag = "NONE"
        rows_out.append((tag, desc, path, a, b, bs, be, bh, len(want)))

    print("=" * 90)
    print("表 A 行号核对  (indexed 区间 vs 全文件最佳匹配窗口)")
    print("=" * 90)
    print(f"{'tag':<7} {'indexed':<20} {'best-match':<16} "
          f"{'tok':>4}  描述")
    print("-" * 90)
    for r in rows_out:
        if r[0] == "NOFILE":
            print(f"{r[0]:<7} {r[2]}:{r[3]}-{r[4]}  FILE MISSING")
            continue
        tag, desc, path, a, b, bs, be, bh, nw = r
        if tag == "HIT":
            continue                      # 对上的不打，减少噪音
        print(f"{tag:<7} {a}-{b:<14} {bs}-{be:<10} "
              f"{bh}/{nw:<4} {desc[:46]}")
    print()
    print(f"  HIT={hit}  DRIFT={drift}  NONE={none}")
    print(f"  => {'PASS (all line refs usable)' if drift + none == 0 else 'ATTENTION'}")
    return 0


if __name__ == "__main__":
    main()
