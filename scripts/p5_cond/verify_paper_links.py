"""出处索引 / 论文 的引用链路核对（第 1 步：路径与行号可解析）。

检查项：
  A. 出处索引里所有反引号路径是否存在（含项目短路径惯例 docs/X -> docs/reports/X）
  B. 形如 `path:123-456` / `path:123` 的行号引用是否在文件长度内
  C. 论文正文里所有反引号路径是否存在
不检查语义，只检查"指得到 / 数得着"。
"""
from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(r"D:\Project\climbing_fall_analysis")
IDX = ROOT / "docs" / "paper" / "论文v2_出处索引.md"
PAPER = ROOT / "docs" / "paper" / "论文v2_数学建模版.md"

SKIP_PARTS = {".venv", ".git", ".pytest_cache", "model", "__pycache__"}

# 短路径回退顺序（项目惯例：写 docs/X 实际在 docs/reports/ 等）
FALLBACK_DIRS = ["docs", "docs/reports", "docs/methods", "docs/notes",
                 "docs/paper", "docs/s5-contact", "docs/archive", "results",
                 "results/opensim_fe", "scripts", "src", "src/climbing",
                 "src/climbing/coupling", ""]


def _prune(p: pathlib.Path) -> bool:
    return bool(set(p.parts) & SKIP_PARTS)


def resolve(ref: str) -> pathlib.Path | None:
    """按短路径惯例解析：先精确，再按 basename 全库回退（含花括号展开）。"""
    ref = ref.split("::")[0].strip()          # pytest 节点法 tests/x.py::Cls
    if ref.startswith("."):                    # 碎片（如 `.csv`）不当路径
        return None
    # 花括号展开 scripts/p5_cond/{a, b}.py -> 多个
    m = re.match(r"^(.*?)\{([^}]*)\}(.*)$", ref)
    if m:
        pre, items, post = m.groups()
        for it in items.split(","):
            c = resolve(f"{pre}{it.strip()}{post}")
            if c is None:
                return None
        return ROOT  # 全部展开项都存在 => 视为有效
    direct = ROOT / ref
    if direct.exists():
        return direct
    name = pathlib.PurePosixPath(ref).name
    if "/" in ref and name == ref.split("/")[-1]:
        pass
    for d in FALLBACK_DIRS:
        c = (ROOT / d / name) if d else (ROOT / name)
        if c.exists():
            return c
    # 最后兜底：全库按 basename 找一次（限制在项目自有目录）
    for pat in ("*.md", "*.py", "*.json", "*.csv", "*.txt"):
        for f in ROOT.rglob(pat):
            if _prune(f) or f.name != name:
                continue
            return f
    return None


def extract_paths(text: str) -> set[str]:
    out = set()
    for m in re.finditer(r"`([^`\n]+?)`", text):
        s = m.group(1).strip()
        if "/" not in s and not re.search(r"\.(md|py|json|csv|txt|feb)$", s):
            continue
        if any(ch in s for ch in (" ", "(", ")", "→", "=")) and not re.search(
                r"\.(md|py|json|csv|txt)$", s):
            continue
        # 去掉行号
        s2 = re.sub(r":\d+(?:-\d+)?$", "", s)
        if s2.startswith(".") and "/" not in s2:
            continue                      # 碎片（`.csv`、`*.py` 后缀片段）
        if re.search(r"\.(md|py|json|csv|txt|feb)$", s2) or "/" in s2:
            out.add(s2)
    return out


def extract_line_refs(text: str) -> list[tuple[str, int, int | None]]:
    out = []
    for m in re.finditer(r"`([^`\n]+?):(\d+)(?:-(\d+))?`", text):
        out.append((m.group(1), int(m.group(2)),
                    int(m.group(3)) if m.group(3) else None))
    return out


def main() -> int:
    fails = []
    print("=" * 74)
    print("A. 出处索引 -> 路径存在性")
    print("=" * 74)
    idx_text = IDX.read_text(encoding="utf-8")
    paths = sorted(extract_paths(idx_text))
    bad_p = 0
    for p in paths:
        r = resolve(p)
        if r is None:
            print(f"  DEAD   {p}")
            bad_p += 1
            fails.append(f"idx path {p}")
    print(f"  checked={len(paths)}  dead={bad_p}  "
          f"{'PASS' if bad_p == 0 else 'FAIL'}")

    print()
    print("=" * 74)
    print("B. 出处索引 -> 行号引用有效")
    print("=" * 74)
    lrefs = extract_line_refs(idx_text)
    bad_l = 0
    for path, a, b in lrefs:
        r = resolve(path)
        if r is None:
            print(f"  DEAD FILE  {path}:{a}" + (f"-{b}" if b else ""))
            bad_l += 1
            fails.append(f"line ref file {path}")
            continue
        try:
            n = len(r.read_text(encoding="utf-8", errors="replace")
                     .splitlines())
        except Exception:
            n = -1
        hi = b if b else a
        if a < 1 or hi > n:
            print(f"  OUT-OF-RANGE  {path}:{a}-{b}  (file has {n} lines)")
            bad_l += 1
            fails.append(f"line range {path}:{a}-{b}")
    print(f"  checked={len(lrefs)}  bad={bad_l}  "
          f"{'PASS' if bad_l == 0 else 'FAIL'}")

    print()
    print("=" * 74)
    print("C. 论文正文 -> 路径存在性")
    print("=" * 74)
    pt = PAPER.read_text(encoding="utf-8")
    ppaths = sorted(extract_paths(pt))
    bad_c = 0
    for p in ppaths:
        if resolve(p) is None:
            print(f"  DEAD   {p}")
            bad_c += 1
            fails.append(f"paper path {p}")
    print(f"  checked={len(ppaths)}  dead={bad_c}  "
          f"{'PASS' if bad_c == 0 else 'FAIL'}")

    print()
    print("=" * 74)
    print("RESULT  " + ("PASS" if not fails else "FAIL: " + repr(fails)))
    print("=" * 74)
    return 0 if not fails else 1


if __name__ == "__main__":
    raise SystemExit(main())
