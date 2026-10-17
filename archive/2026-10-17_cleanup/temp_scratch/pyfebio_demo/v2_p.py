import re
from pathlib import Path

t = Path(r"D:\Project\climbing_fall_analysis\docs/方案更新_v2.md").read_text(encoding="utf-8")
lines = t.split("\n")
# 找所有含 P0..P6 的行，打印其上下文
P = re.compile(r"\bP[0-6]\b")
seen = set()
for i, ln in enumerate(lines):
    if P.search(ln) and i not in seen:
        print(f"{i:4d}| {ln[:130]}")
