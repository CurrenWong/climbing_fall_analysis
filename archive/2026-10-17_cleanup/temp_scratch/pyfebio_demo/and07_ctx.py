"""抽 Anderson 2007 的弹簧/A-P 约束口径。"""
import re
from pathlib import Path

D = Path(r"D:\Project\climbing_fall_analysis\paper")
F = D / "MinerU_markdown_anderson2007_2106695115320401920.md"
t = re.sub(r"\n{3,}", "\n\n", F.read_text(encoding="utf-8", errors="replace"))
for key in ("50 N/mm", "linear spring", "rigid"):
    for m in re.finditer(re.escape(key), t, re.I):
        a, b = max(0, m.start() - 1000), min(len(t), m.end() + 1000)
        print(f"\n{'#'*76}\n### {key!r}\n{'#'*76}")
        print(t[a:b].replace("\n", " "))
