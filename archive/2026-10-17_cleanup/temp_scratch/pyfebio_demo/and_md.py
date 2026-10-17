"""在 paper/*.md（MinerU 抽取版）里搜 Anderson 弹簧/刚性/约束的关键句。"""
import re
from pathlib import Path

D = Path(r"D:\Project\climbing_fall_analysis\paper")
KEYS = ("500 N/mm", "500N/mm", "weakly restrained", "rigid", "linear spring",
        "inversion-eversion", "internal-external", "seating", "pre-loading",
        "preloading", "50 N/mm")
for f in sorted(D.glob("*.md")):
    t = f.read_text(encoding="utf-8", errors="replace")
    hit = {k: t.lower().count(k.lower()) for k in KEYS if k.lower() in t.lower()}
    if hit:
        print(f"\n{'='*76}\n### {f.name}  ({len(t)} chars)\n{'='*76}")
        print("  命中:", hit)
