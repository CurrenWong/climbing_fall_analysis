import re
from pathlib import Path

for f in ("docs/方案更新_v2.md", "阶段总结.md"):
    p = Path(r"D:\Project\climbing_fall_analysis") / f
    if not p.exists():
        print(f"MISS {f}")
        continue
    t = p.read_text(encoding="utf-8")
    print(f"\n{'=' * 70}\n{f}   ({len(t):,} chars)\n{'=' * 70}")
    for ln in t.split("\n"):
        if re.match(r"^#{1,3} ", ln):
            print("  " + ln)
