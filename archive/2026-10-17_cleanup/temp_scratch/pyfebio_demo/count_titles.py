"""精确统计：只数「## [LRN-...]」标题行；并确认攀岩 055 条目的 ID。"""
import re
from pathlib import Path

LT = Path(r"D:\Project\climbing_fall_analysis\.learnings/LEARNINGS.md")
lines = LT.read_text(encoding="utf-8").split("\n")

titles = [(i + 1, l) for i, l in enumerate(lines) if re.match(r"^##\s*\[LRN-", l)]
print(f"文件 {LT.stat().st_size:,} B / {len(lines)} 行；**标题行** {len(titles)} 条\n")
ids = [re.search(r"\[LRN-(\d{8}-\d{3})\]", l).group(1) for _, l in titles]
dup = sorted({x for x in ids if ids.count(x) > 1})
print(f"标题行中的重复 ID（{len(dup)}）：{dup}")
for d in dup:
    print(f"\n  {d} 出现在：")
    for ln, l in titles:
        if d in l:
            print(f"    行{ln:>5}: {l[:100]}")

print("\n" + "=" * 74)
print("攀岩跌落线的 055 条目（搜关键词）")
print("=" * 74)
for i, l in enumerate(lines):
    if "骨折判据" in l and l.startswith("## "):
        print(f"  行{i+1}: {l[:120]}")
    if "Y25" in l and l.startswith("## "):
        print(f"  行{i+1}: {l[:120]}")
print("\n  含 '20261005-055' 的所有行：")
for i, l in enumerate(lines):
    if "20261005-055" in l:
        print(f"    行{i+1}: {l.strip()[:110]}")
print("\n  末尾 3 条标题：")
for ln, l in titles[-3:]:
    print(f"    行{ln}: {l[:110]}")
