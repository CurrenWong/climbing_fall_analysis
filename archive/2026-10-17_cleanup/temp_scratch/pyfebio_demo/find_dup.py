"""定位 LRN-20261004-027/028/029 的重复条目位置与摘要。"""
import io
import re

LEARN = r"D:\Project\climbing_fall_analysis\.learnings\LEARNINGS.md"
lines = io.open(LEARN, encoding="utf-8").read().splitlines()
for tid in ("LRN-20261004-027", "LRN-20261004-028", "LRN-20261004-029"):
    print(f"\n===== {tid} =====")
    for i, l in enumerate(lines):
        if re.match(rf"^#+\s*\[?{tid}\]?", l):
            print(f"  行 {i+1}: {l.strip()[:100]}")
            for j in range(i + 1, min(i + 4, len(lines))):
                s = lines[j].strip()
                if s:
                    print(f"      {s[:110]}")
                    break
# 另外看最大号，便于分配新号
ids = re.findall(r"LRN-(\d{8})-(\d{3})", io.open(LEARN, encoding="utf-8").read())
byday = {}
for d, n in ids:
    byday.setdefault(d, set()).add(int(n))
print("\n各日最大号：")
for d in sorted(byday):
    print(f"  {d}: max={max(byday[d])}  共 {len(byday[d])} 个不同号")
