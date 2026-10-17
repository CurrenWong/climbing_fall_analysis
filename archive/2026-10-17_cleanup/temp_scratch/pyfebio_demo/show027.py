"""打印 027 两处的精确行文本。"""
import io

L = r"D:\Project\climbing_fall_analysis\.learnings\LEARNINGS.md"
for i, l in enumerate(io.open(L, encoding="utf-8").read().splitlines()):
    if "LRN-20261004-027" in l:
        print(i + 1, repr(l))
