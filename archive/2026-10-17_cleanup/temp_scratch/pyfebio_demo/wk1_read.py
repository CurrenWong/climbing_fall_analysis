"""读吴恺论文文本，定位韧带排序与加载条件。"""
import io
import re
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\paper/_text/吴恺_足踝FE.txt")
t = io.open(p, encoding="utf-8", errors="ignore").read()
print(f"chars {len(t)}")
for kw in ("韧带", "排序", "位移", "载荷", "加载", "结论", "压缩", "背屈", "跖屈"):
    print(f"  {kw}: {len(re.findall(kw, t))} 处")

print("\n" + "=" * 78)
print("=== 含「韧带」的段落 ===")
# 按行聚合，打印含韧带的上下文
lines = [l.strip() for l in t.split("\n")]
shown = 0
for i, l in enumerate(lines):
    if "韧带" in l and len(l) > 8:
        ctx = " | ".join(lines[max(0, i - 1):i + 2])
        print(f"\n[{i}] {ctx[:400]}")
        shown += 1
        if shown > 25:
            break

print("\n" + "=" * 78)
print("=== 全文最后 2500 字符（结论常在此）===")
print(t[-2500:])
