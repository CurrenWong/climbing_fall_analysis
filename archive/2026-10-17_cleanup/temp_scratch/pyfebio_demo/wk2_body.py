import io
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\paper/_text/吴恺_足踝FE.txt")
t = io.open(p, encoding="utf-8", errors="ignore").read()
L = t.split("\n")
# 去掉纯空行后的紧凑视图
compact = [(i, l.strip()) for i, l in enumerate(L) if l.strip()]
print(f"非空行 {len(compact)}  总行 {len(L)}")
for n, (i, l) in enumerate(compact):
    if 100 <= n <= 230:
        print(f"{n:4d}[{i:4d}] {l[:150]}")
