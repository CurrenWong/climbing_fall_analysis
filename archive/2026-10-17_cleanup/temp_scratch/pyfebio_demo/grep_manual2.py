"""从已提取的 FEBio 手册文本里查 contact / SurfacePair / 段序 的确切写法。"""
import re
from pathlib import Path

T = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_manual.txt").read_text(
    encoding="utf-8")

for kw in ("<contact ", "sliding-elastic", "surface_pair", "traction Load",
           "3.14", "Contact Section"):
    print(f"\n{'='*20} {kw!r} {'='*20}")
    n = 0
    for m in re.finditer(re.escape(kw), T):
        s = max(0, m.start() - 300)
        print("  ..." + T[s:m.start() + 420].replace("\n", " | ") + "...")
        n += 1
        if n >= 2:
            break
    if n == 0:
        print("  (未找到)")
