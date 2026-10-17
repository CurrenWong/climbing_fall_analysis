"""看 cartilage.py 的偏移生成逻辑。"""
import re
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\cartilage.py")
L = p.read_text(encoding="utf-8").split("\n")
print(f"{len(L)} 行")
for i, l in enumerate(L, 1):
    if re.match(r"^\s*(def |class |# ──|# =)", l) or "normal" in l.lower():
        print(f"{i:5d}| {l[:118]}")
