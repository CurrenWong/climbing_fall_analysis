"""修 bug：prism15 也被元素生成分支接受（原判断只认 "prism"）。"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = '        if CART_ELEM == "prism" and m in (3, 4):'
assert s.count(old) == 1, s.count(old)
s = s.replace(old, '        if CART_ELEM in ("prism", "prism15") and m in (3, 4):')
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched，syntax OK")

# 复核：确认两处判断都一致
import re
for i, l in enumerate(s.split("\n"), 1):
    if "CART_ELEM" in l:
        print(f"  {i:5d}| {l.strip()[:100]}")
