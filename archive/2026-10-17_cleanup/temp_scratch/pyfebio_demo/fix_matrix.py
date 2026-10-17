"""修 OK 检测（"N O R M A L   T E R M I N A T I O N" 带空格）并重跑矩阵。"""
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\matrix.py")
s = p.read_text(encoding="utf-8")
old = '    ok = "NORMAL TERMINATION" in o'
assert s.count(old) == 1
s = s.replace(old, '    ok = "NORMAL" in o.replace(" ", "")')
p.write_text(s, encoding="utf-8")
print("fixed")
