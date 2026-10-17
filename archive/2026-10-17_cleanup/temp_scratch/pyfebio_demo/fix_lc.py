"""去掉未定义的 lc 引用（单元测试没有 <LoadData>）。"""
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\elem_test2.py")
s = p.read_text(encoding="utf-8")
n = s.count('<value lc="1">')
s = s.replace('<value lc="1">', '<value>')
p.write_text(s, encoding="utf-8")
print(f"替换 {n} 处")
