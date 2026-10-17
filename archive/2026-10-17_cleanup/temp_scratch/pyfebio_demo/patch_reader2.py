"""给 read_lig.py 的补丁表补上双引号形式。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\read_lig.py")
s = p.read_text(encoding="utf-8")

old = """for _o, _n in (("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),"""
new = """for _o, _n in (('domain["name"]', 'domain.get("name", str(domain["id"][0]))'),
               ('surface["name"]', 'surface.get("name", str(surface["id"][0]))'),
               (\"{domain['name']}\", \"{domain.get('name', str(domain['id'][0]))}\"),"""
assert s.count(old) == 1, s.count(old)
s = s.replace(old, new)
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
