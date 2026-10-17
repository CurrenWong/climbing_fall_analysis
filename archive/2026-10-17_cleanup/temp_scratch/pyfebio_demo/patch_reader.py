"""给 read_lig.py 加 pyfebio 容错补丁（shell/beam 域没有 name 字段）。"""
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\read_lig.py")
s = p.read_text(encoding="utf-8")

old = "import pyfebio.xplt as X\nX.to_hdf5(str(src), str(dst))"
assert s.count(old) == 1
s = s.replace(old, '''import textwrap
import inspect
import pyfebio.xplt as X

# ── 容错补丁：pyfebio 的 parse_mesh 对 **shell/beam 域**没有 "name" 字段，
#   会 KeyError 'name'（MEMORY 里记为「需容错补丁」）。用 id 兜底。──
_src = textwrap.dedent(inspect.getsource(X.parse_mesh))
for _o, _n in (("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
               ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
               ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
               ("{nodeset['name']}", "{nodeset.get('name', str(nodeset['id'][0]))}"),
               ("{elementset['name']}",
                "{elementset.get('name', str(elementset['id'][0]))}")):
    _src = _src.replace(_o, _n)
_ns = X.__dict__
exec(compile(_src, "<patched_parse_mesh>", "exec"), _ns)
X.parse_mesh = _ns["parse_mesh"]
print("pyfebio parse_mesh 已打容错补丁")

X.to_hdf5(str(src), str(dst))''')

p.write_text(s, encoding="utf-8")
import ast
ast.parse(s)
print("patched + syntax OK")
