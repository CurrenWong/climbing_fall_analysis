"""补 parse_state 的两处 set_id 查找（含 +1 那个）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\read_lig.py")
s = p.read_text(encoding="utf-8")

old = """_ps = _ps.replace(
    'set_name = mesh_dict[set_lut[key]][set_id]',
    'set_name = mesh_dict[set_lut[key]].get(set_id) or '
    'mesh_dict[set_lut[key]].get(int(set_id)) or f"set{set_id}"')"""
assert s.count(old) == 1, s.count(old)
new = """for _a, _b in (
    ('set_name = mesh_dict[set_lut[key]][set_id + 1]',
     'set_name = mesh_dict[set_lut[key]].get(set_id + 1) or '
     'mesh_dict[set_lut[key]].get(int(set_id) + 1) or f"set{set_id}"'),
    ('set_name = mesh_dict[set_lut[key]][set_id]',
     'set_name = mesh_dict[set_lut[key]].get(set_id) or '
     'mesh_dict[set_lut[key]].get(int(set_id)) or f"set{set_id}"'),
):
    _ps = _ps.replace(_a, _b)"""
s = s.replace(old, new)
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
