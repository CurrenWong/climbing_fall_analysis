"""给 parse_state 加兜底：找不到 set_id 时用占位名（pyfebio 对自定义 Surface 会 KeyError）。"""
import textwrap
import inspect
from pathlib import Path
import ast

# 先看源码里出错的几行
import pyfebio.xplt as X
src = textwrap.dedent(inspect.getsource(X.parse_state))
for i, l in enumerate(src.split("\n"), 1):
    if "set_lut" in l or "mesh_dict" in l:
        print(f"{i:4d}| {l}")

p = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\read_lig.py")
s = p.read_text(encoding="utf-8")

old = 'X.parse_mesh = _ns["parse_mesh"]'
assert s.count(old) == 1
s = s.replace(old, old + '''

# ── 第二个容错补丁：parse_state 对**自定义 Surface 集合**会 KeyError
#   （MEMORY 已记录该坑）—— 找不到 set_id 时用占位名兜底。──
_ps = textwrap.dedent(inspect.getsource(X.parse_state))
_ps = _ps.replace(
    'set_name = mesh_dict[set_lut[key]][set_id]',
    'set_name = mesh_dict[set_lut[key]].get(set_id) or '
    'mesh_dict[set_lut[key]].get(int(set_id)) or f"set{set_id}"')
if 'mesh_dict[set_lut[key]].get(set_id)' not in _ps:
    print("!! parse_state 补丁未命中，需人工检查")
_ns2 = X.__dict__
exec(compile(_ps, "<patched_parse_state>", "exec"), _ns2)
X.parse_state = _ns2["parse_state"]
print("pyfebio parse_state 已打容错补丁")''')

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("\npatched + syntax OK")
