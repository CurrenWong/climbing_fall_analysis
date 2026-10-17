"""修：舟骨 NodeSet 必须定义在 <Mesh> 内，不能放 Boundary 段。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

BLOCK = '''if FIX_NAVICULAR and "navicular_spon" in dom_nodes:
    _nav = sorted(dom_nodes["navicular_spon"])
    w(f'    <NodeSet name="navicular_all_ns">'
      f'{",".join(str(x) for x in _nav)}</NodeSet>')
'''
assert s.count(BLOCK) == 1, f"blk {s.count(BLOCK)}"
s = s.replace(BLOCK, "")

# BC 部分不再定义 NodeSet
BC = '''    w('    <bc name="fix_navicular" type="zero displacement" '
      'node_set="navicular_all_ns">')'''
assert s.count(BC) == 1
s = s.replace(BC, '''    w('    <bc name="fix_navicular" type="zero displacement" '
      'node_set="navicular_all_ns">')''')

# NodeSet 插到既有 NodeSet 区（tibia_top_ns 附近）
ANCHOR = '''if LIG_ONLY_NS:
    w(f'    <NodeSet name="lig_only_ns">'
      f'{",".join(str(x) for x in LIG_ONLY_NS)}</NodeSet>')'''
assert s.count(ANCHOR) == 1
s = s.replace(ANCHOR, ANCHOR + "\n" + BLOCK.rstrip())

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
