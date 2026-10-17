"""加 FIX_NAVICULAR 开关（默认开）：舟骨只靠两条软韧带吊着 ⇒ 近奇异。

同时加 FIX_FIBULA_TIP（默认关）备用。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = 'NO_ROTLOCK = os.environ.get("NO_ROTLOCK", "0") == "1"'
assert s.count(old) == 1
s = s.replace(old, old + '''
# ★ 舟骨只靠 胫舟(→胫骨) + 弹簧(→跟骨) 两条软壳约束 ⇒ 近奇异（LRN-016/017 铁律）。
#   真实足部它被大量未建模结构固定，这里固定它是合理的建模选择。
FIX_NAVICULAR = os.environ.get("FIX_NAVICULAR", "1") == "1"''')

old2 = 'print(f"[mode] NO_CONTACT={NO_CONTACT} LIG_SUB={LIG_SUB or \'all\'}")'
assert s.count(old2) == 1
s = s.replace(old2, old2 + '\nprint(f"[mode] FIX_NAVICULAR={FIX_NAVICULAR}")')

# 在 Boundary 里加舟骨全约束
old3 = 'w(\'  </Boundary>\')'
assert s.count(old3) == 1
s = s.replace(old3, """if FIX_NAVICULAR and "navicular_spon" in dom_nodes:
    _nav = sorted(dom_nodes["navicular_spon"])
    w(f'    <NodeSet name="navicular_all_ns">'
      f'{",".join(str(x) for x in _nav)}</NodeSet>')
    w('    <bc name="fix_navicular" type="zero displacement" '
      'node_set="navicular_all_ns">')
    w('      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>')
    w('    </bc>')
w('  </Boundary>')""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
