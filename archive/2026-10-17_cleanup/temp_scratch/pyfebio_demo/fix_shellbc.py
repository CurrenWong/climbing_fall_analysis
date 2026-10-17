"""把无效的 <u_rotation> BC 换成 FEBio 正确写法：zero shell displacement + sx/sy/sz。

手册 §3.11.7：壳的附加自由度叫 back-displacement，dof 名是 sx/sy/sz。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = """    w('    <bc name="lig_rot_lock" type="zero displacement" '
      f'node_set="lig_only_ns">')
    w('      <u_rotation>1</u_rotation>')
    w('      <v_rotation>1</v_rotation>')
    w('      <w_rotation>1</w_rotation>')
    w('    </bc>')"""
assert s.count(old) == 1, f"hit {s.count(old)}"
s = s.replace(old, """    # ★ FEBio 壳的附加自由度叫 back-displacement，dof 名是 sx/sy/sz，
    #   不是 <x_rotation> / <u_rotation>（那两个标签 FEBio 不认）。
    w('    <bc name="lig_back_lock" type="zero shell displacement" '
      f'node_set="lig_only_ns">')
    for _d in ("sx", "sy", "sz"):
        w(f'      <dof>{_d}</dof><value lc="1">0</value><relative>0</relative>')
    w('    </bc>')""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
