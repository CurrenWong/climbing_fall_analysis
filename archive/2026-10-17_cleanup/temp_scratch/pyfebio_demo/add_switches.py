"""加 NO_CONTACT / LIG_SUB 两个二分开关。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = 'ONLY = [x for x in os.environ.get("BONES", "").split(",") if x]'
assert s.count(old) == 1
s = s.replace(old, old + '''
NO_CONTACT = os.environ.get("NO_CONTACT", "0") == "1"
LIG_SUB = [int(x) for x in os.environ.get("LIG_SUB", "").split(",") if x]''')

old2 = 'print(f"[mode] BONES={ONLY or \'all\'}")'
assert s.count(old2) == 1
s = s.replace(old2, old2 + '\nprint(f"[mode] NO_CONTACT={NO_CONTACT} LIG_SUB={LIG_SUB or \'all\'}")')

# 韧带子集过滤
old3 = """    n_lig = int(lg["n_lig"])
    for i in range(n_lig):"""
assert s.count(old3) == 1
s = s.replace(old3, """    n_lig = int(lg["n_lig"])
    for i in range(n_lig):
        if LIG_SUB and (i + 1) not in LIG_SUB:
            continue""")

# 接触开关
old4 = """w('  <Contact>')
w('    <contact name="tibiotalar" surface_pair="tibiotalar_pair" '"""
assert s.count(old4) == 1
s = s.replace(old4, """if not NO_CONTACT:
  w('  <Contact>')
  w('    <contact name="tibiotalar" surface_pair="tibiotalar_pair" '""" + "\n" + """  w('      <laugon>PENALTY</laugon>')
  w('      <penalty>1</penalty>')
  w('      <auto_penalty>1</auto_penalty>')
  w('      <two_pass>0</two_pass>')
  w('      <node_reloc>0</node_reloc>')
  w('      <symmetric_stiffness>0</symmetric_stiffness>')
  w('      <tolerance>0.02</tolerance>')
  w('    </contact>')
  w('  </Contact>')""")
# 删掉原来的接触块正文
old5 = """  w('      <laugon>PENALTY</laugon>')
  w('      <penalty>1</penalty>')
  w('      <auto_penalty>1</auto_penalty>')
  w('      <two_pass>0</two_pass>')
  w('      <node_reloc>0</node_reloc>')
  w('      <symmetric_stiffness>0</symmetric_stiffness>')
  w('      <tolerance>0.02</tolerance>')
  w('    </contact>')
  w('  </Contact>')
w('  <Output>')"""
assert s.count(old5) == 1, f"old5 {s.count(old5)}"
s = s.replace(old5, "w('  <Output>')")
# 清掉重复的 Contact 块（原块的第一行仍在）
old6 = """w('  <Contact>')
  w('    <contact name="tibiotalar" surface_pair="tibiotalar_pair" '
  'type="sliding-elastic">')
  w('      <laugon>PENALTY</laugon>')
  w('      <penalty>1</penalty>')
  w('      <auto_penalty>1</auto_penalty>')
  w('      <two_pass>0</two_pass>')
  w('      <node_reloc>0</node_reloc>')
  w('      <symmetric_stiffness>0</symmetric_stiffness>')
  w('      <tolerance>0.02</tolerance>')
  w('    </contact>')
  w('  </Contact>')
w('  <Output>')"""
if old6 in s:
    s = s.replace(old6, "w('  <Output>')")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
