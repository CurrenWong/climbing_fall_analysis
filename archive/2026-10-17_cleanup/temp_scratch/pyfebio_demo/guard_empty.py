"""健壮性：不写空 <Surface> / 空 <NodeSet>，并跳过引用空面的载荷。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# ① Surface 写出：空的不写
for nm, var in (("tibia_top", "top_tibia"), ("fibula_top", "top_fibula")):
    old = f"""w('    <Surface name="{nm}">')
k4 = [f for f in {var} if len(f) == 4]
k3 = [f for f in {var} if len(f) == 3]
for n, f in enumerate(k4, 1):
    w(f'      <quad4 id="{{n}}">{{",".join(str(c) for c in f)}}</quad4>')
for n, f in enumerate(k3, 1):
    w(f'      <tri3 id="{{n}}">{{",".join(str(c) for c in f)}}</tri3>')
w('    </Surface>')"""
    assert s.count(old) == 1, f"{nm} {s.count(old)}"
    s = s.replace(old, f"""if {var}:
    w('    <Surface name="{nm}">')
    k4 = [f for f in {var} if len(f) == 4]
    k3 = [f for f in {var} if len(f) == 3]
    for n, f in enumerate(k4, 1):
        w(f'      <quad4 id="{{n}}">{{",".join(str(c) for c in f)}}</quad4>')
    for n, f in enumerate(k3, 1):
        w(f'      <tri3 id="{{n}}">{{",".join(str(c) for c in f)}}</tri3>')
    w('    </Surface>')""")

# ② 载荷：只给非空面
old2 = """for surf, F in (("tibia_top", F_t), ("fibula_top", F_f)):
    # traction = 力 / 面积，方向沿 -axis（向下压）
    a = area_proj(top_tibia if surf == "tibia_top" else top_fibula)"""
assert s.count(old2) == 1
s = s.replace(old2, """for surf, F in (("tibia_top", F_t), ("fibula_top", F_f)):
    if surf == "tibia_top" and not top_tibia:
        continue
    if surf == "fibula_top" and not top_fibula:
        continue
    # traction = 力 / 面积，方向沿 -axis（向下压）
    a = area_proj(top_tibia if surf == "tibia_top" else top_fibula)""")

# ③ NodeSet：空的不写
for nm, var in (("tibia_top_ns", "top_tibia"), ("fibula_top_ns", "top_fibula"),
                ("calcaneus_bottom_ns", "bot_calc")):
    old3 = (f"w(f'    <NodeSet name=\"{nm}\">'\n"
            f"  f'{{\",\".join(str(x) for x in ns_of({var}))}}</NodeSet>')")
    assert s.count(old3) == 1, f"ns {nm} {s.count(old3)}"
    s = s.replace(old3, f"if {var}:\n" + "\n".join(
        "    " + l for l in old3.split("\n")))

# ④ 引用空 NodeSet 的 BC 跳过
old4 = """w(f'    <bc name="fix_calcaneus" type="zero displacement" '
  'node_set="calcaneus_bottom_ns">')
w('      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>')
w('    </bc>')"""
assert s.count(old4) == 1
s = s.replace(old4, """if bot_calc:
    w(f'    <bc name="fix_calcaneus" type="zero displacement" '
      'node_set="calcaneus_bottom_ns">')
    w('      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>')
    w('    </bc>')""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
