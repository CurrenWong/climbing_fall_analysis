"""验证"刚体模态"假设：给「只属韧带」的节点加转动零位移约束。

若由此收敛 ⇒ ①假设成立 ②shell 模式即可用（近似但物理上合理：这些节点在韧带
内部/侧边，约束其转动对整体拉伸行为影响很小）。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# 计算 ligament-only 节点（shell 模式）
old = "lig_blocks = {}\nLIG_ETYPE ="
assert s.count(old) == 1
s = s.replace(old, """# 韧带专属节点（仅被壳用、不被任何实体用）—— shell 模式下它们的**转动
# 自由度无约束** ⇒ 刚体模态（实测：首个载荷步 22 个单元翻负）。
LIG_ONLY_NS = []
if LIG_MODE == "shell":
    _solid_nodes = set()
    for _k, _v in blocks.items():
        for _c in _v:
            _solid_nodes.update(_c)
    _lig_nodes = {x for _v in lig_blocks.values() for _c in _v for x in _c}
    LIG_ONLY_NS = sorted(_lig_nodes - _solid_nodes)
    print(f"韧带专属节点（转动自由度原本自由）: {len(LIG_ONLY_NS)}")
    if NO_ROTLOCK:
        LIG_ONLY_NS = []

lig_blocks = {}
LIG_ETYPE =""")

# 开关
old2 = 'NO_CONTACT = os.environ.get("NO_CONTACT", "0") == "1"'
assert s.count(old2) == 1
s = s.replace(old2, old2 + '\nNO_ROTLOCK = os.environ.get("NO_ROTLOCK", "0") == "1"')

# BC：转动约束
old3 = "w('  </Boundary>')"
assert s.count(old3) == 1
s = s.replace(old3, """if LIG_ONLY_NS:
    w('    <bc name="lig_rot_lock" type="zero displacement" '
      f'node_set="lig_only_ns">')
    w('      <u_rotation>1</u_rotation>')
    w('      <v_rotation>1</v_rotation>')
    w('      <w_rotation>1</w_rotation>')
    w('    </bc>')
w('  </Boundary>')""")

# NodeSet 定义（放在 NodeSet 区）
old4 = 'w(\'    <Surface name="tibia_top">\')'
assert s.count(old4) == 1
s = s.replace(old4, """if LIG_ONLY_NS:
    w(f'    <NodeSet name="lig_only_ns">'
      f'{",".join(str(x) for x in LIG_ONLY_NS)}</NodeSet>')
w('    <Surface name="tibia_top">')""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
