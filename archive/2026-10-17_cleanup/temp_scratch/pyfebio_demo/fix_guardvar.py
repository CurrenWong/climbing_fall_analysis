"""修守卫：传真正的列表对象，不是变量名字符串。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = """for ns, nm, _var in (("tibia_top_ns", "tibia", "top_tibia"),
                    ("fibula_top_ns", "fibula", "top_fibula")):
    if not _var:"""
assert s.count(old) == 1, f"hit {s.count(old)}"
s = s.replace(old, """for ns, nm, _faces in (("tibia_top_ns", "tibia", top_tibia),
                       ("fibula_top_ns", "fibula", top_fibula)):
    if not _faces:      # 空 NodeSet 不能写 BC（FEBio: invalid value for node_set）""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
