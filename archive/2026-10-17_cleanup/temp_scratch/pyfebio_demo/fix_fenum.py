"""统一：surface_faces 返回的是 FE 编号，下游不能再过 node_id_of。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")
n = 0


def rep(old, new):
    global s, n
    c = s.count(old)
    assert c >= 1, f"miss: {old[:70]}"
    s = s.replace(old, new)
    n += c


rep('        P = XYZ[[node_id_of[c] - 1 for c in n]]',
    '        P = XYZ[[c - 1 for c in n]]')
rep('    return sorted({node_id_of[c] for f in faces for c in f})',
    '    return sorted({c for f in faces for c in f})')
rep("""for n, f in enumerate(k4, 1):
    w(f'      <quad4 id="{n}">{",".join(str(node_id_of[c]) for c in f)}</quad4>')
for n, f in enumerate(k3, 1):
    w(f'      <tri3 id="{n}">{",".join(str(node_id_of[c]) for c in f)}</tri3>')""",
    """for n, f in enumerate(k4, 1):
    w(f'      <quad4 id="{n}">{",".join(str(c) for c in f)}</quad4>')
for n, f in enumerate(k3, 1):
    w(f'      <tri3 id="{n}">{",".join(str(c) for c in f)}</tri3>')""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print(f"patched {n} 处 + syntax OK")
# 复查残留
import re
left = [l for l in s.split("\n") if "node_id_of[c]" in l or "node_id_of[z]" in l]
print("残留 node_id_of[...] 行:", left)
