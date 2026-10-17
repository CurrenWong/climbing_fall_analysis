"""补守卫：tibia/fibula 横向 BC 跳过空 NodeSet（按行号改，简单直接）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
L = p.read_text(encoding="utf-8").split("\n")

i = L.index('for ns, nm in (("tibia_top_ns", "tibia"), ("fibula_top_ns", "fibula")):')
assert i > 0
L[i:i + 6] = [
    'for ns, nm, _var in (("tibia_top_ns", "tibia", "top_tibia"),',
    '                    ("fibula_top_ns", "fibula", "top_fibula")):',
    '    if not _var:',
    '        continue',
    '    w(f\'    <bc name="{nm}_lateral" type="zero displacement" '
    'node_set="{ns}">\')',
    '    for di, dn in enumerate("xyz"):',
    '        if di != kax:',
    '            w(f\'      <{dn}_dof>1</{dn}_dof>\')',
    '    w(\'    </bc>\')',
]
p.write_text("\n".join(L), encoding="utf-8")
s = p.read_text(encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
