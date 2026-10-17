"""补：tibia/fibula 横向约束 BC 也要守卫空 NodeSet。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
L = p.read_text(encoding="utf-8").split("\n")

i0 = next(i for i, l in enumerate(L) if "_lateral" in l and "for ns, nm in"
          in L[i - 1] if i)
# 找到 for 循环起点
i0 = next(i for i, l in enumerate(L)
          if 'for ns, nm in (("tibia_top_ns"' in l))
# 循环体到下一个顶层语句
j = i0 + 1
while j < len(L) and (L[j].startswith("    ") or not L[j].strip()):
    j += 1
body = ["    " + x if x.strip() else x for x in L[i0:j]]
L[i0:j] = ['for ns, nm, _var in (("tibia_top_ns", "tibia", "top_tibia"),',
           '                    ("fibula_top_ns", "fibula", "top_fibula")):',
           '    if not _var:',
           '        continue'] + body[1:]

p.write_text("\n".join(L), encoding="utf-8")
s = p.read_text(encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
print("\n".join(s.split("\n")[i0 - 2:i0 + 12]))
