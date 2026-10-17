"""加 NO_CONTACT / LIG_SUB 开关（按行缩进包 Contact 块，避开文本匹配的坑）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
L = p.read_text(encoding="utf-8").split("\n")

# ① env 变量
for i, l in enumerate(L):
    if l.startswith('ONLY = [x for x in os.environ.get("BONES"'):
        L.insert(i + 1, 'NO_CONTACT = os.environ.get("NO_CONTACT", "0") == "1"\n'
                        'LIG_SUB = [int(x) for x in os.environ.get("LIG_SUB", "")'
                        '.split(",") if x]')
        break
for i, l in enumerate(L):
    if l.startswith('print(f"[mode] BONES='):
        L.insert(i + 1, 'print(f"[mode] NO_CONTACT={NO_CONTACT} '
                        'LIG_SUB={LIG_SUB or \'all\'}")')
        break

# ② 韧带子集
for i, l in enumerate(L):
    if l.strip() == 'for i in range(n_lig):':
        L.insert(i + 1, '        if LIG_SUB and (i + 1) not in LIG_SUB:')
        L.insert(i + 2, '            continue')
        break

# ③ Contact 块缩进（609-620 附近，按内容定位）
i0 = next(i for i, l in enumerate(L) if l.strip() == "w('  <Contact>')")
i1 = next(i for i, l in enumerate(L) if l.strip() == "w('  </Contact>')")
body = ["  " + x for x in L[i0:i1 + 1]]
L[i0:i1 + 1] = ["if not NO_CONTACT:"] + body

p.write_text("\n".join(L), encoding="utf-8")
s = p.read_text(encoding="utf-8")
ast.parse(s)
print(f"patched: Contact 块缩进 {len(body)} 行 + syntax OK")
