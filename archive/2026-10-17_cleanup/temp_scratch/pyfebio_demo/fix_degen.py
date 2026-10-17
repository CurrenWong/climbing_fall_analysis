"""修退化检查：兼容 tri3(3节点) / quad4(4节点)。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = '''        for ids in lig_blocks[nm]:
            P = XYZ[[i - 1 for i in ids]]
            d = np.cross(P[2] - P[0], P[3] - P[1])
            if np.linalg.norm(d) > 1e-6:'''
assert s.count(old) == 1, f"hit {s.count(old)}"
s = s.replace(old, '''        for ids in lig_blocks[nm]:
            P = XYZ[[i - 1 for i in ids]]
            if len(P) == 3:
                d = np.cross(P[1] - P[0], P[2] - P[0])
            else:
                d = np.cross(P[2] - P[0], P[3] - P[1])
            if np.linalg.norm(d) > 1e-6:''')

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
