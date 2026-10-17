"""修 web 生成的越界：某侧节点 < 2 时跳过。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = '''        if Aset and Bset:
            PA = xyz_all[[idx_all[a] for a in Aset]]'''
assert s.count(old) == 1
s = s.replace(old, '''        if Aset and len(Bset) >= 2:
            PA = xyz_all[[idx_all[a] for a in Aset]]''')

old2 = '''            Dt = D.T
            for j, b in enumerate(Bset):'''
assert s.count(old2) == 1
s = s.replace(old2, '''        if Bset and len(Aset) >= 2:
            if not (Aset and len(Bset) >= 2):
                PB = xyz_all[[idx_all[b] for b in Bset]]
                PA = xyz_all[[idx_all[a] for a in Aset]]
                D = np.linalg.norm(PA[:, None, :] - PB[None, :, :], axis=2)
            Dt = D.T
            for j, b in enumerate(Bset):''')

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
