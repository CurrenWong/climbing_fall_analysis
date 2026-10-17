"""detJ_penta6 改为三点检查（t=0 / 0.5 / 1 的形心），取最小值。

楔形若顶三角与底三角绕向相反（扭转），底面形心 det>0 但顶面 det<0 ——
FEBio 会在 mesh initialization 报 "Negative jacobian"。
"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = '''    dt = (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3.0
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0], dt])))'''
assert s.count(old) == 1
s = s.replace(old, '''    dt = (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3.0
    Jb = np.stack([P[1] - P[0], P[2] - P[0], dt])            # t = 0
    Jt = np.stack([P[4] - P[3], P[5] - P[3], dt])            # t = 1
    Jm = np.stack([(P[1] - P[0] + P[4] - P[3]) / 2.0,        # t = 0.5
                   (P[2] - P[0] + P[5] - P[3]) / 2.0, dt])
    return float(min(np.linalg.det(Jb), np.linalg.det(Jm),
                     np.linalg.det(Jt)))''')
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched，syntax OK")
