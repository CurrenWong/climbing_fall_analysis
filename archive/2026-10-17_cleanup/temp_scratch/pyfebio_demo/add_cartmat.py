"""让软骨材料 E/ν 可由环境变量调（锁死判定实验用）。"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")
n = 0

a = 'CART_ELEM = os.environ.get("CART_ELEM", "tet4").strip() or "tet4"'
assert s.count(a) == 1
s = s.replace(a, a + '''
# 软骨材料（锁死判定：线性 tet4 的体积锁死随 ν→0.5 放大）
CART_E = float(os.environ.get("CART_E", "10.0"))
CART_NU = float(os.environ.get("CART_NU", "0.40"))''')
n += 1

old1 = '"mat_cartilage": ("isotropic elastic", {"density": 1.1e-9, "E": 10.0, "v": 0.40}),'
old2 = '"mat_cartilage": ("isotropic elastic", {"density": 1.1e-9, "E": 260.0, "v": 0.40}),'
for old in (old1, old2):
    if old in s:
        new = old.replace('"E": 10.0', '"E": CART_E') \
                 .replace('"E": 260.0', '"E": CART_E') \
                 .replace('"v": 0.40', '"v": CART_NU')
        s = s.replace(old, new)
        n += 1
        print("  patch", old[:60], "->", new[:60])

assert n == 3, n
p.write_text(s, encoding="utf-8")
ast.parse(s)
print(f"patched {n} 处，syntax OK")
