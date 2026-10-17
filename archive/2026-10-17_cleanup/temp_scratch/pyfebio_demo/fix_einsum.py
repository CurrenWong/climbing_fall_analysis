import io
import ast

p = r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\pt_tri.py"
s = io.open(p, encoding="utf-8").read()
n = 0
for v in "ab ac ab ac ab ac".split():
    pass
for expr, o in (("ab, ap", "d1"), ("ac, ap", "d2"), ("ab, bp", "d3"),
                ("ac, bp", "d4"), ("ab, cp", "d5"), ("ac, cp", "d6")):
    old = f'np.einsum("ij,j->i", {expr})'
    new = f'np.einsum("ij,ij->i", {expr})'
    if old in s:
        s = s.replace(old, new)
        n += 1
        print(f"  [{o}] {old} -> {new}")
io.open(p, "w", encoding="utf-8").write(s)
ast.parse(s)
print(f"修了 {n} 处，syntax OK")
