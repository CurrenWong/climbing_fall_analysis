"""修：坐标是结构化数组 (id,x,y,z)，取 x/y/z 拼成 (n,3)。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\read_lig.py")
s = p.read_text(encoding="utf-8")

old = "X0 = np.asarray(coords)          # (n, 3)"
assert s.count(old) == 1
s = s.replace(old, """_c = np.asarray(coords)
if _c.dtype.names:                       # 结构化数组 (id, x, y, z)
    X0 = np.column_stack([_c["x"], _c["y"], _c["z"]]).astype(float)
else:
    X0 = _c.reshape(-1, 3).astype(float)  # (n, 3)""")

# 位移/支反力可能是 (n,3) 或 (steps,n,3)
old2 = """D = np.asarray(disp)
if D.ndim == 3:
    D = D[-1]
D = np.atleast_2d(D)
if D.shape[0] < coords.shape[0]:
    D = D.T if D.shape[1] >= coords.shape[0] else D
    D = np.atleast_2d(D).reshape(-1, 3) if D.size == coords.size else D
print(f"位移矩阵形状 {D.shape}")"""
assert s.count(old2) == 1
s = s.replace(old2, """D = np.asarray(disp, dtype=float)
if D.ndim == 3:
    D = D[-1]                       # 末帧
D = D.reshape(-1, 3)
print(f"位移矩阵形状 {D.shape}")""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
