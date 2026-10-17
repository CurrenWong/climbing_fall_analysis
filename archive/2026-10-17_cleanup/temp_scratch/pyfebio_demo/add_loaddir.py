"""加 LOADDIR=z：沿全局 Z 加载（验证"垂直"定义对排序的影响）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = 'TRUNC_MM = float(os.environ.get("TRUNC_MM", "0"))'
assert s.count(old) == 1
s = s.replace(old, old + '''
# LOADDIR=axis（默认，沿胫骨长轴）｜z（沿全局 Z）—— 验证"垂直"的定义
LOADDIR = os.environ.get("LOADDIR", "axis").strip() or "axis"''')

# 在 axis 定义之后覆写（保留 axis 做几何选面，用 LAX 做加载方向）
old2 = "print(f\"加载轴（胫骨->距骨）= {np.round(axis,4)}\")"
assert s.count(old2) == 1
s = s.replace(old2, old2 + """\nLAX = axis.copy()
if LOADDIR == "z":
    _z = np.array([0.0, 0.0, 1.0])
    LAX = _z if float(axis @ _z) > 0 else -_z
    print(f"[loaddir] 沿全局 Z 加载，与胫骨轴夹角 "
          f"{np.degrees(np.arccos(abs(float(axis @ _z)))):.1f} deg")""")

# 加载/位移用 LAX
s = s.replace("""    tr = F / a if a > 0 else 0.0""", """    tr = F / a if a > 0 else 0.0""")
s = s.replace("""    w(f'      <traction>{-axis[0]*tr:.6f},{-axis[1]*tr:.6f},{-axis[2]*tr:.6f}'""",
              """    w(f'      <traction>{-LAX[0]*tr:.6f},{-LAX[1]*tr:.6f},{-LAX[2]*tr:.6f}'""")
s = s.replace("    _kax = int(np.argmax(np.abs(axis)))",
              "    _kax = int(np.argmax(np.abs(LAX)))")
s = s.replace("            _v = (-DISP_MM if axis[_di] > 0 else DISP_MM) if _di == _kax else 0.0",
              "            _v = (-DISP_MM if LAX[_di] > 0 else DISP_MM) if _di == _kax else 0.0")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
print("LAX uses:", s.count("LAX["))
