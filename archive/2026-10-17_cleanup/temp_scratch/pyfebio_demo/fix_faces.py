"""修两处：① 韧写出 guard；② 顶/底面按**法向**选（斜轴下投影带会切偏）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# ① guard
old = """LIG_EID0 = eid
if LIG_MODE == "shell":
    for nm, conns in lig_blocks.items():"""
assert s.count(old) == 1
s = s.replace(old, """LIG_EID0 = eid
if LIG_MODE in ("shell", "truss"):
    for nm, conns in lig_blocks.items():""")

# ② 按法向选端面
old2 = '''def extreme_faces(dom, which, band=3.0):
    """该域外表面中沿加载轴最高/最低的一批（阈值同原模型：3 mm 带内）。"""
    fs = []
    for et in ("hex8", "tet4"):
        if (dom, et) in blocks:
            fs += surface_faces(dom, et)
    if not fs:
        return []
    proj = XYZ[[f[0] - 1 for f in fs]] @ axis
    lim = proj.max() if which == "max" else proj.min()
    if which == "max":
        return [f for f, p in zip(fs, proj) if p > lim - band]
    return [f for f, p in zip(fs, proj) if p < lim + band]'''
assert s.count(old2) == 1, f"extreme hit {s.count(old2)}"
new2 = '''def _face_normal(f):
    P = XYZ[[i - 1 for i in f]]
    n = np.zeros(3)
    for i in range(len(P)):
        n += np.cross(P[i], P[(i + 1) % len(P)])
    return n / 2.0


def extreme_faces(dom, which, band=25.0, cos_min=0.6):
    """该域**朝向加载轴**的端面（顶 or 底）。

    ★ 不能只用「沿轴投影带」挑：加载轴是斜的时候（THUMS 里胫骨轴 ≈
      [0.68,0.14,0.72]），水平截面的投影会散开好几 mm ⇒ 3mm 带只切到平台一角
      （实测 56.8 mm²，真值 ~1000 mm²）。**按面法向与轴对齐**才是几何上正确的判据。
    """
    fs = []
    for et in ("hex8", "tet4"):
        if (dom, et) in blocks:
            fs += surface_faces(dom, et)
    if not fs:
        return []
    sgn = 1.0 if which == "max" else -1.0
    aligned = [f for f in fs if float(_face_normal(f) @ axis) * sgn > cos_min]
    if not aligned:
        return []
    proj = XYZ[[f[0] - 1 for f in aligned]] @ axis
    lim = proj.max() if which == "max" else proj.min()
    if which == "max":
        return [f for f, p in zip(aligned, proj) if p > lim - band]
    return [f for f, p in zip(aligned, proj) if p < lim + band]'''
s = s.replace(old2, new2)

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
