"""build_feb.py -> 接触面/骨面过滤退化 facet。
FEBio 会警告 'The surface X has N invalid facets. The model may not run
correctly.' —— 面积接近 0 或节点重复的三角形会让 penalty 接触在 t=0
就产生异常力（实测 2 关节模型第一个载荷步就压翻单元）。"""
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\build_feb.py")
t = P.read_text(encoding="utf-8")

# ① 新增退化 facet 判定
anchor = "TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]\n"
assert t.count(anchor) == 1
t = t.replace(anchor, anchor + '''MIN_FACET_AREA = 1e-8   # mm^2: facets smaller than this are degenerate


def facet_ok(f):
    """Facet must have 3 distinct nodes and non-zero area.
    FEBio warns 'The surface X has N invalid facets' for these, and a penalty
    contact built on them blows up at t=0."""
    if len(set(f)) != len(f):
        return False
    p = XYZ[[i - 1 for i in f]]
    return float(np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))) * 0.5 \\
        > MIN_FACET_AREA

''')

# ② cart_outer 里过滤
old = """            if not all(x in OFFSET_IDS for x in fv):
                continue
            key = tuple(sorted(fv))
            if key in seen:
                continue
            seen.add(key)"""
assert t.count(old) == 1, f"cart_outer x{t.count(old)}"
t = t.replace(old, """            if not all(x in OFFSET_IDS for x in fv):
                continue
            key = tuple(sorted(fv))
            if key in seen:
                continue
            if not facet_ok(fv):
                continue
            seen.add(key)""")

# ③ surface_faces 里过滤
old2 = "            key = tuple(sorted(fv))\n            info[key] = tuple(fv)\n            cnt[key] += 1"
assert t.count(old2) in (0,1), f"surface_faces x{t.count(old2)}"
t = t.replace(old2, """            key = tuple(sorted(fv))
            if not facet_ok(fv):
                continue
            info[key] = tuple(fv)
            cnt[key] += 1""")

P.write_text(t, encoding="utf-8")
import ast
ast.parse(t)
print("OK: 退化 facet 过滤已加入（cart_outer + surface_faces）")
