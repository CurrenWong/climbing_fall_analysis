"""修：截断块内就地算轴（原 axis 在块之后才定义）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = """    _ref = np.mean(_tal_pts, axis=0) if _tal_pts else np.zeros(3)"""
assert s.count(old) == 1
s = s.replace(old, """    _ref = np.mean(_tal_pts, axis=0) if _tal_pts else np.zeros(3)
    _tib_pts = []
    for _c in mesh_elems:
        if _c[0] == "tibia_cort":
            _tib_pts += [xyz_all[idx_all[x]] for x in _c[3]
                         if not isinstance(x, str)]
    _axis = np.mean(_tib_pts, axis=0) - _ref if _tib_pts else np.array(
        [0.0, 0.0, 1.0])
    _axis = _axis / np.linalg.norm(_axis)""")

old2 = "            if float((_pts.mean(0) - _ref) @ axis) > TRUNC_MM:"
assert s.count(old2) == 1
s = s.replace(old2, "            if float((_pts.mean(0) - _ref) @ _axis) > TRUNC_MM:")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
