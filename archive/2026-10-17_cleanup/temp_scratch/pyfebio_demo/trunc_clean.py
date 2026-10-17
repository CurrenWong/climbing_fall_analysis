"""干净地重打 TRUNC 截断补丁（在 mesh_elems 层，节点编号之前），并自校验。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")
print("before: TRUNC_MM count =", s.count("TRUNC_MM"))

# 1) 环境变量
a1 = 'FIX_NAVICULAR = os.environ.get("FIX_NAVICULAR", "1") == "1"'
if "TRUNC_MM" not in s:
    assert s.count(a1) == 1, f"a1 {s.count(a1)}"
    s = s.replace(a1, a1 + '\nTRUNC_MM = float(os.environ.get("TRUNC_MM", "0"))')

a2 = 'print(f"[mode] DRIVER={DRIVER} DISP_MM={DISP_MM}")'
if 'TRUNC_MM={TRUNC_MM}' not in s:
    assert s.count(a2) == 1, f"a2 {s.count(a2)}"
    s = s.replace(a2, a2 + "\nprint(f'[mode] TRUNC_MM={TRUNC_MM}')")

# 2) 截断逻辑（mesh_elems 之后、编号之前）
a3 = 'print(f"骨+软骨：单元 {len(mesh_elems):,}，节点 {len(used):,}")'
if "trunc]" not in s:
    assert s.count(a3) == 1, f"a3 {s.count(a3)}"
    BLK = '''>TRUNC<
'''
    BLK = '''if TRUNC_MM > 0:
    _tal_pts = []
    for _c in mesh_elems:
        if _c[0] == "talus_cort":
            _tal_pts += [xyz_all[idx_all[x]] for x in _c[3]
                         if not isinstance(x, str)]
    _ref = np.mean(_tal_pts, axis=0) if _tal_pts else np.zeros(3)
    _new = []
    _drop = 0
    for _c in mesh_elems:
        _dom = _c[0]
        if _dom.startswith(("tibia_", "fibula_")):
            _pts = np.array([xyz_all[idx_all[x]] for x in _c[3]
                             if not isinstance(x, str)])
            if float((_pts.mean(0) - _ref) @ axis) > TRUNC_MM:
                _drop += 1
                continue
        _new.append(_c)
    mesh_elems = _new
    used = {x for _c in mesh_elems for x in _c[3]
            if not isinstance(x, str)}
    print(f"[trunc] 截断到 {TRUNC_MM} mm 内：剔除 {_drop:,} 单元，"
          f"剩 {len(mesh_elems):,} 单元 / {len(used):,} 节点")

'''
    s = s.replace(a3, BLK + a3)

p.write_text(s, encoding="utf-8")
ast.parse(s)
chk = p.read_text(encoding="utf-8")
print("after : TRUNC_MM count =", chk.count("TRUNC_MM"),
      "| trunc] count =", chk.count("trunc]"))
print("OK" if chk.count("TRUNC_MM") >= 3 and "trunc]" in chk else "FAILED")
