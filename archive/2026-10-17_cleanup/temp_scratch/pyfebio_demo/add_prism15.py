"""加 CART_ELEM=prism15：软骨用 penta15（二次楔形，同网格提高阶次）。

用途：p 收敛研究 —— 若 penta15 的支反力接近 penta6 而不是 tet4，
      说明 penta6 已接近收敛、tet4 偏刚。
实现：先按 prism 生成 penta6，再为 9 条边注册共享中点节点，升级为 penta15。
      （中点必须**全局共享**，否则单元间不连续。）
"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# ── A. penta15 的 Jacobian 代理（用 6 个角点，即 penta6 的采样检查）──
old = "FN = {"
assert s.count(old) == 1
s = s.replace(old, '''def detJ_penta15(P):
    """penta15 的角点落在 penta6 位置上 ⇒ 用其 6 角点做检查（代理）。"""
    return detJ_penta6(P[:6])


FN = {''')
i = s.index("FN = {")
j = s.index("}", i)
s = s[:j] + ', "penta15": detJ_penta15' + s[j:]

# ── B. 接触面表补 penta15（用角点面，与 penta6 同）──
old = 'CART_FACES = {'
assert s.count(old) == 1
s = s.replace(old, '''CART_FACES = {
    # penta15 用角点面即可：中点落在棱上，几何几乎一致（接触只需 tri3）
    "penta15": [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],''')

# ── C. 中点注册 + 升级为 penta15 ──
anchor = 'print(f"合并后节点 {len(coords):,}")'
assert s.count(anchor) == 1
s = s.replace(anchor, anchor + '''

# ── 软骨 p 升级：penta6 -> penta15（9 条边中点，**全局共享**）──
if CART_ELEM == "prism15":
    _id2xyz = {i + 1: np.asarray(coords[i], float)
               for i in range(len(coords))}

    def _xyz(key):
        return _id2xyz[node_id_of[key]]

    # VTK 二次楔形边序：6..14 = mid(0,1) mid(1,2) mid(2,0) mid(3,4) mid(4,5)
    #                             mid(5,3) mid(0,3) mid(1,4) mid(2,5)
    P15_EDGES = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5),
                 (5, 3), (0, 3), (1, 4), (2, 5)]
    _new, _nmid = [], 0
    for _dom, _et, _mat, _conn in mesh_elems:
        if _et == "penta6" and _mat == "mat_cartilage":
            _mids = []
            for _a, _b in P15_EDGES:
                _ka, _kb = _conn[_a], _conn[_b]
                _key = "mid:" + "|".join(sorted([str(_ka), str(_kb)]))
                if _key not in node_id_of:
                    node_id_of[_key] = len(coords) + 1
                    coords.append(tuple((_xyz(_ka) + _xyz(_kb)) / 2.0))
                    _id2xyz[node_id_of[_key]] = np.asarray(coords[-1], float)
                    _nmid += 1
                _mids.append(_key)
            _new.append((_dom, "penta15", _mat, list(_conn) + _mids))
        else:
            _new.append((_dom, _et, _mat, _conn))
    mesh_elems = _new
    XYZ = np.asarray(coords, float)
    print(f"[prism15] 新增边中点 {_nmid:,} 个；软骨改 penta15；"
          f"节点总数 {len(coords):,}")''')

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched，syntax OK")
