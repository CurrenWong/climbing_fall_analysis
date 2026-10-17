"""把 TRUNC 截断移到 node 编号之前的 mesh_elems 层（真正删除节点）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# ① 删掉上一版插在 blocks 之后的截断块
i0 = s.index("# ── 胫/腓骨截断 ──")
i1 = s.index("lig_blocks = {}\nLIG_ETYPE =")
s = s[:i0] + s[i1:]

# ② 插到 mesh_elems 构建完之后（软骨循环之后、打印之前）
anchor = 'print(f"骨+软骨：单元 {len(mesh_elems):,}，节点 {len(used):,}")'
assert s.count(anchor) == 1
NEW = '''# ── 胫/腓骨截断（★ 必须在节点编号之前：否则被删单元会留下孤立节点）──
#   理由：① 速度——整根胫骨带 25,000+ 单元，每轮求解 32k 方程要 20 min，
#         没法迭代；② 物理——**踝关节模型不需要整根胫骨**，吴恺的模型就是
#         "足踝部"，载荷加在截面上正对应"胫腓骨上截面 600 N"。
if TRUNC_MM > 0:
    _tal_pts = []
    for _c in mesh_elems:
        if _c[0] == "talus_cort":
            _tal_pts += [xyz_all[idx_all[x]] for x in _c[3] if not
                         isinstance(x, str)]
    _ref = np.mean(_tal_pts, axis=0) if _tal_pts else np.zeros(3)
    _new, _drop = [], 0
    for _c in mesh_elems:
        _dom, _et, _mat, _conn = _c
        if _dom.startswith(("tibia_", "fibula_")):
            _pts = np.array([xyz_all[idx_all[x]] for x in _conn])
            if float((_pts.mean(0) - _ref) @ axis) > TRUNC_MM:
                _drop += 1
                continue
        _new.append(_c)
    mesh_elems = _new
    used = {x for _c in mesh_elems for x in _c[3] if not isinstance(x, str)}
    for _j in CARTS:
        pass
    print(f"[trunc] 截断到 {TRUNC_MM} mm 内：剔除外侧 {_drop:,} 单元，"
          f"剩 {len(mesh_elems):,} 单元 / {len(used):,} 节点")

''' + anchor
s = s.replace(anchor, NEW)

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
