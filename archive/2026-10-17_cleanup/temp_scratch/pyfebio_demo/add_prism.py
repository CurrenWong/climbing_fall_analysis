"""加 CART_ELEM 开关：软骨直接用 penta6/hex8，不拆成 sliver tet4。

背景：薄软骨(0.703 mm)是 6 节点楔形，现有代码用 WEDGE_T 把它拆成 3 个 tet4
⇒ 3 个必为细长 sliver ⇒ 压缩时 det(F) 掉到 0.52 就翻负（实测 t=0.63 失败）。
FEBio 原生支持 penta6（6 节点楔形），直接用即可。
"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")
n_ok = 0

# ── A. 环境开关 ──
old = 'TOPFIX = os.environ.get("TOPFIX", "lateral").strip() or "lateral"'
assert s.count(old) == 1
s = s.replace(old, old + '''
# 软骨单元类型：tet4（旧，拆楔形）/ prism（penta6+hex8，抗压缩好得多）
CART_ELEM = os.environ.get("CART_ELEM", "tet4").strip() or "tet4"''')
n_ok += 1

# ── B. 单元生成 ──
old = """        full = base + offs_
        for t in (HEX_T if m == 4 else WEDGE_T):
            mesh_elems.append((j, "tet4", "mat_cartilage", [full[z] for z in t]))
        used.update(base)"""
assert s.count(old) == 1, s.count(old)
s = s.replace(old, """        full = base + offs_
        if CART_ELEM == "prism" and m in (3, 4):
            # ★ 直接用楔形/六面体：薄片拆 tet 必为细长 sliver ⇒ 一压就翻负
            #   penta6 = [底三角 0,1,2 | 顶三角 3,4,5]；hex8 = [底 0-3 | 顶 4-7]
            et = "penta6" if m == 3 else "hex8"
            mesh_elems.append((j, et, "mat_cartilage",
                               [full[z] for z in range(2 * m)]))
        else:
            for t in (HEX_T if m == 4 else WEDGE_T):
                mesh_elems.append((j, "tet4", "mat_cartilage",
                                   [full[z] for z in t]))
        used.update(base)""")
n_ok += 1

# ── C. detJ_penta6 ──
old = """def detJ_tet4(P):
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0],
                                         P[3] - P[0]])))


FN = {"hex8": detJ_hex8, "tet4": detJ_tet4}"""
assert s.count(old) == 1
s = s.replace(old, """def detJ_tet4(P):
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0],
                                         P[3] - P[0]])))


def detJ_penta6(P):
    \"\"\"楔形：底三角 0,1,2 + 顶三角 3,4,5；用 P0 出发的三条边。
    符号 = 底三角绕向 x (P3-P0) 方向；为负时按奇置换翻转。\"\"\"
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0],
                                         P[3] - P[0]])))


FN = {"hex8": detJ_hex8, "tet4": detJ_tet4, "penta6": detJ_penta6}""")
n_ok += 1

# ── D. cart_outer：支持 penta6 / hex8 ──
old = """def cart_outer(j):
    out, seen = [], set()
    for ids in blocks[(j, "tet4")]:
        for f in TET_FACES:
            fv = [ids[z] for z in f]
            if not all(x in OFFSET_IDS for x in fv):
                continue
            key = tuple(sorted(fv))
            if key in seen or not facet_ok(fv):
                continue
            seen.add(key)
            opp = ids[[z for z in range(4) if z not in f][0]]
            pts = XYZ[[i - 1 for i in fv]]
            n = np.zeros(3)
            for i in range(3):
                n += np.cross(pts[i], pts[(i + 1) % 3])
            if np.dot(n, pts.mean(0) - XYZ[opp - 1]) < 0:
                fv = fv[::-1]
            out.append(fv)
    return out"""
assert s.count(old) == 1
s = s.replace(old, '''CART_FACES = {
    "tet4": TET_FACES,
    "penta6": [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    "hex8": FACES_HEX,
}


def cart_outer(j):
    """软骨**外表面**（全部由 offset 节点构成的面），朝外定向。

    ★ 支持 tet4 / penta6 / hex8 —— 换单元类型后必须同步换面表，
      否则接触面与网格真实面不重合 ⇒ 接触完全不检测（rc=0 无警告、应力全 0）。
    """
    out, seen = [], set()
    for et in ("tet4", "penta6", "hex8"):
        for ids in blocks.get((j, et), []):
            for f in CART_FACES[et]:
                fv = [ids[z] for z in f]
                if not all(x in OFFSET_IDS for x in fv):
                    continue
                key = tuple(sorted(fv))
                if key in seen or not facet_ok(fv):
                    continue
                seen.add(key)
                rest = [ids[z] for z in range(len(ids)) if z not in f]
                cen = XYZ[[i - 1 for i in rest]].mean(0)
                pts = XYZ[[i - 1 for i in fv]]
                n = np.zeros(3)
                for i in range(len(pts)):
                    n += np.cross(pts[i], pts[(i + 1) % len(pts)])
                if np.dot(n, pts.mean(0) - cen) < 0:
                    fv = fv[::-1]
                out.append(fv)
    return out''')
n_ok += 1

p.write_text(s, encoding="utf-8")
ast.parse(s)
print(f"patched {n_ok}/4 处，syntax OK")
