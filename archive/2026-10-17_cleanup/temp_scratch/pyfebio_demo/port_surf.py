"""用「只保留出现 1 次的面」替换 extreme_faces（真外表面 + 正确朝向）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

a1 = 'def extreme_faces(dom, which):'
i0 = s.index(a1)
i1 = s.index('top_tibia = extreme_faces("tibia_cort", "max")')
old = s[i0:i1]
assert "return sel" in old, old[-300:]

new = '''def surface_faces(dom, etype):
    """该域的**真外表面**：朝外定向后，只保留出现 1 次的面。

    ★ 这两点都不能省：
      ① 既遍历单元的**所有**面再按「出现次数==1」筛 —— 只看极值会把
         **内部面**当成载荷面（实测胫骨 2415 个六面体却选出 1551 个"顶面"）。
      ② 朝向必须朝外（法向背离单元质心），否则 contact/load 的正负号会反。
    """
    cnt, info = {}, {}
    for ids in blocks[(dom, etype)]:
        P = XYZ[[i - 1 for i in ids]]
        cen = P.mean(0)
        for f in (FACES_TET if etype == "tet4" else FACES_HEX):
            fv = [ids[z] for z in f]
            pts = XYZ[[i - 1 for i in fv]]
            n = np.zeros(3)
            for i in range(len(pts)):
                n += np.cross(pts[i], pts[(i + 1) % len(pts)])
            if np.dot(n, pts.mean(0) - cen) < 0:
                n = -n
                fv = fv[::-1]
            key = tuple(sorted(fv))
            if not facet_ok(fv):
                continue
            info[key] = tuple(fv)
            cnt[key] = cnt.get(key, 0) + 1
    return [info[k] for k, c in cnt.items() if c == 1]


def extreme_faces(dom, which, band=3.0):
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
    return [f for f, p in zip(fs, proj) if p < lim + band]


'''
s = s[:i0] + new + s[i1:]

# facet_ok 现在必须先于 surface_faces 定义 —— 提到 axis 之后
a2 = 'def surface_faces(dom, etype):'
assert s.count(a2) == 1
s = s.replace(a2, '''TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
FACES_TET = TET_FACES
MIN_FACET_AREA = 1e-8


def facet_ok(f):
    if len(set(f)) != len(f):
        return False
    P = XYZ[[i - 1 for i in f]]
    n = np.zeros(3)
    for i in range(len(P)):
        n += np.cross(P[i], P[(i + 1) % len(P)])
    return np.linalg.norm(n) / 2.0 >= MIN_FACET_AREA


def surface_faces(dom, etype):''')

# 去掉接触端口里重复的 facet_ok / TET_FACES 定义
a3 = '''TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
MIN_FACET_AREA = 1e-8


def facet_ok(f):
    if len(set(f)) != len(f):
        return False
    P = XYZ[[i - 1 for i in f]]
    n = np.zeros(3)
    for i in range(len(P)):
        n += np.cross(P[i], P[(i + 1) % len(P)])
    return np.linalg.norm(n) / 2.0 >= MIN_FACET_AREA


def cart_outer(j):'''
assert s.count(a3) == 1, f"cart block hit {s.count(a3)}"
s = s.replace(a3, 'def cart_outer(j):')

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
