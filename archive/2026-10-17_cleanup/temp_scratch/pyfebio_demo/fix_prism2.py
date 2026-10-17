"""prism 模式改为：四边形底也拆成 2 个 penta6（避免歪斜 hex8 被剔除）。"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = """        if CART_ELEM == "prism" and m in (3, 4):
            # ★ 直接用楔形/六面体：薄片拆 tet 必为细长 sliver ⇒ 一压就翻负
            #   penta6 = [底三角 0,1,2 | 顶三角 3,4,5]；hex8 = [底 0-3 | 顶 4-7]
            et = "penta6" if m == 3 else "hex8"
            mesh_elems.append((j, et, "mat_cartilage",
                               [full[z] for z in range(2 * m)]))
        else:"""
assert s.count(old) == 1, s.count(old)
s = s.replace(old, """        if CART_ELEM == "prism" and m in (3, 4):
            # ★ 全用 penta6（6 节点楔形）：薄片拆 tet 必为细长 sliver ⇒ 一压就翻负。
            #   四边形底再拆成 2 个三角 —— 实测直接用 hex8 会因偏移方向与面不垂直
            #   而歪斜（det(J)≈0.5），被 DET_MIN 闸门剔除 42/61 个 ⇒ 软骨破洞、
            #   接触面从 150/197 掉到 51/64。拆楔形后 0 剔除。
            tris = ((0, 1, 2),) if m == 3 else ((0, 1, 2), (0, 2, 3))
            for tri in tris:
                conn = ([full[z] for z in tri]
                        + [full[z + m] for z in tri])
                mesh_elems.append((j, "penta6", "mat_cartilage", conn))
        else:""")
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched，syntax OK")
