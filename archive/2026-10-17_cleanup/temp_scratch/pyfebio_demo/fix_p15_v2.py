"""修 prism15 的两个 bug（按文件实际内容）：
① cart_outer 的单元类型遍历漏 penta15 ⇒ 接触面 0 面
   ⇒ FEBio "Missing contact primary surface"
② 朝向修正的奇置换对 15 节点无效 ⇒ 多剔 36 个。
   penta15 直接跳过修正（penta6 构造已保证朝向）。
"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# ---- ① ----
old = '    for et in ("tet4", "penta6", "hex8"):'
assert s.count(old) == 1, ("iter list", s.count(old))
s = s.replace(old, '    for et in ("tet4", "penta6", "penta15", "hex8"):')

# ---- ② ----
old2 = """        alt = ([ids[1], ids[0]] + ids[2:]) if et == "tet4" else \\
              (ids[len(ids) // 2:] + ids[:len(ids) // 2])
        d2 = fn(XYZ[[i - 1 for i in alt]])"""
assert s.count(old2) == 1, ("alt block", s.count(old2))
s = s.replace(old2, """        if et == "tet4":
            alt = [ids[1], ids[0]] + ids[2:]                 # 1 次对换 = 奇
        elif et == "penta15":
            # ★ 15 节点的奇置换不好构造：
            #   交换上下三角 = 角 3 对换 + 边 3 对换 = 6 次 = **偶**，符号不变；
            #   再加对换又会把中点节点挪到别的棱上（几何就变了）。
            #   penta6 的构造本身已保证楔形朝向正确 ⇒ 继承角点朝向，直接剔除。
            alt = None
        else:
            alt = ids[len(ids) // 2:] + ids[:len(ids) // 2]   # swap halves = 奇
        d2 = (fn(XYZ[[i - 1 for i in alt]])
              if alt is not None else -1.0)""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched 2 处，syntax OK")
