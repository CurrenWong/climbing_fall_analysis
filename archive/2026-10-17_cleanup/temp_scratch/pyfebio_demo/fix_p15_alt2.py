"""把 penta15 的朝向修正改为「不修正、只剔除」。

理由：15 节点的奇置换不好构造 ——
  · 交换上下三角（角 3 对换 + 边 3 对换）= 6 次 = **偶** ⇒ 行列式不变号
  · 想再加一次对换又会把中点节点挪到别的棱上（几何就变了）
而 penta6 的构造本身已保证楔形朝向正确（实测剔除仅 14+29，都是"太扁"而非"翻反"）
⇒ penta15 继承其角点朝向，直接剔除即可。
"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

start = s.index('        if et == "tet4":')
end = s.index("        d2 = fn(XYZ[[i - 1 for i in alt]])")
new = '''        if et == "tet4":
            alt = [ids[1], ids[0]] + ids[2:]                # 1 次对换 = 奇
        elif et == "penta15":
            alt = None      # 见文件头说明：不修正，直接剔除
        else:
            alt = ids[len(ids) // 2:] + ids[:len(ids) // 2]  # swap halves = 奇
'''
s = s[:start] + new + s[end:]

old = "        d2 = fn(XYZ[[i - 1 for i in alt]])"
assert s.count(old) == 1
s = s.replace(old, "        d2 = (fn(XYZ[[i - 1 for i in alt]])\n"
                   "              if alt is not None else -1.0)")
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched，syntax OK")
