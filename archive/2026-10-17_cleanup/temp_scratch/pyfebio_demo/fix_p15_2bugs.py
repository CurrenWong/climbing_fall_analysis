"""修 prism15 的两个 bug：
① cart_outer 的单元类型遍历列表漏了 penta15 ⇒ 接触面 0 面 ⇒
   FEBio "Missing contact primary surface"
② 朝向修正的奇置换对 15 节点无效（ids[7:]+ids[:7] 会打乱）⇒ 多剔 36 个。
   penta15 正解：交换上下三角（角 0-2 ↔ 3-5，边 6-8 ↔ 9-11，12-14 不动）
   —— 等价于 3 次对换 = 奇置换 ⇒ 行列式变号。
"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = '    for et in ("tet4", "penta6", "hex8"):'
assert s.count(old) == 1, s.count(old)
s = s.replace(old, '    for et in ("tet4", "penta6", "penta15", "hex8"):')

old2 = """        alt = ([ids[1], ids[0]] + ids[2:]) if et == "tet4" else \\
              (ids[len(ids) // 2:] + ids[:len(ids) // 2])"""
assert s.count(old2) == 1, s.count(old2)
s = s.replace(old2, """        if et == "tet4":
            alt = [ids[1], ids[0]] + ids[2:]                # 1 次对换
        elif et == "penta15":
            # 交换上下三角：角 0-2 <-> 3-5；边 6-8 <-> 9-11；边 12-14 不动
            # = (0 3)(1 4)(2 5)(6 9)(7 10)(8 11)，3+3 = 6 次对换 —— 偶！
            # ★ 必须取**奇**置换。实测：直接用上面的 6 对换是偶置换、号不变。
            #   正解 = 再加上一次对换（交换前两个角），共 7 次 = 奇。
            alt = [ids[4], ids[3], ids[5], ids[1], ids[0], ids[2],
                   ids[9], ids[10], ids[11], ids[6], ids[7], ids[8],
                   ids[12], ids[13], ids[14]]
        else:
            alt = ids[len(ids) // 2:] + ids[:len(ids) // 2]""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched，syntax OK")
