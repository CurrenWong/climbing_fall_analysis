"""移动 LIG_ONLY_NS 计算到 lig_blocks 构建之后。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

i0 = s.index("# 韧带专属节点（仅被壳用")
i1 = s.index("lig_blocks = {}\nLIG_ETYPE =")
block = s[i0:i1]
s = s[:i0] + s[i1:]

# 找到 lig_blocks 构建结束的位置（下一个空行后的注释/代码）
anchor = '            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]\n                                     for t in r["tris"]]'
assert s.count(anchor) == 1, s.count(anchor)
s = s.replace(anchor, anchor + "\n\n" + block.rstrip() + "\n")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("moved + syntax OK")
