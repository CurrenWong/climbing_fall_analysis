"""查 blocks 的构建方式（判断新格式能否用字符串键）。"""
import io
import re
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
L = io.open(P, encoding="utf-8").read().split("\n")
print("=== blocks / node_id_of 解析相关 ===")
for i, l in enumerate(L, 1):
    if re.search(r"blocks\s*[=\[]|blocks\.setdefault|node_id_of\[|isinstance\(.*str|_resolve|idx_all|xyz_all", l):
        print(f"{i:>4}|{l[:112]}")
print("\n=== mesh_elems 被消费的位置 ===")
for i, l in enumerate(L, 1):
    if re.search(r"for .* in mesh_elems|mesh_elems\b", l):
        print(f"{i:>4}|{l[:112]}")
