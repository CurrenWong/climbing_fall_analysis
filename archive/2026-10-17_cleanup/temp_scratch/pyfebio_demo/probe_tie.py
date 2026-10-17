"""摸 ankle_lig_feb.py 里可复用的助手（决定 TIE 移植量）。"""
import io
import re
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
L = io.open(P, encoding="utf-8").read().split("\n")
print("=== 顶层的 def / 常量 ===")
for i, l in enumerate(L, 1):
    if re.match(r"^(def |[A-Z][A-Z0-9_]* *=)", l):
        print(f"{i:>4}|{l[:104]}")
print("\n=== 软骨加载段（125~205）关键行 ===")
for i, l in enumerate(L[124:205], 125):
    if re.search(r"cg\[|CARTS|conn|base|offs_|full|mesh_elems|CART_|HEX_T|WEDGE", l):
        print(f"{i:>4}|{l[:112]}")
print("\n=== mesh_elems 的构建方式 ===")
for i, l in enumerate(L, 1):
    if re.search(r"mesh_elems\.append|mesh_elems *=", l):
        print(f"{i:>4}|{l[:112]}")
print("\n=== 是否已有 facet / surface_faces / node_id_of ===")
for k in ("facet_ok", "surface_faces", "TET_FACES", "node_id_of", "_cname", "_cflag",
          "elem_blocks", "TIE"):
    hits = [i for i, l in enumerate(L, 1) if k in l]
    print(f"  {k:16} {len(hits):>3} 处  {hits[:6]}")
