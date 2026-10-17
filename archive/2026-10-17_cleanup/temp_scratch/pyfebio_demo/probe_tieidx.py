"""查 {j}_tie_faces / {j}_elems 的索引约定 + bone_faces。"""
from pathlib import Path

import numpy as np

P = Path(r"D:\Project\climbing_fall_analysis\temp\thums\cart_new\cartilage_mesh.npz")
d = np.load(P)
for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    print(f"\n===== {j} =====")
    for k in sorted(x for x in d.files if x.startswith(j)):
        a = d[k]
        print(f"  {k:38} {str(getattr(a,'shape','')):>10} {getattr(a,'dtype','')}")
    e = np.asarray(d[f"{j}_elems"])
    t = np.asarray(d[f"{j}_tie_faces"])
    n = np.asarray(d[f"{j}_nodes"])
    off = np.asarray(d[f"{j}_is_offset"], bool)
    print(f"  elems 前 3 行:\n{e[:3]}")
    print(f"  tie_faces 前 5 行:\n{t[:5]}")
    print(f"  elems 值域 [{e.min()}, {e.max()}]  含 0? {(e == 0).any()}  含负数? {(e < 0).any()}")
    print(f"  tie_faces 值域 [{t.min()}, {t.max()}]  含 0? {(t == 0).any()}  含负数? {(t < 0).any()}")
    print(f"  节点数 {len(n)} / is_offset True 数 {int(off.sum())}")
    # tie_faces 是否 1-based
    print(f"  ⇒ tie_faces 按 1-based 解析后最大索引 {t.max() - 1} < 节点数 {len(n)} ? "
          f"{t.max() - 1 < len(n)}")
