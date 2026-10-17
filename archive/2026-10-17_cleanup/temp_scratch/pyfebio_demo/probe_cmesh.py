"""只读探针：cartilage_mesh.npz 结构 + 编码规则确认。"""
import numpy as np
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\temp\thums\cartilage_mesh.npz")
d = np.load(P, allow_pickle=True)
print("keys:", d.files)
for k in d.files:
    a = d[k]
    print(f"  {k}: shape={getattr(a,'shape',None)} dtype={a.dtype}")
    if a.dtype == object or a.size and a.ndim <= 1 and a.size:
        try:
            v = a[0] if a.ndim == 1 else a
            print(f"      样例: {str(v)[:150]}")
        except Exception as e:
            print("      ", e)

print("\n--- 若有 conn/cells，看编码 ---")
for key in ("conn", "cells", "elems", "cart_conn"):
    if key in d.files:
        c = d[key]
        print(f"{key}: {c.shape}")
        print(c[:4])
        break
