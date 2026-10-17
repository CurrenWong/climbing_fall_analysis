"""直接列出 npz 每个 key 的节点数分布（别再猜）。"""
from collections import Counter
from pathlib import Path

import numpy as np

for f in ("thums_lowerlimb.npz", "cartilage_mesh_t0d.npz"):
    d = np.load(Path(r"D:\Project\climbing_fall_analysis\temp\thums") / f,
                allow_pickle=True)
    print(f"\n=== {f} ===")
    for k in d.files:
        a = d[k]
        if a.ndim == 2 and a.dtype.kind in "iu" and k.endswith(("solid", "shell")):
            c = Counter(len([x for x in row if int(x) != 0]) for row in a)
            print(f"  {k:26s} {str(a.shape):>14s}  节点数分布 {dict(sorted(c.items()))}")
        elif k.endswith(("elems",)):
            c = Counter(len([x for x in row if int(x) != 0]) for row in a)
            print(f"  {k:26s} {str(a.shape):>14s}  宽度 {dict(sorted(c.items()))}")
