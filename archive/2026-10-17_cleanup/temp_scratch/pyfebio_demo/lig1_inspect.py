"""L1: 摸清 npz 里 8 组踝韧带的实际形态。"""
import numpy as np
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\temp\thums\thums_lowerlimb.npz")
z = np.load(P, allow_pickle=True)

print("=== npz 全部 key ===")
ks = list(z.keys())
print(f"共 {len(ks)} 个")
for k in ks:
    a = z[k]
    print(f"  {k:42s} {str(a.shape):>16s}  {a.dtype}")

print("\n=== 含 LIG / 81xx 的 key ===")
lig = [k for k in ks if "LIG" in k.upper() or k.startswith("811")]
for k in lig:
    a = z[k]
    print(f"  {k:42s} {str(a.shape):>16s}  {a.dtype}")
    if a.ndim == 1 and a.size <= 40:
        print(f"      {a}")
