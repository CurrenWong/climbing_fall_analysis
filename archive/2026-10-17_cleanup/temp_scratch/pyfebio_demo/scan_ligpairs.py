"""看 lig_read npz 的键与 lig_pairs 的编号空间。"""
import numpy as np
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
for nm in ("lig_read_WK.npz",):
    p = W / nm
    if not p.exists():
        print(f"✗ 无 {p}")
        continue
    z = np.load(p, allow_pickle=True)
    print(f"=== {p.name} 键 {list(z.keys())}")
    for k in z.keys():
        a = z[k]
        print(f"  {k:<18} shape={getattr(a,'shape',None)} dtype={a.dtype}")
    tr = z["lig_pairs"]
    for i in range(min(3, len(tr))):
        print(f"  lig_pairs[{i}] = {tr[i][:6]}")
    print(f"  n_truss = {z['lig_npair']}")
# ligaments.npz 对照
q = Path(r"D:\Project\climbing_fall_analysis\temp\thums\ligaments.npz")
z2 = np.load(q, allow_pickle=True)
print(f"\n=== ligaments.npz 名称 = {list(z2['names'])}")
print(f"  truss[0][:5] = {z2['truss'][0][:5]}   (THUMS id)")
print(f"  length = {np.round(z2['length'],2)}")
