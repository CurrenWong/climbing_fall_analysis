"""A5: .xplt -> HDF5 -> 提取应力，与吴恺 3.74 MPa 量级对照。"""
from pathlib import Path
import numpy as np
import h5py
import pyfebio.xplt as X

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
src = W / "ankle_tibiotalar.xplt"
dst = W / "ankle_tibiotalar.h5"

X.to_hdf5(str(src), str(dst))
print(f"转换完成 {dst.stat().st_size/1024/1024:.1f} MB")


def walk(g, pre="", depth=0):
    for k in g:
        o = g[k]
        if isinstance(o, h5py.Group):
            print("  " * depth + f"[G] {pre}{k}")
            if depth < 2:
                walk(o, "", depth + 1)
        else:
            try:
                print("  " * depth + f"[D] {pre}{k}  shape={o.shape} dtype={o.dtype}")
            except Exception:
                print("  " * depth + f"[D] {pre}{k}")


print("\n=== HDF5 结构 ===")
with h5py.File(dst, "r") as f:
    walk(f)
