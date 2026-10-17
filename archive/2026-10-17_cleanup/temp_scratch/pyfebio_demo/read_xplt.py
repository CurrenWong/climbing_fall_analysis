"""用 pyfebio 读 .xplt，核验刚体位移 = F/k。"""
import sys
from pathlib import Path

import numpy as np

X = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\rt.xplt")
try:
    import pyfebio
    from pyfebio import xplt
    print("[ok] pyfebio", getattr(pyfebio, "__version__", "?"))
except Exception as e:
    print("[no pyfebio]", e)
    sys.exit(1)

h5 = X.with_suffix(".h5")
try:
    xplt.to_hdf5(str(X), str(h5))
    print("[ok] xplt -> h5:", h5)
except Exception as e:
    print("[to_hdf5 failed]", type(e).__name__, e)

import h5py
with h5py.File(h5, "r") as f:
    def walk(g, pre=""):
        for k in g:
            o = g[k]
            if isinstance(o, h5py.Group):
                print(f"  {pre}{k}/")
                walk(o, pre + "  ")
            else:
                print(f"  {pre}{k}  {o.shape} {o.dtype}")
    walk(f)
