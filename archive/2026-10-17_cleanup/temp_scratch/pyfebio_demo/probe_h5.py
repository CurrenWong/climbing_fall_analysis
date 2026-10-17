"""探查 xplt->h5 里各 dataset 的真实 dtype/shape/内容布局。"""

from pathlib import Path

import h5py

HERE = Path(__file__).parent
dst = HERE / "my_model_h5.h5"

with h5py.File(dst, "r") as f:
    for path in [
        "meshes/0/nodes",
        "states/0/node_data/displacement/nodes",
        "states/0/element_data/stress/box",
    ]:
        d = f[path]
        arr = d[:]
        print(f"\n=== {path} ===")
        print(f"  dtype={d.dtype}  shape={d.shape}")
        print(f"  arr type={type(arr).__name__}")
        print(f"  arr = {arr}")
