"""核对 cart_new 里的网格到底是哪一份 + 确认 CART_NEW 判定。"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
for f in ["temp/thums/cartilage_mesh.npz", "temp/thums/cartilage_mesh_patch_m2.npz",
          "temp/thums/cart_new/cartilage_mesh.npz"]:
    p = ROOT / f
    if not p.exists():
        print(f"(缺) {f}")
        continue
    d = np.load(p)
    ks = sorted(d.files)
    has = "tibiotalar_tibia_nodes" in d.files
    print(f"\n=== {f}  ({p.stat().st_size/1e6:.2f} MB, {len(ks)} keys) ===")
    print(f"  含 tibiotalar_tibia_nodes ? {'YES  -> 新格式' if has else 'NO   -> 旧格式'}")
    for k in ks:
        if k.startswith("tibiotalar_tibia"):
            print(f"    {k:36} {getattr(d[k], 'shape', '')}")
