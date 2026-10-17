"""CARTFIX 诊断 2：源 npz 定位（cart_patch.py 的输入长什么样）。"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
for p in sorted((ROOT / "temp/thums").glob("*.npz")):
    try:
        d = np.load(p, allow_pickle=True)
    except Exception as e:                                  # noqa: BLE001
        print(f"{p.name}: 读失败 {e}")
        continue
    ks = list(d.files)
    has_off = [k for k in ks if k.endswith("_offsets")]
    has_el = [k for k in ks if k.endswith("_elems")]
    print(f"\n### {p}  ({p.stat().st_size/1e6:.2f} MB, {len(ks)} 键)")
    print(f"   *_offsets: {has_off}")
    print(f"   *_elems  : {has_el}")
    if has_el:
        E = np.asarray(d[has_el[0]])
        print(f"   {has_el[0]} shape={E.shape} 采样行={E[0]}")
