"""对比 cart_patch 输出与旧软骨网格的 key 集合（判断能否直接替换）。"""
import io
import re
from pathlib import Path

import numpy as np

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
FILES = ["cartilage_mesh.npz", "cartilage_mesh_patch.npz", "cartilage_mesh_patch_m2.npz"]

sets = {}
for f in FILES:
    p = THUMS / f
    if not p.exists():
        print(f"(缺) {f}")
        continue
    d = np.load(p)
    sets[f] = set(d.files)
    print(f"\n=== {f}  ({p.stat().st_size/1e6:.2f} MB, {len(d.files)} keys) ===")
    for k in sorted(d.files):
        a = d[k]
        sh = getattr(a, "shape", ())
        print(f"   {k:34} {str(sh):>18}  {getattr(a,'dtype','')}")

if "cartilage_mesh.npz" in sets:
    base = sets["cartilage_mesh.npz"]
    for f in FILES[1:]:
        if f not in sets:
            continue
        only_new = sorted(sets[f] - base)
        only_old = sorted(base - sets[f])
        print(f"\n--- {f} vs cartilage_mesh.npz ---")
        print(f"   {f} 独有: {only_new}")
        print(f"   旧网格独有: {only_old}")
        print(f"   交集 {len(sets[f] & base)} keys")

# build_feb.py 如何取软骨
bf = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\build_feb.py")
if bf.exists():
    s = io.open(bf, encoding="utf-8").read()
    print("\n=== build_feb.py 中 CART / TIE / 网格相关行 ===")
    for i, l in enumerate(s.split("\n"), 1):
        if re.search(r"cartilage_mesh|CART_DIR|CART_ELEM|cart_patch|_patch", l):
            print(f"  {i:>4}|{l.strip()[:110]}")
