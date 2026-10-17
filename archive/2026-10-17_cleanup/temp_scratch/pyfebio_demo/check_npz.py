"""核对 cartilage_mesh.npz 是否与 cartilage.py 的防穿透逻辑一致（是否 stale）。"""
import datetime
import json
from pathlib import Path

import numpy as np

T = Path(r"D:\Project\climbing_fall_analysis\temp\thums")


def ts(p):
    return f"{datetime.datetime.fromtimestamp(p.stat().st_mtime):%Y-%m-%d %H:%M}"


print("=" * 74)
print("文件时间戳")
print("=" * 74)
for p in (T / "cartilage.py", T.parent.parent / "scripts/ankle_fe/cartilage.py",
          T / "cartilage_mesh.npz", T / "cartilage.json", T / "cartilage_mesh_hex.npz",
          T / "thums_lowerlimb.npz"):
    if p.exists():
        print(f"  {ts(p)}  {p.stat().st_size:>12,}  {p}")
    else:
        print(f"  {'—':16}  {'—':>12}  {p}  (不存在)")

j = T / "cartilage.json"
if j.exists():
    d = json.loads(j.read_text(encoding="utf-8"))
    print("\n" + "=" * 74)
    print("cartilage.json 记录的参数与检查")
    print("=" * 74)
    print(f"  params: {d.get('params')}")
    for jn, v in d.get("joints", {}).items():
        th = v.get("thickness_mm") or {}
        pv = v.get("penetration_check") or {}
        print(f"  {jn}:")
        print(f"    厚度 min {th.get('min')} / 中位 {th.get('median')} / max {th.get('max')}"
              f"  target {th.get('target')}")
        print(f"    穿透检查: {pv}")

cg = np.load(T / "cartilage_mesh.npz", allow_pickle=True)
print("\n" + "=" * 74)
print("cartilage_mesh.npz 内容")
print("=" * 74)
print(f"  数组: {sorted(cg.files)}")
for k in sorted(cg.files):
    a = cg[k]
    if hasattr(a, "shape"):
        print(f"    {k:26s} {str(a.shape):>16s} {a.dtype}")
