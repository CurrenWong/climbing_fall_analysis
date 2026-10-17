"""回归检查（纯 ASCII 输出）：原 26 个部件与备份逐位一致？"""
import numpy as np
from pathlib import Path

T = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
new = np.load(T / "thums_lowerlimb.npz", allow_pickle=True)
old = np.load(T / "thums_lowerlimb.npz.bak-preLIG", allow_pickle=True)

bad = []
for k in old.files:
    if k in ("node_ids", "node_xyz"):
        continue
    a, b = old[k], new[k]
    if a.shape != b.shape or not np.array_equal(a, b):
        bad.append((k, a.shape, b.shape))
print(f"old parts checked: {len(old.files) - 2}")
print(f"mismatches: {len(bad)}")
for k, s1, s2 in bad:
    print(f"  XX {k}: {s1} -> {s2}")

# 节点表：原 21,608 个必须仍是前缀/同值
oi = {int(v): i for i, v in enumerate(old["node_ids"])}
ni = {int(v): i for i, v in enumerate(new["node_ids"])}
miss = [k for k in oi if k not in ni]
diff = []
for k, io in oi.items():
    if k in ni and not np.allclose(old["node_xyz"][io], new["node_xyz"][ni[k]]):
        diff.append(k)
print(f"\nold node ids still present: {len(oi) - len(miss)}/{len(oi)}")
print(f"coord mismatches: {len(diff)}")
print("node_ids/xyz identical prefix:",
      np.array_equal(old["node_ids"], new["node_ids"][:len(old["node_ids"])])
      and np.allclose(old["node_xyz"], new["node_xyz"][:len(old["node_xyz"])]))
print("\nRESULT:", "PASS - existing parts untouched"
      if not bad and not miss and not diff else "FAIL")
