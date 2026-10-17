"""核对各部件包围盒, 确认是否同一坐标系 (装配前提)。"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from msh_info import parse  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
MESH = ROOT / "temp" / "ankle_mesh"

print(f"{'部件':<24}{'x (cm)':>20}{'y (cm)':>20}{'z (cm)':>20}")
print("-" * 84)
info = {}
for f in sorted(MESH.glob("*.msh")):
    nodes, tets = parse(f)
    if not nodes:
        continue
    a = np.array(list(nodes.values()))
    lo, hi = a.min(0), a.max(0)
    info[f.name] = (lo, hi)
    print(
        f"{f.name:<24}"
        f"{lo[0]*100:>9.2f}~{hi[0]*100:>9.2f}"
        f"{lo[1]*100:>10.2f}~{hi[1]*100:>9.2f}"
        f"{lo[2]*100:>10.2f}~{hi[2]*100:>9.2f}"
    )

print("\n=== 相对位置检查 (踝关节应在距骨/跟骨/胫骨重叠处) ===")
for k, (lo, hi) in info.items():
    c = (lo + hi) / 2
    print(f"  {k:<24} 中心 = [{c[0]*100:>7.2f}, {c[1]*100:>7.2f}, {c[2]*100:>7.2f}] cm")
