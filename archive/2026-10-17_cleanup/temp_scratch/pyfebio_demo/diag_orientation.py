"""针对单个软骨 hex 元素，实测各节点排序下的 det(J)，找出 FEBio 认可的朝向。
同时给出 FEBio hex8 的节点顺序约定验证：用已知合法的 THUMS 骨 hex 做对照。
"""
import math
from pathlib import Path
import numpy as np

NPZ = Path(r"D:\Project\climbing_fall_analysis\temp\thums\thums_lowerlimb.npz")
CG = Path(r"D:\Project\climbing_fall_analysis\temp\thums\cartilage_mesh.npz")
db = np.load(NPZ)
nid, xyz_all = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}
cg = np.load(CG)

GP = 1 / math.sqrt(3)
REF = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)


def detJ_hex(P, strict=False):
    pts = [np.array([[-GP, -GP, -GP], [GP, -GP, -GP], [GP, GP, -GP], [-GP, GP, -GP],
                     [-GP, -GP, GP], [GP, -GP, GP], [GP, GP, GP], [-GP, GP, GP]])]
    samples = [np.array([[-GP, -GP, -GP], [GP, -GP, -GP], [GP, GP, -GP], [-GP, GP, -GP],
                         [-GP, -GP, GP], [GP, -GP, GP], [GP, GP, GP], [-GP, GP, GP]])]
    if strict:
        samples.append(REF)          # 也检查 8 个节点
    worst = 1e30
    for S in samples:
        for s in S:
            dN = np.zeros((8, 3))
            for i in range(8):
                a, b, c = REF[i]
                dN[i, 0] = a * (1 + b * s[1]) * (1 + c * s[2]) / 8
                dN[i, 1] = b * (1 + a * s[0]) * (1 + c * s[2]) / 8
                dN[i, 2] = c * (1 + a * s[0]) * (1 + b * s[1]) / 8
            worst = min(worst, float(np.linalg.det(P.T @ dN)))
    return worst


# 对照：THUMS 骨 hex（已知 FEBio 接受）
print("=== 对照：THUMS tibia_cort 的 hex（FEBio 接受，det 应为正）===")
arr = db["81000700_solid"]
n = 0
for row in arr:
    conn = []
    for x in row:
        xi = int(x)
        if xi != -1 and xi not in conn:
            conn.append(xi)
    if len(conn) != 8:
        continue
    P = xyz_all[[idx[v] for v in conn]]
    print(f"  elem: det(J) gauss={detJ_hex(P):.4f}  strict(含节点)={detJ_hex(P, True):.4f}")
    n += 1
    if n >= 3:
        break

# 软骨 hex
print("\n=== 软骨 tibiotalar_tibia hex ===")
j = "tibiotalar_tibia"
elems, offs = cg[f"{j}_elems"], cg[f"{j}_offsets"]
shown = 0
for row in elems:
    vals = [int(x) for x in row if int(x) != 0]
    if not vals or vals[0] <= 0:
        continue
    m = len(vals) // 2
    if m != 4:
        continue
    base_id = vals[:m]
    off_id = [-x - 1 for x in vals[m:]]
    Pb = xyz_all[[idx[v] for v in base_id]]
    Po = offs[off_id]
    combos = {
        "base+, off+": np.vstack([Pb, Po]),
        "base-, off-": np.vstack([Pb[::-1], Po[::-1]]),
        "base-, off+": np.vstack([Pb[::-1], Po]),
        "base+, off-": np.vstack([Pb, Po[::-1]]),
        "off+ then base+": np.vstack([Po, Pb]),
        "off- then base-": np.vstack([Po[::-1], Pb[::-1]]),
    }
    nn = np.zeros(3)
    for i in range(4):
        nn += np.cross(Pb[i], Pb[(i + 1) % 4])
    oc, bc = Po.mean(0), Pb.mean(0)
    print(f"\n  base node ids {base_id}  off ids {off_id}")
    print(f"  base 面 Newell 法向 (原序) = {np.round(nn,3)}")
    print(f"  (off心 - base心)          = {np.round(oc-bc,4)}")
    print(f"  两者点积 = {np.dot(nn, oc-bc):+.4f}  (正=offset 朝法向侧)")
    for nm, P in combos.items():
        print(f"    {nm:16s} det(J): gauss={detJ_hex(P):11.4f}  strict={detJ_hex(P,True):11.4f}")
    shown += 1
    if shown >= 2:
        break
