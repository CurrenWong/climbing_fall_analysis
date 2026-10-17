"""ATFL 几何真值：距骨/腓骨表面最近距离在哪？A/B 簇到底在解剖何处？"""
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
z = np.load(ROOT / "temp/thums/thums_lowerlimb.npz", allow_pickle=True)
nid, xyz = z["node_ids"], z["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}

BONES = {
    "tibia": [81000700, 81000600, 81000601],
    "fibula": [81000900, 81000800],
    "talus": [81001100, 81001000],
    "calcaneus": [81001300, 81001200],
    "navicular": [81001400, 81001401],
}
B = {}
for b, pids in BONES.items():
    s = set()
    for pid in pids:
        for kind in ("solid", "shell"):
            k = f"{pid:08d}_{kind}"
            if k in z:
                s |= {int(v) for v in z[k].ravel()}
    B[b] = np.array([xyz[i_of[n]] for n in sorted(s)])

print("=== 各骨包围盒 / 质心（THUMS 原始坐标）===")
for b, Q in B.items():
    print(f"  {b:<10} n={len(Q):>5}  min={np.round(Q.min(0),1)}  "
          f"max={np.round(Q.max(0),1)}  对角={np.round(Q.max(0)-Q.min(0),1)}  "
          f"质心={np.round(Q.mean(0),1)}")

print("\n=== 骨-骨表面最近距离 ===")
def mind(A, Bx):
    best, arg = 1e9, None
    for i in range(0, len(A), 1):
        d = np.sqrt(((Bx - A[i]) ** 2).sum(1))
        j = int(d.argmin())
        if d[j] < best:
            best, arg = float(d[j]), (A[i], Bx[j])
    return best, arg

pairs = [("talus", "fibula"), ("tibia", "talus"), ("fibula", "calcaneus"),
         ("tibia", "fibula"), ("talus", "calcaneus")]
for a, b in pairs:
    d, arg = mind(B[a], B[b])
    print(f"  {a:<10} ↔ {b:<10} 最近 {d:>7.2f} mm   在 "
          f"{np.round(arg[0],1)} / {np.round(arg[1],1)}")

# ATFL 的 A/B 簇
ATFL_PIDS = (81100901,)
e = z["81100901_shell"]
nodes = sorted({int(v) for v in e.ravel()})
BONE_IDS = {b: {int(v) for pid in pids for kind in ("solid", "shell")
                if f"{pid:08d}_{kind}" in z
                for v in z[f"{pid:08d}_{kind}"].ravel()}
            for b, pids in BONES.items()}
A = [n for n in nodes if n in BONE_IDS["talus"]]
Bb = [n for n in nodes if n in BONE_IDS["fibula"]]
Pa = xyz[[i_of[n] for n in A]]
Pb = xyz[[i_of[n] for n in Bb]]
print(f"\n=== ATFL: A(距骨) {len(A)} 个 / B(腓骨) {len(Bb)} 个")
print(f"  A 簇质心 {np.round(Pa.mean(0),1)}")
print(f"  B 簇质心 {np.round(Pb.mean(0),1)}")
print(f"  A 簇到腓骨面最近 = {mind(Pa, B['fibula'])[0]:.2f} mm")
print(f"  B 簇到距骨面最近 = {mind(Pb, B['talus'])[0]:.2f} mm")
print(f"  A 簇到胫骨面最近 = {mind(Pa, B['tibia'])[0]:.2f} mm")
# 各簇在各自骨包围盒中的相对位置（0=min, 1=max）
for lbl, Pc, bone in (("A(距骨)", Pa, "talus"), ("B(腓骨)", Pb, "fibula")):
    Q = B[bone]
    rel = (Pc.mean(0) - Q.min(0)) / (Q.max(0) - Q.min(0) + 1e-12)
    print(f"  {lbl} 在 {bone} 包围盒中的相对位置 = {np.round(rel,3)}  (0=min边 1=max边)")
