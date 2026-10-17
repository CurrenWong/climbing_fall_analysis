"""验证绕向假设：底面 quad/tri 的 Newell 法向 vs 挤出方向，符号一致吗？"""
import numpy as np
from pathlib import Path

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}


def newell(pts):
    n = np.zeros(3)
    for i in range(len(pts)):
        n += np.cross(pts[i], pts[(i + 1) % len(pts)])
    return n


print("=" * 82)
print("底面面片：Newell 法向 · 挤出方向 的符号（+ = 绕向正确，− = 需要反转）")
print("=" * 82)
print(f"{'关节':22s} {'面数':>5s} {'+ 一致':>8s} {'− 反绕':>8s} {'反绕占比':>9s}")
for j in ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]:
    E = np.asarray(src[f"{j}_elems"], np.int64)
    off = src[f"{j}_offsets"]
    npos = nneg = 0
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        pb = np.array([bpos[r] for r in v[:m]])
        po = np.array([off[-r - 1] for r in v[m:]])
        n = newell(pb)
        ext = po.mean(axis=0) - pb.mean(axis=0)
        if np.dot(n, ext) > 0:
            npos += 1
        else:
            nneg += 1
    tot = npos + nneg
    print(f"{j:22s} {tot:5d} {npos:8d} {nneg:8d} {100*nneg/tot:8.1f}%")

print()
print("=" * 82)
print("★ 若反绕占比 ≈ 40% ⇒ 与 Gauss det(J)<0 的 40% 完全吻合 ⇒ 假设成立")
print("=" * 82)

# 再验：把绕向修正后，hex 的 det(J) 是否全部转正？
G = 1 / np.sqrt(3)
GP = [(-G, -G, -G), (G, -G, -G), (G, G, -G), (-G, G, -G),
      (-G, -G, G), (G, -G, G), (G, G, G), (-G, G, G)]
REF = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
       (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]


def detJ_min(P):
    w = 1e30
    for s in GP:
        dN = np.zeros((8, 3))
        for i in range(8):
            a, b, c = REF[i]
            dN[i, 0] = a * (1 + b * s[1]) * (1 + c * s[2]) / 8
            dN[i, 1] = b * (1 + a * s[0]) * (1 + c * s[2]) / 8
            dN[i, 2] = c * (1 + a * s[0]) * (1 + b * s[1]) / 8
        w = min(w, float(np.linalg.det(P.T @ dN)))
    return w


new = np.load(THUMS / "cartilage_mesh_hex.npz", allow_pickle=True)
print("\n修正绕向后（底面倒序）重测新网格 hex 的 Gauss det(J)：")
print(f"{'关节':22s} {'hex':>6s} {'修正前<0':>10s} {'修正后<0':>10s}")
for j in ["tibiotalar_tibia", "subtalar_calcaneus"]:
    nd = new[f"{j}_nodes"]
    E = new[f"{j}_elems"]

    def P_of(r):
        return nd[-r - 1] if r < 0 else bpos[r]

    nbad_b, nbad_a, nhex = 0, 0, 0
    for row in E:
        blk = [int(x) for x in row if int(x) != 0]
        m = len(blk) // 2
        if m != 4:
            continue
        nhex += 1
        P = np.array([P_of(r) for r in blk])
        if detJ_min(P) <= 0:
            nbad_b += 1
        # 正确做法：先按挤出方向把底面调成 CCW
        pb = P[:4]
        ext = P[4:].mean(axis=0) - pb.mean(axis=0)
        if np.dot(newell(pb), ext) < 0:
            P = P[[3, 2, 1, 0, 7, 6, 5, 4]]
        if detJ_min(P) <= 0:
            nbad_a += 1
    print(f"{j:22s} {nhex:6d} {nbad_b:10d} {nbad_a:10d}")
