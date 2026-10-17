"""CARTFIX 诊断 3：disp 场为什么跳变 + 受控实验（决定修复方案）。

A. 源数据：逐面 offset 向量 = old_off[-o-1] - bone_pos。
   |offset| 是否恒为 1.5mm？方向是否 = 该面法向？（若是 ⇒ 逐面法向 ⇒ 折边上中位数变短）
B. 关键骨节点：列出其所有相邻面的 offset 向量/模长/与面法向夹角
C. 受控实验：坏单元在 (i) 原 disp (ii) disp 取常量 (iii) 逐顶点替换为另两点均值
   (iv) disp 场 Laplacian 平滑 —— 各自 worst detJ
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
SRC = np.load(ROOT / "temp/thums/cartilage_mesh.npz", allow_pickle=True)
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
bids, bxyz = SRC["bone_node_ids"], SRC["bone_xyz"]
bpos = {int(v): bxyz[i] for i, v in enumerate(bids)}

RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)


def det_pt(P, r, s, t):
    G = np.array([
        [-(1 - t), -(1 - t), -(1 - r - s)],
        [1 - t, 0.0, -r],
        [0.0, 1 - t, -s],
        [-t, -t, (1 - r - s)],
        [t, 0.0, r],
        [0.0, t, s],
    ])
    return float(np.linalg.det(P.T @ G))


def worst(P):
    w, arg = None, None
    for r in RS:
        for s in RS:
            if r + s > 1 + 1e-9:
                continue
            for t in TS:
                v = det_pt(P, r, s, t)
                if w is None or v < w:
                    w, arg = v, (r, s, t)
    return w, arg


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


# ══ A. 源数据：逐面 offset ══
print("=" * 80)
print("A. 源数据逐面 offset（= old_off[-o-1] - bone_pos）")
print("=" * 80)
FACES = {}
for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    E = np.asarray(SRC[f"{j}_elems"], np.int64)
    OFF = np.asarray(SRC[f"{j}_offsets"], float)
    faces = []
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        faces.append((v[:m], v[m:]))
    FACES[j] = (faces, OFF)
    mags, angs, osa = [], [], []
    for base, outer in faces:
        P = np.array([bpos[r] for r in base])
        nf = np.cross(P[1] - P[0], P[2] - P[0])
        nf = nf / (np.linalg.norm(nf) + 1e-15)
        for r, o in zip(base, outer):
            v = OFF[-o - 1] - bpos[r]
            mags.append(np.linalg.norm(v))
            c = np.dot(v, nf) / (np.linalg.norm(v) + 1e-15)
            angs.append(np.degrees(np.arccos(np.clip(abs(c), -1, 1))))
            osa.append(v / (np.linalg.norm(v) + 1e-15))
    mags = np.array(mags)
    angs = np.array(angs)
    osa = np.array(osa)
    print(f"\n[{j}] 面 {len(faces)}  offset 样本 {len(mags)}")
    print(f"  |offset| : min={mags.min():.6f} max={mags.max():.6f} 中位={np.median(mags):.6f}")
    print(f"  |offset| == 1.5 (1e-6) 的比例: {(np.abs(mags-1.5) < 1e-6).mean()*100:.1f}%")
    print(f"  与所在面法向夹角: 中位={np.median(angs):.3f}° max={angs.max():.3f}°  "
          f"(*注：4 节点面用前三点的法向，故有偏差)")
    # 同一骨节点上的多个 offset 方向是否发散
    acc = {}
    for base, outer in faces:
        for r, o in zip(base, outer):
            acc.setdefault(r, []).append(OFF[-o - 1] - bpos[r])
    spread = []
    for r, vs in acc.items():
        if len(vs) < 2:
            continue
        V = np.array(vs)
        sp = max(np.degrees(np.arccos(np.clip(np.dot(a / np.linalg.norm(a), b / np.linalg.norm(b)), -1, 1)))
                 for i, a in enumerate(V) for b in V[i + 1:])
        spread.append((r, len(vs), sp, np.linalg.norm(np.median(V, axis=0)),
                       np.mean([np.linalg.norm(v) for v in V])))
    spread.sort(key=lambda z: -z[2])
    print("  ── 同节点上 offset 方向发散最大的 8 个（⇒ 中位数变短）──")
    print(f"     {'骨id':>11} {'面数':>4} {'最大夹角':>8} {'|中位|':>8} {'平均|offset|':>12}")
    for r, n, sp, mmed, mmean in spread[:8]:
        print(f"     {r:>11} {n:>4} {sp:8.2f} {mmed:8.4f} {mmean:12.4f}")

# ══ B. 关键节点的 disp 全貌 ══
print()
print("=" * 80)
print("B. 关键骨节点：其相邻面的 offset 明细 vs 中位数 disp")
print("=" * 80)
KEY = {"tibiotalar_tibia": [81001898, 81001749, 81001662, 81002570, 81002571,
                            81001752, 81002572, 81001754, 81002565, 81002580],
       "tibiotalar_talus": [81000070, 81013288, 81000667, 81013283, 81000539,
                            81000090, 81013282, 81000557, 81000580]}
for j, keys in KEY.items():
    faces, OFF = FACES[j]
    acc = {}
    for base, outer in faces:
        for r, o in zip(base, outer):
            acc.setdefault(r, []).append(OFF[-o - 1] - bpos[r])
    print(f"\n[{j}]")
    for r in keys:
        if r not in acc:
            continue
        V = np.array(acc[r])
        med = np.median(V, axis=0)
        print(f"  骨id {r}: 面数 {len(V)}  |中位|={np.linalg.norm(med):.4f}  "
              f"各面 |offset|={[round(float(np.linalg.norm(v)),4) for v in V]}")
        # 面法向夹角
        nms = []
        for base, outer in faces:
            if r not in base:
                continue
            P = np.array([bpos[q] for q in base])
            nf = np.cross(P[1] - P[0], P[2] - P[0])
            nms.append(nf / (np.linalg.norm(nf) + 1e-15))
        if len(nms) > 1:
            nms = np.array(nms)
            sp = max(np.degrees(np.arccos(np.clip(np.dot(a, b), -1, 1)))
                     for i, a in enumerate(nms) for b in nms[i + 1:])
            print(f"        所在面法向之间最大夹角 = {sp:.2f}°")

# ══ C. 受控实验 ══
print()
print("=" * 80)
print("C. 受控实验：坏单元 detJ 对 disp 扰动的敏感性")
print("=" * 80)
for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    off = np.asarray(CM1[f"{j}_is_offset"], bool)
    bt = [i for i in range(len(N)) if not off[i]]
    pair = {}
    for row in E:
        c = [ix(x) for x in row]
        for b, t in zip(c[:3], c[3:]):
            pair[b] = t
    disp = {b: N[pair[b]] - N[b] for b in pair}
    adj = {b: set() for b in bt}
    for row in E:
        c = [ix(x) for x in row]
        for a in range(3):
            for bb in range(3):
                if a != bb:
                    adj[c[a]].add(c[bb])
    print(f"\n[{j}]")
    for k, row in enumerate(E):
        c = [ix(x) for x in row]
        B, T = c[:3], c[3:]
        w0, arg0 = worst(N[c])
        if w0 > 1e-2:
            continue
        D = np.array([disp[b] for b in B])
        Dm = D.mean(axis=0)
        res = {}
        for q in range(3):
            D2 = D.copy()
            D2[q] = (D[(q + 1) % 3] + D[(q + 2) % 3]) / 2
            P2 = np.array([N[B[q2]] + (Dm if q2 == q and False else D2[q2]) for q2 in range(3)])
            P2 = np.concatenate([N[B], N[B] + D2], axis=0)
            res[f"仅改顶点{q}"] = worst(P2)[0]
        Pc = np.concatenate([N[B], N[B] + Dm], axis=0)
        # 局部平滑：只用「本次三角形 3 顶点的邻居中位数」替换各自 disp（1 次 Jacobi）
        Ds = []
        for b in B:
            nb = list(adj[b])
            loc = np.array([disp[x] for x in nb] + [disp[b]])
            Ds.append(np.median(loc, axis=0))
        Ds = np.array(Ds)
        Ps = np.concatenate([N[B], N[B] + Ds], axis=0)
        print(f"   #{k:>4} 原 detJ={w0:+9.3f} @{arg0}  |  常量 disp -> {worst(Pc)[0]:+9.3f}  "
              f"| 仅改单个顶点 -> " +
              ", ".join(f"{v:+.1f}" for v in res.values()) +
              f"  | 1×邻域中位 -> {worst(Ps)[0]:+9.3f}")
