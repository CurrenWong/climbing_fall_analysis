"""CARTFIX：对比「脚本内 tris/disp」与「cart_m1 反推的 tris/disp」，并分别解缠。"""
import io
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
SRC = io.open(ROOT / "scripts/ankle_fe/cart_patch.py", encoding="utf-8").read()
ns = {"__name__": "cp_head"}
exec(compile(SRC.split("\nsrc = np.load(")[0], "head", "exec"), ns)
untangle_disp = ns["untangle_disp"]

S0 = np.load(ROOT / "temp/thums/cartilage_mesh.npz", allow_pickle=True)
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
bids, bxyz = S0["bone_node_ids"], S0["bone_xyz"]
bpos = {int(v): bxyz[i] for i, v in enumerate(bids)}


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


from scipy.spatial import Delaunay   # noqa: E402


def build_A(j):
    """cart_patch 的做法：源面 -> 最优对角线三角化 + 逐节点中位 disp。"""
    E = np.asarray(S0[f"{j}_elems"], np.int64)
    OFF = np.asarray(S0[f"{j}_offsets"], float)
    faces = []
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        faces.append((v[:m], v[m:]))
    nodes = sorted({r for b, _ in faces for r in b})
    dacc = {r: [] for r in nodes}
    for base, outer in faces:
        for r, o in zip(base, outer):
            dacc[r].append(OFF[-o - 1] - bpos[r])
    disp = {r: np.median(np.array(v), axis=0) for r, v in dacc.items() if v}

    def _ma(P):
        a = []
        for k in range(3):
            u = P[(k - 1) % 3] - P[k]
            v = P[(k + 1) % 3] - P[k]
            c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-15)
            a.append(np.arccos(np.clip(c, -1, 1)))
        return min(a)
    tris = []
    for base, _ in faces:
        if len(base) == 3:
            tris.append(tuple(base))
            continue
        Q = [np.array(bpos[r]) for r in base]
        m02 = min(_ma(np.array([Q[0], Q[1], Q[2]])), _ma(np.array([Q[0], Q[2], Q[3]])))
        m13 = min(_ma(np.array([Q[1], Q[2], Q[3]])), _ma(np.array([Q[1], Q[3], Q[0]])))
        if m13 > m02:
            tris.append((base[1], base[2], base[3]))
            tris.append((base[1], base[3], base[0]))
        else:
            tris.append((base[0], base[1], base[2]))
            tris.append((base[0], base[2], base[3]))
    return tris, disp


def build_B(j):
    """从 cart_m1 反推：tri = 存下来的底面 3 点对应的骨 id；disp = 顶-底。"""
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    C = np.array([[ix(x) for x in row] for row in E])
    key = {tuple(np.round(bxyz[i], 9)): int(v) for i, v in enumerate(bids)}
    tris, disp = [], {}
    for e in range(len(E)):
        t = []
        for a in range(3):
            p = N[int(C[e, a])]
            bid = key.get(tuple(np.round(p, 9)))
            if bid is None:
                d = np.linalg.norm(bxyz - p, axis=1)
                bid = int(bids[int(np.argmin(d))])
            t.append(bid)
            disp[bid] = N[int(C[e, a + 3])] - N[int(C[e, a])]
        tris.append(tuple(t))
    return tris, disp


for j in ["tibiotalar_talus"]:
    tA, dA = build_A(j)
    tB, dB = build_B(j)
    print(f"{j}: A tri={len(tA)} disp={len(dA)}   B tri={len(tB)} disp={len(dB)}")
    sA = {frozenset(t) for t in tA}
    sB = {frozenset(t) for t in tB}
    print(f"  三角集合相同? {sA == sB}  仅 A {len(sA - sB)}  仅 B {len(sB - sA)}")
    if sA != sB:
        print("  样例（仅 A）:", [tuple(x) for x in list(sA - sB)[:3]])
    diff = [(r, np.round(dA[r] - dB[r], 9).tolist()) for r in dA if r in dB
            and not np.allclose(dA[r], dB[r], atol=1e-9)]
    print(f"  disp 不同的节点数: {len(diff)} / {len(dA)}  {diff[:3]}")
    onlyA = set(dA) - set(dB)
    onlyB = set(dB) - set(dA)
    print(f"  disp 只在 A: {len(onlyA)}  只在 B: {len(onlyB)}")
    for nm, (tt, dd) in (("A(脚本内)", (tA, dA)), ("B(cart_m1)", (tB, dB))):
        _, st = untangle_disp(dict(dd), tt, bpos)
        print(f"   {nm}: {st['before']:+.4f} → {st['after']:+.4f}  ≤0 {st['nneg_before']}"
              f"→{st['nneg_after']}  Δd {st['dmove']:.3f}", flush=True)
