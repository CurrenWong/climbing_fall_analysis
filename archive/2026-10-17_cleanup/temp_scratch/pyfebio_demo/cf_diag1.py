"""CARTFIX 诊断 1：结构 + 抬回机制验证 + detJ 完整采样 + 扭转拆解。

输出：
  A. 源 npz 的面/三角化计数（复算 cart_patch 的 tri 生成）
  B. cart_m1 的楔形逐单元诊断（detJ 完整采样 / 底顶法向 / 扭角 / 面积 / 厚度）
  C. 抬回机制验证：顶面 == 底面 + 逐节点 disp 的重心插值？（用 cart_m1 反推）
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
SRC = np.load(ROOT / "temp/thums/cartilage_mesh.npz", allow_pickle=True)
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
bone_ids, bone_xyz = SRC["bone_node_ids"], SRC["bone_xyz"]
bpos = {int(v): bone_xyz[i] for i, v in enumerate(bone_ids)}

RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)
DEN = np.zeros((6, 3))


def det_penta_full(P):
    """完整采样（含 (r,s)=(1,0)/(0,1)）——与 ankle_lig_feb.detJ_penta6_full 同源。"""
    worst = None
    for r in RS:
        for sq in RS:
            if r + sq > 1.0 + 1e-9:
                continue
            for t in TS:
                G = np.array([
                    [-(1 - t), -(1 - t), -(1 - r - sq)],
                    [1 - t, 0.0, -r],
                    [0.0, 1 - t, -sq],
                    [-t, -t, (1 - r - sq)],
                    [t, 0.0, r],
                    [0.0, t, sq],
                ])
                det = float(np.linalg.det(P.T @ G))
                if worst is None or det < worst:
                    worst = det
    return worst


def det_penta_n(P, n):
    """细采样（n=5 -> 105 点）。"""
    rs = tuple(i / n for i in range(n + 1))
    ts = tuple(i / (n - 1) for i in range(n)) if n > 1 else (0.0,)
    ts = tuple(sorted(set(ts + (0.0,))))
    worst = None
    for r in rs:
        for sq in rs:
            if r + sq > 1.0 + 1e-9:
                continue
            for t in ts:
                G = np.array([
                    [-(1 - t), -(1 - t), -(1 - r - sq)],
                    [1 - t, 0.0, -r],
                    [0.0, 1 - t, -sq],
                    [-t, -t, (1 - r - sq)],
                    [t, 0.0, r],
                    [0.0, t, sq],
                ])
                det = float(np.linalg.det(P.T @ G))
                if worst is None or det < worst:
                    worst = det
    return worst


print("=" * 78)
print("A. 源 npz 复算 cart_patch 的 tri 生成")
print("=" * 78)
for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    E = np.asarray(SRC[f"{j}_elems"], np.int64)
    faces = []
    for row in E:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        faces.append((v[:m], v[m:]))
    nl = [len(b) for b, _ in faces]
    print(f"\n[{j}] 源面 {len(faces)}  底节点数分布 {dict((x, nl.count(x)) for x in sorted(set(nl)))}")
    print(f"       行宽 {E.shape[1]}  首行={E[0]}")
    # 复算 tris（DIAG=best）
    def _ma(P):
        a = []
        for k in range(3):
            u = P[(k - 1) % 3] - P[k]
            v = P[(k + 1) % 3] - P[k]
            c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-15)
            a.append(np.arccos(np.clip(c, -1, 1)))
        return min(a)
    tris, nsw = [], 0
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
            nsw += 1
        else:
            tris.append((base[0], base[1], base[2]))
            tris.append((base[0], base[2], base[3]))
    ncm = len(np.asarray(CM1[f"{j}_elems"]))
    print(f"       复算 tri = {len(tris)} (换对角 {nsw}) ；but cart_m1 ncell = {ncm}")

print()
print("=" * 78)
print("B. cart_m1 楔形逐单元诊断（完整采样）")
print("=" * 78)
store = {}
for j in ["tibiotalar_tibia", "tibiotalar_talus", "subtalar_talus", "subtalar_calcaneus"]:
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    off = np.asarray(CM1[f"{j}_is_offset"], bool)

    def ix(x):
        x = int(x)
        return (-x - 1) if x < 0 else (x - 1)

    W, W2 = [], []
    rec = []
    for k, row in enumerate(E):
        c = [ix(x) for x in row if int(x) != 0]
        if len(c) != 6:
            continue
        B, T = np.array(c[:3]), np.array(c[3:])
        P = N[c]
        w = det_penta_full(P)
        w2 = det_penta_n(P, 5)
        W.append(w)
        W2.append(w2)
        nb = np.cross(N[B[1]] - N[B[0]], N[B[2]] - N[B[0]])
        nt = np.cross(N[T[1]] - N[T[0]], N[T[2]] - N[T[0]])
        ab, at = np.linalg.norm(nb), np.linalg.norm(nt)
        cb, ct = N[B].mean(0), N[T].mean(0)
        d = N[T] - N[B]
        rec.append(dict(k=k, det=w, det105=w2, B=B, T=T,
                        area_b=ab / 2, area_t=at / 2,
                        twist=float(np.degrees(np.arccos(np.clip(nb @ nt / (ab * at + 1e-15), -1, 1)))),
                        dotnn=float(nb @ nt / (ab * at + 1e-15)),
                        disp=d, thick=float(np.linalg.norm(d, axis=1).mean()),
                        ctr=cb, offB=bool(off[B].any()), offT=bool(off[T].all())))
    W = np.array(W)
    W2 = np.array(W2)
    store[j] = (N, E, rec)
    print(f"\n[{j}]  {len(W)} 个楔形")
    print(f"   detJ 完整采样 (6x6x5=...): min={W.min():+.4e}  ≤0 的 {int((W<=0).sum())} 个   "
          f"≤1e-2 的 {int((W<=1e-2).sum())} 个")
    print(f"   detJ 细采样(105 点)       : min={W2.min():+.4e}  ≤0 的 {int((W2<=0).sum())} 个   "
          f"≤1e-2 的 {int((W2<=1e-2).sum())} 个")
    bad = [r for r in rec if r["det"] <= 1e-2]
    gd = [r for r in rec if r["det"] > 1e-2]
    if bad:
        print(f"   坏单元 idx (detJ≤1e-2): {[r['k'] for r in bad]}")
        print(f"   坏单元最差 5: {[(r['k'], round(r['det'],4)) for r in sorted(bad, key=lambda z: z['det'])[:5]]}")
    if bad and gd:
        def st(key, lst):
            v = np.array([r[key] for r in lst], float)
            return f"{np.median(v):.4f} (min {v.min():.4f} / max {v.max():.4f})"
        print("   ── 坏 vs 好 对比 ──")
        for key in ("twist", "area_b", "area_t", "thick"):
            print(f"     {key:8s} 坏: {st(key, bad):38s} 好: {st(key, gd)}")
        print(f"     坏单元 detJ 与 twist 相关: " +
              " ".join(f"#{r['k']}:det={r['det']:.2f},tw={r['twist']:.1f}°" for r in
                       sorted(bad, key=lambda z: z["det"])[:8]))

np.savez(ROOT / "temp/pyfebio_demo/cf_diag1.npz",
         **{f"{j}_det": np.array([r["det"] for r in store[j][2]]) for j in store},
         **{f"{j}_tw": np.array([r["twist"] for r in store[j][2]]) for j in store},
         **{f"{j}_ab": np.array([r["area_b"] for r in store[j][2]]) for j in store})
print("\n[saved] temp/pyfebio_demo/cf_diag1.npz")
