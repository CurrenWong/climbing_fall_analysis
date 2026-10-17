"""独立验收 subagent 产物：对比 cart_m1（基线）与 cart_m1_fix 的 worst-det(J) + 厚度。"""
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")


def dNw(r, s, t):
    return np.array([-(1-t), -(1-t), -(1-r-s), (1-t), 0, -r, 0, (1-t), -s,
                     -t, -t, (1-r-s), t, 0, r, 0, t, s]).reshape(6, 3)


# 完整采样：RS 含 1.0（覆盖三角远端顶点），TS 含端点
RS = (0.0, 1/6, 1/3, 1/2, 2/3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)
PTS = [(r, sq, t) for r in RS for sq in RS if r + sq <= 1 + 1e-9 for t in TS]
print(f"采样点数 = {len(PTS)}")


def worst_det(P):
    return min(float(np.linalg.det(dNw(*p).T @ P)) for p in PTS)


def thick(P):
    return float(np.linalg.norm(P[[3, 4, 5]].mean(0) - P[[0, 1, 2]].mean(0)))


for tag in ("cart_m1", "cart_m1_fix"):
    d = R / "temp/thums" / tag
    f = d / "cartilage_mesh.npz"
    if not f.exists():
        print(f"✗ {tag} 缺 npz")
        continue
    z = np.load(f, allow_pickle=True)
    joints = sorted({k.split("_nodes")[0] for k in z.files if k.endswith("_nodes")})
    print("\n" + "=" * 92)
    print(f"{tag}   ({f.stat().st_size:,} B)   键: {len(z.files)}   关节: {joints}")
    print("=" * 92)
    gneg = gsm = gt = 0
    allt = []
    for j in joints:
        N = np.asarray(z[f"{j}_nodes"], float)
        E = np.asarray(z[f"{j}_elems"], int)
        w, tt = [], []
        for e in E:
            try:
                P = N[np.asarray(e, int)]
            except Exception:
                continue
            if P.shape[0] < 6:
                continue
            w.append(worst_det(P))
            tt.append(thick(P))
        w = np.array(w)
        tt = np.array(tt)
        allt.append(tt)
        gneg += int((w <= 0).sum())
        gsm += int((w <= 1e-2).sum())
        gt += len(w)
        print(f"  {j:<24} 单元 {len(w):>5}  worst-det(J) min {w.min():>12.5e}  "
              f"≤0 {int((w<=0).sum()):>3}  ≤1e-2 {int((w<=1e-2).sum()):>3}  "
              f"厚度 min {tt.min():.4f} 中位 {np.median(tt):.4f} "
              f"max {tt.max():.4f}")
    allt = np.concatenate(allt) if allt else np.array([])
    print(f"  ⇒ 合计 {gt} 单元：**worst-det(J) ≤ 0 的 {gneg} 个**、≤1e-2 的 {gsm} 个")
    if len(allt):
        print(f"  ⇒ 厚度：min {allt.min():.4f} / 中位 {np.median(allt):.4f} / "
              f"max {allt.max():.4f} mm；<0.3 的 {int((allt<0.3).sum())} 个、"
              f">1.6 的 {int((allt>1.6).sum())} 个")
