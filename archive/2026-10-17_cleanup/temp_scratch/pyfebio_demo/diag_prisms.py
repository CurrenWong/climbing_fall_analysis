"""按**实际发射的 penta6 楔形**逐个算 det(J)，定位退化者。

发射规则（与 ankle_lig_feb.py 一致）：
  三角底(m=3) → 1 个 penta6: [b0,b1,b2, o0,o1,o2]
  四边底(m=4) → 拆 2 个: (0,1,2)/(4,5,6) 与 (0,2,3)/(4,6,7)
det(J) 用采样版（与闸门同一公式）。
"""
import numpy as np
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
cg = np.load(W / "cartilage_mesh.npz", allow_pickle=True)
nid = cg["bone_node_ids"]
xyz = np.asarray(cg["bone_xyz"], float)
idx = {int(v): i for i, v in enumerate(nid)}

RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)


def detJ_penta6(P):
    """与 ankle_lig_feb.py 的采样版完全一致。"""
    worst = None
    for r in RS:
        for sq in RS:
            if r + sq > 1.0 + 1e-9:
                continue
            for t in TS:
                d = np.array([
                    [-(1 - t), -(1 - t), -(1 - r - sq)],
                    [1 - t, 0.0, -r],
                    [0.0, 1 - t, -sq],
                    [-t, -t, (1 - r - sq)],
                    [t, 0.0, r],
                    [0.0, t, sq],
                ])
                det = float(np.linalg.det(P.T @ d))
                if worst is None or det < worst:
                    worst = det
    return worst


for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    E = cg[f"{j}_elems"]
    O = np.asarray(cg[f"{j}_offsets"], float)
    pris = []          # (m, base_ids, off_idx, tri)
    for row in E:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        m = len(vals) // 2
        base = vals[:m]
        oi = [-x - 1 for x in vals[m:]]
        if m == 3:
            pris.append((3, base, oi, (0, 1, 2)))
        else:
            pris.append((4, base, oi, (0, 1, 2)))
            pris.append((4, base, oi, (0, 2, 3)))
    dets, diag_used = [], []
    bad = []
    for m, base, oi, tri in pris:
        if any(b not in idx for b in base):
            continue
        B = np.array([xyz[idx[b]] for b in base])
        T = O[oi]
        P = np.array([B[z] for z in tri] + [T[z] for z in tri])
        d = detJ_penta6(P)
        dets.append(d)
        if d <= 0.5:
            # 这个楔形是哪个底？面积多少？
            Ab = 0.5 * np.linalg.norm(np.cross(B[tri[1]] - B[tri[0]],
                                               B[tri[2]] - B[tri[0]]))
            bad.append((d, m, Ab, tri))
    a = np.array(dets)
    print(f"\n{'='*74}\n{j}: 发射 {len(pris)} 个 penta6，可算 {len(a)}")
    print(f"  det(J)：中位 {np.median(a):9.3f}  p5 {np.percentile(a,5):9.3f}  "
          f"p1 {np.percentile(a,1):9.3f}  min {a.min():9.3f}")
    print(f"  det(J) <= 0.5（被剔除）: {(a<=0.5).sum()} 个 "
          f"({(a<=0.5).mean()*100:.1f}%)")
    print(f"  det(J) <= 0（会翻负）: {(a<=0).sum()}")
    if bad:
        bad.sort()
        print("  最差 8 个：")
        for d, m, Ab, tri in bad[:8]:
            print(f"    det={d:10.4f}  底 m={m} 底面积 {Ab:7.3f} mm²  三角 {tri}")
    # 按底类型分
    for mm in (3, 4):
        sub = [d for (m, _, _, _), d in zip(pris, dets) if m == mm]
        if sub:
            s = np.array(sub)
            print(f"  底 m={mm}: {len(s)} 个, det 中位 {np.median(s):9.3f}, "
                  f"min {s.min():9.4f}, <=0.5 有 {(s<=0.5).sum()}")
