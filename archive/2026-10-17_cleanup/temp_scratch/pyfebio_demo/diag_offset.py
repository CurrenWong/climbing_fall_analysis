"""诊断软骨偏移质量：逐楔形算
  · 底面/顶面三角形法向夹角（扭转角）
  · 3 条偏移方向 (o_i - b_i) 的两两夹角（偏移场一致性）
  · 归一化体积（六面体/楔形的"有效性"）
并统计 det(J) 分布，找出退化楔形集中在哪。
"""
import numpy as np
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
cg = np.load(W / "cartilage_mesh.npz", allow_pickle=True)
print("npz keys:", list(cg.keys())[:20])

db = np.load(W / "thums_bones.npz", allow_pickle=True)
print("bones keys:", list(db.keys())[:10])


def unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.maximum(n, 1e-12)


for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    E = cg[f"{j}_elems"]
    offs = cg[f"{j}_offsets"]
    print(f"\n{'='*72}\n{j}: 元素 {len(E)} 行, offsets {len(offs)}")
    # 建立 offset 索引：负值 -k-1 -> off:k
    # 行 = [base..., offset 反索引...]
    stats = []
    for row in E:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        m = len(vals) // 2
        base = vals[:m]
        oidx = [-x - 1 for x in vals[m:]]
        stats.append((m, base, oidx))
    n3 = sum(1 for s in stats if s[0] == 3)
    n4 = sum(1 for s in stats if s[0] == 4)
    print(f"  三角形底 {n3} 个，四边形底 {n4} 个")

    # 用 offset 坐标算偏移方向一致性
    O = np.asarray(offs, float)
    print(f"  offsets 形状 {O.shape}")
    # 需要底面坐标 —— 从 thums_bones 里取
    nid = db["node_ids"]
    xyz = db["node_xyz"]
    idx = {int(v): i for i, v in enumerate(nid)}

    ang_list, twist_list, voln_list = [], [], []
    for m, base, oidx in stats:
        if any(b not in idx for b in base):
            continue
        B = np.array([xyz[idx[b]] for b in base])
        T = O[oidx]
        # 偏移方向
        d = T - B
        dn = unit(d)
        # 两两夹角
        if m == 3:
            c = [float(dn[0] @ dn[1]), float(dn[1] @ dn[2]), float(dn[2] @ dn[0])]
            ang = np.degrees(np.arccos(np.clip(c, -1, 1)))
            ang_list.append(float(ang.max()))
        # 底/顶面法向夹角（扭转）
        if m == 3:
            nb = np.cross(B[1] - B[0], B[2] - B[0])
            nt = np.cross(T[1] - T[0], T[2] - T[0])
        else:
            nb = np.cross(B[1] - B[0], B[3] - B[0])
            nt = np.cross(T[1] - T[0], T[3] - T[0])
        nb, nt = nb / max(np.linalg.norm(nb), 1e-12), nt / max(np.linalg.norm(nt), 1e-12)
        twist_list.append(float(np.degrees(np.arccos(
            np.clip(nb @ nt, -1, 1)))))
        # 归一化体积：V / (A_b * h_mean)
        if m == 3:
            Ab = 0.5 * np.linalg.norm(np.cross(B[1] - B[0], B[2] - B[0]))
            h = np.linalg.norm(d, axis=1).mean()
            V = abs(np.linalg.det(np.stack(
                [B[1] - B[0], B[2] - B[0],
                 (T[0] + T[1] + T[2] - B[0] - B[1] - B[2]) / 3.0]))) / 2.0
            voln_list.append(V / max(Ab * h, 1e-12))

    if ang_list:
        a = np.array(ang_list)
        print(f"  偏移方向两两最大夹角（度）：中位 {np.median(a):.1f}  "
              f"p90 {np.percentile(a,90):.1f}  max {a.max():.1f}")
        print(f"    >20° 的有 {(a>20).sum()} 个 ({(a>20).mean()*100:.1f}%)")
        print(f"    >40° 的有 {(a>40).sum()} 个 ({(a>40).mean()*100:.1f}%)")
    if twist_list:
        t = np.array(twist_list)
        print(f"  底/顶面法向夹角（度）：中位 {np.median(t):.1f}  "
              f"p90 {np.percentile(t,90):.1f}  max {t.max():.1f}")
        print(f"    >20° 的有 {(t>20).sum()} 个 ({(t>20).mean()*100:.1f}%)")
    if voln_list:
        v = np.array(voln_list)
        print(f"  归一化体积 V/(A·h)：中位 {np.median(v):.3f}  "
              f"p1 {np.percentile(v,1):.3f}  min {v.min():.3f}")
        print(f"    <0.5 的有 {(v<0.5).sum()} 个 ({(v<0.5).mean()*100:.1f}%)")
