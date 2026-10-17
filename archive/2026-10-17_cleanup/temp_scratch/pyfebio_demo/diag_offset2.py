"""诊断软骨偏移质量（修正：骨坐标在同一 npz）。"""
import numpy as np
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
cg = np.load(W / "cartilage_mesh.npz", allow_pickle=True)
nid = cg["bone_node_ids"]
xyz = np.asarray(cg["bone_xyz"], float)
idx = {int(v): i for i, v in enumerate(nid)}
print(f"骨节点 {len(nid):,}")


def unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.maximum(n, 1e-12)


for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    E = cg[f"{j}_elems"]
    O = np.asarray(cg[f"{j}_offsets"], float)
    rows = []
    for row in E:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        m = len(vals) // 2
        rows.append((m, vals[:m], [-x - 1 for x in vals[m:]]))
    n3 = sum(1 for r in rows if r[0] == 3)
    n4 = sum(1 for r in rows if r[0] == 4)
    print(f"\n{'='*74}\n{j}: 共 {len(rows)} 楔形/六面体（三角底 {n3}，四边底 {n4}）"
          f"  offsets {len(O)}")

    angs, twists, vols, depths = [], [], [], []
    for m, base, oidx in rows:
        if any(b not in idx for b in base):
            continue
        B = np.array([xyz[idx[b]] for b in base])
        T = O[oidx]
        d = T - B
        h = np.linalg.norm(d, axis=1)
        depths.append(h.mean())
        dn = unit(d)
        if m == 3:
            angs.append(float(np.degrees(np.arccos(np.clip(
                [dn[0] @ dn[1], dn[1] @ dn[2], dn[2] @ dn[0]], -1, 1))).max()))
            nb = np.cross(B[1] - B[0], B[2] - B[0])
            nt = np.cross(T[1] - T[0], T[2] - T[0])
            Ab = 0.5 * np.linalg.norm(nb)
            V = abs(np.linalg.det(np.stack(
                [B[1] - B[0], B[2] - B[0],
                 (T[0] + T[1] + T[2] - B[0] - B[1] - B[2]) / 3.0]))) / 2.0
        else:
            nb = np.cross(B[1] - B[0], B[3] - B[0])
            nt = np.cross(T[1] - T[0], T[3] - T[0])
            Ab = np.linalg.norm(nb)
            # 六面体近似体积：底面积 × 平均厚
            V = Ab * h.mean()
        nb, nt = unit(nb), unit(nt)
        twists.append(float(np.degrees(np.arccos(np.clip(nb @ nt, -1, 1)))))
        vols.append(V / max(Ab * h.mean(), 1e-12))

    for nm, arr, hi in (("偏移方向两两最大夹角", angs, 20),
                        ("底/顶面法向夹角(扭转)", twists, 20)):
        if arr:
            a = np.array(arr)
            print(f"  {nm}（度）：中位 {np.median(a):5.1f}  p90 "
                  f"{np.percentile(a, 90):5.1f}  max {a.max():6.1f}   "
                  f">{hi}° 占 {(a > hi).mean()*100:4.1f}%")
    if vols:
        v = np.array(vols)
        print(f"  体积/（底面积×平均厚）：中位 {np.median(v):.3f}  "
              f"p1 {np.percentile(v,1):.3f}  min {v.min():.3f}   "
              f"<0.5 占 {(v<0.5).mean()*100:4.1f}%")
    if depths:
        hh = np.array(depths)
        print(f"  厚度（mm）：中位 {np.median(hh):.3f}  p1 "
              f"{np.percentile(hh,1):.3f}  min {hh.min():.4f}  "
              f"max {hh.max():.3f}")
