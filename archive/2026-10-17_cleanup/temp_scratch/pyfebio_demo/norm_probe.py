"""探针：关节面节点的「面积加权法向」质量 vs「指向对侧」方向。

目的：验证代码注释里那句「粗网格上相邻节点的面积加权法向严重不一致（相加几乎抵消）」。
若真抵消 ⇒ 节点法向不可用，必须按 Anderson 用**四边形底 + 顶点法向**重建；
若不抵消 ⇒ 之前是别的原因，可以直接改偏移方向试 hex8。
"""
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(r"D:\Project\climbing_fall_analysis")
sys.path.insert(0, str(ROOT / "scripts" / "ankle_fe"))
d = np.load(ROOT / "temp/thums/thums_lowerlimb.npz")
nid, xyz = d["node_ids"], d["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}

TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
WEDGE_FACES = [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
GROUPS = {"tibia": ["81000700_solid", "81000600_solid", "81000601_solid"],
          "talus": ["81001100_solid", "81001000_solid"]}


def oriented_boundary(keys):
    acc = {}
    for k in keys:
        for row in d[k]:
            nds = [int(x) for x in row]
            uniq = sorted(set(nds))
            if len(uniq) == 4:
                tmpl = TET_FACES
            elif len(uniq) == 8:
                tmpl = HEX_FACES
            elif len(uniq) == 6:
                tmpl = WEDGE_FACES
            else:
                continue
            P = np.array([xyz[idx[v]] for v in uniq], float)
            cen = P.mean(0)
            for f in tmpl:
                if max(f) >= len(uniq):
                    continue
                verts = [uniq[z] for z in f]
                pts = np.array([xyz[idx[v]] for v in verts], float)
                n = np.zeros(3)
                m = len(pts)
                for i in range(m):
                    a, b = pts[i], pts[(i + 1) % m]
                    n += np.cross(a, b)
                fc = pts.mean(0)
                if np.dot(n, fc - cen) < 0:
                    n, verts, pts = -n, verts[::-1], pts[::-1]
                key = tuple(sorted(verts))
                rec = acc.get(key)
                if rec is None:
                    acc[key] = {"verts": tuple(verts), "n": n, "cnt": 1, "pts": pts}
                else:
                    rec["cnt"] += 1
    return [v for v in acc.values() if v["cnt"] == 1]


print("=== 提取外表面 ===")
surf = {}
for g, keys in GROUPS.items():
    f = oriented_boundary(keys)
    surf[g] = f
    nt = sum(1 for x in f if len(x["verts"]) == 3)
    nq = sum(1 for x in f if len(x["verts"]) == 4)
    print(f"  {g:8s} {len(f):>6,} 面  (tri {nt:,} / quad {nq:,})  "
          f"quad 占比 {nq/max(len(f),1)*100:.1f}%")


def node_normal_detail(faces, nodes):
    """返回 (单位法向, 未归一化的合成模长, 贡献面数)。模长归一化前大小 = 抵消程度。"""
    acc = {v: np.zeros(3) for v in nodes}
    cnt = {v: 0 for v in nodes}
    for f in faces:
        pts = f["pts"]
        a, b = pts[1] - pts[0], pts[2] - pts[0]
        fn = np.cross(a, b) * 0.5          # 面积向量（三角形面）
        nrm = np.linalg.norm(fn)
        unit = fn / nrm if nrm > 1e-12 else f["n"] / max(np.linalg.norm(f["n"]), 1e-12)
        for v in f["verts"]:
            if v in acc:
                acc[v] += unit * max(nrm, 1e-12)
                cnt[v] += 1
    out = {}
    for v, s in acc.items():
        ln = np.linalg.norm(s)
        out[v] = (s / ln if ln > 1e-12 else np.zeros(3), ln, cnt[v])
    return out


print("\n=== 关节面法向质量 ===")
fA = surf["tibia"]
nA = set().union(*[set(f["verts"]) for f in fA])
ndet = node_normal_detail(fA, nA)

fB = surf["talus"]
nB = set().union(*[set(f["verts"]) for f in fB])
PB = xyz[np.fromiter((idx[int(v)] for v in sorted(nB)), np.int64, len(nB))]
extra = np.array([f["pts"].mean(0) for f in fB])
QB = np.vstack([PB, extra])
treeQ = cKDTree(QB)

idsA = np.array(sorted(nA), dtype=np.int64)
PA = xyz[np.fromiter((idx[int(v)] for v in idsA), np.int64, len(idsA))]
gap, nq = treeQ.query(PA, k=1)
dirA = QB[nq] - PA
dirA /= np.maximum(np.linalg.norm(dirA, axis=1, keepdims=True), 1e-12)

art = gap < 6.0
print(f"  tibia 外表面节点 {len(nA):,}；关节面节点(gap<6mm) {art.sum():,}")
if art.sum():
    ids = idsA[art]
    Narr = np.array([ndet[int(v)][0] for v in ids])
    mag = np.array([ndet[int(v)][1] for v in ids])
    nfc = np.array([ndet[int(v)][2] for v in ids])
    dA = dirA[art]
    cos = np.einsum("ij,ij->i", Narr, dA)
    typ_area = np.median(mag) / max(np.median(nfc), 1)
    print(f"  贡献面数/节点: 中位 {np.median(nfc):.0f} / min {nfc.min()} / max {nfc.max()}")
    print(f"  合成模长(未归一化) 中位 {np.median(mag):.4f} / min {mag.min():.4f}"
          f"   ⇒ 单面平均面积 {typ_area:.4f} mm²")
    print(f"  抵消判据 |Σ|/(Σ面积) : 中位 {np.median(mag/np.maximum(nfc,1)):.3f}")
    print(f"  法向·指向对侧 cos : 中位 {np.median(cos):+.3f} / "
          f"p10 {np.percentile(cos,10):+.3f} / p90 {np.percentile(cos,90):+.3f}")
    print(f"  夹角 >60° 的节点: {int((cos < 0.5).sum())}/{art.sum()} "
          f"({(cos<0.5).mean()*100:.0f}%)")
    print(f"  夹角 >90° (反向) 的节点: {int((cos < 0).sum())}/{art.sum()} "
          f"({(cos<0).mean()*100:.0f}%)")
    # 关键对比：用节点法向 vs 用 dirA 挤出，偏移面的位置差多少
    th = 1.5
    off_n = PA[art] + Narr * th
    off_d = PA[art] + dA * th
    dd = np.linalg.norm(off_n - off_d, axis=1)
    print(f"\n  两种偏移方向（法向 vs 指向对侧）在同一节点上偏移 1.5mm 后距离差：")
    print(f"    中位 {np.median(dd):.3f} mm / p90 {np.percentile(dd,90):.3f} / "
          f"max {dd.max():.3f} mm")
    # 偏移面是否比 base 面更小（退化迹象）
    print(f"  （若该值接近 1.5mm×2 ⇒ 两个方向几乎相反，法向不可用）")
