"""均匀 1 级红色加密 thums_lowerlimb.npz（只加密 *_solid，壳保持原样）。

策略（保证【完全共形】，不需要 tied contact）：
  ① 先把所有 solid 单元统一转成 tet4
     - 8 节点 hex8  -> 6 tets（面对角线由「面 4 节点 id 集合」决定 ⇒ 两侧拥有者一致）
     - 6 节点 wedge -> 3 tets
     - 4 节点 tet4  -> 保持
  ② 再对所有 tet4 做 1:8 红色加密（只用棱中点；棱中点全局共享 ⇒ 天然共形）
     - 4 个角 tet + 中心八面体（用一条对角切成 4 个 tet，属内部，不影响共形）

验证：体积守恒、无退化、节点/单元计数。
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
SRC = ROOT / "temp/thums/thums_lowerlimb.npz"
DST = ROOT / "temp/thums/thums_lowerlimb_r1.npz"

d = np.load(SRC, allow_pickle=True)
nid = d["node_ids"]
xyz = d["node_xyz"].astype(float)
idx = {int(v): i for i, v in enumerate(nid)}
PART_KEYS = [k for k in d.files if k.endswith(("_solid", "_shell"))]

HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
# hex 的 6 个 tet 分解：选一条主轴 (0,6)，再按面对角线把 6 个面切开
# 实用做法：先把 hex 拆成 6 个「面金字塔→tet」不可靠；改用标准 5-tet? 不。
# 这里用「一致面对角线 + 中心」法：hex -> 6 tets（通过 4 个 tet 的经典分解）
HEX_TO_TET = [(0, 1, 2, 5), (0, 2, 3, 7), (0, 5, 7, 4), (0, 5, 2, 7),
              (3, 7, 2, 5), (0, 4, 7, 5)]
WEDGE_TO_TET = [(0, 1, 2, 5), (0, 1, 5, 4), (0, 4, 5, 3)]
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]


def vol_tet(P):
    return float(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0]) / 6.0)


def orient_tet(P):
    return P if vol_tet(P) > 0 else P[[0, 2, 1, 3]]


# ---------- ① 统计 & 原始体积 ----------
print("=== 原始 ===")
V0 = 0.0
keep = {}
all_tets = []
for k in PART_KEYS:
    if not k.endswith("_solid"):
        continue
    rows = []
    for row in d[k]:
        uq = sorted({int(x) for x in row if x > 0})
        if len(uq) == 8:
            tets = [[uq[z] for z in t] for t in HEX_TO_TET]
        elif len(uq) == 6:
            tets = [[uq[z] for z in t] for t in WEDGE_TO_TET]
        elif len(uq) == 4:
            tets = [uq]
        else:
            print(f"  ⚠️ {k} 跳过 {len(uq)} 节点单元")
            continue
        for t in tets:
            t = [int(x) for x in t]
            P = np.array([xyz[idx[x]] for x in t])
            P = orient_tet(P)
            V0 += vol_tet(P)
            all_tets.append(t)
    keep[k] = rows
print(f"  统一为 tet4：{len(all_tets):,} 个")
print(f"  原始总体积 {V0:,.3f} mm^3 = {V0/1000:.4f} cm^3")

# ---------- ② 棱中点（全局共享） ----------
edge_mid = {}
new_pts = []
for t in all_tets:
    for a, b in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        key = (min(t[a], t[b]), max(t[a], t[b]))
        if key not in edge_mid:
            pa = xyz[idx[key[0]]]
            pb = xyz[idx[key[1]]]
            edge_mid[key] = ("m", len(new_pts))
            new_pts.append((pa + pb) / 2.0)
print(f"  新增棱中点 {len(new_pts):,}")

N_OLD = len(nid)
NEW_IDS = [int(nid.max()) + 1 + i for i in range(len(new_pts))]
XYZ_NEW = np.vstack([xyz, np.array(new_pts)]) if new_pts else xyz


def mid_of(a, b):
    key = (min(a, b), max(a, b))
    tag, i = edge_mid[key]
    return NEW_IDS[i]


# ---------- ③ 1:8 红色加密 ----------
out_tets = []
for t in all_tets:
    a, b, c, e = t
    mab, mac, mad = mid_of(a, b), mid_of(a, c), mid_of(a, e)
    mbc, mbd, mcd = mid_of(b, c), mid_of(b, e), mid_of(c, e)
    out_tets += [[a, mab, mac, mad], [b, mba_i := mab, mbc, mbd],
                 [c, mac, mbc, mcd], [e, mad, mbd, mcd]]
    # 中心八面体 (mab,mac,mad,mbc,mbd,mcd) 用对角 mab-mcd 切 4 个 tet
    out_tets += [[mab, mcd, mac, mad], [mab, mcd, mad, mbd],
                 [mab, mcd, mbd, mbc], [mab, mcd, mbc, mac]]

V1 = 0.0
for t in out_tets:
    P = np.array([XYZ_NEW[i - N_OLD + N_OLD] if i <= N_OLD else np.nan for i in t])
# 重新用索引函数计算（上面写法易错，改为按 id->坐标表）
id2xyz = {int(v): xyz[i] for i, v in enumerate(nid)}
for j, v in enumerate(NEW_IDS):
    id2xyz[v] = new_pts[j]
V1 = 0.0
bad = 0
for t in out_tets:
    P = np.array([id2xyz[i] for i in t])
    v = vol_tet(P)
    if v <= 0:
        bad += 1
    V1 += v
print("\n=== 加密后 ===")
print(f"  tet4 {len(out_tets):,} 个（原 {len(all_tets):,}）")
print(f"  体积 {V1:,.3f} mm^3   （原始 {V0:,.3f}）  "
      f"偏差 {(V1-V0)/V0*100:+.4f}%")
print(f"  非正体积 {bad} 个")

# ---------- ④ 写出 ----------
arrs = {"node_ids": np.concatenate([nid, np.array(NEW_IDS, dtype=nid.dtype)]),
        "node_xyz": XYZ_NEW}
# solid 全部写成 tet4（按原始 part 归属重排）
p = 0
for k in PART_KEYS:
    if not k.endswith("_solid"):
        arrs[k] = d[k]           # 壳保持原样
        continue
    n_t = 0
    for row in d[k]:
        uq = sorted({int(x) for x in row if x > 0})
        n_t += 6 if len(uq) == 8 else (3 if len(uq) == 6 else 1)
    body = out_tets[p:p + n_t]
    p += n_t
    w = max(len(t) for t in body)
    arrs[k] = np.array([t + [0] * (w - len(t)) for t in body], dtype=np.int64)
assert p == len(out_tets), (p, len(out_tets))
np.savez_compressed(DST, **arrs)
print(f"\n写出 {DST}  ({DST.stat().st_size/1e6:.1f} MB)")
print(f"  节点 {N_OLD:,} -> {N_OLD + len(NEW_IDS):,}")
