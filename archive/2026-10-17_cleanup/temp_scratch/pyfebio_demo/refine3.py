"""均匀 1 级加密（干净版）：hex8->8 hex8；tet4->8 tet4；quad4 壳->4 quad4。

新增节点全部全局共享（棱中点/面心/体心）⇒ hex-hex、tet-tet 面匹配自动保持共形。
红色加密保留全部原节点 ⇒ 未加密的壳仍挂在原节点上。
"""
from pathlib import Path
import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
SRC = ROOT / "temp/thums/thums_lowerlimb.npz"
DST = ROOT / "temp/thums/thums_lowerlimb_r1.npz"

d = np.load(SRC, allow_pickle=True)
nid = np.asarray(d["node_ids"], dtype=np.int64)
xyz = d["node_xyz"].astype(float)
src = {int(v): xyz[i] for i, v in enumerate(nid)}

# FEBio hex 局部顺序：底 0,1,2,3(CCW) / 顶 4,5,6,7（4 在 0 上方）
HEX_UVW = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
           (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
HEX_E = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
         (0, 4), (1, 5), (2, 6), (3, 7)]
HEX_F = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
         (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
TET_E = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
QUAD_E = [(0, 1), (1, 2), (2, 3), (3, 0)]

_new = []


def alloc(p):
    _new.append(np.asarray(p, float))
    return int(nid.max()) + len(_new)


def _coord(x):
    """旧 id -> 原坐标；新 id -> _new。"""
    return src[x] if x in src else _new[x - int(nid.max()) - 1]


_em, _fm, _cm = {}, {}, {}


def E(a, b):
    k = (min(a, b), max(a, b))
    if k not in _em:
        _em[k] = alloc((src[k[0]] + src[k[1]]) / 2.0)
    return _em[k]


def F(ids):
    k = tuple(sorted(ids))
    if k not in _fm:
        _fm[k] = alloc(np.mean([src[x] for x in k], axis=0))
    return _fm[k]


def C(ids):
    k = tuple(sorted(ids))
    if k not in _cm:
        _cm[k] = alloc(np.mean([src[x] for x in k], axis=0))
    return _cm[k]


def detJ_hex_center(P):
    """中心 detJ 符号 = 六面体朝向。u,v,w 定义在 [0,1]。"""
    u = v = w = 0.5
    du = np.array([-(1-v)*(1-w), (1-v)*(1-w), v*(1-w), -v*(1-w),
                   -(1-v)*w, (1-v)*w, v*w, -v*w], float)
    dv = np.array([-(1-u)*(1-w), -u*(1-w), u*(1-w), (1-u)*(1-w),
                   -(1-u)*w, -u*w, u*w, (1-u)*w], float)
    dw = np.array([-(1-u)*(1-v), -u*(1-v), -u*v, -(1-u)*v,
                   (1-u)*(1-v), u*(1-v), u*v, (1-u)*v], float)
    J = np.stack([du @ P, dv @ P, dw @ P], axis=0)
    return float(np.linalg.det(J))


def fix_hex(conn, coord):
    """朝向不对就交换底面/顶面（4,5,6,7,0,1,2,3）——合法翻转，不改形状。"""
    P = np.array([coord(x) for x in conn])
    if detJ_hex_center(P) < 0:
        return conn[4:] + conn[:4]
    return conn


def hvol(P):
    """hex8 体积 = 6 个 tet 的带符号和。
    ★ 分解必须【不退化且无缝铺满】：(0,4,7,5) 这类是共面的（体积 0），不能用。
    已验证正确：绕体对角线 0-6 的扇形分解。"""
    v = 0.0
    for t in [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
              (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]:
        A, B, Cc, D = P[t[0]], P[t[1]], P[t[2]], P[t[3]]
        v += np.dot(np.cross(B - A, Cc - A), D - A) / 6.0
    return abs(v)


def tvol(P):
    return abs(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0]) / 6.0)


out, stat = {}, []
V0 = V1 = 0.0
for k in d.files:
    if not k.endswith(("_solid", "_shell")):
        continue
    rows_out = []
    n0 = n1 = 0
    for row in d[k]:
        uq = sorted({int(x) for x in row if x > 0})
        if len(uq) == 8:                       # ---- hex8 -> 8 hex8
            # ★ 先规范化母单元朝向：npz 里的顺序未必符合 FEBio 约定
            #   （不符合 ⇒ 棱/面识别错 ⇒ 子单元铺不满母单元 ⇒ 体积丢 ~9%）
            uq = fix_hex(uq, _coord)
            P = np.array([src[x] for x in uq])
            V0 += hvol(P)
            pmap = {HEX_UVW[i]: uq[i] for i in range(8)}
            for a, b in HEX_E:
                pa, pb = HEX_UVW[a], HEX_UVW[b]
                mid = tuple((pa[i] + pb[i]) / 2.0 for i in range(3))
                pmap[mid] = E(uq[a], uq[b])
            for f in HEX_F:
                cen = tuple(np.mean([HEX_UVW[z] for z in f], axis=0))
                pmap[(round(cen[0], 6), round(cen[1], 6), round(cen[2], 6))] = \
                    F([uq[z] for z in f])
            pmap[(0.5, 0.5, 0.5)] = C(uq)
            for i in (0.0, 0.5):
                for j in (0.0, 0.5):
                    for m in (0.0, 0.5):
                        vs = [(i + a * 0.5, j + b * 0.5, m + c * 0.5)
                              for (a, b, c) in HEX_UVW]
                        conn = [pmap[(round(u, 6), round(v, 6), round(w, 6))]
                                for (u, v, w) in vs]
                        rows_out.append(fix_hex(conn, _coord))   # ★ 朝向修正
                        n1 += 1
            n0 += 1
        elif len(uq) == 4:                     # ---- tet4 -> 8 tet4
            P = np.array([src[x] for x in uq])
            V0 += tvol(P)
            mab, mac, mad = E(uq[0], uq[1]), E(uq[0], uq[2]), E(uq[0], uq[3])
            mbc, mbd, mcd = E(uq[1], uq[2]), E(uq[1], uq[3]), E(uq[2], uq[3])
            rows_out += [[uq[0], mab, mac, mad], [uq[1], mab, mbc, mbd],
                         [uq[2], mac, mbc, mcd], [uq[3], mad, mbd, mcd],
                         [mab, mcd, mac, mad], [mab, mcd, mad, mbd],
                         [mab, mcd, mbd, mbc], [mab, mcd, mbc, mac]]
            n0 += 1; n1 += 8
        elif len(uq) == 6:                     # ---- wedge6：保持原样（仅 ~10 个）
            rows_out.append(uq); n0 += 1; n1 += 1
        elif len(uq) == 4 and k.endswith("_shell"):
            pass
        else:
            rows_out.append(uq); n0 += 1; n1 += 1
    # quad4 壳（4 节点，闭环）
    if k.endswith("_shell"):
        q_out = []
        for row in d[k]:
            uq = [int(x) for x in row if int(x) > 0]
            if len(uq) != 4:
                q_out.append(uq); continue
            cen = F(list(uq))
            e = [E(uq[a], uq[b]) for a, b in QUAD_E]
            for i in range(4):
                q_out.append([uq[i], e[i], cen, e[(i + 3) % 4]])
        rows_out = q_out
        n1 = len(q_out)
    out[k] = rows_out
    stat.append((k, n0, n1))

# 体积校验
def _coord(x):
    """旧 id -> src；新 id -> _new（★ 不能用 dict.get 的默认参数，它总会被求值）。"""
    if x in src:
        return src[x]
    return _new[x - int(nid.max()) - 1]


V1 = 0.0
for k in d.files:
    if not k.endswith("_solid") or k not in out:
        continue
    for row in out[k]:
        uq = sorted({int(x) for x in row if x > 0})
        if len(uq) == 4:
            V1 += tvol(np.array([_coord(x) for x in uq]))
        elif len(uq) == 8:
            V1 += hvol(np.array([_coord(x) for x in uq]))

print("=== 加密 ===")
for k, a, b in stat:
    print(f"  {k:22s} {a:>6,} -> {b:>8,}")
print(f"\n体积 {V0:,.3f} -> {V1:,.3f} mm^3   偏差 {(V1-V0)/V0*100:+.4f}%")
print(f"节点 {len(nid):,} -> {len(nid)+len(_new):,}"
      f"（新增 棱{len(_em):,} 面{len(_fm):,} 体{len(_cm):,}）")

arrs = {"node_ids": np.concatenate([nid, np.array([int(nid.max()) + 1 + i
                                                   for i in range(len(_new))])]),
        "node_xyz": np.vstack([xyz, np.array(_new)]) if _new else xyz}
for k in d.files:
    if k in ("node_ids", "node_xyz"):
        continue          # ★ 必须跳过：否则会把上面拼好的新节点表覆盖回原始值
    if k in out:
        w = max(len(r) for r in out[k])
        arrs[k] = np.array([list(r) + [0] * (w - len(r)) for r in out[k]],
                           dtype=np.int64)
    else:
        arrs[k] = d[k]
np.savez_compressed(DST, **arrs)
print(f"写出 {DST}  ({DST.stat().st_size/1e6:.1f} MB)")
