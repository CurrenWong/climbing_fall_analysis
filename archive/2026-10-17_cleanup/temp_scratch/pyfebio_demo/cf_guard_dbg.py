"""调试：对侧骨面护栏的有符号距（含内部最近点测试）在**原始**顶节点上是什么。"""
import io
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(r"D:\Project\climbing_fall_analysis")
SRC = io.open(ROOT / "scripts/ankle_fe/cart_patch.py", encoding="utf-8").read()
ns = {"__name__": "cp_head"}
exec(compile(SRC.split("\nsrc = np.load(")[0], "head", "exec"), ns)
guard_of = ns["guard_of"]
newell = ns["newell"]
_tri_closest = ns["_tri_closest"]

raw = np.load(ROOT / "temp/thums/cartilage_mesh.npz", allow_pickle=True)
bids, bxyz = raw["bone_node_ids"], raw["bone_xyz"]
bpos = {int(v): bxyz[i] for i, v in enumerate(bids)}
ns["src"] = raw
ns["bpos"] = bpos

CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


for A, B in (("tibiotalar_tibia", "tibiotalar_talus"),
             ("tibiotalar_talus", "tibiotalar_tibia")):
    pts, nrm, tri, trin = guard_of(B, bpos, newell)
    print(f"\n{ A } 的护栏 = {B} 的底面: 三角 {len(tri)}")
    # 护栏三角法向与 A 的原始抬升方向的关系（应大致反向/同向？）
    N = np.asarray(CM1[f"{A}_nodes"], float)
    E = np.asarray(CM1[f"{A}_elems"], np.int64)
    off = np.asarray(CM1[f"{A}_is_offset"], bool)
    C = np.array([[ix(x) for x in row] for row in E])
    Bn, Tn = C[:, :3], C[:, 3:]
    d0 = N[Tn] - N[Bn]
    print(f"   A 原始抬升方向均值 {np.round(d0.reshape(-1,3).mean(0),3)}  "
          f"|d| 中位 {np.median(np.linalg.norm(d0,axis=2)):.3f}")
    print(f"   护栏法向均值 {np.round(trin.mean(0),3)}")
    # 原始顶节点 -> 护栏的有符号距（内部最近点）
    tops = N[off]
    gt = cKDTree(tri.mean(axis=1))
    k = min(6, len(tri))
    _, iq = gt.query(tops, k=k)
    R = np.broadcast_to(tops[:, None, :], (len(tops), k, 3)).reshape(-1, 3)
    T3 = tri[iq.reshape(-1)]
    q, inner = _tri_closest(R, T3[:, 0], T3[:, 1], T3[:, 2])
    nr = trin[iq.reshape(-1)]
    dd = np.einsum("ij,ij->i", nr, R - q)
    dd = dd.reshape(len(tops), k)
    inn = inner.reshape(len(tops), k)
    # 取「内部且最小」的那个
    use = np.where(inn, dd, np.inf)
    dsel = use.min(axis=1)
    n_int = int(inn.any(axis=1).sum())
    print(f"   顶节点 {len(tops)} 个：有内部最近三角的 {n_int} 个；"
          f"其有符号距 min {np.nanmin(dsel):+.4f} 中位 {np.nanmedian(dsel):+.4f} "
          f"max {np.nanmax(dsel):+.4f}")
    print(f"      <0 的 {int((dsel<0).sum())} 个  <0.01 的 {int((dsel<0.01).sum())} 个")
