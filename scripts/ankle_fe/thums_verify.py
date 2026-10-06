"""验证 v2：修正 v1 的 bug，并对体积做**两种独立算法交叉验证**。

v1 的 bug：
  ① 退化统计把 (n,8,3) reshape 成 (n,24) -> 唯一值数变成 12/24（物理上不可能，8 节点最多 8 个唯一）
  ② json 落盘时 numpy.int64 作 key -> TypeError

本脚本新增：
  ③ 用「边界面积分（散度定理）」独立算体积，与「hex 6-四面体分解」比对
  ④ 合成单位立方体自检，先证明分解公式本身是对的
  ⑤ CORT 与 SPON 共享节点数 -> 判断两者是「共形贴合」还是「各自独立可能重叠」
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

NPZ = Path(r"D:\Project\climbing_fall_analysis\temp\thums\thums_lowerlimb.npz")

HEX_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
            (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
# hex 的 6 个面（外法向一致的绕向）：底 0-1-2-3 反，顶 4-5-6-7 正，四个侧面
HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def tet_signed_vol(p: np.ndarray) -> np.ndarray:
    return np.einsum("ij,ij->i",
                     np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                     p[:, 3] - p[:, 0]) / 6.0


def self_test() -> bool:
    """单位立方体 (0,0,0)-(1,1,1)：分解法应给出 1.0"""
    h = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)[None, ...]
    s = sum(tet_signed_vol(h[:, list(t), :]).sum() for t in HEX_TETS)
    # 边界面积分
    b = boundary_volume(h)
    print(f"  自检 单位立方体：分解法 {s:.10f}  边界积分 {b:.10f}  (期望 1.0)")
    return abs(s - 1.0) < 1e-9 and abs(b - 1.0) < 1e-9


def boundary_volume(coords: np.ndarray) -> float:
    """coords: (n,8,3) 所有 hex。对边界三角面做散度定理积分 V = (1/6)Σ (a×b)·c"""
    faces: dict[tuple[int, ...], int] = {}
    tri_list = []
    for hexc in coords:
        for f in HEX_FACES:
            quad = [tuple(np.round(hexc[i], 9)) for i in f]
            key = tuple(sorted(quad))
            if key in faces:
                faces[key] += 1
            else:
                faces[key] = 1
                tri_list.append((hexc[f[0]], hexc[f[1]], hexc[f[2]]))
                tri_list.append((hexc[f[0]], hexc[f[2]], hexc[f[3]]))
    # 重新只对「只出现一次」的面做积分
    total = 0.0
    seen: dict[tuple[int, ...], int] = {}
    for hexc in coords:
        for f in HEX_FACES:
            quad = [tuple(np.round(hexc[i], 9)) for i in f]
            key = tuple(sorted(quad))
            seen[key] = seen.get(key, 0) + 1
    for hexc in coords:
        for f in HEX_FACES:
            quad = [tuple(np.round(hexc[i], 9)) for i in f]
            if seen[tuple(sorted(quad))] != 1:
                continue
            for a, b, c in ((f[0], f[1], f[2]), (f[0], f[2], f[3])):
                p0, p1, p2 = hexc[a], hexc[b], hexc[c]
                total += float(np.dot(np.cross(p0, p1), p2))
    return total / 6.0


def main() -> None:
    print("=== 算法自检 ===")
    ok = self_test()
    print(f"  自检通过: {ok}\n")

    d = np.load(NPZ)
    nid = d["node_ids"]
    xyz = d["node_xyz"]
    idx = {int(v): i for i, v in enumerate(nid)}

    parts = {k: (int(k.rsplit("_", 1)[0]), k.rsplit("_", 1)[1])
             for k in d.files if k not in ("node_ids", "node_xyz")}

    # 部件 -> 唯一节点集合（用于 CORT/SPON 共享检查）
    node_sets: dict[str, set[int]] = {}

    print("① 单元形状 / 体积（两种算法交叉验证）")
    print(f"{'pid':>10s} {'kind':6s} {'elem':>7s} {'uniqN':>6s} "
          f"{'degen':>6s} {'nuniq_hist':>28s} {'vol_decomp':>12s} {'vol_bnd':>12s} {'diff%':>7s}")
    results: dict[str, dict] = {}
    for key, (pid, kind) in sorted(parts.items(), key=lambda kv: kv[1][0]):
        arr = d[key]
        row2i = np.vectorize(lambda x: idx.get(int(x), -1))
        flat = row2i(arr)
        coords = xyz[flat]
        node_sets[key] = set(np.unique(arr).tolist())

        if kind == "solid" and arr.shape[1] == 8:
            v = coords
            dec = float(sum(tet_signed_vol(v[:, list(t), :]).sum() for t in HEX_TETS))
            bnd = boundary_volume(v)
            ndist = np.array([len(set(r)) for r in arr.tolist()])
            hist = {str(int(k)): int(v2) for k, v2 in zip(*np.unique(ndist, return_counts=True))}
            n_deg = int((ndist < 8).sum())
            diff = abs(dec - bnd) / max(abs(bnd), 1e-12) * 100
        else:
            dec = bnd = math.nan if False else float("nan")
            n_deg = -1
            hist = {}
            diff = float("nan")

        print(f"{pid:>10d} {kind:6s} {len(arr):>7,} {len(node_sets[key]):>6,} "
              f"{n_deg:>6,} {str(hist):>28s} {dec:>12.1f} {bnd:>12.1f} {diff:>7.2f}")
        results[key] = {
            "pid": pid, "kind": kind, "n_elem": int(len(arr)),
            "n_uniq_nodes": len(node_sets[key]),
            "n_degenerate": n_deg, "uniq_nodes_per_elem_hist": hist,
            "vol_decomp_raw": dec, "vol_boundary_raw": bnd,
            "vol_diff_pct": diff,
            "bbox_min": coords.reshape(-1, 3).min(0).tolist(),
            "bbox_max": coords.reshape(-1, 3).max(0).tolist(),
        }

    # ── CORT/SPON 共享节点 ──
    print("\n② CORT 与 SPON 是否共形贴合（共享节点数）")
    cfg = [("81000900_solid", ["81000801_solid"], "腓骨"),
           ("81000700_solid", ["81000600_solid", "81000601_solid"], "胫骨"),
           ("81001100_solid", ["81001000_solid"], "距骨"),
           ("81001300_solid", ["81001200_solid"], "跟骨")]
    for cort, spon_list, name in cfg:
        if cort not in node_sets:
            continue
        c = node_sets[cort]
        for s in spon_list:
            if s not in node_sets:
                continue
            sh = len(c & node_sets[s])
            print(f"  {name:4s} {cort:16s} ∩ {s:16s} = {sh:5d} 个共享节点 "
                  f"(cort {len(c)}, spon {len(node_sets[s])})"
                  f"  {'✅ 共形' if sh > 0 else '⚠️ 不共形/独立'}")

    # ── 单位判定 ──
    print("\n③ 单位判定（按解剖长度）")
    for key, mid_cm in (("81000100_solid", "股骨 ~45"), ("81000700_solid", "胫骨 ~33"),
                        ("81000900_solid", "腓骨 ~34"), ("81001100_solid", "距骨 ~5"),
                        ("81001300_solid", "跟骨 ~7")):
        if key not in results:
            continue
        lo = np.array(results[key]["bbox_min"])
        hi = np.array(results[key]["bbox_max"])
        L = float((hi - lo).max())
        print(f"  {key:18s} 最长轴 {L:8.2f} raw -> mm 解释 {L/10:6.2f} cm "
              f"(预期 {mid_cm})   {'✅' if 0.5 < L/10 < 60 else '❌'}")

    (NPZ.parent / "verify.json").write_text(
        json.dumps({"self_test_ok": ok, "parts": results},
                   ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n写出 {NPZ.parent/'verify.json'}")


import math  # noqa: E402
main()
