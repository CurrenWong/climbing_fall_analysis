"""软骨**面内**红色加密（骨网格保持粗网格），为 tied-elastic 绑定做准备。

── 为什么需要它 ────────────────────────────────────────────────────────
之前的 P1d 路线是「整根骨一起加密」：关节面节点 73 -> 692（9.5x），
证明接触区确实欠分辨；但骨 hex8 加密后出现**负 Jacobian**（rc=1），
且加密后软骨体积涨 2.3 倍 ⇒ 加密骨不可用。

本脚本把两件事解耦：
  * **骨**：完全不动，保持已验证的粗网格（rc=0 / 17 s / 436.8 N）
  * **软骨**：面内任意加密，底面脱离骨节点

脱离开的底面用 FEBio 的 ``tied-elastic`` 接触（手册 §3.14.6：「连接非共形
网格，强制界面两侧位移连续」）绑回粗骨。primary = 细软骨底面，
secondary = 粗骨关节面（primary 节点投影到 secondary 三角/四边形面上）。

── 编码（对 L=0 向后兼容）─────────────────────────────────────────────
``{j}_nodes``  (n,3)   软骨专属节点（旧文件里叫 ``_offsets``）
``{j}_elems``  (M,8)   >0 = 骨节点 ID；<0 = ``-(k+1)`` -> ``_nodes[k]``
``{j}_is_offset`` (n,)  该节点是否落在**外表面**（旧文件里全是 True）
``{j}_tie_faces``  (F,4) 细软骨底面（绑定的 primary）
``{j}_bone_faces`` (F,4) 粗骨关节面（绑定的 secondary，全为骨节点 ID）

L=0 且 DUPLICATE=0 时与旧 ``cartilage_mesh.npz`` 语义完全一致。
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
THUMS = ROOT / "temp/thums"

LEVEL = int(os.environ.get("CART_L", "1"))          # 加密级数
# 1 = 底面节点全部复制成软骨专属节点（纯 mortar 界面，L>0 时必须开）
DUPLICATE = os.environ.get("CART_DUP", "1") == "1"

JOINTS = ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]

src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bone_ids = src["bone_node_ids"]
bone_xyz = src["bone_xyz"]
bpos = {int(v): bone_xyz[i] for i, v in enumerate(bone_ids)}

out = {"bone_node_ids": bone_ids, "bone_xyz": bone_xyz}
report: dict[str, dict] = {}

for j in JOINTS:
    elems = np.asarray(src[f"{j}_elems"], dtype=np.int64)
    old_offsets = np.asarray(src[f"{j}_offsets"], float)

    # ---- 解码成 (base_refs, outer_refs)，并建位置表 ----
    cells = []
    for row in elems:
        v = [int(x) for x in row if int(x) != 0]
        m = len(v) // 2
        base, outer = v[:m], v[m:]
        cells.append((base, outer))
        assert m in (3, 4), f"{j}: 只处理 tri3/quad4 底面，实际 {len(v)} 节点"

    nodes: list[np.ndarray] = []       # 软骨专属节点
    node_tag: list[bool] = []          # True = 外表面

    def pos_of(ref):
        if ref > 0:
            return bpos[ref]
        return nodes[-ref - 1]

    def add_node(p, is_off):
        nodes.append(np.asarray(p, float))
        node_tag.append(bool(is_off))
        return -(len(nodes))            # 负引用 = -(k+1)

    def pos_map(refs):
        return {r: (bpos[r] if r > 0 else nodes[-r - 1]) for r in refs}

    # ---- 把原来的 offset 节点搬进 nodes（保持它们的局部下标）----
    for k in range(len(old_offsets)):
        nodes.append(old_offsets[k])
        node_tag.append(True)

    # ---- DUPLICATE：底面节点复制一份（变成纯 mortar 界面）----
    # ★ 必须按骨节点 ID 去重：同一骨节点被多个面片共用，若每个面片各建一个
    #   本地副本，相邻单元的底面就**裂开**（内部自由面），网格非共形。
    dup: dict[int, int] = {}
    base_of: list[list[int]] = []
    for base, outer in cells:
        nb = []
        for r in base:
            if r > 0:
                if r not in dup:
                    dup[r] = add_node(bpos[r], False)
                nb.append(dup[r])
            else:
                nb.append(r)
        base_of.append(nb)

    # ---- 记录粗骨关节面（tie 的 secondary 面）----
    bone_faces = [list(b) for b, _ in cells]

    n0_cell, n0_base = len(cells), len({r for b in base_of for r in b})
    # ★ 必须在这里就把 cells 换成「复制后的底面」，否则 L=0 时 cells 仍是
    #   原始 cells（底面直接共用骨节点），导出的 npz 里软骨与骨**仍然共形**，
    #   而 _tie_faces 用的却是复制后的节点 —— tie 会绑在一组不属于软骨的节点上，
    #   模型退化成「原共形模型 + 死 tie」，测出来的结果毫无意义。
    cells = [(list(b), list(o)) for b, (_, o) in zip(base_of, cells)]
    coarse_cells = [tuple(c) for c in cells]

    # ---- 逐级红色加密 ----
    emap: dict[tuple, int] = {}
    for _lv in range(LEVEL):
        new_cells = []
        for base, outer in zip(base_of, [c[1] for c in cells]):
            bp, op = pos_map(base), pos_map(outer)

            def mid(refs, pts, i, j, is_off):
                k = (min(refs[i], refs[j]), max(refs[i], refs[j]))
                if k not in emap:
                    emap[k] = add_node((pts[refs[i]] + pts[refs[j]]) / 2.0, is_off)
                return emap[k]

            def cen(refs, pts, is_off):
                k = tuple(sorted(refs))
                if k not in emap:
                    emap[k] = add_node(np.mean([pts[r] for r in refs], axis=0), is_off)
                return emap[k]

            nb = len(base)
            bm = [mid(base, bp, i, (i + 1) % nb, False) for i in range(nb)]
            bcen = cen(base, bp, False) if nb == 4 else None
            om = [mid(outer, op, i, (i + 1) % nb, True) for i in range(nb)]
            ocen = cen(outer, op, True) if nb == 4 else None

            b = base + bm
            o = outer + om
            if len(base) == 4:
                # quad -> 4 quad（b[4..7] = 4 条棱中点，bcen = 面心）
                new_cells += [
                    ([b[0], b[4], bcen, b[7]], [o[0], o[4], ocen, o[7]]),
                    ([b[4], b[1], b[5], bcen], [o[4], o[1], o[5], ocen]),
                    ([bcen, b[5], b[2], b[6]], [ocen, o[5], o[2], o[6]]),
                    ([b[7], bcen, b[6], b[3]], [o[7], ocen, o[6], o[3]]),
                ]
            else:
                # tri -> 4 tri（b[3..5] = 3 条棱中点；tri 没有面心）
                # ★★ 必须是 1→4，**不能只写 3 个角三角形**！
                #   3 个角各占 1/4 面积，剩下的 1/4 是中央那个 (m01,m12,m20)；
                #   漏掉它 = 静默丢掉 25% 的体积，而且**不会报负体积**
                #   （3 个子楔形各自都是合法正体积）—— 我第一版就漏了，
                #   还把 −14.1% 的体积偏差错归因成"非凸单元 fan 分解偏差"。
                new_cells += [
                    ([b[0], b[3], b[5]], [o[0], o[3], o[5]]),
                    ([b[3], b[1], b[4]], [o[3], o[1], o[4]]),
                    ([b[5], b[4], b[2]], [o[5], o[4], o[2]]),
                    ([b[3], b[4], b[5]], [o[3], o[4], o[5]]),   # ★ 中央，别漏
                ]
        cells = new_cells
        base_of = [c[0] for c in cells]

    # ---- 体积校验（加密是保体积的重划分）----
    def vol(c):
        b, o = c[0], c[1]
        P = np.array([pos_of(r) for r in b + o])
        if len(b) == 4:
            ts = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
                  (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]      # hex8
        else:
            ts = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]      # wedge6
        t = 0.0
        for q in ts:
            a, b2, c2, e = (P[i] for i in q)
            t += np.dot(np.cross(b2 - a, c2 - a), e - a) / 6.0
        return t

    v_new = sum(vol(c) for c in cells)
    v_base = sum(vol(c) for c in coarse_cells)
    neg_base = sum(1 for c in coarse_cells if vol(c) <= 0)
    neg_new = sum(1 for c in cells if vol(c) <= 0)
    print(f"    负体积单元 {neg_base:,} -> {neg_new:,}")
    print(f"    体积(同分解法) {v_base:,.2f} -> {v_new:,.2f} mm^3   "
          f"偏差 {(v_new - v_base) / v_base * 100:+.4f}%")

    # 旧体积直接用 cartilage.json 的报告值更省事
    v_ref = None
    try:
        import json
        v_ref = json.loads((THUMS / "cartilage.json").read_text(encoding="utf-8")) \
            ["joints"][j]["volume_mm3"]
    except Exception:
        pass

    tie_faces = [list(b) for b in base_of]
    # 混合 hex8 / wedge6 => 补 0 成一维等长（0 不是合法 ref，读方按 !=0 过滤）
    W = max(len(b) + len(o) for b, o in cells)
    pad = np.zeros((len(cells), W), dtype=np.int64)
    for i, (b, o) in enumerate(cells):
        row = b + o
        pad[i, :len(row)] = row
    def _rect(rows):
        w = max(len(r) for r in rows)
        a = np.zeros((len(rows), w), dtype=np.int64)
        for i, r in enumerate(rows):
            a[i, :len(r)] = r
        return a

    out[f"{j}_nodes"] = np.asarray(nodes, float)
    out[f"{j}_is_offset"] = np.asarray(node_tag, bool)
    out[f"{j}_elems"] = pad
    out[f"{j}_tie_faces"] = _rect(tie_faces)
    out[f"{j}_bone_faces"] = _rect(bone_faces)

    n_off = sum(node_tag)
    n_bnd = sum(1 for t in node_tag if not t)
    report[j] = {
        "cells": [n0_cell, len(cells)],
        "base_nodes": [n0_base, len({r for b in base_of for r in b})],
        "cart_nodes": [n_off + n_bnd, len(nodes)],
        "volume": v_new,
        "volume_ref": v_ref,
        "dup": DUPLICATE,
    }
    print(f"\n[{j}]")
    print(f"  单元       {n0_cell:>7,} -> {len(cells):>8,}")
    print(f"  底面节点   {n0_base:>7,} -> "
          f"{len({r for b in base_of for r in b}):>8,}")
    print(f"  软骨节点   {n_off + n_bnd:>7,} -> {len(nodes):>8,}"
          f"   (底面 {n_bnd:,} / 外表面 {n_off:,})")
    print(f"  体积       {v_new:,.2f} mm^3"
          + (f"   参考(旧) {v_ref:,.2f}  偏差 {(v_new - v_ref) / v_ref * 100:+.4f}%"
             if v_ref else ""))
    neg = sum(1 for c in cells if vol(c) <= 0)
    print(f"  负体积单元 {neg:,}")

dst = THUMS / f"cartilage_mesh_t{LEVEL}{'d' if DUPLICATE else 's'}.npz"
np.savez_compressed(dst, **out)
print(f"\n写出 {dst}  ({dst.stat().st_size / 1e6:.2f} MB)")
print(f"LEVEL={LEVEL} DUPLICATE={DUPLICATE}")
