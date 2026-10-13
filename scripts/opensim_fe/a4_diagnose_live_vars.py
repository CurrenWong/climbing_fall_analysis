# -*- coding: utf-8 -*-
"""A4 可测性诊断：PCA 到底在几个「活变量」上做？

背景：`s4_alignment_bc2.py:257-282` 的 a4_from_matrix 里，
  - 长骨：sig[cn] = ends_stat * (loads/fe_solve_n)  ⇒ 随 h 变
  - 非长骨：sig[cn] = fm["sigma_vm_mpa"]["median"][bone]  ⇒ **与 h 无关（常数）**
⇒ 若多数部件是常数，中心化后整列为 0 ⇒ 对协方差/PCA 零贡献，
   而并行分析 null 仍按 9 变量算 ⇒ **门槛与数据不匹配**。

本脚本**只读冻结件**，回答三个问题：
  Q1 9 个部件里几个有非零方差（活变量）？
  Q2 只在活变量上做 PCA，PC2 比是多少？
  Q3 并行分析 null 对「活变量个数」重算后是多少？PC2 还超不超 null？
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import paper_alignment as pa  # noqa: E402

FM = ROOT / "results" / "opensim_fe" / "fracture_matrix_bc2.json"
YI = ROOT / "results" / "opensim_fe" / "s4_alignment_bc2.json"

REGION_NAMES = list(pa.REGION_NAMES)
LONG = set(getattr(pa, "LONG_BONES", []) or [])


def sig_matrix(yi: dict, heights: np.ndarray, fm: dict):
    """照抄 a4_from_matrix 的构造（只读，不改）。"""
    sig = {}
    for cn, bone, joint, subs in pa.REGIONS:
        d = yi["detail"][cn]
        if bone in pa.LONG_BONES:
            loads = np.asarray(fm["loads_n"][joint], float)
            scale = loads / d["fe_solve_n"]
            sig[cn] = np.asarray(d["ends_stat_mpa"], float) * scale
        else:
            sig[cn] = np.asarray(fm["sigma_vm_mpa"]["median"][bone], float)
    X = np.column_stack([np.broadcast_to(sig[cn], (len(heights),)) for cn in REGION_NAMES])
    return X


def parallel_analysis_null(X: np.ndarray, n_iter: int = 3000, seed: int = 0) -> np.ndarray:
    """并行分析：同等形状随机数据的各主成分解释方差比 95 分位（数据驱动门槛）。"""
    rng = np.random.default_rng(seed)
    n, p = X.shape
    ratios = np.empty((n_iter, min(p, 3)))
    for i in range(n_iter):
        R = rng.standard_normal((n, p))
        Rc = R - R.mean(axis=0, keepdims=True)
        w = np.linalg.eigvalsh(np.cov(Rc, rowvar=False))[::-1]
        ratios[i] = (w[: ratios.shape[1]] / w.sum())
    return np.percentile(ratios, 95, axis=0)


def main() -> int:
    fm = json.loads(FM.read_text(encoding="utf-8"))
    yi_raw = json.loads(YI.read_text(encoding="utf-8"))
    yi = {"detail": yi_raw["detail_yi"]}   # 存盘键 = detail_yi
    heights = np.asarray(fm["heights_m"], float)

    X = sig_matrix(yi, heights, fm)
    Xc = X - X.mean(axis=0, keepdims=True)

    print(f"矩阵: {X.shape[0]} 高度 x {X.shape[1]} 部件    长骨集合 = {sorted(LONG)}")
    print()
    print("=== Q1 逐部件 σ(h) 变异性 ===")
    live = []
    for i, cn in enumerate(REGION_NAMES):
        col = Xc[:, i]
        rng_ = float(col.max() - col.min())
        sd = float(col.std())
        mean_abs = float(np.abs(X[:, i]).mean())
        is_live = rng_ > 0.0
        if is_live:
            live.append(cn)
        tag = "活" if is_live else "**死（常数）**"
        rel = (rng_ / mean_abs) if mean_abs > 0 else float("nan")
        print(f"  {cn:<6} σ 均值={mean_abs:>12.4g} MPa  σ 极差={rng_:>12.4g}  "
              f"相对极差={rel:>7.2%}  -> {tag}")
    print(f"\n  活变量 {len(live)}/{len(REGION_NAMES)}：{live}")
    print(f"  死变量 {len(REGION_NAMES)-len(live)}/{len(REGION_NAMES)}："
          f"{[c for c in REGION_NAMES if c not in live]}")

    # ---- 全部 9 变量（现行口径）----
    C = np.cov(Xc, rowvar=False)
    w = np.linalg.eigvalsh(C)[::-1]
    ratio_all = w / w.sum()
    null_all = parallel_analysis_null(X)

    print("\n=== Q2/Q3 现行口径（9 变量全用）===")
    print(f"  协方差秩（数值） = {int(np.sum(w > 1e-9 * w[0]))}")
    print(f"  PC1/PC2/PC3 解释比 = {ratio_all[0]:.5%} / {ratio_all[1]:.5%} / {ratio_all[2]:.5%}")
    print(f"  并行分析 null(95%) = {null_all[0]:.5%} / {null_all[1]:.5%} / {null_all[2]:.5%}")
    print(f"  PC2 {'>' if ratio_all[1] > null_all[1] else '<='} null2 ⇒ rank-2 判据 "
          f"{'PASS' if ratio_all[1] > null_all[1] else 'FAIL'}")

    # ---- 只在活变量上（口径诊断，不代表应改判据）----
    if live and len(live) < len(REGION_NAMES):
        Xi = Xc[:, [REGION_NAMES.index(c) for c in live]]
        Ci = np.cov(Xi, rowvar=False)
        wi = np.linalg.eigvalsh(Ci)[::-1]
        ratio_i = wi / wi.sum()
        rng_ = np.random.default_rng(0)
        n, p = Xi.shape
        it = 3000
        rr = np.empty((it, min(p - 1, 3)))
        for k in range(it):
            R = rng_.standard_normal((n, p))
            Rc = R - R.mean(axis=0, keepdims=True)
            wk = np.linalg.eigvalsh(np.cov(Rc, rowvar=False))[::-1]
            rr[k] = wk[: rr.shape[1]] / wk.sum()
        null_i = np.percentile(rr, 95, axis=0)
        print(f"\n=== 口径诊断：只用 {len(live)} 个活变量 ===")
        print(f"  PC1/PC2 解释比 = {ratio_i[0]:.5%} / {ratio_i[1]:.5%}")
        print(f"  并行分析 null(95%, {len(live)} 变量) = {null_i[0]:.5%} / {null_i[1]:.5%}")
        print(f"  PC2 {'>' if ratio_i[1] > null_i[1] else '<='} null2 ⇒ "
              f"{'PASS' if ratio_i[1] > null_i[1] else 'FAIL'}")

    # ---- 决定性实验：标准化（相关阵）vs 原始协方差 ----
    Z = Xc / Xc.std(axis=0, keepdims=True)
    R = np.corrcoef(Z, rowvar=False)
    rz = np.linalg.eigvalsh(R)[::-1]
    ratio_z = rz / rz.sum()
    null_z = parallel_analysis_null(Z)
    print("\n=== 决定性实验：口径对比（同一份数据，只换标准化）===")
    print(f"  原始协方差  : PC1={ratio_all[0]:.4%}  PC2={ratio_all[1]:.4%}")
    print(f"  标准化(相关阵): PC1={ratio_z[0]:.4%}  PC2={ratio_z[1]:.4%}   "
          f"（PC2 放大 {ratio_z[1]/ratio_all[1]:.1f}×）")
    print(f"  并行分析 null : 协方差口径 {null_all[1]:.4%}   标准化口径 {null_z[1]:.4%}")
    print(f"  PC2 vs null   : 协方差 {'PASS' if ratio_all[1] > null_all[1] else 'FAIL'}   "
          f"标准化 {'PASS' if ratio_z[1] > null_z[1] else 'FAIL'}")

    # 量级离散度（解释为什么）
    sds = Xc.std(axis=0)
    print(f"\n  各部件 σ 标准差跨度: {sds.min():.2f} … {sds.max():.2f} MPa "
          f"⇒ 方差比 {(sds.max()/sds.min())**2:.0f}×  (原始协方差由此被最大部件主导)")

    # ---- 能力探针（存在性测试，非调参）----
    # 问：σ(h) 后处理层**有没有能力**造出 PC2 > null？
    # 做法：把全部部件沿卡内最强软化端（E_t=0，完美平台）映射，扫统一屈服阈值 σ_y。
    # 若连「全平台」都做不到 ⇒ 架构上不可能 ⇒ A4 = 模型类别不足（非本构强度问题）。
    print("\n=== 能力探针：σ(h) 后处理层能否造出 rank-2？（E_t=0 完美平台，扫 σ_y）===")
    print(f"  {'σ_y(MPa)':>9} {'PC1':>10} {'PC2':>10} {'null2':>10}  判定")
    best = -1.0
    for sy in [10.0, 25.0, 50.0, 100.0, 150.0, 200.0, 300.0, 500.0, 1e9]:
        Y = np.where(X > sy, sy, X)          # E_t = 0 平台
        Yc = Y - Y.mean(axis=0, keepdims=True)
        live_j = Yc.std(axis=0) > 0
        Yl = Yc[:, live_j]
        if Yl.shape[1] < 2:
            continue
        Cl = np.cov(Yl, rowvar=False)
        wl = np.linalg.eigvalsh(Cl)[::-1]
        rl = wl / wl.sum()
        nll = parallel_analysis_null(Yl, n_iter=800)
        ok = rl[1] > nll[1]
        best = max(best, rl[1] / nll[1])
        tag = "PASS" if ok else "FAIL"
        note = " (= 无损伤基线)" if sy >= 1e9 else ""
        print(f"  {sy:>9.0f} {rl[0]:>9.4%} {rl[1]:>9.4%} {nll[1]:>9.4%}  {tag}{note}")
    print(f"\n  ⇒ 全平台扫描下 PC2/null2 最大值 = {best:.3f}  "
          f"({'存在可达阈值' if best > 1 else '**任何阈值都不可达** ⇒ 架构上不可能'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
