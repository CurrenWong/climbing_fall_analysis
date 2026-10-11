# G8 判定补录 —— 腓骨渐进 OR（1D 载荷分享路径 · 路线 A）

> 生成 2026-10-11 12:22；仓库 `D:\Project\climbing_fall_analysis`。耗时 2s。
> 预注册设计：`docs/methods/G8_腓骨OR_1D载荷分享设计.md`（含 §3 勘误：`σ_fib = risk × σc_row(160)`，非 ×70）。**纯增量**：不改任何既有文件、不重跑 FE。
> 口径纪律：OR 一律 **per-metre + 主设定分担比 + 混合口径**（腓骨 1D 行 + 其余 FE 行）。

---

## 0. 一句话结论

把腓骨一行换成 **1D 载荷分享名义应力**（`σ_fib(h)=risk(h)×160`，5% 分担）接入 G3 后：**腓骨 status=`ok`，OR=1.152（靶 1.682 ±30%）**。
⇒ **DoD-1 FAIL**：主设定腓骨未落靶（status=`ok`，OR=1.152）。
A3 = **2/3（PASS）**。

---

## 1. 预注册参数快照（回显）

| 项 | 值 |
|---|---|
| 主设定分担比 | 腓骨 **0.05** / 胫骨 0.95（`rows_with_share_lit05`）|
| 敏感性分担比 | 腓骨 0.10 / 胫骨 0.90（`rows_with_share_lit10`）|
| n_subj / seeds | 400 × [20261003, 20261004, 20261005]（官方 20261003）|
| 主设定 COV | COV_B=0.15, COV_E=0.10, COV_A=0.15（G3 冻结值）|
| 腓骨 σc（弱区，G3） | 70 MPa |
| 1D 行 σc_row | 160 MPa（骨干）|
| 高度栅格 | 1…50 m（1D 行与 G3 同名：True）|
| [Y25] 腓骨靶 | 1.682，±30% ⇒ [1.177, 2.187] |

---

## 2. 腓骨接入后的结果（主设定 vs 敏感性）

| 分担比 | 腓骨 status | 腓骨 OR | Wald 95% CI | n_fracture | util@5m | util_max | 确定性首骨折 h |
|---|---|---:|---|---:|---:|---:|---:|
| 0.05 | ok | 1.152 | 1.136–1.167 | 433 | 0.177 | 0.791 | >50 |
| 0.10 | ok | 1.147 | 1.143–1.151 | 8142 | 0.354 | 1.583 | 24 |

> `util = σ_fib(h)/σc_weak(70)`；`util_max<1` ⇒ 确定性不骨折，概率层的混合来自被试变异。

**主设定 A3 全表（per-metre + 主设定 COV）：**

| 部位 | 状态 | OR | [Y25] 靶 | ±30% 内 |
|---|---|---:|---:|:--:|
| 腓骨 | ok | 1.152 | 1.682 | ❌ |
| 骨盆 | ok | 1.147 | 1.236 | ✅ |
| 颅骨 | ok | 1.417 | 1.576 | ✅ |

**主设定 A3 = 2/3 ⇒ PASS**

---

## 3. 腓骨概率矩阵 p(h)（主设定，节选）

| 部位 | h=1 | h=5 | h=10 | h=20 | h=30 | h=40 | h=50 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 腓骨 | 0.000 | 0.000 | 0.000 | 0.000 | 0.007 | 0.045 | 0.102 |

> p(h) 单调非降；主设定腓骨不再是 `all_one`（全 1）而是 `0→…→p<1` 的渐进曲线。

---

## 4. 其余 8 部位不变性（DoD-2）

对照 `g3_or_estimates.json` 主设定 per_site：**全部逐值相同 = True**（max|ΔOR| = 0.00e+00）。

| 部位 | 与 G3 相同 |
|---|---|
| 足部 | ✅ |
| 胫骨 | ✅ |
| 股骨 | ✅ |
| 骨盆 | ✅ |
| 腰椎 | ✅ |
| 胸椎 | ✅ |
| 颈椎 | ✅ |
| 颅骨 | ✅ |

> 原理：只替换腓骨行 ⇒ 其余行的 σ_nom 未变，被试因子（逐被试 `default_rng([seed,s])`）与高度无关、完全一致 ⇒ 其余行 y 不变 ⇒ OR 逐值相同。

---

## 5. DoD 逐条

- **DoD-1 腓骨渐进 + OR∈[1.18,2.19]** — **FAIL**（主设定 status=`ok`，OR=1.152）。
- **DoD-2 其余部位不变** — **PASS**。
- **DoD-3 A3 ≥2/3** — **PASS**（2/3）。
- **DoD-4 极限退化** — **PASS**（COV=0.001 且质量不入载荷，多数判定复现含 1D 腓骨行的确定性基准；max|Δ|=3.58e-01）。
- **DoD-5 零回归** — 纯增量（见 §7）；回归基线 ≥687 passed / 8 skipped；新增 `tests/test_g8_fibula_or.py`（见 §9）。
- **DoD-6 口径纪律** — OR 均标 per-metre + 主设定分担比 + 混合口径；未出现禁用比值。

---

## 6. 诚实结论与风险

1. **混合口径**：腓骨行 = 1D `F/A` 均匀轴向名义应力；其余 8 行 = FE 几何解 —— 两者**不同源**，本件为**显式建模选择**。
2. **数据源勘误**：初稿误用 `×70`；本件按 `×σc_row(160)`（真实名义应力）实现，并已回写设计稿 §3。此更正不改变任何预注册决策。
3. **OR 是否落靶**：未落靶（如实报，未回解）。
4. **敏感性**：分担比 0.10 的结果一并报告（未挑）。

---

## 7. 纯增量核对

G8 新增（应只见以下条目，无既有文件被改）：

- `scripts/opensim_fe/g8_fibula_or_share.py`
- `results/opensim_fe/g8_fibula_or.json`
- `results/opensim_fe/G8_FIBULA_OR_REPORT.md`
- `tests/test_g8_fibula_or.py`
- `docs/methods/G8_腓骨OR_1D载荷分享设计.md`

`git status --porcelain` 当前条目数 52（既有未提交改动 + 本件新增；详细清单见仓库）。

---

## 8. 复现命令与产物

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\g8_fibula_or_share.py
# 只读复用：
#   results/opensim_fe/risk_1d_loadshare.json        (腓骨 1D 行)
#   results/opensim_fe/fracture_matrix_bc2.json       (其余 8 部位 median)
#   scripts/opensim_fe/g3_subject_variability.py      (被试变异机制, 零改动)
# 产物（纯增量）：
#   results/opensim_fe/g8_fibula_or.json
#   results/opensim_fe/G8_FIBULA_OR_REPORT.md
```

**弱区 σc（route-乙，G3 同款）：**

| 部位 | 足部 | 胫骨 | 腓骨 | 股骨 | 骨盆 | 腰椎 | 胸椎 | 颈椎 | 颅骨 |
|---|---|---|---|---|---|---|---|---|---|
| σc (MPa) | 150 | 70 | 70 | 80 | 180 | 150 | 150 | 150 | 160 |

---

## 9. 实测记录

- G8 脚本：`scripts/opensim_fe/g8_fibula_or_share.py`（官方 seed 20261003，n_subj=400，2 档分担比 × 3 seeds）。
- 测试：`tests/test_g8_fibula_or.py`；全测试套件基线 687 passed / 8 skipped（G8 为纯增量，实测见仓库运行日志）。

