---
title: G8 · 腓骨渐进 OR —— 1D 载荷分享路径接入（路线 A）
date: 2026-10-11
type: 设计（待审批）
status: 预注册（pre-registration）· 待批准即锁
tags:
  - G8
  - 腓骨
  - 渐进 OR
  - 载荷分享
  - A3
  - 预注册
---

# G8 · 腓骨渐进 OR —— 1D 载荷分享路径接入（路线 A · 预注册）

> **一句话**：把**腓骨**的 σ(h) 从**既有 1D 胫-腓载荷分享路径**（已是生理量级、不饱和）接入 G3 概率机制，使腓骨从 `all_one` 变为**渐进**、Logistic OR 可定义，并对标 [Y25] 靶 **1.682**（per-metre，±30%）。
> **性质 = 预注册**：数据源、主设定分担比、判定规则在跑数**前**锁定；跑后**禁止**回解分担比/σc 凑 OR（红线：不调参凑靶）。
> **上游件**：`results/opensim_fe/FIBULA_IDENTIFIABILITY_REPORT.md`、`results/opensim_fe/FIBULA_FE_AUDIT_REPORT.md`（根因）、`docs/methods/G3_被试级变异_统计层设计.md`（G3 机制）、`docs/methods/目标重构_刚性验证与垫上推断.md` §三（G3 行）。
> **路线定位**：三路线（A/B/C）中的 **A（便宜·零 FE）**；B（修腓骨 FE）/C（改 σc 口径）为备选，见 §8。

---

## 0. 目标与范围

| | 内容 |
|---|---|
| **直接解锁** | **A3·腓骨**：腓骨 Logistic OR 从 `all_one`（无定义）变为 **渐进 + OR ∈ [1.18, 2.19]**（= [Y25] 1.682 ±30%） |
| **副产品** | 完整三渐进部位 OR（腓骨 / 骨盆 1.147 / 颅骨 1.417）全落靶 ⇒ A3 有望 **3/3** |
| **方法** | 复用 G3 被试级变异概率机制（`logistic_fit` 原封）；**只换腓骨一行的 σ(h) 来源** |
| **零 FE** | 腓骨用 1D 路径；其余 8 部位用**已冻结**的 FE σ(h)；**无任何 FE 重跑** |
| **不解锁** | 胫骨两端"概率层形状"（见 §8 路线 C）；A4；V8b；绝对骨折高度 |

---

## 1. 背景与缺口

- **[Y25]**：腓骨 ∈ **渐进部位**（`paper_alignment.PAPER_GRADUAL = {腓骨, 骨盆, 颅骨}`），per-metre OR = **1.682**（`A3_TARGETS`），验收 ±30% ⇒ **[1.18, 2.19]**。
- **现状**（`S4_A3_G3_REPORT.md`）：G3 主设定行腓骨 = **`all_one`**（OR 无定义）⇒ A3 = **2/3**（骨盆 1.147 ✅、颅骨 1.417 ✅、腓骨 ❌）。
- **根因**（`FIBULA_FE_AUDIT_REPORT.md`）：FE 把**整条 `ankle_r` 反力**同时喂给胫骨与腓骨（**未做胫-腓载荷分配**），叠加**孤立腓骨 FE 在 ~550–700 N 屈曲**（骨间膜支撑被忽略）⇒ 腓骨 σ_nom/σc ≈ **6–44×** ⇒ 确定性列全 1 ⇒ 被试变异（各因子 ≈1±0.15）无法 lift 成概率。

---

## 2. 关键事实（使路线 A 便宜且合理）

| 事实 | 出处 | 用法 |
|---|---|---|
| **1D 胫-腓载荷分享路径已存在**：`build_rows` 对每骨算 `risk(h) = mult·loads[joint](h) / (A_section·σc)`，`mult` = 分担比；含 `share_lit05={胫0.95,腓0.05}` / `share_lit10={胫0.90,腓0.10}` | `scripts/opensim_fe/risk_1d_loadshare.py`（:195–242、:304–306） | **腓骨 σ(h) 的全部来源** |
| **既有 1D 锚**：`fibula_ends.load_fraction = 0.05`（单侧腓骨 ~5%） | `src/climbing/bone.py::fibula_ends` | 主设定分担比（**非新值**） |
| **1D 腓骨不饱和/不屈曲**：`risk@5m` = **0.077（5%）/ 0.155（10%）**（`risk=F/(A·σc_row)`，`σc_row=160`）；全高度 `risk<1` ⇒ **确定性不骨折**，概率层的 0/1 混合**由被试变异产生**（非确定性） | `RISK_1D_LOADSHARE_REPORT.md` / `fibula_identifiability_analysis` §6 | 保证 OR 可定义 |
| **冻结的腓骨行输入**：`rows_with_share_lit05` / `rows_with_share_lit10`（含 `risk_by_height` 50 档 + **`sigma_c_mpa=160`**〔骨干强度，**非** 70〕） | `results/opensim_fe/risk_1d_loadshare.json` | 只读抽取 |
| **G3 机制本身有效**：骨盆/颅骨 OR 落靶，极限退化/单调/种子稳定全 PASS | `S4_A3_G3_REPORT.md` §6 | 只换腓骨输入，不动机制 |
| **FE σ(h) 其余 8 行**：`sigma_vm_mpa["median"][bone]`（只读冻结） | `results/opensim_fe/fracture_matrix_bc2.json` | 混合矩阵的非腓骨部分 |

---

## 3. 数据来源与口径（⭐ 冻结）

| 项 | 取值 | 性质 |
|---|---|---|
| **腓骨 σ_fib(h)** | `σ_fib(h) = risk_by_height(h) · σc_row`，**`σc_row=160 MPa`**（骨干强度）= **真实名义轴向应力 `F/A`**。⚠️ **易错点**：行 `σc` 是 **160**（骨干），**不是**弱区 70；若误用 `×70` 则 `util=risk<1` ⇒ 全高度恒不骨折（all_zero），**已勘误**。取自 `risk_1d_loadshare.json::rows_with_share_lit05`（腓骨行） | 1D，`F/A` 均匀轴向 |
| **主设定分担比** | **腓骨 0.05 / 胫骨 0.95**（`share_lit05`） | = `bone.py::fibula_ends` 既有锚，**非调参** |
| **敏感性分担比** | 腓骨 0.10 / 胫骨 0.90（`share_lit10`） | 全量报告，不挑 |
| **其余 8 部位 σ(h)** | FE `fracture_matrix_bc2.json::sigma_vm_mpa["median"][bone]` | 冻结，只读 |
| **σc** | route-乙 弱区口径（足150/胫70/腓70/股80/盆180/腰胸颈150/颅160），**不动** | 与 G3 一致 |
| **高度栅格** | 1…50 m（与 FE 矩阵、G3 同栅格；接入前**校验**同名高度） | 一致性检查 |

> **勘误（2026-10-11，实现前核对数据源发现）**：初稿 §2/§3 写 `σ_fib = risk_by_height × 70 MPa`，**错把 1D 行的归一化强度（`σc_row=160`，骨干）当成了弱区 70**。正确：`σ_fib(h) = risk_by_height(h) × 160`（= 真实名义应力 F/A），再与 G3 的弱区阈值 `σc=70` 相比 ⇒ `util = risk × 160/70`。此更正**只修数据源换算**，不改任何预注册决策（F1–F5：来源=1D 路径、主设定 0.05、敏感性 0.10、判定规则、混合口径均不变）。核对事实：`rows_with_share_lit05` 腓骨行 `sigma_c_mpa=160.0`、`risk@5m=0.0775`、`risk(50m)=0.346`。

---

## 4. 接入方法（零 FE）

1. **构造混合 σ_nom 矩阵**：9×50，其中**腓骨行 ← §3 的 1D σ_fib(h)**，其余 8 行 ← FE（§3）；
2. **跑 G3 被试变异**：COV_B=0.15 / COV_E=0.10 / COV_A=0.15；n_subj=400 × seeds {20261003, 20261004, 20261005}；`MASS_MODEL` 原样；`logistic_fit` **原封复用**；
3. **算腓骨 OR / status**（per-metre，Wald CI）；
4. **判定**：`s4_alignment_bc2.a3_verdict` **原封**（腓骨靶 1.682，±30%，A3 ≥2/3）。

---

## 5. 预注册决策点（已建议冻结）

| # | 决策 | 冻结值 | 依据 |
|---|---|---|---|
| **F1** | 腓骨 σ(h) 来源 | **1D 载荷分享路径**（`rows_with_share_lit05`） | 不饱和/不屈曲；零 FE |
| **F2** | 主设定分担比 | **腓 0.05 / 胫 0.95** | `bone.py::fibula_ends.load_fraction` 既有锚（非调参） |
| **F3** | 口径 | **混合**：腓骨 1D 行 + 其余 FE 行；本件**明示该混用** | 承认口径差，如实标 |
| **F4** | 判定规则 | `logistic_fit` / `a3_verdict` **零改动** | 冻结 |
| **F5** | 敏感性 | 分担比 **0.10** 一档（grid 全量报告，**不挑**） | 防挑最好看 |

> **冻结即锁**：跑后**禁止**回解分担比 / σc 去使 OR 落进 ±30%。主设定 FAIL 就是 FAIL。

---

## 6. 验收判据（DoD）

| # | 判据 | 线 |
|---|---|---|
| 1 | **腓骨** | `status = ok`（渐进）**且 OR ∈ [1.18, 2.19]**（主设定 0.05） |
| 2 | **其余部位不变** | 骨盆 1.147、颅骨 1.417、及其余 6 部位 OR/status 与 G3 主设定**逐值相同** |
| 3 | **A3** | ≥2/3（不降；若腓骨入靶则 3/3） |
| 4 | **极限退化** | COV→0 时混合矩阵收敛回其确定性基准（含 1D 腓骨行） |
| 5 | **零回归** | 既有产物逐字节不动；全测试套件 **≥687 passed / 8 skipped**；**纯增量**（新脚本+新结果+新测试） |
| 6 | **口径纪律** | OR 标 "per-metre + 主设定分担比 + 混合口径"；不出现 σ/σ_law 比值类禁用数字 |

---

## 7. 诚实风险（必须记录）

1. **混合口径**：腓骨行是 1D `F/A` 均匀轴向，其余行是 FE 几何解 —— 两者 σ 量级**可比性须声明**；本件把"腓骨用 1D、其余用 FE"作为**显式建模选择**登记，不得声称两行同源。
2. **OR 可能仍出 ±30%**：1D 腓骨 OR 未必等于 1.682 ⇒ **如实 FAIL**，不回解。
3. **栅格一致性**：需核对 1D `heights` 与 FE 矩阵 `heights_m` 同名；不一致则按 FE 栅格插值（登记为口径）。
4. **1D 载荷链与 FE 载荷链同源**（均 OpenSim `ankle_r`），但 1D 不含 FE 的降载折算 —— 登记。

---

## 8. 范围外 / 备选

- **路线 B（修腓骨 FE）**：载荷分配（`fracture_matrix.py::load_fraction`，已实现）+ 骨间膜侧向支撑抑屈曲 + 重解 —— FE 级，另立项。
- **路线 C（胫骨两端概率层形状）**：命中 A2 的"阶跃"已成立（deterministic `complete_separation`）；若要把概率层 sigmoid 收成真阶跃，需动 σc 口径 —— 单独预注册，本件**不做**。
- 修 A2 其余 miss（**股骨 all_zero / 胸椎 all_one**）—— 与本件无关，另议。

---

## 9. 依赖与复现

- **复用（只读）**：`results/opensim_fe/risk_1d_loadshare.json`、`results/opensim_fe/fracture_matrix_bc2.json`、`scripts/opensim_fe/risk_1d_loadshare.py`、`scripts/opensim_fe/g3_subject_variability.py`、`paper_alignment.py::{logistic_fit, REGIONS}`、`s4_alignment_bc2.py::{A3_TARGETS, a3_verdict}`、`scenarios.py::MASS_MODEL`。
- **命令**：`$env:PYTHONPATH="src"; & .venv\Scripts\python.exe scripts\opensim_fe\g8_fibula_or_share.py`（RMS 纪律：必用 `.venv`）。
- **依赖**：`numpy`、`scipy`。**零 FE 重跑**。

---

## 10. 交付物

| 件 | 内容 |
|---|---|
| `scripts/opensim_fe/g8_fibula_or_share.py` | 腓骨 1D σ(h) 抽取 + 混合矩阵 + G3 复算（增量，不改既有文件） |
| `results/opensim_fe/g8_fibula_or.json` | 腓骨 OR/CI/status × {主设定 0.05, 敏感性 0.10} × 3 seeds + 混合矩阵 meta（预注册参数快照） |
| `results/opensim_fe/G8_FIBULA_OR_REPORT.md` | 判定补录（含 DoD 1–6 逐条 + 诚实风险） |
| `tests/test_g8_fibula_or.py` | 混合矩阵构造 + 非腓骨行不变 + 极限退化 + `logistic_fit` 接口锁 |
| 文档同步 | `最终结论_2026-10-10.md` A3 行、`目标重构` §三 G3/G8 行、`总体进度.md` A3 行、`docs/README.md` 索引按结果更新 |

---

*本件为 G8 预注册设计稿（路线 A）。批准即冻结 §3/§5；实现须严格照此执行，跑后如实判定。*
