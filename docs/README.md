# docs/

抱石落地骨与韧带损伤建模项目的设计文档、阶段报告、论文稿与文献笔记。

> **本页自 2026-10-10 起是「权威索引」** —— 解决"文档太多、不知道信哪份"。
> **规则：一个主题只有一份权威。** 其余是过程件或归档件，**不要从里面取结论**。

---

## 🎯 权威索引（一个主题 → 只看这一份）

| 我想知道… | 只看这一份 |
|---|---|
| **现在整体到哪了 / 下一步做什么** | [`reports/总体进度.md`](reports/总体进度.md) —— 一页总览，**§四 = 下一步**，§五 = V/A 靶子总表 |
| **目标是什么 / 做到什么算够 / 推断怎么合法** | [`methods/目标重构_刚性验证与垫上推断.md`](methods/目标重构_刚性验证与垫上推断.md) —— **目标层权威件**（两阶段：刚性验证 → 垫上推断） |
| **论文写什么 / 数字从哪来** | [`paper/论文v2_数学建模版.md`](paper/论文v2_数学建模版.md) ＋ [`paper/论文v2_出处索引.md`](paper/论文v2_出处索引.md)（逐条溯源） |
| **当前结论与可信边界** | [`reports/最终结论_2026-10-10.md`](reports/最终结论_2026-10-10.md)（重述版；旧件 `最终结论_2026-10-07.md` 保留为历史） |
| **线 1 路线 P0–P5** | [`methods/方案更新_v2.md`](methods/方案更新_v2.md) |
| **结构件 G1 多段链（= P3 二切）** | [`methods/G1_多段链设计.md`](methods/G1_多段链设计.md)（预注册，已冻结）＋ 判定 [`../results/G1_CHAIN2D_REPORT.md`](../results/G1_CHAIN2D_REPORT.md)、[`../results/G1_WAVE_ROD_REPORT.md`](../results/G1_WAVE_ROD_REPORT.md) |
| **结构件 G3 被试级变异（统计层）** | [`methods/G3_被试级变异_统计层设计.md`](methods/G3_被试级变异_统计层设计.md)（预注册，已执行）＋ 判定 [`../results/opensim_fe/S4_A3_G3_REPORT.md`](../results/opensim_fe/S4_A3_G3_REPORT.md) |
| **结构件 G4 损伤演化** | [`methods/G4_损伤演化设计.md`](methods/G4_损伤演化设计.md)（预注册，已冻结）＋ 判定 [`../results/opensim_fe/G4_DAMAGE_REPORT.md`](../results/opensim_fe/G4_DAMAGE_REPORT.md)（最小实现已跑 10-11，A4 如实 FAIL） |
| **结构件 G4 载荷链顺应性杠杆（第二杠杆）** | 判定 [`../results/opensim_fe/G4_LOADCHAIN_REPORT.md`](../results/opensim_fe/G4_LOADCHAIN_REPORT.md)（已跑 10-11；唯一冻结预设 D6 重跑顺应载荷链；A4 rank-2/LP 如实 FAIL；**A2 股骨获阶跃、胸椎 all_zero**；§7 触发未成立） |
| **结构件 G4 · 全量 FE 损伤材料（② FE 级深水区）** | 计划 [`methods/G4_FE全量_工作计划.md`](methods/G4_FE全量_工作计划.md)（**v2.5**）—— 把 σ(h) 后处理损伤升级为 **FEBio 原生 `elastic damage` + 逐高度 FE 重解**；目标 A4 rank-2（统计 1/4→2/4）；**腓骨 OR 补不上（结构上限）** |
| **G4 ② 五步实证链（P0→P1 试点，均本机实测）** | P0 [`G4_FE_SPIKE_REPORT.md`](../results/opensim_fe/G4_FE_SPIKE_REPORT.md)（材料 ✅，1.0 s/解）· P1-a [`G4_FE_STAB_REPORT.md`](../results/opensim_fe/G4_FE_STAB_REPORT.md)（stock 杠杆 ❌，差 23.2×）· 可视化 [`g4_fe_viz/README.md`](../results/opensim_fe/g4_fe_viz/README.md)（损伤**萌生于人工加载带**、失效=屈曲折叠）· erosion [`G4_FE_EROSION_REPORT.md`](../results/opensim_fe/G4_FE_EROSION_REPORT.md)（❌ 且更早发散）· 改 BC [`G4_FE_P1B_REPORT.md`](../results/opensim_fe/G4_FE_P1B_REPORT.md)（✅ 峰值 1066→147 MPa，7.2×）· P1 试点 [`G4_FE_PILOT_REPORT.md`](../results/opensim_fe/G4_FE_PILOT_REPORT.md)（σ(h) 非线性 ✅／逐高度重解 143/150 发散 ❌）。**总判定：材料可用，瓶颈 = FE 载荷 BC；须完整 P1-b（真接触/压板）方可交付** |
| **结构件 G8 腓骨 OR（1D 载荷分享，路线 A）** | [`methods/G8_腓骨OR_1D载荷分享设计.md`](methods/G8_腓骨OR_1D载荷分享设计.md)（预注册）＋ 判定 [`../results/opensim_fe/G8_FIBULA_OR_REPORT.md`](../results/opensim_fe/G8_FIBULA_OR_REPORT.md)（已跑 10-11；腓骨变渐进 status=ok 但 OR=1.152 差 2% ⇒ DoD-1 如实 FAIL） |
| **结构件 G9 腓骨 FE 侧向支撑（路线 B）** | [`methods/G9_腓骨FE侧向支撑设计.md`](methods/G9_腓骨FE侧向支撑设计.md)（预注册）＋ 判定 [`../results/opensim_fe/G9_FIBULA_FE_SUPPORT_REPORT.md`](../results/opensim_fe/G9_FIBULA_FE_SUPPORT_REPORT.md)（已跑 10-11；单带 FE 分叉非唯一 + 三段带稳定档 5%=all_zero/10% OR=1.130 ⇒ DoD FAIL） |
| **腓骨 all_one 根因 + FE 审计** | [`../results/opensim_fe/FIBULA_IDENTIFIABILITY_REPORT.md`](../results/opensim_fe/FIBULA_IDENTIFIABILITY_REPORT.md)（FE 量级：无胫-腓载荷分配）＋ [`../results/opensim_fe/FIBULA_FE_AUDIT_REPORT.md`](../results/opensim_fe/FIBULA_FE_AUDIT_REPORT.md)（10-11：分配必要但不充分，孤立腓骨 FE 屈曲） |
| **[Y25] OR 口径重审（我们 vs 论文是不是同一个量）** | [`reports/重审_Y25_OR口径.md`](reports/重审_Y25_OR口径.md) —— 我三部位 OR 均在 [Y25] 各自 95% CI 内、log-OR 差异均不显著；A3 的 ±30% 带**比 [Y25] 自身精度更严**（腓骨 CI 下界 1.106<1.177）⇒ 腓骨"未落靶"含判据口径效应；**OR 次序腓/颅互换为真实差异** |
| **OR（优势比）是什么 / 怎么读 / 怎么用** | [`notes/OR_优势比_说明.md`](notes/OR_优势比_说明.md) —— 概率→几率→OR、读法、用途、与 RR 的区别、本项目"变异口径"坑（写论文/答辩可直接引） |
| **阶段二推断 · P-1/P-2/P-3 输出清单** | [`reports/阶段二推断_输出清单.md`](reports/阶段二推断_输出清单.md) —— 三条垫侧推断的三件套（点估计 + 区间 + 工况标签【垫】）+ 可信度链 + 后验检查单；含 P-1 能量闭环 BREAKS 归属（= **口径差**非线1 误差）与 "接触窗 34 ms = **力峰时刻 tF**" 的口径澄清 |
| **P-4 踝骨折率 · 偏差入账** | [`reports/P4_踝骨折率_偏差入账.md`](reports/P4_踝骨折率_偏差入账.md) —— 踝骨折率**上界 22.8%**（Wilson [20.31,25.50]%）+ 系统偏差 **2.6–3.1×**（三旋钮全不敏感，结构性）+ 引用规约（**禁校正点估计**）；闭合阶段二 DoD 第 2 条 |
| **主线 P2 踝（解析模型）结论** | [`reports/P2_踝结论_2026-10-07.md`](reports/P2_踝结论_2026-10-07.md) |
| **工况错配到底占多少（P5 前置侦察）** | [`notes/P5_CONDITION_ALIGN.md`](notes/P5_CONDITION_ALIGN.md) |
| **踝 FE 子线怎么用 / 到哪了** | [`notes/ANKLE_FIX_TODO.md`](notes/ANKLE_FIX_TODO.md) **顶部「🧭 接手入口」** |
| **踝 FE 接触根因 + 罚标定** | [`notes/CASEB_CONTACT_ROOTCAUSE.md`](notes/CASEB_CONTACT_ROOTCAUSE.md) §1–§14 |
| **踝子线的定位（供不供论文数）** | [`reports/踝子线_scope决策_2026-10-07.md`](reports/踝子线_scope决策_2026-10-07.md)（§八 = 在论文里的角色） |
| **线 2 OpenSim → FE 怎么做** | [`reports/OpenSim_FE交接.md`](reports/OpenSim_FE交接.md) |
| **S5 唯一可引用数字** | [`reports/S5_QUOTABLE.md`](reports/S5_QUOTABLE.md)（正文结论必须引它） |
| **S5 脚-垫接触契约 / 门禁** | [`s5-contact/S5_contact_contract.md`](s5-contact/S5_contact_contract.md) |
| **踝 FE 运行入口（命令）** | [`methods/踝FE可用模型使用说明.md`](methods/踝FE可用模型使用说明.md) |

---

## 🗂 目录

| 目录 | 作用 | 权威件 |
| --- | --- | --- |
| `paper/` | **论文成品**（v2 稿 ＋ 出处索引 ＋ 评审 ＋ 汇总报告） | `论文v2_*.md` |
| `reports/` | **阶段结论与进度**（这是取结论的地方） | `总体进度.md` |
| `methods/` | 仿真方案、方法学说明、建模交接 | `方案更新_v2.md` |
| `notes/` | 阶段性技术记录（根因、研究、结论） | `CASEB_CONTACT_ROOTCAUSE.md`、`P5_CONDITION_ALIGN.md` |
| `s5-contact/` | S5 脚-垫接触收敛系列 17 篇，**成套阅读** | `S5_contact_contract.md` |
| `archive/` | **已归档（只读历史，不再是权威）** | 见 [`archive/README.md`](archive/README.md) |
| （根目录） | `README.md`、`rigid_vtarget_2026-10-08.md` | — |

---

## ⚠️ 三层文档，别混（这是"搞不清楚"的主因）

| 层 | 特征 | 数量 | 能不能取结论 |
| --- | --- | --- | --- |
| **① 持久权威** | `docs/**` 且**不以 `_TODO.md` 结尾**，进了 git | **120** | ✅ **唯一可以取结论的层** |
| **② 会话草稿** | **`*_TODO.md`**（根目录 **13** ＋ `docs/notes/` **23**） | **36** | ❌ **被 `.gitignore` 排除，只在本机**，换机器/清库即消失 |
| **③ 归档** | `docs/archive/`（5 份文档 ＋ 1 份 README） | 5 | ❌ 只读历史 |

> 口径：项目自有 markdown **156** 份（不含 `.venv/`、根 `paper/` 第三方文献、`model/`）。
> 归档/待判/绝不动的处置判断 → [`reports/文档处置说明.md`](reports/文档处置说明.md)。

⇒ **看到 `*_TODO.md` 就知道那是"当时的待办/接手笔记"，不是结论。**
⇒ 真正的结论只在 ① 里；② 过期后应把结论**搬进 ①** 再归档或删除。

> 另：项目文档普遍用**短路径**引用（写 `docs/X.md`、实际在 `docs/reports/X.md`），
> 这是全库惯例，**不是死链**。

---

## 📖 推荐阅读顺序

1. [`reports/总体进度.md`](reports/总体进度.md) —— 先看一页总览与 §四 下一步
2. [`reports/最终结论_2026-10-07.md`](reports/最终结论_2026-10-07.md) —— 结论与可信边界
3. [`paper/论文v2_数学建模版.md`](paper/论文v2_数学建模版.md) —— 论文本体（数字溯源走 `论文v2_出处索引.md`）
4. [`methods/方案更新_v2.md`](methods/方案更新_v2.md) —— 路线与 P0–P5 定义
5. 需要细节再进对应子目录

## 📦 已归档

[`archive/`](archive/) 存放已被取代的稿件（当前：**论文初稿** → 由 v2 取代）。
**归档 ≠ 删除**，内容完整可查；⚠️ 但**章节号不可与 v2 互换**，详见 [`archive/README.md`](archive/README.md)。

## 第三方模型与文献

`model/`（AM50、THUMS、足踝 FE 模型等）与 `paper/` 下的 `MinerU_markdown_*.md`
（吴恺论文、Anderson 双论文、FEBio 手册等第三方资料）**不随本仓库分发**，因各自许可约束。
复现所需的获取方式、版本与放置约定见仓库根目录 `README.md`。
