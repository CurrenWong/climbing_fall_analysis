# 论文图 1–9 图注（FIGURE CAPTIONS）

> 生成脚本：`scripts/opensim_fe/make_paper_figures.py`（仅新增，只读复用既有 JSON 产物；不调用 FEBio、不重跑 OpenSim、不改任何既有文件）。
> 输出目录：`results/opensim_fe/figures/`；分辨率 **300 dpi**（>=200 dpi 要求已满足）。
> 用途：供 `docs/论文初稿_抱石落地骨与韧带损伤建模.md` 图表清单引用。
> **字体路由**：本机检测到 CJK 字体 **Microsoft YaHei**（备选 Microsoft YaHei UI / SimHei / SimSun / Noto Sans CJK SC），故图中轴标签/图例/注记均以**中文渲染**；若在无 CJK 字体的环境重跑，脚本会自动**回退英文轴标签**（判定逻辑见脚本 `_pick_cjk_font()` / `FONT_ROUTE`），中文图注仍保留于本文件。

**标签约定**：
- `measured` = OpenSim/FE/THUMS 实测或从实测读出的量；
- `modeled` = 由模型公式/几何代理导出的量（绝对值不可直接引用，只看相对与是否越阈）；
- `assumed` = 无实测、取文献或量级估计的参数。

---

## 图 2  `fig2_alignment.png` — 纯轴向+弱子区域 1D 首骨折高度与 [Y25] 对齐

**图注（中文）**：纯轴向 + [Y25] 弱子区域 1D 口径下 9 个部位的首骨折高度。红条 `腓骨 = 7 m` 落在 [Y25] 腓骨两端 7–9 m 区间（红色阴影带）内，是本路线**唯一**对齐论文的单点；`足部/胫骨 = 8 m`、`股骨 = 19 m`、`胸椎 = 47 m`；`腰椎/颈椎/骨盆/颅骨` 在 1–50 m 内**从不越阈**（灰斜纹条，标 `>50`，为外推值）。柱色：腓骨红、足蓝、胫绿、股灰橙；其余灰。**结论**：腓骨首骨折高度对齐论文区间，但整体排序与论文“颅骨>骨盆渐进”不一致。

- **数据源**：`results/opensim_fe/axial_subregion.json`（`route_a_axial_ends`；`paper_facts.fibula_ends_fracture_m=[7,9]`）。
- **标签**：
  - 首骨折高度 — **modeled**（纯轴向名义压缩 σ(h) 判据，1D）；
  - σ(h) — derived（`risk_1d_loadshare.json::risk_by_height × σ_c`，载荷链 measured）；
  - σ_c（腓骨两端 70 / 股骨颈 80 / 整域 150–180 MPa）— **assumed**（`src/climbing/bone.py::MATERIAL_STRENGTH_MPA`，[Y25] Table 1 压缩材料强度）；
  - `>50` — **modeled 外推**（1–50 m 内未越阈）。
- **对齐检验**：腓骨 7 m ∈ [7,9] m ✓。（对应 `AXIAL_SUBREGION_REPORT.md` §1.2 / `PAPER_ALIGNMENT_REPORT.md`。）

---

## 图 3  `fig3_single_foot.png` — 单脚落地对承力侧关节反力的放大

**图注（中文）**：(a) 单脚（全部冲量经一足，`split=(0,1)`）与对称（`(0.5,0.5)`）在 5 m 处的承力侧关节反力纵向分量 `F_y`（N）对比；(b) 载荷因子（单脚/对称）。下肢承力链 `subtalar_r ×2.014`、`ankle_r ×2.015`、`hip_r ×2.448`（≈ 翻倍，非严格 2：刚性锁死腿 + 偏心加载驱动骨盆转动）；中轴 `lumbar ×0.959`（≈ 不变）。**结论**：单脚使承力侧局部关节反力约翻倍，中轴近乎不变。

- **数据源**：`results/opensim_fe/nonvertical_s2.json`（`joint_load_factors_5m`）。
- **标签**：
  - 关节反力 sym/single @5 m — **measured**（OpenSim 逐高度载荷链）；
  - 载荷因子 — **derived**（measured 反力之比）；
  - 注：单脚**首骨折高度投影**为 **modeled**（载荷因子沿高度线性插值；见 `NONVERTICAL_S2_REPORT.md` §5），本图未使用投影值，仅用 5 m 实测反力。
- **对应报告**：`NONVERTICAL_S2_REPORT.md` §3 / `nonvertical_s2.json`。

---

## 图 4  `fig4_supination.png` — 踝旋后使腓骨远端在真实高度内越阈

**图注（中文）**：刚性地面下，叠加旋后角 β 后腓骨远端侧向 risk 随 β 的变化。(a) 侧向 risk（log 轴）：β=0 时 risk=0（逐位复现纯轴向），β≥15° 时在全部测点（2/3/3.5/4.5 m）**越阈**（红色虚线 `risk=1`）；模型临界角 β_c ≈ 0.82°（**上界式，不可直接引用**）。(b) 冠状面横向力 `F_lat = F_y·sinβ`（N）随 β。**结论**：纯轴向下腓骨在 ≤4.5 m 安全（首骨折 ≈7 m），一旦叠加旋后即在同一真实高度内越阈 —— 这正是纯轴向模型漏掉的 #1 临床机制（[Heck2024] §5.2.3）。

- **数据源**：`results/opensim_fe/nonvertical_s5.json`（`scenario_hard`，12 格 = 4 高度 × 3 角）。
- **标签**：
  - `F_vert`（单足峰值 GRF）— **measured**（OpenSim）；
  - `F_lat / M_inv / σ_lat / risk` — **modeled**（上界式 1D 梁弯曲+剪切线性叠加）；
  - `A_section = 101.94 mm²` — **measured**（THUMS CORT，`bc_robust_metric_route1.json`）；
  - `I = 826.9 mm⁴`、`c = 5.696 mm` — **modeled**（等效圆截面代理远端/外踝）；
  - `d = 30 mm`（CoP→距下轴力臂）— **assumed**；
  - `σ_c = 70 MPa` — **assumed**（[Y25] fibula_ends 压缩强度，用作弯曲极限；拉伸侧仅 30 MPa）。
- **诚实边界**：绝对值不可引用，只看相对与是否越阈（`NONVERTICAL_S5_SUPINATION_REPORT.md` §8）。

---

## 图 5  `fig5_pad_tradeoff.png` — 垫子刚度 × 旋后两难：无内部最优、k*→刚性极限

**图注（中文）**：(a) 4 个高度下 `total_risk(k)` 随垫子刚度 k（N/mm，log 轴）的变化，均为**单调递减**、**无内部极小**，`k*` 落在**刚性极限 k→∞**（右端 `∞(刚性)` 标记；黑色箭头“无内部最优 k*→刚性极限”）。红色虚线为 `total_risk=1` 阈值；刚性极限下 total_risk 由 h=2 m 的 0.659 升到 h=4.5 m 的 0.779（安全余量随高度收缩）。(b) h=4.5 m 的**分量分解**：轴向 axial（蓝）、侧向 lateral（红）、合计 total（黑）；红点线为阈值 k_th ≈ 2.37e+05 N/mm。**结论**：软垫省下的**轴向**风险（S1 实测最大折减 23.6%）远小于它因脚沉入垫（δ=F/k → β↑）引入的**侧向弯曲**风险 → 软端必输、无内部最优。

- **数据源**：`results/opensim_fe/pad_supination_tradeoff.json`（`sweep.cells`，中央 L_roll=40 mm、等效圆截面；`sweep.threshold_by_height`）。
- **标签**：
  - `R_FE(k)`（刚性比折减，FE 探针）— **measured**（S1 `pad_stiffness.json`）；
  - `F_peak(k)=F_peak_rigid·R_FE(k)`、`δ(k)=F/k`、`β(k)`、`axial/lateral/total risk` — **modeled**；
  - `L_roll = 40 mm`（β 几何杠杆）— **modeled**（灵敏度 25/60 mm）；
  - `d = 30 mm` — **assumed**；`σ_c = 70 MPa` — **assumed**。
- **稳健性**：L_roll∈{25,40,60} mm 三档、4 高度下均为“无内部极小、k* 落硬端”（`robust_no_interior_optimum=True`）。

---

## 图 6  `fig6_ligament.png` — 踝韧带判据：ATFL 关节级扭伤 risk 随 β（两口径）

**图注（中文）**：ATFL 关节级扭伤 risk（`risk = M_inv / M_fail`，log 轴）随旋后角 β 的变化，两口径并排。左：**保守口径**（无预载，`M_fail = 21 N·m`，Funk 2002）下 β≥5° 全部越阈；模型临界角 ≈ 2.61°（21 kN 冲击 GRF 下的输入旋后角）。右：**高承载口径**（2 kN 预载，`M_fail = 77 N·m`）下 β≥10° 越阈；临界角 ≈ 9.59°。红色虚线为扭伤阈值 `risk=1`。**结论**：保守口径下模型在真实抱石高度（2–4.5 m）内即预测 ATFL 扭伤，定性对齐 [Heck2024] 踝伤**扭伤 71.3%** / 骨折 27.4%，补上骨模型漏掉的 #1 主线。

- **数据源**：`results/opensim_fe/ankle_ligament_criterion.json`（`scenario.rows`，40 行 = 4 高度 × 5 角 × 2 口径；`critical_thresholds`）。
- **标签**：
  - `F_vert`（单足峰值 GRF）— **measured**（OpenSim）；`M_inv = F_lat·d` — **modeled**（承接 R1）；
  - ATFL 关节级 risk — **modeled**（`M_inv/M_fail` 主口径，避开 r_lig 力臂假设）；
  - `M_fail` 锚点（0 N→21 N·m、2 kN→77 N·m）— **文献 assumed**（Funk 2002；>2 kN 外推不确定，故并报两口径）；
  - 韧带极限载荷 ATFL 200 N / 应变 14% — **文献 measured**（Siegler 1988 等，本图未用该力/应变判据）；
  - `r_lig = 22 mm` — **assumed**（仅用于力/应变判据，本图未用）；
  - 模型临界角 2.61°/9.59° 是 **21 kN 冲击 GRF 下的输入旋后角**，**不等于** Funk/Parenteau 尸体测试失效角（33–40°），两者不可混用。
- **诚实边界**：力/应变判据为上界式（全部 M_inv 归单条韧带）；关节级为主口径（`ANKLE_LIGAMENT_CRITERION_REPORT.md` §9）。

---

## 图 7  `fig7_external_rotation.png` — 踝外旋使腓骨远端 / 外踝在真实抱石高度内越阈（SER 外踝臂）

**图注（中文）**：刚性地面下，腓骨远端 / 外踝的外旋 risk（`σ_vm / σ_c`，log 轴）随高度 h 的变化，5 条曲线对应 5 个外旋角 θ ∈ {0°, 5°, 10°, 15°, 30°}。**θ = 0°**（纯轴向）risk 严格为 0，曲线贴近图底（逐位复现纯轴向基线），真实抱石高度内安全；**θ = 5°** 在 4 个高度已全面越阈（外旋 risk ≈ 4.7 / 4.9 / 5.1 / 5.6，对应 h = 2 / 3 / 3.5 / 4.5 m）；**θ ≥ 10°** 在全部测点进一步加深（10° ≈ 9.4–11.1、15° ≈ 14.0–16.6、30° ≈ 27.0–32.1），4 个真实高度内**全部越阈**。橙色淡填充覆盖 h ∈ [2, 4.5] m = 真实抱石高度带；红色虚线 `risk = 1` 为失效阈值；**临界外旋角 θ_c ≈ 0.89–1.06°**（随 h 升高而单调减小，上界式 1D，绝对值不可引用）以小注记给出（褐色 `5a1f02`，区别于图 4 旋后）。**结论**：纯轴向 + 内翻（R1 / 图 4）口径下腓骨在 ≤ 4.5 m 内安全（首骨折 ≈ 7 m），**一旦叠加横截面外旋**（即 Lauge-Hansen trimalleolar 的 SER 机制，临床上占踝骨折 57–85%）即在同一真实高度内即越阈 —— 这正是纯轴向 + 旋后模型漏掉的 #2 临床机制，与图 4（旋后 / 冠状面）互补：S6 = 腓骨远端**横截面扭转**（螺旋 / 斜形骨折），R1 = 腓骨远端**冠状面弯曲 + 剪切**。

- **数据源**：`results/opensim_fe/ankle_external_rotation.json`（`scenario_hard`，24 行 = 4 高度 × 6 外旋角；`fibula_external_rotation_risk`、`critical_external_rotation_deg_this_h`）。
- **标签**：
  - `F_vert`（单足峰值 GRF）— **measured**（OpenSim `ground_reaction.f_right_n`）；
  - `T_ext = F_y · d · sinθ`（绕腓骨长轴的扭矩）— **modeled**（横截面扭转 1D）；
  - `J = 2·I = 1653.8 mm⁴`（等效圆极惯性矩）、`σ_vm = √3 · T_ext · c / J`（纯扭转 von Mises）— **modeled**；
  - `A_section = 101.94 mm²` — **measured**（THUMS CORT hex，`bc_robust_metric_route1.json`）；
  - `I = 826.9 mm⁴`、`c = 5.696 mm` — **modeled**（等效圆截面代理远端 / 外踝，无独立远端截面）；
  - `d = 30 mm`（CoP → 胫骨长轴水平偏距）— **assumed**（足部半宽 50–70 mm 中点偏内；灵敏度 15/30/45 mm）；
  - `σ_c = 70 MPa` — **assumed**（[Y25] fibula_ends 压缩材料强度，用作扭转极限；扭转 / 弯曲拉伸侧 ~30 MPa，本口径略偏乐观 → 结论只会更强）；
  - **θ_c ≈ 0.89–1.06°** 是 **modeled 上界式 1D 的输入临界外旋角**，**不等于** Funk/Parenteau 尸体测试真实外伤角（5–15°），两者不可混用。
- **与图 4（R1 旋后）对照**：两者**同为 1D 上界**、**都把同一条腓骨**推到真实抱石高度内越阈，但机制互补：R1 = 内翻（冠状面）弯曲 + 剪切；S6 = 外旋（横截面）扭转 + von Mises。临床 trimalleolar 三踝骨折的**外踝臂**（lateral malleolus）由本图 7 直接覆盖。
- **诚实边界**：上界式 / 1D 筛选（材料强度 + 中段截面 + 扭转 + 线性叠加）—— 忽略剪切修正、应力集中、韧带 / 腱卸载、接触摩擦与能量吸收；**绝对值不可引用**，只看相对与是否越阈（`ANKLE_EXTERNAL_ROTATION_REPORT.md` §10）。

---

## 图 8  `fig8_trimalleolar_ser.png` — 三踝骨折机制（Lauge-Hansen 旋后-外旋 SER；真实骨几何 + 受伤体位渲染）

**图注（中文）**：**J 后方观**（后侧）。基于真实下肢体网格（THUMS AM50 右下肢：胫骨 / 腓骨 / 距骨 / 跟骨 / 足骨）渲染的三踝骨折（trimalleolar）机制图。足部被摆在 **SER 受伤体位 = 内翻（旋后，脚心向内）28° + 外旋 20°**；胫腓骨保持不动（论文初稿 §3.6 / §4.2）。红色 **骨折线**由骨网格平面切片生成、贴合骨面：**后踝骨折**（胫骨远端后缘，PITFL 撕脱 / 距骨楔入）、**内踝骨折**（内踝基部）、**外踝骨折**（腓骨远端斜形，SER-II 位）。机制箭头以真 3D 弧线表示：**内翻**（绕足长轴、足跟下方大弧，表示脚心向内翻转）+ **外旋**（脚掌下方水平弧，箭头朝外，表示足绕小腿长轴外旋）+ **轴向载荷**（蓝色向下直箭头，软垫主要卸载的分量）。

- **数据源 / 几何**：`temp/opensim_fe/parts_*/…/*_surface.stl`（THUMS AM50 真实骨网格；胫骨用 gmsh 流形面 + Taubin 平滑）；骨折线 = 网格平面切片（`_fracture_v9.py`）；**无仿真实测数值**。
- **受伤体位依据（关键）**：三踝 / 后踝骨折最常见的机制是 **旋后-外旋 SER**——`supination = inversion（脚心内翻）`、`pronation = eversion（脚心外翻）` [Lauge-Hansen 1950]；队列数据：三踝骨折 SER **66.7%** vs PER 30.7%（1000 例 CT 队列）；后踝骨折中三踝-SER **71.0%**（472 例）；SER 占全部踝骨折 40–70%。故本图取**内翻 + 外旋**位，而非外翻。
- **标签**：距骨 / 胫骨远端 / 跟骨（黑色，引线指骨）；内翻 / 外旋 / 轴向载荷（机制箭头）；三条踝骨折红线（骨折部位）。
- **诚实边界**：本图为**机制示意**（真实几何 + 姿态），不携带量化"是否越阈"判断；定量见 **图 4（R1 旋后/内翻，冠状面弯曲+剪切）**与 **图 7（S6 外旋，横截面扭转 `risk=σ_vm/σ_c`）**；后踝臂 / 内踝臂由文献 [26,27] 支撑、未建模（论文 §5.2 第 3 条）。
- **与图 4 / 图 7 对照**：图 7 = 三踝**外踝臂**定量（外旋 / S6）；图 4 = 腓骨远端**冠状面**机制（旋后 / R1）；本图把两者收束到 SER 的**真实受伤体位**上，是图 7 的**概念伴随图**。

---

## 图 9  `fig9_talus_mechanism.png` — 距骨骨折机制（真实骨几何 + 背屈体位渲染）

**图注（中文）**：**G 内侧观**。与图 8 同源的真实骨几何渲染；足部被摆在 **强制背屈位（绕踝背屈轴 +ml 转动 24°）**，用于说明抱石落地（强制背屈 + 轴向撞击）下**距骨骨折**（距骨颈 / 体）的力学机制（论文 §3.6）。红色 **骨折线**为距骨网格**冠状面**切片生成、贴合骨面（A/G 视角呈竖线、D/J 呈环）。机制箭头：**轴向撞击**（蓝色向下直箭头，胫→距→跟载荷柱）+ **背屈**（绕踝横轴的真实 3D 弧箭头，趾端向上）。

- **数据源 / 几何**：`temp/opensim_fe/parts/talus_r/*` 等真实骨网格；骨折线 = 网格切片；背屈 = 距骨/跟骨/足绕踝横轴（+ml）旋转（`_fracture_v9.py`）；**无仿真实测数值**。
- **与其他图对照**：图 4 / 6 / 7 / 8 走**腓骨远端**链（旋后 / 内翻 / 外旋 → 弯曲 + 剪切 + 扭转）；本图走**距骨**的**胫–距–跟轴向载荷柱**链。**机制互补**：图 4 / 7 解释单脚踝为何绕开轴向模型；本图解释距骨骨折为何绕开旋后模型。
- **文献**：[28] Hawkins 1970 / [29] Canale & Kelly 1978（距骨颈骨折机制"飞行员距骨"）；[30] Peterson 1976（距骨 = 载荷柱唯一悬臂）；[31] Wong 2016（阈值 ≈ 6–10 kN）。
- **诚实边界**：本图为**机制示意**，不携带任何量化结论；箭头方向即力学作用方向。凡未在正文引用的数字 / 阈值，均不应被引用。

---

## 附：复现与文件清单

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\make_paper_figures.py
```

| 图 | 文件 | 数据源（只读） |
|---|---|---|
| 图 1 | `figures/fig1_pipeline.png` | （无数据，纯文献示意图，见论文初稿 §2.1 / §2.5） |
| 图 2 | `figures/fig2_alignment.png` | `axial_subregion.json` |
| 图 3 | `figures/fig3_single_foot.png` | `nonvertical_s2.json` |
| 图 4 | `figures/fig4_supination.png` | `nonvertical_s5.json` |
| 图 5 | `figures/fig5_pad_tradeoff.png` | `pad_supination_tradeoff.json` |
| 图 6 | `figures/fig6_ligament.png` | `ankle_ligament_criterion.json` |
| 图 7 | `figures/fig7_external_rotation.png` | `ankle_external_rotation.json` |
| 图 8 | `figures/fig8_trimalleolar_ser.png` | （无仿真数据；真实骨几何 + SER 受伤体位渲染，见论文 §3.6 / §4.2；外踝臂定量见 `ankle_external_rotation.json`） |
| 图 9 | `figures/fig9_talus_mechanism.png` | （无仿真数据；真实骨几何 + 背屈体位渲染，见论文 §3.6） |
| 动画 | `figures/fall_simulation.mp4` / `.gif`（+ `_passive` 对比版） | OpenSim 正动力学 states 轨迹（`run_dead_drop`，见附节） |
| 图 10 | `figures/fig_pad_paradox.png` | 「垫子悖论」：旋后着地时 total_risk(k) 单调递减 | `pad_supination_tradeoff.json`（R4，见附 2） |
| 动画 A | `figures/fall_scenario_ser.mp4` / `.gif` | 三踝 / SER 姿势演示（单脚 + 内翻30° + 外旋25°） | 见附 2 |
| 动画 B | `figures/fall_scenario_talus.mp4` / `.gif` | 距骨 / 背屈楔入姿势演示（单脚 + 踝背屈25°） | 见附 2 |

**通用诚实边界**：本组图仅复现既有 1D/FE 筛选口径的**相对结论与阈值跨越**，绝对值不可直接引用；`measured` 量来自 OpenSim/FE/THUMS，`modeled/assumed` 量的敏感性见对应 `*_REPORT.md` 的“假设与诚实边界”节。

---

## 图 1  `fig1_pipeline.png` — 管线总览（四段水平流 + 末端二判据分支）

**图注（中文）**：本文方法由四段串联构成（论文初稿 §2.1 + §2.5）：
**地/垫反力 GRF**（`pad.py`：1D 双质点 → 每足 3D 力向量）→ **OpenSim 多体正动力学**（关节 3D 反力，`subtalar_reaction`）→ **载荷传递**（关节 wrench → FE 边界，`load_transfer`）→ **FEBio**（单骨 von Mises 应力）。末端分支到两判据盒：**骨骼判据**（按子区域材料强度，论文初稿 §2.3） / **踝韧带判据**（ATFL/CFL 关节级内翻力矩，论文初稿 §2.4）。
**顶部横幅**：非垂直落地 6 自由度（S1–S4：①方向 ②地面倾角 ③姿势 ④单脚 ⑤力矩；⑥未做）。
**底部注记**：两端 1D 接口 → 升为 3D（中段本已 3D） — 即上游 `pad.py` / `opensim_fall._add_foot_force` 与下游 `load_transfer` 三处被写死成 1D 标量/无矩，本文扩展集中于把这两端升级为 3D 并保持物理一致。

- **数据源**：无 JSON 输入；本图为文献示意图（pipeline schematic），文字严格对齐论文初稿 §2.1（管线总览）与 §2.5（表 2 自由度）及 `docs/非垂直落地扩展方案.md` §0（关键洞察）。
- **标签**：
  - 阶段 1 `地/垫反力 GRF`：`pad.py` 1D 双质点 → 每足 3D 力向量 — **架构前提**（论文初稿 §2.2 地面/垫反力 + `docs/非垂直落地扩展方案.md` §1 A2/A3）；
  - 阶段 2 `OpenSim 多体正动力学`：AM50 全身模型坠落解算 → 关节 3D 反力与力矩（距下/踝/膝/髋/脊柱）— **架构前提**（论文初稿 §2.2 + `OpenSim_FE交接.md`）；
  - 阶段 3 `载荷传递`：关节 3D wrench → FE 边界（力+力矩）— **架构前提**（论文初稿 §2.2 + `docs/非垂直落地扩展方案.md` §1 A6）；
  - 阶段 4 `FEBio`：跟骨/胫骨/腓骨/股骨/骨盆/L3/T6/C5/颅骨单骨 FE（双材料：松质+皮质壳）→ von Mises — **架构前提**（论文初稿 §2.2 有限元段）；
  - 判据盒 1 `骨骼判据`：按子区域材料强度（MPa），子区域表见论文初稿 §2.3 表 1（[Y25] Table 1）；
  - 判据盒 2 `踝韧带判据`：ATFL/CFL 关节级内翻力矩 `M_inv=F_lat·d vs M_fail(预载)`（Funk 等），保守口径 21 N·m / 高承载口径 77 N·m。
- **诚实边界**：本图为**架构示意**，不携带任何数值；箭头方向即数据流方向，无物理假设。

---

## 附：坠落仿真动画  `fall_simulation.mp4` / `.gif`　·　对比版 `fall_simulation_passive.mp4` / `.gif`

**内容**：把**全身多体坠落仿真过程**渲染成动画（左：Rajagopal2015 全身骨模型 + 垫 + 地面；右：同步的 GRF(t) 脉冲曲线 + 时间游标）。两段：
1. **自由落体段**（24 帧）：身体从释放高度以 v = g·t 下落（纯运动学，非 FD），`h` 读数递减到 0；
2. **冲击段**（60 帧，慢放 ≈ 80×）：`run_dead_drop` 的**正动力学**积分结果（每帧从 states 轨迹取一个状态，`getTransformInGround` 逐刚体摆姿），GRF 曲线上的游标同步推进，`F_GRF` 读数实时更新。

**第二版渲染修正（2026-10-05，三点）**：
1. **垫子在地面上被压缩，而不是整块下沉**。垫厚取 `pad.CrashPad.thickness_m = 0.20 m`，**底面固定在 y = −0.20 m 的地面**、只有**顶面**随足下沉（画面另绘**未压缩参考线框**以便读出压缩量）+ 固定地面板作视觉锚点。此前版本把顶面画在 y=0 而骨穿板（"沉入垫子"）、以及顶面整块平移（"整块垫子下沉"）两个误读，均已消除。
2. **角标标注关节模式**：`joints: LOCKED (S1 rigid leg, [Y25])` 或 `UNLOCKED + passive stiffness (S3b)`，并实时打印 `knee_r / ankle_r` 屈曲角 —— 免得观众期待模型并不求解的屈膝。
3. **新增 `--passive` 对比版**：S3 path (b)（解锁 `LEG_FLEX_JOINTS` + `PassiveStiffness.nominal_landing()`）；实测冲击中踝/膝只弯 **≈1.3° / 0.7°**（画面角标即为实时读数），与 `NONVERTICAL_S3_POSTURE_REPORT.md`「被动刚度拦住了'飞'，但没换来'弯'」一致。

**参数（本次渲染）**：释放高度 **2.0 m**、质量 **77 kg（AM50/[Y25] 口径）**、**落在软垫上**（`on_pad=True`）；撞击速度 **v₀ = 6.26 m/s**，**峰值 GRF ≈ 34.8 kN**，冲击窗 **≈ 34 ms**。对照：同一工况落**刚性地面**峰值 **42.4 kN**、窗 **20 ms** —— 即垫子把峰值削低 ≈ 18%、并把冲击时间拉长 ≈ 70%（与论文 §3「垫子削峰」的定性结论一致）。

**顺带得到一个结论性数字**：FD 里身体在冲击窗内**下沉 ≈ 16.2 cm**（`probe_sink.py`：pelvis 0.940 → 0.781 m）；而垫子总厚只有 **20 cm** —— 即 **2 m 坠落已把垫子压到约 81%（接近触底）**。这与论文"垫子在中高坠落里会触底、不再非线性吸能"的论点是同一件事，可直接用于演示。

> ⚠️ **量级/口径保留（2026-10-06，依 S5 更新）**：上述 **16.2 cm 是 FD 自由行程**（无接触约束），其量级由 **1D `pad.py` 律**决定（该律自标"文献量级、非实测标定"）。**3D FE 脚-垫接触（S5）尚不能给出可信的绝对压入量**：接触解只在**名义应变 ξ ≤ 0.25**（下沉 ≤ 50 mm）稳健，该区间内 3D 相对 1D 单轴应力律的增强为 **≈1.6**（穿透 ≈ 0）；而 2 m 坠落的工况需 ξ≈0.85，**超出 tet4 泡沫稳健域**（同一配置 3 次给出 3 个不同结果、全部失败 ⇒ **不可复现**）。此前试出的 "15.14 cm / 2.70×" **已作废，不得引用**（`docs/S5_contact_convergence.md`）。
> ⇒ **"接近触底（约 81%）"只应作量级陈述**，依据是 **1D 律（~18 cm）**；**不要**给单一确定值，**也不要**引用 3D 接触的 15 cm / 2.70×。

> **⚠️ 引用限定（2026-10-06，针对本节提及的 ≈1.6、15.14 cm、2.70×、16.2 cm）**
>
> **(i) 配方**：「≈1.6」是 **OLD AUGLAG pen=1** 配方在 tet4 (3,3,2) 单网格上、ξ ≤ 0.25
> 区间的平滑值；15.14 cm / 2.70× 来自 **NEW PENALTY 0.1** 配方在 sink=170 mm tet4 (6,6,4)
> 单次测量（事后证不可复现）；16.2 cm 来自 **FD 自由行程**（`pad.py` 1D 双质点律），
> 与 FE 接触无关。
>
> **(ii) 网格未收敛**：NEW PENALTY 在 3 网格 (3,3,2)→(6,6,4)→(12,12,8) 上 σ/σ_law
> 漂移 **+30.62 %**（`S5_g7_verify_mesh.md` §3），固定 sink 时 `F_n ≈ 34.4 kN` 与
> `σ_contact ≈ 0.574 MPa` 漂移 **<0.6 %**（当前唯一在 2 m 落地可引用且网格收敛的量）。
> OLD 配方未做 mesh-convergence 扫描。
>
> **(iii) 不可跨配方对比**：NEW vs OLD σ/σ_law 在 4 个 anchor sinks 上单调偏移
> **−5.85 %…−9.38 %**（`S5_g7_verify_anchors.md` §4.1），两条曲线不可叠加；
> OLD 配方装上 NEW stepper 在 closure `negJac=55` 立即失败 ⇒ 两配方**不可互换**。
>
> **(iv) 权威裁决**：`docs\S5_QUOTABLE.md` —— 2 m 落地**唯一可引用且网格收敛**的量是
> **`F_n ≈ 34.4 kN`** 与 **`σ_contact ≈ 0.574 MPa`**；任何 σ/σ_law 比值（含 2.205）
> 一律不得作为「收敛」单一值引用；15.14 cm / 2.70× 已撤销；本节保留为**事件时间线记录**，
> 实际引用数字一律以 `S5_QUOTABLE.md` 为准。

**数据源 / 复现**：
- 地面反力 `pad.py` 1D 双质点 → `ground_reaction(height_m, mass_kg, on_pad)`（`src/climbing/coupling/opensim_grf.py`）；
- 正动力学 `run_dead_drop(grf, scale_mass_kg=77)`（`src/climbing/coupling/opensim_fall.py`）：把 GRF 以 `PrescribedForce` 注入，RungeKuttaMerson 积分，**下肢关节锁死**（[Y25] 的"直立绷直腿"刚性柱近似）；
- 渲染脚本 `scripts/opensim_fe/make_fall_animation.py`（= `temp/opensim_fe/_fall_anim.py`）。

```powershell
$env:PYTHONPATH="src"
.venv\Scripts\python.exe scripts\opensim_fe\make_fall_animation.py --height 2.0 --on-pad --fps 25
.venv\Scripts\python.exe scripts\opensim_fe\make_fall_animation.py --height 2.0 --on-pad --passive --fps 25
.venv\Scripts\python.exe temp\opensim_fe\probe_sink.py      # 诊断：下沉量 / 冲量-动量
```

**诚实边界**：
- 这是 **S1「竖直 + 对称双足 + 绷直腿」**工况（`docs/OpenSim_FE交接.md` 已声明的适用范围），**不是**真实抱石坠落的完整运动学（无实测动作捕捉 → 见 `docs/OpenSim_FE复现方案.md` §A「凭空给轨迹 = 伪造物理」）；自由落体段为纯运动学插值，仅冲击段是 FD 解算。
- ⚠️ **FD 里没有接触约束**：GRF 是 `PrescribedForce`（力），不是"挡脚的垫子"。所以身体在冲击窗内的 **16.2 cm 下沉是 FD 自由行程**，不是垫子本构算出的压缩量 —— 两者在本动画里被当作同一件事来画（因为腿锁死时二者相等），**这是渲染约定，不是求解结果**。真实的"垫子压缩 vs 身体下沉"分离需要 S5（FE 脚-垫接触）。
- **OpenSim 刚体不会弯**：骨的变形/应力属 FEBio 单骨 FE 阶段，本动画只体现刚体运动；`--passive` 也只让关节弯 ~1–2°。
- 模型质量已缩放到 77 kg（`scale_model_mass`）；GRF 是 1D 等效模型的**载荷边界近似**，非实测。
- 峰值 GRF 数量级很大（数十 kN），属**上界式**口径，与论文 §3 的"上界 / 相对结论"一致；**绝对值不可直接引用**。

---

## 附 2：危险姿势演示（2026-10-05）

**背景**：论文 §3 的「软垫盲区」论点 —— 垫子只吸**轴向**能量，不吸**旋转 / 楔入**。下面一图两动画把它坐实。

| 产物 | 内容 | 关键数字 |
|---|---|---|
| `fig_pad_paradox.png` | **R4 数据**：踝总风险 `total_risk(k)` vs 垫子等效刚度 k（4 个高度） | 刚性极限 risk **0.66–0.78**（=纯轴向）；软垫（k ≤ 1000 N/mm）risk **26–45**；安全阈值 **k_th ≈ 1.0–2.4×10⁵ N/mm** → **越软越危险、无内部最优** |
| `fall_scenario_ser.mp4` | **三踝（SER）姿势**：单脚先落 + 内翻 β=30° + 外旋 θ=25°，h=2 m，垫 vs 刚性 | **单脚口径**：M_inv 峰值 ≈ **518 N·m = 24.7 × 21 N·m**；削峰 40.0 → 35.2 kN |
| `fall_scenario_ser_twofoot.mp4` | 同上，但 **S5 报告口径**（`--two-foot`：split=(0.5,0.5) + 75.337 kg） | M_inv 峰值 ≈ **264 N·m**（≈ 报告 256.9）→ 供与 R1/S5 数字逐条对齐 |
| `fall_scenario_talus.mp4` | **距骨（Hawkins 背屈楔入）姿势**：单脚 + 踝背屈 +25°，h=2 m | 削峰 **42.4 → 34.8 kN**、窗 20 → 33 ms；**无量化判据**（机制演示） |

### 与计算的一致性审计（必读）

| 动画元素 | 一致？ | 说明 |
|---|---|---|
| 冲击段姿态 | ✅ **真 FD 解** | `run_dead_drop` 的 `states`，逐帧 `getTransformInGround` |
| GRF 曲线 / `F_GRF` 读数 | ✅ | 就是 `ground_reaction`（pad.py 1D）那条力，也正是驱动 FD 的 `PrescribedForce` |
| 峰值 / 时长 | ✅ | 垫 35.2 kN / 34 ms；刚性 40.0–42.4 kN / 20 ms |
| `--passive` 关节读数 | ✅ | 真解出的 `knee_angle_r` / `ankle_angle_r` |
| SER 的 `M_inv` | ⚠️ **公式一致、载荷口径由我选** | 用同项目 `supination_load(F_right, 30°)`（= S5 的 `max(grf.f_right_n)` 口径）。**默认单脚** `split=(1,0)` + 77 kg → 518 N·m；**`--two-foot`** `split=(0.5,0.5)` + 75.337 kg → 264 N·m ≈ 报告的 256.9。**两者差 ≈2×，画面角标已写明口径** |
| 自由落体段 | ❌ **非计算** | 纯运动学 `h−½gt²` + 初始位形插值；画面已标「〔非 FD〕」 |
| 垫子压缩 15–18 cm | ❌ **非计算** | FD 无接触约束 → 是**自由行程**；画面已标「〔渲染约定〕」。真值需 S5（FE 脚-垫接触） |
| SER 的「内翻 / 外旋」姿态 | ❌ **未渲染** | `subtalar_angle` / `hip_rotation` 锁死，画面足部为中立位；姿态见 `fig8` |
| 右侧「刚性」虚线 | ⚠️ **仅载荷对比** | 曲线是刚性 GRF，但**姿态仍是垫工况解出的** → 只能读作载荷对比，不是两套姿态 |

**版面**：主视图（全身坠落 + 垫压缩）+ **踝部特写内嵌**（镜头跟右脚距骨）+ 右侧 **垫 vs 刚性 GRF 曲线**（灰虚线=刚性 / 红=垫）。角标实时给出 `M_inv / risk`（SER）或 `踝背屈角`（距骨）。

**诚实边界**：
- **SER 动画里的内翻 / 外旋不是模型自由度**（`subtalar_angle` / `hip_rotation` 锁死）→ 画面里的足部**没有**内翻姿态；`M_inv` 是用 `ankle_supination.supination_load(F, 30°)` 由**载荷时程**解析算出的（项目 R1 / S5 口径）。内翻足姿态见 `fig8`。
- 距骨**无量化判据**、`subtalar` 未解锁 → 只能机制演示。
- 两者都是 **S1「竖直 + 绷直腿」**工况的扩展（单脚 + 锁定姿势），不是真实抱石坠落运动学。
- `posture="one-leg-awkward"` 与 `split=(1,0)` 只改**载荷边界**（接触面积 / 单足分配），不改 OpenSim 关节姿态。

**复现**：
```powershell
$env:PYTHONPATH="src"
.venv\Scripts\python.exe scripts\opensim_fe\make_pad_paradox_figure.py
.venv\Scripts\python.exe scripts\opensim_fe\make_scenario_animations.py --scenario both
.venv\Scripts\python.exe scripts\opensim_fe\make_scenario_animations.py --scenario ser --two-foot  # S5 口径
```
> 性能：渲染器**复用同一个 `pv.Plotter`、逐帧只重写点坐标**（不重建 81 个 actor）；两个场景合计 **≈ 93 s**。早期每帧重建 Plotter 的版本需要 **≈ 15 min**。`--quick` 可只出 7 帧做冒烟测试。
