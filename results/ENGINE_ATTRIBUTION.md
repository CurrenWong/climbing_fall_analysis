# 引擎归属表：**哪些是有限元算的、哪些是 1D 算的**（2026-10-16）

> **触发**：用户要求「说明哪些是有限元算的，哪些是 1d 刚度模型算的」。
> **性质**：**清点 + 归属**。**零新求解、不改任何数字与判定。**
> **目的**：堵住本项目最大的口径混淆源 —— **同一结论在不同引擎下产出，却常被并排引用。**

---

## §0 · 一句话

**本项目有**三个**引擎（不是一个），外加一个第三方模型：**

| 引擎 | 是什么 | 谁写的 | 产出什么 |
|---|---|---|---|
| **① 1D / 自研多体** | 质量-弹簧、平面多刚体、截面梁 + 强度判据 | **本项目自研**（纯 Python） | **时刻 / 高度 / 概率 / 分类**（V1–V9 的绝大多数） |
| **② OpenSim 多体动力学** | Rajagopal2015 全身肌骨模型 | 第三方 | **载荷**（GRF、距下关节反力 `F_sub(h)`）—— **只供载荷，不出应力** |
| **③ FEBio 有限元** | 单骨/多骨连续体 FE | 第三方 | **应力场 σ_vm / 接触力**（V1 的 FE 支、A1–A4、S5、踝子线） |

**⇒ 最容易混的一句：**「载荷从哪来」（②）和「应力谁算」（③）是**两个引擎**，中间由 `joint_loads.py` 桥接。

---

## §1 · ★ 验证靶 V1–V9 的引擎（**表 5 逐个标注**）

| 靶 | 判据内容 | **引擎** | 依据 |
|---|---|---|---|
| **V1** | 足部骨折临界高度 6–9 m | **③ FE ＋ ① 1D（并列，且互不重合）** | 「1D 纯轴向+子区域强度：足 8 m、腓 7 m」＝**①**；「其余口径 2–26.5 m」含**③**（单骨 FE `gauge_max`=2 m、踝 FE 26.5 m）⇒ **跨口径 ~9×** |
| **V2** | 「低于 7 m 无骨折」 | **① 1D 截面 × 组 A（唯一依据）＋ ③（并列反证）** | ✅ 来自 `validate_vs_fem.py`（**1D**，组 A＝胫骨干/胫远端/股骨颈）；**③** 矩阵 `median` **3/8 违反** ⇒ 「FE 线未通过」 |
| **V3** | 头部应力峰值时刻 ≤12.5 ms | **① 1D 垫模型** | `pad.RIGID_STANDING_CHAIN`、`k_neck`、`alpha_pad` 全在 `src/climbing/pad.py`；**无 FE** |
| **V4** | 两端先断 vs 骨干 | **① 1D** | 1D 截面 + 子区域强度；含载荷分配系数 |
| **V5** | 加权部位分布（下肢 67%、踝第一） | **① 1D** | [4] 为抱石=垫队列；`ankle_injury.py` 出踝损伤率 |
| **V6** | 垫刚度 vs 韧带应变 非单调 | **① 1D 垫模型 ＋ 韧带解析** | `pad.py` + `ankle_ligament.py`；9/9 |
| **V7** | 阶跃 vs 渐进分类 | **① 1D** | 分类量在 1D 矩阵上算 |
| **V8** | 峰值时刻随高度提前（@20 m ≤5 ms） | **① 1D 链** | `chain2d.py`（6 段软关节链）+ `pad.py`；**V8b 的 14.0 ms 是 1D 链结构性极限** |
| **V9** | 内翻轻于外翻 | **① 1D** | `ankle_rotation.py` 力臂编码；6/6 |

**⇒ ★ 结论：V1–V9 里**只有 V1 和 V2 涉及 FE，而且两者都是「并列」而非唯一依据；V3–V9 全部是 ① 1D。****

---

## §2 · **纯 FE（③）的独立产线**（这些数字**只能**来自 FE）

| 产线 | 链上的件 | 关键产物 |
|---|---|---|
| **单骨跟骨线（V1 的 FE 支）** | `thums_feb.py` / `s2_thums_sweep.py` / `thums_feb_run.py` | `gauge_max`、首折、S2/S3 族、**本会话的 S10 / 跖面 BC 族 / 界面锐度** |
| **A1–A4 骨折矩阵** | `fracture_matrix.py` / `rerun_matrix_bc2.py` / `s4_alignment_bc2.py` | `fracture_matrix_bc2.json`、A2 5/5、A3 3/3 |
| **端部 BC（P2 大件②）** | `g4_fe_plate.py` / `g4_fe_p1b.py` / `p2_t6_*.py` / `fe_bc_*.py` | 压板 z 注册 bug、T6 能力门 |
| **P5 式接触 / 稳定化** | `g4_fe_stab.py` / `g4_fe_pilot*.py` / `g4_fe_area*` | G4②（BC 屈曲 `P_cr≈252 N`） |
| **S5 脚-垫接触** | `nonvertical_s5_contact*.py` / `proof_s5_facet_fix.py` | `F_n≈34.4 kN`、`σ_contact≈0.574 MPa`（**唯一可引用且网格收敛的 FE 量**） |
| **踝关节 FE 子线** | `ankle_ligament_criterion.py` / `ankle_*` 相关 | 双靶 2.96 MPa / 287.4 mm²（**相对可用、绝对不可用**） |
| **G5/G6/G8/G9** | `g4_damage_evolution.py` / `g8_fibula_or_share.py` / `g9_*` | 粘弹性窗、腓骨载荷分享 |
| **网格/分层/接触基础设施** | `meshing.py` / `part_meshing.py` / `plantar_bc.py` / `joint_contact.py` / `pad_contact.py` / `fe_post.py` / `febio_run.py` / `xplt_shim.py` | — |

## §3 · **纯 ① 1D / 自研**（这些数字**只能**来自 1D）

| 模块 | 位置 | 产出 |
|---|---|---|
| 1D 垫-人（单自由度） | `src/climbing/pad.py` | V3、V6、V8a/b、冲击窗、`alpha_pad` |
| 2D 平面（2 质量 + 1 转动） | `src/climbing/pad2d.py` | **P3 首切**（4 不变量） |
| 2D 多刚体链 | `src/climbing/chain2d*.py` | **P3 二切**（G1 6 段链，DOF=12） |
| 1D 骨强度/截面 | `src/climbing/bone.py` + `scripts/opensim_fe/risk_1d_bending.py` / `risk_1d_loadshare.py` / `axial_subregion.py` / `paper_alignment.py` | **V1 的 1D 支（足 8 m/腓 7 m）**、V4、V7 |
| 1D 截面 vs FE 对照 | **`scripts/paper/validate_vs_fem.py`** | **V2 的 ✅ 唯一依据**（组 A 2/3） |
| 踝解析式 | `src/climbing/coupling/ankle_ligament.py` / `ankle_rotation.py` / `ankle_supination.py` | V6、V9、R1（旋后）、S6（外旋） |
| 损伤/统计 | `src/climbing/injury.py` / `ankle_injury.py` / `metrics.py` / `scenarios.py` | V5、踝骨折率 22.8% |
| 头颈分支 | `src/climbing/head_branch.py` | V3 相关 |
| 主控/批处理 | `src/climbing/batch.py` / `bodies.py` / `rope.py` / `belay.py` / `wave_rod.py` | 场景与回归 |

## §4 · **桥（② → ③）**：载荷链

| 件 | 作用 |
|---|---|
| `coupling/opensim_grf.py` | ② OpenSim → **地面反力 GRF**（`ground_reaction(height_m=)`） |
| `coupling/opensim_fall.py` | ② 坠落仿真 → 关节反力时程 |
| `coupling/joint_loads.py` | **②→③ 桥**：取距下关节**峰值竖向反力** `F_sub(h)` 交给 FE |
| `coupling/loadchain.py` + `scripts/opensim_fe/loadchain_audit.py` | **闸门 10「基线同链」**：审计某产物是否与当前载荷链同链 |

**⇒ 引用任何 FE 数字，必须同时声明**载荷链口径**（本会话已发现 2 份异链产物，形状跨度 17.49%）。**

---

## §5 · 三个**同名不同物**、最易混的量（务必分清）

| 概念 | ① 1D 里的名字 | ③ FE 里的名字 | 说明 |
|---|---|---|---|
| **垫的刚度** | `alpha_pad`（把泡沫三应力项同乘一个因子） | `pad_e_mpa`（对偶面模量）/ `k`（Winkler 接触刚度） | **不是同一个量**，不可互算 |
| **首折高度** | 逐高度二值（σ>σ_c 的首个高度） | 同一形式，但 σ 来自 FE 场 | 口径同为「首个越阈高度」，但**应力来源不同** |
| **判据阈值** | `bone.py` 的 `*_ENDS=70` 等**子区域**强度 | 皮质域 `gauge_max > 150 MPa`（`THRESHOLD_MPA`） | **同一 ✓ 同骨 σ_c 也分口径**（如 `fibula` 矩阵里 160 MPa，`bone.py` 另有 `FIBULA_ENDS=70`） |

---

## §6 · 使用规则（建议提升为项目级）

1. **每个数字必须带引擎标签**：`[FE]` / `[1D]` / `[OpenSim 载荷]`，并带**载荷链口径**与**网格成员**。
2. **同一结论若两引擎都算过且不一致 ⇒ 必须并列书写**（例：V2 的 1D ✅ ＋ FE 线未通过）。
3. **不得跨引擎换算**（1D 的 `alpha_pad` ↔ FE 的 `pad_e_mpa`；1D 首折 ↔ FE 首折）。
4. **不得用 ① 的结果充当 ③ 的依据**（例：不得用 1D 足 8 m 去证明 FE 单骨线的首折）。
5. **`F_sub(h)` 属 ②**，不是 ① 也不是 ③；它本身的不确定度会同时进入两条线。

---

## §7 · 溯源

| 项 | 位置 |
|---|---|
| 验证靶表（V1–V9） | `docs/paper/论文v2_数学建模版.md` 表 5 |
| V2 的 1D 依据 | `scripts/paper/validate_vs_fem.py`（组 A/组 B 定义在 `:206-215`）、`results/vs_fem_thresholds.txt` |
| V1/V2 口径裁定 | `results/opensim_fe/V1V2_CALIBER_RULING.md` |
| 载荷链审计 | `src/climbing/coupling/loadchain.py`、`scripts/opensim_fe/loadchain_audit.py`、`results/opensim_fe/LOADCHAIN_RESIDUAL_CLOSURE.md` |
| FE 基础设施 | `src/climbing/coupling/{febio_model,febio_run,fe_post,xplt_shim,meshing,plantar_bc,joint_contact}.py` |
| 1D 主体 | `src/climbing/{pad,pad2d,chain2d,chain2d_soft,bone,injury,ankle_injury,head_branch}.py` |
| 阈值定义 | `scripts/opensim_fe/s2_thums_sweep.py:66`（FE 侧）、`src/climbing/bone.py`（1D 侧） |
