# S5-contact 波5 —— 真实泡沫垫矩阵：G7 配方 + hold 弹簧在 real THUMS 网格上的验证

> 建立于 2026-10-06（波5）。本文件是**波5 的耐久计划 + 逐步日志**（prior subagent 在
> 2h14m 撞配额，只有落盘内容幸存 —— 故本文件先写、每步追写）。
> 只读/追加约束见任务书；**不覆盖任何既有 `results/opensim_fe/*`**。

---

## 0. 一句话目标

波4 的 7 条真实 `contact` 行全部 `rc=1 / end_t=None / negJac=6 / invalid=4`，**对 μ、载荷、
pad 场景不变量** —— 是 deck 构造缺陷（tri3 facet + FEBio 默认 `<Control>` + 缺 hold 弹簧），
不是物理。commit `315c373` 在 **rigid-plate proof case** 上修了三处（quad4 + G7 配方 +
opt-in hold 弹簧）。本波：**把同一修复接进真实泡沫垫矩阵**，逐 case 判定 contact 是否
CONVERGED 且 CARRIES LOAD。

## 1. 依赖的接口事实（已核实，勿重查）

| 事实 | 值 / 出处 |
|---|---|
| `build_thums_feb` 签名 | `scripts/opensim_fe/thums_feb.py:641-685`。**无 `penalty` kwarg**（故 penalty 只能对 emitted `.feb` 做文本 `re.sub`，与 `proof_s5_facet_fix.py:73-75` 同法） |
| G7 守卫 | `thums_feb.py:777`：`g7_solver_recipe=True` 仅允许 `analysis=="STATIC" and int(time_steps)==1`。`analysis` 默认 `"STATIC"`（L659），`_solve_case` 传 `time_steps=1` ⇒ 只需显式加 `analysis="STATIC"` 保险 |
| G7 emit | `thums_feb.py:911-933`：`time_steps=2400` / `step_size=dtmax=1/2400` / `cutback=0.125` / `max_retries=20` / `opt_iter=15` |
| quad4 自动重建 | `thums_feb.py:984-985`：`plantar_bc=="contact" and "elements_cort" in mesh and "boundary_polys" in mesh` ⇒ `_build_thums_plantar_quads`；单测 `test_thums_feb_real_plantar_is_quad4_in_deck` 断言 **146 quad4 / 0 tri3** |
| hold 弹簧 kwarg | `contact_hold_spring_k`（`thums_feb.py:684`, → `plantar_bc._apply_spring`，Winkler `1 N/mm`/plantar节点，方向 normal+x+z） |
| load ramp | `load_ramp_time_s`（静态 ramp 0→1，`thums_feb.py:966-977`），proof 用 `1.0` |
| `F_n` 反应求和**不可信** | FEBio 4.13 本 deck 不报固定 DOF 反力（`proof_s5_facet_fix.py:31-35` 实测；即使收敛的 `fixed` 参考 deck 也读 ~0）。波4 脚本 `nonvertical_s5_contact.py:310-340` 的 `F_n` 提取是**已知缺陷路径** ⇒ 本波只报 **INFERRED**（力平衡：静力下 `F_n = 施加载荷`），显式标注 |
| QUOTABLE 锚 | `F_n ≈ 34.4 kN` @ 2m（tet4 pad / 固定 sink / NEW PENALTY，网格收敛 <0.6%）。本矩阵载荷 = 每高度距下关节力（h=2.0 → **21086.2 N**）。**不同 mesh / 不同 load path ⇒ 只作量级锚，不并列** |
| 无接触 sink | 21086 N / 191 N/mm ≈ **110 mm**（任务书）。判定"承载"用它：penetration ≪ 110 mm |

## 2. 已落盘的前置证据（波5 开工前的免费证据）

`temp/opensim_fe/s5_contact/proof_h2.0_contact_rigid_mu0.6_pad_pen*.log`（real THUMS 网格，
pad 场景，G7 + hold k=1 + quad4）：

| penalty | end_t | negJac 警告 | invalid facets |
|---|---:|---:|---:|
| 0.1 | **1.0** | 73 | 0 |
| 1 | **1.0** | 6 | 0 |
| 100 | **1.0** | **0** | 0 |
| 10000 | **1.0** | **0** | 0 |

⇒ 真实网格上 G7+hold+quad4 **收敛到 end_t=1**；`invalid facets` 归零；penalty≥100 negJac 归零。
波5 待补：**penetration 数值 + 承载判定 + hold 弹簧必要性**。

## 3. 波5 case 计划

接触子集（先前失败的 7 条）+ hold 弹簧必要性对照（2 条）：

| id | h(m) | μ | pad_scenario | penalty | hold k (N/mm) | 目的 |
|---|---:|---:|---|---:|---:|---|
| h2.0_contact_rigid_mu0.0_hard_surface | 2.0 | 0.0 | hard_surface | 100 | 1.0 | 硬面 + μ=0 |
| h2.0_contact_rigid_mu0_pad | 2.0 | 0.0 | on_pad | 0.1 | 1.0 | 泡沫 pad + μ=0 |
| h2.0_contact_rigid_mu0.3_pad | 2.0 | 0.3 | on_pad | 0.1 | 1.0 | 泡沫 pad + μ=0.3 |
| h2.0_contact_rigid_mu0.6_pad | 2.0 | 0.6 | on_pad | 0.1 | 1.0 | 主 case（proof 已证收敛） |
| h3_contact_rigid_mu0.6_pad | 3.0 | 0.6 | on_pad | 0.1 | 1.0 | 高度扫描 |
| h3.5_contact_rigid_mu0.6_pad | 3.5 | 0.6 | on_pad | 0.1 | 1.0 | 高度扫描 |
| h4.5_contact_rigid_mu0.6_pad | 4.5 | 0.6 | on_pad | 0.1 | 1.0 | 高度扫描 |
| **h2.0_contact_rigid_mu0.6_pad_nohold** | 2.0 | 0.6 | on_pad | 0.1 | **None** | hold 是否必需 |
| h2.0_contact_rigid_mu0.6_pad_pen100 | 2.0 | 0.6 | on_pad | 100 | 1.0 | penalty 敏感性 |

**penalty 政策（场景相关）**：on_pad（泡沫垫）→ `0.1`（G7 配方为泡沫标定，软罚允许较大
penetration）；hard_surface（硬面）→ `100`（proof 证硬面 0.1 太软、penetration ~39mm；
`penalty>=100` 使 negJac=0）。

## 4. 判定标准（acceptance = 承载，不是 rc=0）

每个 contact case：
- **CONVERGED**：`febio_rc==0` 且 `end_t==1.0` 且 `NORMAL TERMINATION`。
- **CARRYING**：`penetration_mm` ≪ 无接触 sink 深度（≈110 mm @ h2.0）。承载时 contact
  刚度随 penalty 上升，penetration 有界在几何 gap 量级；不承载时骨会一直沉到 hold 弹簧
  单独平衡载荷（≈ 110 mm 量级）。
- **UNLOADED**：收敛但 penetration ≈ 无接触 sink（≈110 mm）。
- **DIVERGED**：`rc!=0` / `end_t<1`。
- `F_n`：**INFERRED**（力平衡，静力 = 施加载荷），绝不当 measured 报。`F_n` 的
  FEBio 反应求和路径已知失效，仅记录原始读出以便对照，不作为证据。

## 5. 耐久日志（每步追写）

- [T0] 2026-10-06 建立本文件；核实 §1 接口事实；发现 §2 前置证据（proof 4 penalty 全 end_t=1）。
- [T1] 读取 `proof_s5_facet_fix.py` / `nonvertical_s5_contact.py` / `thums_feb.py` /
  `plantar_bc.py` / 5 套门禁测试；确认 `_solve_case` 未传 `g7_solver_recipe`。
- [T2] **代码改动**：`nonvertical_s5_contact.py::_solve_case` 加 `penalty` /
  `contact_hold_spring_k` / `g7_contact` 参数（contact 分支注入 G7+ramp+hold，落盘后
  `re.sub` penalty）；新增 `_override_penalty` / `_contact_metrics`；contact 的 `F_n` 改
  标 **INFERRED**（力平衡），reaction-sum 单列 `unreliable_diagnostic`；新增 `plantar_drop_mm`
  / `penetration_mm` / `contact_area_mm2` / `no_contact_sink_mm` / `carries_load` /
  `contact_verdict`。**新增** `scripts/opensim_fe/nonvertical_s5_contact_wave5.py`（驱动器）。
- [T2a] **deck 核验（anatframe 矩阵网格，无 FEBio）**：跖面 Surface = **137 quad4 / 0 tri3**
  （⚠️ 任务书写的 146 是 proof 用的 `calcnframe` 网格；矩阵用 `anatframe`，故 137。
  两者都是 quad4 / 0 tri3，缺陷 A 修复生效）。`<Control>` = `time_steps=2400` +
  `dtmax=1/2400` + `cutback=0.125` ✔；`<penalty>0.1</penalty>` 覆盖 ✔；hold 弹簧 block ✔。
- [T3] **单 case 验证（3 次实跑）——关键发现：加载路径是剩余阻塞**：
  | 试验 | 网格 | use_rigid | hold | rc | end_t | negJac | 结论 |
  |---|---|---:|---:|---:|---:|---:|---|
  | `h2.0_contact_rigid_mu0.6_pad`（波4 原配置：pressure） | anatframe | False | 1.0 | 1 | 0.419 | 0 | **DIVERGED**（闭合处） |
  | `anat_rigid_hold` | anatframe | **True** | 1.0 | **0** | **1.0** | 29 | **CONVERGED+CARRYING**（pen 40.5mm，drop 41.8mm，are 2953.8mm²） |
  | `anat_rigid_nohold` | anatframe | True | **None** | 1 | None | 21 | **DIVERGED**（2.7s，闭合前刚体模态） |

  ⇒ (1) **hold 弹簧必需**（无它立即发散，与 proof 一致）；(2) **波4 接触行用
  `use_rigid=False`（PressureLoad）是剩余阻塞**：G7+hold+quad4 已修复 facet/solver，
  但压力加载路径在接触闭合（t≈0.419）处 FEBio 收敛判据退化为「零残差」→
  `Max nr of reformations reached` → ERROR（1101 步已收敛，第 1102 步失败）。
  proof（calcnframe）用的是 `use_rigid=True`（幽灵刚体 RigidForceLoad），矩阵用
  pressure ⇒ **proof 与矩阵的真实差异是加载路径，不只是「rigid plate vs foam pad」**。
  对照：`preact_contact_pressure_hold.log`（prior agent，pressure 路径）也是 end_t=0.285 失败；
  `final_default_hold.log` / `final_swap_hold.log`（rigid 路径）end_t=1。
  ⇒ 波5 contact 行统一改用 `use_rigid=True`（= proof 与 S4 主口径），并保留 1 条
  pressure 行（`_press`）作阴性对照。
- [T4] **全矩阵实跑完成（17 行）** → `results/opensim_fe/nonvertical_s5_contact_wave5.json`
  + `NONVERTICAL_S5_CONTACT_WAVE5_REPORT.md`（**未触碰** `nonvertical_s5_contact.json`）。
  结果：
  - 6 条非接触参照（fixed/spring）rc=0 end_t=1.0（不变）。
  - **7 条先前失败 contact 行全部 CONVERGED+CARRYING**（rc=0/end_t=1.0）：
    pad h2.0 μ0/0.3/0.6 → penetration 40.45 mm；h3/3.5/4.5 → 42.3/44.0/48.1 mm；
    hard_surface μ0（pen100）→ 2.32 mm。contact_area 2953.8 mm²（pen0.1）/130–146 mm²（pen100）。
  - `_nohold` → **DIVERGED**（rc=1，3 s）⇒ **hold 弹簧必需**。
  - `_pen100` → CONVERGED+CARRYING，penetration 2.19 mm，negJac=0。
  - `_press`（波4 原 `use_rigid=False`）→ **DIVERGED**（rc=1，end_t=0.419）。
  - 承载判据：无接触 sink = 施加载荷/191 N·mm⁻¹（h2.0 → 110 mm）；所有承载行 penetration ≪ sink。
- [T5] **门禁全绿**：S5 `test_pad_*` **93 passed**；G0 `test_opensim_fe+test_febio_run`
  **19 passed / 1 skipped**；`test_nonvertical_s4.py` **11**；`test_plantar_contact.py` **11**；
  `test_nonvertical_s5_contact.py` **13**；`test_s5_contact_facet_fix.py` **15**。
- [T6] 清理：删除 `temp/opensim_fe/s5_contact/` 下 21 个 >50 MB 的 `.xplt`/`.hdf5`
  （可再生；`.feb` + `.log` 证据保留），目录 8.0 GB → 0.20 GB。

## 6. 结论（per case）
- 7 条先前失败 contact 行：**CONVERGED+CARRYING**（需 `use_rigid=True` + hold k=1 + G7）。
- `_nohold`：**DIVERGED**（hold 必需）。
- `_press`（PressureLoad 路径）：**DIVERGED**（加载路径是剩余阻塞）。
- 无一例 CONVERGED+UNLOADED。

## 7. 给出 proof 与矩阵的真实差异
proof（calcnframe）与矩阵（anatframe）的差异**不止**「rigid plate vs foam pad + hold」：
1. mesh 帧不同（calcnframe 146 quads vs anatframe 137 quads）；
2. **加载路径不同**（proof `use_rigid=True`；波4 矩阵 contact 行 `use_rigid=False`）。
第 2 点是矩阵仍发散的根因；hold 弹簧在两种网格上都必需。



