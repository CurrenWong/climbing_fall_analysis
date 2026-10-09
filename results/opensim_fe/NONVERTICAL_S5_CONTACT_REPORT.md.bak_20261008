# Phase-S5-contact 波4 · 跟骨跖面真实接触 + Coulomb 摩擦（FE 场景矩阵）报告

> 生成时间 2026-10-06 17:18:22；仓库 `D:\Project\climbing_fall_analysis`；单位 mm-N-MPa-s；FEBio 4.13。
> 设计 `docs/S5_wave4_design.md`（APPROVED）；契约 `docs/S5_contact_contract.md` §4（G0–G7）/§5/§6。
> **本脚本只报绝对量，任何 σ/σ_law 比值一律不得作为「收敛」单一值引用**（`docs/S5_QUOTABLE.md` §4 #1、#2）。
> 引用限定：仅 `docs/S5_QUOTABLE.md` 的 QUOTABLE 行；跨配方对比禁止。

---

## 0. 一句话结论

1. **默认回归（硬门槛 G0）: PASS** —— `axial_regression` vs `s1_h5_opensim.json` 最大相对误差 0.00e+00；接口逐位检查 default_unchanged=True。
2. **范围**：波4 = 最小版（刚性面对偶 + Coulomb μ，无泡沫），只在 §5.1 限定的可信口径内；G1 → `N/A_minimum_version`，G7 → `RETIRED_by_contract_条6-8`（均不省略）。
3. **A7 消失（G5）未能判定**：7 条 contact 行**全部发散**（wave-3 接触面在真实 hex8 网格上被 FEBio 判为 invalid facets / negJac）。详见 §6/§9 与各 run 的 `diagnostics`。
4. **摩擦锥（G4）**：无收敛 contact 行 ⇒ 无可判定数字（FAIL）。
5. 诚实边界见 §10。

> ⚠️ **本次矩阵的真实结果**：fixed / spring 行收敛并给出应力；**contact 行全部 rc≠0**（FEBio `negative jacobians` / `invalid facets`）。这是 wave-3 接触面在真实 THUMS hex8 皮质网格上的首次真实 FE 验证结论 —— 详见 §6/§9/§10。

---

## 1. API（opt-in contact kind）

| 层 | 文件 / 符号 | 新增（默认） | 作用 |
|---|---|---|---|
| BC | `plantar_bc.apply_plantar_bc` | `kind="contact"`（默认 `"fixed"`） | 跖面 ⇄ 固定薄板（数值刚性壁）+ Coulomb μ |
| 构建 | `thums_feb.build_thums_feb` | `contact_mu=0.6` | 透传 fric_coeff（opt-in，默认不变） |
| 调度 | `scripts/opensim_fe/thums_feb_run.py` | `--plantar-bc contact` | CLI opt-in（wave-3） |

**默认 `plantar_bc="fixed"` 一个字节不变**（接口回归逐位检查 = True）。

---

## 2. 默认回归（硬门槛）

| 检查 | 结果 |
|---|---|
| 默认 `plantar_bc` 逐位不变 | **True** |
| contact opt-in 才发射 `<Contact>` | **True** |

默认链路重跑（`ground_reaction`→`run_dead_drop`→`subtalar_reaction`）vs 缓存：

| 量 | 本次默认 | 缓存 | 相对误差 |
|---|---:|---:|---:|
| `impulse_ns` | 789.328187 | 789.328187 | 0.00e+00 |
| `grf_peak_n` | 53064.8332 | 53064.8332 | 0.00e+00 |
| `subtalar_peak_vertical_n` | 26383.1338 | 26383.1338 | 0.00e+00 |
| `subtalar_peak_force_n` | 26383.134 | 26383.134 | 0.00e+00 |
| `subtalar_peak_moment_nm` | 20.0250922 | 20.0250922 | 0.00e+00 |

- 最大相对误差 **0.00e+00**（阈值 1e-06）→ **PASS**

---

## 3. 接触配方（FEBio recipe block）

- recipe id：`minimum_penalty0.1_rigid_plate`；段序 `Loads → Boundary → Contact`（契约 §5）。
- contact：`sliding-elastic`，`laugon=PENALTY`，`penalty=0.1`，`search_radius=20.0`，`node_reloc=0`，`tolerance=0.005`，`fric_coeff=μ`。
- 主面 = `plantar top facet set`；从面 = `pad_floor top facet set (fixed hex8 plate)`。
- 刚性对偶面：`rigid plate (fully fixed hex8; no rigid_wall in FEBio 4.13)`，法向 [0, -1, 0]，gap 0.5 mm，锁 ['tx', 'ty', 'tz']。
- 出处：`docs/S5_contact_udg.md` §7（NEW PENALTY 0.1）+ `docs/S5_QUOTABLE.md` §3 #1/#5。
- **FOAM 配方（G1/G7）不在范围内**（scope_limit.recipe_constraint）。

---

## 4. FE 接触输出（per-run 摘要）

| id | BC | pad | F_n | F_t | fric_max | σ_avg | σ_max | σ_p95 | gauge_max | rc | end_t | verdict |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| g0_fixed_axial_ref | fixed | on_pad | 26383 | 0 | — | 0.4397 | 8.9320 | 8.9320 | 186.76 | 0 | 1 | PASS |
| h2_fixed_pad | fixed | on_pad | 21086 | 0 | — | 0.3514 | 7.1387 | 7.1387 | 149.88 | 0 | 1 | PASS |
| h3_fixed_pad | fixed | on_pad | 22051 | 0 | — | 0.3675 | 7.4652 | 7.4652 | 156.61 | 0 | 1 | PASS |
| h3.5_fixed_pad | fixed | on_pad | 22913 | 0 | — | 0.3819 | 7.7572 | 7.7572 | 162.62 | 0 | 1 | PASS |
| h4.5_fixed_pad | fixed | on_pad | 25079 | 0 | — | 0.4180 | 8.4906 | 8.4906 | 177.67 | 0 | 1 | PASS |
| h2.0_spring_k100_pad | spring | on_pad | 21086 | 0 | — | 0.3514 | 7.1387 | 7.1387 | 122.66 | 0 | 1 | PASS |
| h2.0_contact_rigid_mu0_pad | contact | on_pad | — | — | — | — | — | — | — | 1 | — | FAIL |
| h2.0_contact_rigid_mu0.3_pad | contact | on_pad | — | — | — | — | — | — | — | 1 | — | FAIL |
| h2.0_contact_rigid_mu0.6_pad | contact | on_pad | — | — | — | — | — | — | — | 1 | — | FAIL |
| h3_contact_rigid_mu0.6_pad | contact | on_pad | — | — | — | — | — | — | — | 1 | — | FAIL |
| h3.5_contact_rigid_mu0.6_pad | contact | on_pad | — | — | — | — | — | — | — | 1 | — | FAIL |
| h4.5_contact_rigid_mu0.6_pad | contact | on_pad | — | — | — | — | — | — | — | 1 | — | FAIL |
| h2.0_fixed_hard_surface | fixed | hard_surface | 26021 | 0 | — | 0.4337 | 8.8093 | 8.8093 | 184.22 | 0 | 1 | PASS |
| h2.0_contact_rigid_mu0.0_hard_surface | contact | hard_surface | — | — | — | — | — | — | — | 1 | — | FAIL |

> 单位：F_n/F_t = N；σ = MPa；tolerance 见 §3。`—` = 未收敛/无结果（见 §10）。

---

## 5. 场景表（on_pad / hard_surface）

### on_pad

| id | h (m) | BC | μ | spring k | F_n (N) | verdict |
|---|---:|---|---:|---:|---:|---|
| g0_fixed_axial_ref | 5 | fixed | — | — | 26383 | PASS |
| h2_fixed_pad | 2 | fixed | — | — | 21086 | PASS |
| h3_fixed_pad | 3 | fixed | — | — | 22051 | PASS |
| h3.5_fixed_pad | 3.5 | fixed | — | — | 22913 | PASS |
| h4.5_fixed_pad | 4.5 | fixed | — | — | 25079 | PASS |
| h2.0_spring_k100_pad | 2 | spring | — | 100 | 21086 | PASS |
| h2.0_contact_rigid_mu0_pad | 2 | contact | 0 | — | — | FAIL |
| h2.0_contact_rigid_mu0.3_pad | 2 | contact | 0.3 | — | — | FAIL |
| h2.0_contact_rigid_mu0.6_pad | 2 | contact | 0.6 | — | — | FAIL |
| h3_contact_rigid_mu0.6_pad | 3 | contact | 0.6 | — | — | FAIL |
| h3.5_contact_rigid_mu0.6_pad | 3.5 | contact | 0.6 | — | — | FAIL |
| h4.5_contact_rigid_mu0.6_pad | 4.5 | contact | 0.6 | — | — | FAIL |

### hard_surface

| id | h (m) | BC | μ | spring k | F_n (N) | verdict |
|---|---:|---|---:|---:|---:|---|
| h2.0_fixed_hard_surface | 2 | fixed | — | — | 26021 | PASS |
| h2.0_contact_rigid_mu0.0_hard_surface | 2 | contact | 0 | — | — | FAIL |

---

## 6. 验收 G0..G6 逐条

| Gate | verdict | max/rel | evidence |
|---|---|---:|---|
| g0_default_unchanged | **PASS** | 0.000e+00 | axial_regression, baseline_ref |
| g1_pad_stress_reproduced | **N/A_minimum_version** | — | 最小版无泡沫实体 ⇒ 无 ξ σ 曲线可复现（scope_limit.recipe_constraint） |
| g2_limit_self_check | **FAIL** | — | contact μ=0 hard_surface 行未收敛（见 runs.diagnostics.neg_jacobians） |
| g3_contact_conservation | **FAIL** | — | 所有 contact 行发散；无接触守恒数字 |
| g4_friction_cone | **FAIL** | — | 所有 contact 行发散；无摩擦锥数字 |
| g5_a7_relief | **FAIL** | — | contact 行未收敛；无法对比 gauge_max（fixed 行有效） |
| g6_sensitivity_monotone_bounded | **FAIL** | — | contact 扫描行发散；μ/k 单调性不可判定 |
| g7_penetration_window | **RETIRED_by_contract_条6-8** | — | 契约 §4 条 6–8：最小版不在 G7 阻塞内；最小版不给出可信压入量 |

> G1 = `N/A_minimum_version`（最小版无泡沫 ⇒ 无 ξ σ 曲线可复现）；G7 = `RETIRED_by_contract_条6-8`（最小版不在 G7 阻塞内）。二者**均不省略**。

---

## 7. 摩擦锥 + 接触守恒

| id | force_balance_rel | friction_cone_ratio_max |
|---|---:|---:|
| g0_fixed_axial_ref | 0.000e+00 | — |
| h2_fixed_pad | 0.000e+00 | — |
| h3_fixed_pad | 0.000e+00 | — |
| h3.5_fixed_pad | 0.000e+00 | — |
| h4.5_fixed_pad | 0.000e+00 | — |
| h2.0_spring_k100_pad | 0.000e+00 | — |
| h2.0_contact_rigid_mu0_pad | — | — |
| h2.0_contact_rigid_mu0.3_pad | — | — |
| h2.0_contact_rigid_mu0.6_pad | — | — |
| h3_contact_rigid_mu0.6_pad | — | — |
| h3.5_contact_rigid_mu0.6_pad | — | — |
| h4.5_contact_rigid_mu0.6_pad | — | — |
| h2.0_fixed_hard_surface | 0.000e+00 | — |
| h2.0_contact_rigid_mu0.0_hard_surface | — | — |

> G3（守恒 ≤2%）/ G4（锥 ≤1）只在**收敛**的 contact 行上有意义；本次 contact 行全部发散，故 G3/G4 判 FAIL（见 §6）。fixed/spring 行的 force_balance_rel 是支撑反力 vs 施加载荷的残差。

---

## 8. 敏感性（μ / penalty / k）

| μ | friction_cone_holds | evidence |
|---|---|---|
| 0 | False | h2.0_contact_rigid_mu0_pad |
| 0.3 | False | h2.0_contact_rigid_mu0.3_pad |
| 0.6 | False | h2.0_contact_rigid_mu0.6_pad |

> penalty_sweep / k_sweep：penalty 固定 0.1（§5.1 #1 单配方）；k_sweep 仅 spring 参照。本次 contact 发散 ⇒ G6 判 FAIL。

---

## 9. A7 消失验证（fixed vs contact gauge 对照）

- `plantar_bc="fixed"`（h=2.0）gauge_max = —；对照 `nonvertical_s4.json::cases[axial_baseline_fixed]` = 213.50/186.76 MPa（不同载荷口径，仅形状对照）。
- `plantar_bc="contact"`（h=2.0, μ=0.6）gauge_max = 未收敛。
- **G5 = FAIL**（contact 行未收敛；无法对比 gauge_max（fixed 行有效））。
- **不掺入**任何 FOAM 版或 σ/σ_law 比值。

---

## 10. 假设（ASSUMPTIONS）与诚实边界

1. [assumed] 最小版（刚性面对偶 + μ，无泡沫）为波4 唯一配置（docs/S5接触方案.md §3.1 + docs/S5_contact_contract.md §4 条 8）。G1/G7/G7' 超范围。
2. [assumed] 配方 = NEW PENALTY 0.1（laugon=PENALTY, penalty=0.1, search_radius=20, node_reloc=0, tolerance=0.005），取自 docs/S5_contact_udg.md §7。禁止 OLD AUGLAG 与跨配方对比。
3. [measured] μ=0.6 来自 docs/S5_QUOTABLE.md §3 #1 的设定（引用，非重新标定）。
4. [measured] F_n ≈ 34.4 kN / σ_contact ≈ 0.574 MPa 是 2 m 落地（固定 sink 170 mm、NEW PENALTY 配方）网格收敛 <0.6% 的锚点；本波只作量级/口径锚，不复用为结论。
5. [assumed] 不引用任何 σ/σ_law 比值；只报绝对量 σ_contact（MPa）。
6. [assumed] 只引 docs/S5_QUOTABLE.md 的 QUOTABLE 行；禁止跨配方对比（QUOTABLE §4 #3）。
7. [assumed] 对偶面 = 全节点固定 hex8 薄板（FEBio 4.13 无 rigid_wall）；不建泡沫实体。
8. [measured] 骨 + 跖面面片来自 THUMS AM50 跟骨网格（temp/opensim_fe/thums_calcaneus/）。
9. [measured] contact 行 7 条，收敛 0 条。未收敛行的诊断记于 runs[*].diagnostics（FEBio negJac / invalid facets 计数）。
10. [measured] 所有 FEBio 调用经 climbing.coupling.febio_run.run_febio（契约 §5）；.feb 段序 Loads→Boundary→Contact；跖面对偶为真实单元面。
11. [assumed] 硬门槛 G0 是「默认路径逐位不变」的唯一数值证据。
12. [assumed] 每高度施加载荷取 s2_thums_height_sweep.json（OpenSim 距下关节峰值力）线性插值；hard_surface 载荷按 1D GRF 比值缩放（modeled）。

### 未解决的风险与下一步

- **wave-3 接触面在真实 hex8 皮质网格上非法**：`mesh['surfaces']['plantar']` 是 **292 tri3**，FEBio 判为 **292 invalid facets**（tri3 不是 hex8 的面；hex8 需 quad4）。接触因而**不承载**，骨在软罚下穿透 → 元素反演（negJac）。
- **换 quad4 也不收敛**：把跖面重建为真实 hex8 quad4 面（340/232 quads）后，invalid facets 归零，但 deck 仍在首个时间步 `negative jacobians`（负载/penalty ∈ {0.1,1,10,1000}、载荷 2 kN–21 kN、rigid/pressure 两路均如此）——疑似**皮质 hex 网格过粗 + 曲面跖面 + 平面刚性对偶**的接触力集中。
- **下一步**：(a) 把接触面建到 **SPON tet4** 域（真实 tet facet）或细分皮质壳；(b) 用 `node_reloc=1` + `auto_penalty` + 分段加载；(c) 若需可信 σ_contact 分布，改用 wave-3 的固定薄板 + 更细的接触网格。在这些完成前，波4 **不给出 contact 的 F_n/σ_contact/friction 数字**，只给出真实失败诊断。
- **structural boundary**：波4 是**最小版**；G1/G7/G7' 超范围；脚本只提供 FE-contact 能力；σ_contact 的绝对量继承 QUOTABLE §3 #1/#2 的网格收敛界（未进一步细化）。

---

## 11. 复现命令

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
# 全矩阵（含 OpenSim 硬门槛 G0）：
& .venv\Scripts\python.exe scripts\opensim_fe\nonvertical_s5_contact.py
# 只跑 FE（跳过 OpenSim 轴向回归重跑）：
& .venv\Scripts\python.exe scripts\opensim_fe\nonvertical_s5_contact.py --fast
```

## 12. 产物

- `results\opensim_fe\nonvertical_s5_contact.json`（机器可读）
- `results\opensim_fe\NONVERTICAL_S5_CONTACT_REPORT.md`（本报告）
- `temp\opensim_fe\s5_contact\*.feb/.log/.xplt`（每次求解输入/结果，可再生）

