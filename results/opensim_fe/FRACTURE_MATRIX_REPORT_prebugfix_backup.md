# 1–50 m × 9 部位骨折矩阵骨架（[Y25] 对标用）

> **本报告只做与 [Y25] 的定性比较，不声称已对齐 [Y25]。** 所有数字均为本机实测；映射与近似逐条写明；失败骨以原文记录。

## 0. 一句话结论

- 9×50 二值矩阵**成形**（9 行 × 50 列）。主用统计量：**主承载域 σ_vm 的 p95**（正则化热点；`max`/`mean`/`median` 见 §4 敏感性）。
- 在 ≤50 m 内，9 个代表骨中 **9/9** 个至少在某高度判为骨折（判据 σ_vm(h) ≥ σ_c）。
- 其中 **8/9** 在 **1 m** 即判骨折 —— 材料强度判据在「点载荷单骨 FE」上**严重过判**（见 §7 边界 3）；`median` 口径给出更有区分度的梯度（§4）。
- 载荷口径：**逐高度实测**（h=1..50 m 每档一次 `ground_reaction → run_dead_drop → joint_reaction`），**未**用 GRF 峰值比缩放。
- FE 参考载荷：6/9 骨在 @5 m 关节峰值下直接求解；3 骨因 FEBio 发散**降载**后线性折算（见 §3）。

## 1. 方法链与口径

| 环节 | 模块 | 说明 |
|---|---|---|
| 载荷端 | `climbing.coupling.joint_reactions.joint_reaction` | 关节远端刚体子树自由体；修正状态时间口径（`set_state_time=True`） |
| 通用 FE | `scripts/opensim_fe/bone_feb.py` | `extract_mesh` → `build_bone_feb` → `run_bone_feb` → `post_process`，与 `bone_feb_validate.py` 同写法 |
| 判据 | `climbing.bone.MATERIAL_STRENGTH_MPA` | [Y25] 表 1 **压缩**强度 σ_c |

- 每骨抽网格一次；在参考载荷下 build+solve+post 一次（失败则降载，见 §3）。
- 参考高度 `h_ref = 5 m`；`load_ref` = 该关节 @5 m 纵向峰值。
- 线弹性严格线性：`σ_vm(h) = σ_vm_ref · load(h)/load_ref`；判据 `fracture(h) = σ_vm(h) ≥ σ_c`。
- 求解器：`D:\Program\FEBioStudio\bin\febio4.exe`（probe=4.13.0）；`time_steps=1`。

## 2. 逐高度关节纵向峰值 load(h)（实测）

| 关节 | @5 m load_ref (kN) | @50 m (kN) | 50 m/5 m 比 | GRF 50 m/5 m 比 |
|---|---|---|---|---|
| `subtalar_r` | 25.354 | 113.700 | 4.484 | 4.680 |
| `ankle_r` | 25.266 | 112.934 | 4.470 | 4.680 |
| `hip_r` | 15.228 | 57.941 | 3.805 | 4.680 |
| `lumbar` | 22.088 | 86.203 | 3.903 | 4.680 |

> 关节载荷的 50 m/5 m 缩放比**不等于** GRF 峰值比 4.68 —— 这正是采用**逐高度实测**而非 GRF 比例外推的理由。完整表见 JSON `loads_n`。

## 3. 每骨规模 / 材料 / 参考 σ_vm / 求解载荷

| 部位 | 代表骨 | 关节 | 主承载域 | 规模 (nodes/tets+hex/shell) | E_max (MPa) | 求解载荷 (N) | 降载因子 | σ_vm_ref p95 (MPa) | σ_vm_ref max (MPa) | 线性自检比 | σ_c (MPa) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 足部 | `calcaneus_r` | `subtalar_r` | solid_calc_cort (solid) | 1544/3957/0 | 1.5e+04 | 25354 | 1.00× | 257.15 | 405.49 | 2.062 vs 2.000 | 150 |
| 胫骨 | `tibia_r` | `ankle_r` | solid_tibia_cort (solid) | 7245/23550/0 | 1.8e+04 | 5000 | 5.05× | 871.28 | 1352.90 | 2.182 vs 2.000 | 200 |
| 腓骨 | `fibula_r` | `ankle_r` | solid_fib_cort (solid) | 1837/3820/0 | 1.85e+04 | 1000 | 25.27× | 656.87 | 1061.45 | 1.842 vs 2.000 | 160 |
| 股骨 | `femur_r` | `hip_r` | solid_fem_cort (solid) | 8536/25909/0 | 1.73e+04 | 15228 | 1.00× | 1168.84 | 1417.62 | 2.617 vs 2.000 | 220 |
| 骨盆 | `R_HIPBONE` | `hip_r` | shell_hip_cort_83500201 (shell) | 3175/11767/2232 | 1.73e+04 | 15228 | 1.00× | 145.68 | 258.39 | 2.210 vs 2.000 | 180 |
| 腰椎 | `L3` | `lumbar` | shell_l3_cort (shell) | 2203/8792/1340 | 1.3e+04 | 22088 | 1.00× | 303.14 | 878.79 | 2.154 vs 2.000 | 150 |
| 胸椎 | `T6` | `lumbar` | shell_t6_cort (shell) | 933/3186/712 | 1.3e+04 | 10000 | 2.21× | 777.24 | 1419.72 | 2.552 vs 2.000 | 150 |
| 颈椎 | `C5` | `lumbar` | shell_c5_cort (shell) | 550/1660/466 | 1.3e+04 | 22088 | 1.00× | 641.14 | 1008.96 | 2.271 vs 2.000 | 150 |
| 颅骨 | `parietal_r` | `lumbar` | shell_par_ext (shell) | 2952/1430/2860 | 1.49e+04 | 22088 | 1.00× | 423.41 | 902.71 | 2.666 vs 2.000 | 160 |

> 规模列 = 全局节点 / 实体单元(tet+hex) / 壳单元。主承载域：**有壳→壳域**（多壳取该统计量最大者），**无壳→实体 CORT 域**。E_max = 该骨最硬域的材料 E。

## 4. 每部位首次骨折高度（1–50 m 或 >50）

主口径 = **σ_vm p95**（加粗）；同时给出 `max`/`mean`/`median` 敏感性。

| # | 部位 | 代表骨 | 关节 | σ_c (MPa) | σ_vm@5m p95 (MPa) | 利用率@5m | **首次骨折 (p95)** | 首次骨折 (max) | 首次骨折 (mean) | 首次骨折 (median) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 足部 | `calcaneus_r` | `subtalar_r` | 150 | 257.15 | 1.714 | **1** | 1 | 11 | 17 |
| 2 | 胫骨 | `tibia_r` | `ankle_r` | 200 | 871.28 | 4.356 | **1** | 1 | 7 | >50 |
| 3 | 腓骨 | `fibula_r` | `ankle_r` | 160 | 656.87 | 4.105 | **1** | 1 | 3 | 14 |
| 4 | 股骨 | `femur_r` | `hip_r` | 220 | 1168.84 | 5.313 | **1** | 1 | 4 | >50 |
| 5 | 骨盆 | `R_HIPBONE` | `hip_r` | 180 | 145.68 | 0.809 | **8** | 2 | 44 | >50 |
| 6 | 腰椎 | `L3` | `lumbar` | 150 | 303.14 | 2.021 | **1** | 1 | 10 | 16 |
| 7 | 胸椎 | `T6` | `lumbar` | 150 | 777.24 | 5.182 | **1** | 1 | 1 | 1 |
| 8 | 颈椎 | `C5` | `lumbar` | 150 | 641.14 | 4.274 | **1** | 1 | 2 | 7 |
| 9 | 颅骨 | `parietal_r` | `lumbar` | 160 | 423.41 | 2.646 | **1** | 1 | 3 | 4 |

> 判定用主承载域的 σ_vm 统计量；`max`（局部极值）对网格/点载荷应力集中最敏感，`median` 最稳健但丢失热点。四者差异本身即是「判据口径」不确定性的量化。

## 5. 与 [Y25] 的定性比较（仅定性）

| [Y25] 定性特征 | 本骨架观察（主口径） | 方向 |
|---|---|---|
| 足部最早骨折（约 7–9 m） | 足部(calcaneus) 首次骨折 = **1 m** | ⚠️ 足部与另外 7 个部位**并列最早**（1 m，胫骨、腓骨、股骨、腰椎、胸椎、颈椎、颅骨）——未复现 [Y25]「足部最早」的独特性 |
| 部位间先后排序（OR 序列） | 足部(1m)、胫骨(1m)、腓骨(1m)、股骨(1m)、腰椎(1m)、胸椎(1m)、颈椎(1m)、颅骨(1m)、骨盆(8m) | 见上；仅方向性参考 |
| 阶跃/渐进（临界高度） | 判据 σ_vm(h) ≥ σ_c 天然给出阶跃 | 形式一致，数值不可对标 |
| 双路径（接触 + 轴向） | 只用关节反力**单路径**传入单骨 FE | 只覆盖一条路径 |

**不声称对齐 [Y25]**：本骨架是单骨截断 + 刚性腿伪影 + 材料强度判据，可用于比较的是**排序与单调性**，不是临界高度。

## 6. 映射假设与近似

| # | 部位 | 代表骨 | 关节 | σ_c 键 | 映射理由 |
|---|---|---|---|---|---|
| 1 | 足部 | `calcaneus_r` | `subtalar_r` | `calcaneus` | 跟骨是 [Y25] 足部加载的代表骨；无壳 → 取实体 CORT 域。σ_c=150（calcaneus）。 |
| 2 | 胫骨 | `tibia_r` | `ankle_r` | `tibia_shaft` | 合并 3 PID 成一根骨（PARTS_MAPPING §1）；THUMS 胫 CORT 为实体 hex。无壳 → 取实体 CORT 域。σ_c=200（tibia_shaft，按任务默认；胫骨两端 70 见诚实边界）。 |
| 3 | 腓骨 | `fibula_r` | `ankle_r` | `fibula_shaft` | 合并 3 PID；无壳 → 实体 CORT 域。σ_c=160（fibula_shaft）。 |
| 4 | 股骨 | `femur_r` | `hip_r` | `femur_shaft` | 合并 3 PID；无壳 → 实体 CORT 域。σ_c=220（femur_shaft；股骨颈 80 见诚实边界）。 |
| 5 | 骨盆 | `R_HIPBONE` | `hip_r` | `pelvis` | 取髋骨代表骨（骶/耻骨联合未并，PARTS_MAPPING §3）；有壳 → 取最受载壳域（每统计量取跨 9 壳的最大）。σ_c=180（pelvis）。 |
| 6 | 腰椎 | `L3` | `lumbar` | `spine` | 取单代表椎体 L3（左右 SPON 两半合并，PARTS_MAPPING §3）；有壳 → 壳域。σ_c=150（spine）。 |
| 7 | 胸椎 | `T6` | `lumbar` | `spine` | 取单代表椎体 T6；有壳 → 壳域。σ_c=150（spine）。 |
| 8 | 颈椎 | `C5` | `lumbar` | `spine` | 取单代表椎体 C5；有壳 → 壳域。σ_c=150（spine）。 |
| 9 | 颅骨 | `parietal_r` | `lumbar` | `skull` | 取单块代表颅骨 parietal_r（整颅是薄壳穹顶，不可当实心骨，PARTS_MAPPING §4）；有壳 → 壳域。σ_c=160（skull）。载荷沿用 lumbar（默认映射，非物理结论，见诚实边界）。 |

- **参考载荷**：每部位取其关节 @5 m 纵向峰值（实测）。
- **主承载域**：有壳→壳域（多壳取统计量最大），无壳→实体 CORT 域。
- **FE 边界**：沿用 `bone_feb_validate.py` 的 bbox 最长轴轴向压缩 + 顶/底 30% 带节点力 + 底带三向固定；`time_steps=1`（模块默认）——实测步数 >1 会让 FEBio 对本 NodalForce 设置误加载而发散，故不采用。
- **降载**：发散骨取候选载荷中最大可解者，σ_vm 按线性折回参考载荷。
- **线性外推**：`σ_vm(h)=σ_vm_ref·load(h)/load_ref`（线弹性）。

## 7. 诚实边界

1. **单骨截断 vs 全身**：每个 FE 是一根代表骨在两条节点带间受压，无关节/韧带/肌肉/相邻骨，不体现全身载荷重分配。
2. **刚性腿伪影**：dead-drop 锁定下肢关节 + 禁用肌肉，实测近端 |F| **偏低**（髋 15.2 kN < 距下 25.4 kN @5 m），与活体趋势可能相反。
3. **判据用材料强度×正则化**：[Y25] 的 σ_c 是**材料强度**，本骨架直接与单骨 FE 的 σ_vm 比较 → 点载荷/薄截面处存在**人为应力集中**，`max` 对网格敏感（本骨架主口径改取 **p95** 以正则化；`max`/`median` 一并给出）。整体骨失效还涉及弯曲、小梁骨、应力集中，材料强度会**高估**承载能力。
4. **长骨 σ_c 取骨干值**：任务默认胫/腓/股用 200/160/220（shaft）。两端（tibia/fibula_ends=70、femoral_neck=80）约弱 3 倍，若改取骨端值首次骨折高度会**显著提前**——最大的一处口径分歧。
5. **载荷方向**：关节反力是全局纵向，FE 却沿骨 bbox 最长轴施加；对非竖直骨是方向近似。
6. **载荷端近似**：GRF 来自 `pad.py` 1D 等效双质点模型、施加于 calcn 原点；关节反力为自由体结果；仅 GRF 窗口内有效。
7. **颅骨**：[Y25] 明确颅骨失效由 HIC/接触主导、不在轴向通路；本骨架按任务默认映射到 `lumbar`，该行**仅为骨架占位**，无物理意义。
8. **代表骨 ≠ 整个部位**：足部=跟骨（非 26 块足骨）；脊柱各取单椎体；骨盆只取髋骨（骶/耻骨联合未并）；颅骨取 parietal_r 一块。
9. **线性损伤**：不建模塑性/损伤/失效模式转换/应变率；材料卡中的塑性/损伤部分被 `bone_feb.read_materials` 取**弹性近似**（模块既有约定）。
10. **FE 降载与线性假设的实测偏差**：3 骨在 @5 m 载荷下发散（近不可压松质骨/薄壳点载荷），在较低载荷下求解后线性折回；双载荷自检实测比值 σ(2L)/σ(L) 为 **1.84–2.67**（严格线性应为 2.0），即 FE 响应在本载荷区间**并非严格线性**，`σ_vm_ref` 与 `σ_vm(h)` 的线性外推可带 **约 ±8–33%** 的量级误差（见 §3「线性自检比」列）。这是本骨架最需注意的定量边界之一。
11. **不声称对齐 [Y25]**：仅定性比较排序/单调性/阶跃形式。

## 8. 失败骨原文记录

本次运行 **9/9 骨全部成功**（`run_bone_feb` rc=0，`post_process` 取到主域）。无失败需记录。

> 以下为**降载搜索中被 FEBio 拒绝的尝试**（原文），供复现与审计：
> - `tibia_r` load=25266 N: BoneFebError: FEBio 失败 FEBio 运行失败（非 0 退出码）。 命令: D:\Program\FEBioStudio\bin\febio4.exe -i D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\tibia_r\tibia_r.feb 工作目录: D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\tibia_r 退出码: 1  ---- 输出（尾部）----  	Total calls to linear solver ........ : 72  	Avg iterations per solve ............ : 1  	Time in linear solver: 0:00:07
> - `tibia_r` load=10000 N: BoneFebError: FEBio 失败 FEBio 运行失败（非 0 退出码）。 命令: D:\Program\FEBioStudio\bin\febio4.exe -i D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\tibia_r\tibia_r.feb 工作目录: D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\tibia_r 退出码: 1  ---- 输出（尾部）----  	Total calls to linear solver ........ : 204  	Avg iterations per solve ............ : 1  	Time in linear solver: 0:00:1
> - `fibula_r` load=25266 N: BoneFebError: FEBio 失败 FEBio 运行失败（非 0 退出码）。 命令: D:\Program\FEBioStudio\bin\febio4.exe -i D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\fibula_r\fibula_r.feb 工作目录: D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\fibula_r 退出码: 1  ---- 输出（尾部）----  	Total calls to linear solver ........ : 36  	Avg iterations per solve ............ : 1  	Time in linear solver: 0:00
> - `fibula_r` load=10000 N: BoneFebError: FEBio 失败 FEBio 运行失败（非 0 退出码）。 命令: D:\Program\FEBioStudio\bin\febio4.exe -i D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\fibula_r\fibula_r.feb 工作目录: D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\fibula_r 退出码: 1  ---- 输出（尾部）----  	Total calls to linear solver ........ : 72  	Avg iterations per solve ............ : 1  	Time in linear solver: 0:00
> - `fibula_r` load=5000 N: BoneFebError: FEBio 失败 FEBio 运行失败（非 0 退出码）。 命令: D:\Program\FEBioStudio\bin\febio4.exe -i D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\fibula_r\fibula_r.feb 工作目录: D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\fibula_r 退出码: 1  ---- 输出（尾部）----  	Total calls to linear solver ........ : 72  	Avg iterations per solve ............ : 1  	Time in linear solver: 0:00
> - `fibula_r` load=2500 N: BoneFebError: FEBio 失败 FEBio 运行失败（非 0 退出码）。 命令: D:\Program\FEBioStudio\bin\febio4.exe -i D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\fibula_r\fibula_r.feb 工作目录: D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\fibula_r 退出码: 1  ---- 输出（尾部）----  	Total calls to linear solver ........ : 162  	Avg iterations per solve ............ : 1  	Time in linear solver: 0:0
> - `fibula_r` load=2000 N: BoneFebError: FEBio 失败 FEBio 运行失败（非 0 退出码）。 命令: D:\Program\FEBioStudio\bin\febio4.exe -i D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\fibula_r\fibula_r.lincheck.feb 工作目录: D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\fibula_r 退出码: 1  ---- 输出（尾部）----  	Total calls to linear solver ........ : 84  	Avg iterations per solve ............ : 1  	Time in linear sol
> - `T6` load=22088 N: BoneFebError: FEBio 失败 FEBio 运行失败（非 0 退出码）。 命令: D:\Program\FEBioStudio\bin\febio4.exe -i D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\T6\T6.feb 工作目录: D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\T6 退出码: 1  ---- 输出（尾部）----  	Total calls to linear solver ........ : 126  	Avg iterations per solve ............ : 1  	Time in linear solver: 0:00:01   Peak memory
> - `T6` load=20000 N: BoneFebError: FEBio 失败 FEBio 运行失败（非 0 退出码）。 命令: D:\Program\FEBioStudio\bin\febio4.exe -i D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\T6\T6.lincheck.feb 工作目录: D:\Project\climbing_fall_analysis\temp\opensim_fe\fracture_matrix\T6 退出码: 1  ---- 输出（尾部）----  	Total calls to linear solver ........ : 138  	Avg iterations per solve ............ : 1  	Time in linear solver: 0:00:02   Pe

## 9. 复现命令

```powershell
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\fracture_matrix.py
# 只重跑载荷链 / 只重跑 FE：
& .venv\Scripts\python.exe scripts\opensim_fe\fracture_matrix.py --loads-only
& .venv\Scripts\python.exe scripts\opensim_fe\fracture_matrix.py --fe-only
# 换主统计量 / 忽略缓存强制重算：--stat max|p95|mean|median --force-loads --force-fe
```

回归（必须保持 19 passed / 1 skipped）：

```powershell
& .venv\Scripts\python.exe -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q
```

## 10. 产物

- `results/opensim_fe/fracture_matrix.json`（9×50 矩阵 + 逐高度载荷、σ_vm、四统计量）
- `results/opensim_fe/fracture_matrix.csv`（宽表：行=部位，列=h1..h50，主口径）
- `results/opensim_fe/fracture_matrix.png`（左：二值骨折矩阵；右：利用率+等值线）
- `results/opensim_fe/FRACTURE_MATRIX_REPORT.md`（本文件）
- `temp/opensim_fe/fracture_matrix/`（每骨 .feb/.xplt/.hdf5/.log + per-bone JSON；`load_chain.npz`；`run.log`）

_生成时间：2026-10-04 21:23:03；总耗时 0 s。_
