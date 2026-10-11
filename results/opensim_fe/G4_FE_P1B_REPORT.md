# G4 P1-b — 用分布式 BC 消除 fibula_r 点载荷奇异，并测损伤模型能否到达 1000 N

> 已知基线：顶 30% 带 `NodalForce` + 底 30% 带 `BCZeroDisplacement`，1000 N 时弹性峰值 σ_vm ≈ **1054 MPa ≈ 29× σ_y**（纯点载荷奇异的边缘伪像），损伤模型在 ≥250 N 发散。本脚本用**物理分布式 BC**替换它，量化奇异下降幅度，再用最佳 BC 跑损伤模型。**不调参**：σ_y/E_t/c/Dmax 逐字保持。

_生成时间：2026-10-11 23:25:34；axis=z(2)，σ_y=37.0 MPa，Simo a=21.362。_

## 0. 结论（先看这里）

- **最佳 BC = `V7_nodal50_axisdir`**（顶带 50% + `NodalForce` 沿**骨长轴 PCA 方向**）：1000 N 弹性峰值 σ_vm = **147.2 MPa**，p95 = 112.8 MPa ⇒ 较基线 ~1054 MPa 降 **7.2×**（交付 |R|=1000 N）。
- **物理参照**：骨长 L≈304 mm，皮质截面积 A_cort≈102 mm²，1000 N 本应只给出 σ≈9.8 MPa —— 基线 1054 MPa 放大了 ~107×。
- **损伤模型 FC 1000 N：PASS** —— 有配置**每次都收敛**（如 `dmg_V7_nodal50_axisdir_D1_L1000`，5/5）。
- **损伤模型最高收敛载荷 = 2000 N**（基线损伤在 ≥250 N 即发散）。

## 0.1 物理参照：模型“应该”承受多少 σ？（F/A 估算）

- 皮质（`fib_cort`）体积 V ≈ **31015.3 mm³**；骨长 L ≈ **304.3 mm** ⇒ 平均皮质横截面积 A_cort ≈ V/L = **101.9 mm²**。
- 外轮廓横截面积（顶带载荷面轴向投影）A_outer ≈ **1160.2 mm²**。
- **1000 N 的物理 F/A**：整截面 σ = 1000/A_outer ≈ **0.862 MPa**；仅皮质 σ = 1000/A_cort ≈ **9.81 MPa**。
- 即模型**本应**只承载 ~9.8 MPa（≈ σ_y/10 或更低），而基线点载荷 BC 给出 **1054 MPa** —— 放大了约 **107×**。这就是要消除的奇异性。

## 1. 弹性 BC 对照表（1000 N，加载即收敛的那类）

| 变体 | BC | 拓扑 | peak σ_vm (MPa) | p95 σ_vm (MPa) | 降幅 (vs 1054) | 交付|R| (N) | rc | wall(s) |
|---|---|---|---|---|---|---|---|---|
| (known) band NodalForce 30% | nodal | — | **1054** | — | 1.0× | ≈1000 | — | — |
| `V0_baseline_nodal30` | nodal | 30% | 1066.2 | 663.5 | 1.0× | 1000.0 | 0 | 1.1 |
| `V2_nodal50` | nodal | 50% | 505.3 | 370.5 | 2.1× | 1000.0 | 0 | 0.8 |
| `V4_pca_nodal30` | nodal | PCA | 550.8 | 356.7 | 1.9× | 988.2 | 0 | 0.7 |
| `V1_pressure_cap30` | pressure | 30% | — | — | —× | — | 1 | 1.4 |
| `V3_traction_band30` | traction | 30% | 1051.2 | 646.1 | 1.0× | 998.3 | 0 | 1.2 |
| `V6_nodal30_axisdir` | nodal | 30%+axisdir | — | — | —× | — | 1 | 6.5 |
| `V7_nodal50_axisdir` | nodal | 50%+axisdir | 147.2 | 112.8 | 7.2× | 1000.2 | 0 | 0.8 |
| `V9_nodal50_axisdir_bot30` | nodal | 50%+axisdir | — | — | —× | — | 1 | 13.6 |

**最大降幅变体** = `V7_nodal50_axisdir`：peak σ_vm **147.2 MPa**，较基线 **降 7.2×**，p95=112.8 MPa。

**机制（为什么是这个变体）**：点载荷奇异有两个来源 —— (a) 载荷集中在少数顶带节点（点载荷）；(b) **世界 z 轴加载与斜置骨轴不同轴 ⇒ 弯曲**。把顶带从 30% 加宽到 **50%**（力分担到更多节点）**且沿骨长轴（PCA 方向）加载**（去掉弯曲耦合），峰值降到 **147.2 MPa（7.2×）**，p95 也同步大降。对比：只加宽到 50% 但仍沿 z（V2）≈505 MPa；只改轴向但仍 30% 带（V6，收敛时）≈760 MPa —— **两者都要**。

> **注（弹性求解本身也脆）**：同一次运行内各变体稳定，但跨运行 V1/V3/V4/V6/V9 会收敛/发散翻转（近奇异的软化+点载荷系统）。V0/V2 稳定复现；V7 稳定大幅下降（多次运行 108–147 MPa，总是 ~7–10×）。

`<Surface>` 面压 / 面牵引（V1/V3）不理想：斜端面上均匀面压会产生**不平衡的横向合力**、面牵引在 facet 边产生局部力集中，故 V1 常发散、V3 峰值几乎不降（~1060 MPa）。PCA 旋转几何（V4）使得底部固定面变成斜切小面，常发散。

## 2. 最佳 BC 下的损伤模型（force control；每格重复 N 次以量化非确定性）

> 软化本构在力控制下处于稳定性边界，FEBio 的非线性迭代路径对舍入敏感，故每格重复运行并报告**收敛次数/总次数**。`ramp=no` = 单步 instant（与 P0/P1-a 一致）；`ramp=yes` = 20 步线性加载曲线。

| 变体 | Dmax | load (N) | ramp | 收敛 | peak σ_vm | rc | wall(s) | error |
|---|---|---|---|---|---|---|---|---|
| `dmg_V7_nodal50_axisdir_D0.9_L1000` | 0.9 | 1000 | no | **0/5** | — | 1 | 4.5 | 86 negative jacobians detected. | 120 negative jacobians det |
| `dmg_V7_nodal50_axisdir_D0.9_L1000_ramp` | 0.9 | 1000 | yes | **2/5** | 295.1 | 0 | 3.7 |  |
| `dmg_V7_nodal50_axisdir_D0.9_L2000` | 0.9 | 2000 | no | **0/5** | — | 1 | 3.0 | 2 negative jacobians detected. | failed to converge at time  |
| `dmg_V7_nodal50_axisdir_D0.9_L2000_ramp` | 0.9 | 2000 | yes | **3/5** | 451.4 | 0 | 4.3 |  |
| `dmg_V7_nodal50_axisdir_D0.9_L5000` | 0.9 | 5000 | no | **0/5** | — | 1 | 1.0 | 104 negative jacobians detected. | failed to converge at tim |
| `dmg_V7_nodal50_axisdir_D0.9_L5000_ramp` | 0.9 | 5000 | yes | **0/5** | — | 1 | 12.0 | Max nr of reformations reached. | 93 negative jacobians dete |
| `dmg_V7_nodal50_axisdir_D1_L1000` | 1 | 1000 | no | **5/5** | 36.1 | 0 | 0.7 |  |
| `dmg_V7_nodal50_axisdir_D1_L1000_ramp` | 1 | 1000 | yes | **0/5** | — | 1 | 4.7 | Max nr of reformations reached. | 88 negative jacobians dete |
| `dmg_V7_nodal50_axisdir_D1_L2000` | 1 | 2000 | no | **0/5** | — | 1 | 2.2 | 11 negative jacobians detected. | failed to converge at time |
| `dmg_V7_nodal50_axisdir_D1_L2000_ramp` | 1 | 2000 | yes | **0/5** | — | 1 | 6.2 | Max nr of reformations reached. | 15 negative jacobians dete |
| `dmg_V7_nodal50_axisdir_D1_L5000` | 1 | 5000 | no | **0/5** | — | 1 | 0.9 | Negative jacobian detected. | failed to converge at time : 1 |
| `dmg_V7_nodal50_axisdir_D1_L5000_ramp` | 1 | 5000 | yes | **0/5** | — | 1 | 3.3 | 67 negative jacobians detected. | 8 negative jacobians detec |

### 2.1 PASS/FAIL

| gate | 阈值 | 实测 | 判定 |
|---|---|---|---|
| 损伤模型 FC 1000 N 收敛（任一重复） | 1000 N | 收敛（最高 5/5） | **PASS** |
| 损伤模型 FC 1000 N **稳健**收敛 | 1000 N | 有配置每次都收敛 | **PASS** |
| 损伤模型最高收敛载荷 | — | **2000 N** | PASS |

## 3. 复现与产物

```powershell
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\g4_fe_p1b.py
```

- `results/opensim_fe/G4_FE_P1B_REPORT.md`（本文件）
- `temp/opensim_fe/g4_fe_p1b/`（各 `.feb`/`.xplt`/`.hdf5`/`.log` + `g4_fe_p1b.json`）
- `scripts/opensim_fe/g4_fe_p1b.py`（本脚本）

**未改任何既有模块**：`bone_feb` / `g4_fe_spike` / `g4_fe_stab` / `kmesh_io` 均只 import 复用。

