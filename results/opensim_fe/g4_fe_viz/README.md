# g4_fe_viz — fibula_r 损伤模型失败可视化

_生成时间：2026-10-11 22:05:28。脚本：`scripts/opensim_fe/g4_fe_viz.py`（只新增，不改任何既有模块）。_

## 主导出（点载荷奇异 / 单元反转）

| 项 | 值 |
|---|---|
| 运行配置 | **位移控制** `disp_d2.0x1000` = P1-a 的 `disp@-397.858`, `Dmax=1.0`（ramped prescribed z-displacement，∂ 输入位移 −397.858 mm） |
| `.feb` | `D:\Project\climbing_fall_analysis\temp\opensim_fe\g4_fe_stab\disp_d2.0x1000.feb` |
| `.xplt` | `D:\Project\climbing_fall_analysis\temp\opensim_fe\g4_fe_stab\disp_d2.0x1000.xplt` |
| **导出 state 索引** | **92**（共 93 个写入 state；FEBio 在**收敛增量**才写 state，故 state 92 = **发散前最后一个收敛态**） |
| 节点数 / 单元块 | 1837 / fib_spon_end(tet4,762)、fib_spon_center(tet4,2216)、fib_cort(hex8,842) |
| **max\|u\|** | **217.846 mm** |
| **负 Jacobian 单元数** | **0** / 3820 |
| 最大损伤 D | 0.9968 |

### 单元 Jacobian（末收敛态；位置按**初始**沿 z 归一化，0=底 1=顶）

| 域 | 单元数 | 负 J 数 | 顶带(z≥0.70) | 底带(z≤0.30) | 中部 | min detJ/\|detJ0\| | 挤到<10% |
|---|---|---|---|---|---|---|---|
| `fib_spon_end` | 762 | 0 | 0 | 0 | 0 | 1.0000 | 0 |
| `fib_spon_center` | 2216 | 0 | 0 | 0 | 0 | 0.0054 | 2 |
| `fib_cort` | 842 | 0 | 0 | 0 | 0 | 0.5293 | 0 |

## 副导出（力控制 1000 N——经典 `98 negative jacobians`）

- **该运行没有写出任何收敛 state**（`re-run also wrote no converged state`）。力控制 1000 N 在**第一个增量**就发散（`98 negative jacobians detected. | failed to converge at time : 1`），故 `.xplt` 只有 state 0（未变形、位移=0），无法用于「变形网格」可视化——这正是需要位移控制来**逐增量**逼近失败的原因。

## 产物

- `initial_mesh.vtk` — 未变形网格（hex8/tet4 连通性保留；point_data 含 `top_load`/`bottom_fix` 掩码）。
- `deformed_first_failure.vtk` — 末收敛态网格；point_data `displacement`(向量)、`abs_u`(标量)；cell_data `damage`、`von_mises_mpa`。
- `overlay_initial_vs_deformed.png` — 初始（灰线框）vs 变形（按 |u| 着色），X–Z 与 Y–Z 两视图。
- `damage_field.png` — 变形网格按单元损伤 D 着色。

## 观察

FEBio 日志原文（节选）：

```
Number of time steps completed .................... : 92
* 32 negative jacobians detected.                                       *
* Negative jacobian detected.                                           *
* Negative jacobian detected.                                           *
* Negative jacobian detected.                                           *
* Negative jacobian detected.                                           *
```

- 导出的 state **92** 是**最后一个收敛态**（FEBio 只在收敛增量写 state）——顶带被下压到 `max|u|=217.8 mm`。
- 该态**负 Jacobian 单元数 = 0**。**单元反转发生在下一个（发散的、未写入的）增量**：FEBio 在试探更新处检测到反转即 `ERROR TERMINATION`，该构型**不进 `.xplt`** —— 所以任何「已写入 state」按定义都是未反转的。这是「文件里的最后状态」与「崩溃瞬间构型」的差别。
- 压溃热点（最接近反转）：`min detJ/|detJ0| = 0.0054`，域 `fib_spon_center`，初始 z 归一化位置 `0.006` （**bottom_fix 带（底固定带）**）→ 反转由该带起始。
- 损伤 D 达 **0.997**（`fib_cort` 近饱和），集中在压溃区。
- VTK 为 legacy ASCII（`meshio` 写出并**回读校验**：hex8=8 节点、tet4=4 节点）。


## 忠实首失败（viz2：最小过冲）+ 损伤起始

_生成时间：2026-10-11 22:12:10。脚本：`scripts/opensim_fe/g4_fe_viz2.py`。_

- **δ\*** = **211.362 mm**（在弹性 δ=198.929 mm 与已知发散 δ=397.9 mm 之间二分得到的最小发散 preset δ；**过冲 = 1.062×**，即仅比弹性-1000 N 位移大 6.2%）。
- **忠实首失败态**：preset δ*=211.36 mm 运行的**末收敛 state 63**（共 64；t=0.9949）→ `deformed_first_failure_faithful.vtk`。
  - `max|u| = 217.69 mm`；**负 Jacobian 单元数 = 0**（min detJ/|detJ0| = 0.00252，即压到 0.25% 体积、即将反转）。
- **损伤起始态**：state 1（`fib_cort` 首个 D>0；`max|u| = 3.14 mm`，`Dmax = 0.579`，727 个单元）→ `damage_onset.vtk`。

### 新增图（各 1–2 行说明）

- `overlay_first_failure_faithful.png`：δ* 运行**末收敛态**——初始灰线框 vs 变形网格（按 |u| 着色）；折叠幅度最小、且含最接近反转的压溃单元（min detJ→0）。
- `damage_onset.png`：**损伤起始**态（`fib_cort` 首次 D>0）按单元 D 填充着色；此时损伤仅出现在点载荷热点，网格几乎未折叠。

### δ* 二分记录

| preset δ (mm) | rc | last t | 结果 | wall (s) |
|---|---|---|---|---|
| 298.393 | 1 | — | 发散 | 13.8 |
| 248.661 | 1 | — | 发散 | 13.3 |
| 223.795 | 1 | — | 发散 | 17.1 |
| 211.362 | 1 | — | 发散 | 9.6 |
| 205.145 | 0 | 1.0000 | 收敛 | 9.6 |
| 208.254 | 0 | 1.0000 | 收敛 | 13.3 |

