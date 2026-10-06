# 距下关节载荷表示 × 分析类型对皮质 gauge 的影响（THUMS AM50 跟骨）

**问题**：现有方向 C anat 单骨跟骨 FE 用跖面三向全固定 + 距下**均匀面压**
`PressureLoad p=F/A` + **STATIC**。跖面实验已证明跖面**法向刚度**是主因之一；
本实验查**另一半**：距下**均匀压力**（vs 偏心梯度 / 逐节点力）与 **STATIC**
（vs 短时瞬态 DYNAMIC）各自把 gauge 抬高/压低多少。

- 网格：`calcaneus_r_anatframe.npz`（AM50 THUMS 右跟骨，双域 CORT hex8 + SPON tet4）
- 载荷：h=5 m, a=0（无跟腱）→ F_subt = 26.383 kN（取自 anat 扫描真跑；T1 冲量误差 +5.8%）
- 判据：皮质域 `calcaneus` 过程区正则化 gauge_max（体积平均 R=4 mm）；屈服参考 [Y25] = 150 MPa（**后处理判据，无单元删除**）
- 默认行为未变：`joint_load='pressure'` + `analysis='STATIC'` 与旧路径**逐字节一致**

## 0. 基线复现（硬门槛）

- anat 扫描 h=5, a=0 基线：gauge_max = **186.756 MPa**（raw max = 213.50 MPa）
- 本实验复现：gauge_max = **186.756 MPa**（rc=0，2 state）
- 相对误差：+0.00%（门槛 ±1%）→ **通过**

## 1–2. 各工况 gauge

| tag | joint_load | analysis | rc | gauge_max (MPa) | vs 基线 | gauge_p95 | raw max | p95 | 峰值位置 |
|---|---|---|---|---|---|---|---|---|---|
| baseline_pressure_static | pressure | STATIC | 0 | 186.756 | +0.00% | 136.928 | 213.503 | 136.297 | 跖面边缘 (plantar rim) |
| gradient_c030_static | gradient | STATIC | 0 | 191.849 | +2.73% | 141.898 | 233.451 | 141.331 | 跖面边缘 (plantar rim) |
| gradient_c040_static | gradient | STATIC | 0 | 194.679 | +4.24% | 143.079 | 240.186 | 143.651 | 跖面边缘 (plantar rim) |
| gradient_c060_static | gradient | STATIC | 0 | 205.231 | +9.89% | 150.144 | 253.782 | 150.144 | 跖面边缘 (plantar rim) |
| gradient_c040_dirZ_static | gradient | STATIC | 0 | 167.537 | -10.29% | 119.345 | 167.537 | 120.101 | 跖面边缘 (plantar rim) |
| nodal_static | nodal | STATIC | 0 | 199.519 | +6.83% | 158.087 | 221.145 | 157.898 | 跖面边缘 (plantar rim) |
| pressure_dynamic | pressure | DYNAMIC | 0 | 188.789 | +1.09% | 138.587 | 217.851 | 138.222 | 跖面边缘 (plantar rim) |
| gradient_c040_dynamic | gradient | DYNAMIC | 0 | 198.591 | +6.34% | 145.609 | 244.866 | 145.578 | 跖面边缘 (plantar rim) |

> 静态工况 headline = 末态（满载荷）；动态工况 headline = **全状态峰值**（含峰值时刻）。

### 峰值位置明细

| tag | 质心 (mm) | 最近具名面 | 距离 (mm) | d(joint) | d(plantar) | 分类 |
|---|---|---|---|---|---|---|
| baseline_pressure_static | (33.4, -0.9, -16.7) | plantar | 2.99 | 13.5 | 3.0 | 跖面边缘 (plantar rim) |
| gradient_c030_static | (33.4, -0.9, -16.7) | plantar | 2.99 | 13.5 | 3.0 | 跖面边缘 (plantar rim) |
| gradient_c040_static | (33.4, -0.9, -16.7) | plantar | 2.99 | 13.5 | 3.0 | 跖面边缘 (plantar rim) |
| gradient_c060_static | (33.4, -0.9, -16.7) | plantar | 2.99 | 13.5 | 3.0 | 跖面边缘 (plantar rim) |
| gradient_c040_dirZ_static | (27.3, -16.8, 20.5) | plantar | 3.53 | 20.2 | 3.5 | 跖面边缘 (plantar rim) |
| nodal_static | (2.9, -25.0, -14.8) | plantar | 2.54 | 5.4 | 2.5 | 跖面边缘 (plantar rim) |
| pressure_dynamic | (33.4, -0.9, -16.7) | plantar | 2.99 | 13.5 | 3.0 | 跖面边缘 (plantar rim) |
| gradient_c040_dynamic | (33.4, -0.9, -16.7) | plantar | 2.99 | 13.5 | 3.0 | 跖面边缘 (plantar rim) |

## 3. DYNAMIC 峰值时刻 + T8 不变量

T8（方案 §七）：瞬态峰值 ≤ 2× static 峰值。

| 动态工况 | 峰值 t (s) | peak gauge | final gauge | 对应 static gauge | 比值 peak/static | T8 (≤2×) |
|---|---|---|---|---|---|---|
| pressure_dynamic | 0.0050 | 188.789 | 188.789 | 186.756 | 1.0109 | **通过** |
| gradient_c040_dynamic | 0.0050 | 198.591 | 198.591 | 194.679 | 1.0201 | **通过** |

## 4. 结论

- **距下压力分布（偏心/梯度，vs 均匀）**：`gradient_c030_static` +2.73%、`gradient_c040_static` +4.24%、`gradient_c060_static` +9.89%、`gradient_c040_dirZ_static` -10.29%；最大绝对偏离 **10.29%** → **明显 (>10%)**。（方向敏感：沿 +X 偏心提高 gauge，沿 +Z 偏心降低。）
- **逐节点力（vs 均匀面压）**：gauge = 199.52 MPa（+6.83%）→ **不显著 (≤10%)**。
- **DYNAMIC vs STATIC（均匀面压）**：峰值 gauge = 188.79 MPa @ t=5.00 ms（+1.09%）→ **不显著 (≤10%)**；末态 188.79 MPa。
- **DYNAMIC + 梯度（coeff 0.40, dir X）**：峰值 gauge = 198.59 MPa @ t=5.00 ms（+6.34%）；T8 peak/static=1.020。

### 「均匀压力 / STATIC 各自把 gauge 抬高/压低多少？是否明显 (>10%)？」

- **均匀压力（vs 偏心梯度 / 逐节点力）**：最大绝对偏离 **10.29%** → **明显 (>10%)**。逐项：沿 +X 偏心 coeff 0.30→0.60 单调抬高 +2.7%→+9.9%，逐节点力 +6.8%，**沿 +Z 偏心 coeff 0.40 反而 -10.3%**（唯一越过 10% 门槛者）。即：一般载荷重分布把 gauge 改变 ≤10%（不显著），只有**强偏心且方向不利**时才>10%。
- **STATIC（vs 短时 DYNAMIC）**：瞬态峰值相对 static 偏离 **1.09%** → **不显著 (≤10%)**。即：本 **2 ms ramp + 隐式动力学**下，惯性放大可忽略，STATIC 是充分的（这是**负结论**）。

> **诚实边界**：所有结论限于**本几何 / 本网格 / 本线弹性材料 / 本固定跖面**下的
> 筛选级对比；动态用隐式积分 + 自适应时间步（`<solver/>` 默认），
> 峰值可能仍受离散 / 时间步分辨率影响；绝对应力不是临床预测。
> 若某因素影响很小，同样是有效负结论。
