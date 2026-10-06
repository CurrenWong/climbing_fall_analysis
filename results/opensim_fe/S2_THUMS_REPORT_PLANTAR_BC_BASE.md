# S2-THUMS · AM50 跟骨刚性地面全高度扫描（Option C） · plantar_bc_base

**模型**：`temp\opensim_fe\thums_calcaneus\calcaneus_r_anatframe.npz`（AM50 THUMS 右跟骨，mm，OpenSim `calcn_r` 帧；原点 = 距下关节面质心）。
双材料域：皮质 `calcaneus`（hex8，E=15000 MPa，ν=0.30）+ 松质 `trabecular`（tet4，E=73.4 MPa，ν=0.45）。

- 单元：CORT hex8 = 634，SPON tet4 = 3323；节点 1544。
- 具名面：subtalar_joint 119 面 / plantar 137 面 / achilles 105 面。
- **距下关节载荷面**：dist<25 mm 且 n·(+Y)>0.2 → 2495.3 mm²（旧 foot-frame 配准 = 222.0 mm²，14 面）。
- 求解：FEBio 4.13 STATIC × 1 次，**全部退出码 0**；距下关节载荷用 `PressureLoad p=F/A`（`use_rigid=False`，与旧管线同口径），跟腱用 `TractionLoad`。
- 判据：皮质域 `gauge_max`（体平均正则化，半径 R=4 mm = ℓ）> **150 MPa**（[Y25] 皮质压缩强度）。**后处理判据，无单元删除。**

## 载荷链（复用 OpenSim 管线原样）

```
h -> ground_reaction(h) -> run_dead_drop(grf) -> subtalar_reaction(fall, grf, side='r')
  -> peak_vertical_n = F_subt(h)      # 真跑正动力学
  -> load_transfer.transfer(F_subt, achilles_activation=a)
       subtalar_dir=(0,-1,0), achilles_dir=(0.15,0.985,0), F_MAX=10885 N
```

## 首次骨折高度（gauge_max > 150 MPa）

| 激活 a | 跟腱力 (kN) | 首次骨折高度 | 触发点 gauge_max (MPa) |
|---|---|---|---|
| 0 | 0.00 | 5 m | 186.8 |

**与 [Y25] 足部骨折簇（7–9 m）对比**：

- **a=0**：首次骨折高度 5 m → ≈5.00 m，**早于** [Y25] 7–9 m。

- 旧整体足网格（均匀 E=15000，平滑网格，gauge ℓ=4 mm）首次骨折高度 ≈ **26.5 m**（文档记录，`docs/OpenSim_FE交接.md`）。新模型改用真实 AM50 跟骨 + CORT/SPON 双域，绝对高度**不可直接移植**，只作量级/方向核对。

## 逐点结果

| h (m) | a | F_subt (kN) | F_ach (kN) | max | p99 | p95 | mean | gauge_max |
|---|---|---|---|---|---|---|---|---|
| 5 | 0 | 26.38 | 0.00 | 213.50 | 177.22 | 136.30 | 61.279 | 186.76 |

## 诚实说明 / 既有 caveat（必读）

- **目标几何替换**：旧 S1–S3 的 `"calcaneus"` 域实际是 `r_foot.vtp` 的**整足**网格（121.7 cm³），不是单一跟骨；本模型的 `"calcaneus"` 是真实 AM50 跟骨（CORT 18279 + SPON 77794 mm³）。两者数值**不可直接比较**，只能看趋势/方向。
- **骨位置**：AM50 跟骨按其**解剖标志**置于 `calcn_r` 帧 —— 距下关节面质心映射到原点，跖面朝 −Y，后结节朝 −X；并非随意移动。
- **距下关节载荷块**：本配准把**距下关节面质心**放到原点，故 `subtalar_joint` 面（dist<25 mm 且 n·(+Y)>0.2）覆盖整个距下关节区，形成正常的载荷块（面积 2495.3 mm²，旧 foot-frame 仅 222.0 mm²）。
- **峰值受边界奇异支配**：载荷面 / 跖面固定棱边处 σ_vm 随网格发散，故报告以 `gauge_max`（过程区 ℓ=4 mm 正则化）为强度判据，裸 `max` 仅作诊断。
- 跟腱 `a=1` 施加 10.885 kN（Rajagopal2015 三头肌 F_max 之和的等长上界），属**参数化锚定**，非具体激活时程。
- 本扫描是**筛选性**研究：绝对应力不是临床预测；σ_vm≥150 MPa 判据仅在后处理中应用，不做单元删除 / 不把断裂作为求解动作。

## 产物

- `s2_thums_height_sweep_plantar_bc_base.json` / `.csv`
- `s2_thums_height_sweep_plantar_bc_base.png`
- `S2_THUMS_REPORT_PLANTAR_BC_BASE.md`（本文件）
