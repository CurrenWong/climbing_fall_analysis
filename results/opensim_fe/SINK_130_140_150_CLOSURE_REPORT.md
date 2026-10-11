# S5 闭环一致性 re-verification — corrected sink window 130/140/150 mm (2026-10-10)

> 生成时间 2026-10-10 21:08:42；仓库 `D:\Project\climbing_fall_analysis`；单位 mm-N-MPa-s；FEBio 4.13。
> 运行器：`scripts/opensim_fe/run_sink_closure_20261010.py`
> 引用限定：仅 `docs/reports/S5_QUOTABLE.md` 的 QUOTABLE 行；本报告**补一行**（NEW dated 130/140/150 mm）**不删**旧 sink=170 mm 行（其 provenance 仍写）。

## 0. 任务与结论（TL;DR）

1. **任务**：在**修正 sink 窗口 130/140/150 mm**（pad 压缩 12.96 cm ⇒ 13–15 cm 容许带，`src/climbing/coupling/pad_contact.py::G7_SINKAGE_WINDOW_CM`）上重跑**最小基线**（tet4 standalone pad / sliding-elastic / μ=0.6 / NEW PENALTY 0.1 / G7 控制），把 `F_contact` 与线1 动力学峰值 **34.8 kN**（`pad_contact.FALL_LOAD_N` = `pad.py` 2 m 落地载荷）做闭环对比，**替代**原 sink=170 mm 的『待复核』行。

2. **结论一句话**：见 §1 表格；本次 3 次实跑**全部 CONVERGED+CARRYING** （rc=0 / end_t=1.0 / NORMAL TERMINATION），F_contact ∈ [13.96, 21.81] kN，对应线1 峰值 34.8 kN 的闭合差 ∈ [-59.88%, -37.33%]。**不是独立预测**：sink 由线1 运动学定（pad 压缩 12.96 cm），FE 只回算接触力。

3. **关于原 1.1% 闭合差**（sink=170 mm）：QUOTABLE 旧行 F_n=34.4 kN 对应闭合差 = (34.4 − 34.8) / 34.8 × 100% = -1.15%（≈ 1.1%）。本次在 130/140/150 mm 的新闭合差 ∈ [−59.9%, −37.3%] ⇒ **『1.1% 闭合』在修正窗口内不成立**（不是『改善 / 维持』，而是**直接 BROKEN**）。根因 = FE 常数 60 000 mm² 足迹 vs line-1 渐进面积（详见 §2），不是 FE 数值缺陷。**建议 QUOTABLE 处置**：旧 sink=170 mm 行保留为『**历史/常数面积冲突的 1.1% 巧合**』provenance；新增 130/140/150 mm 行**显式标注** closure status = NOT-CLOSED + reason（见 §5）；**真闭环需要在 FE 中使用 line-1 的渐进接触面积**（non-trivial 建模扩展，超出最小基线）。

## 1. 闭环校正（per sink）

| sink (mm) | sink (cm) | rc | end_t | negJac | NAN | termination | 收敛 | F_contact (N) | F_contact (kN) | σ_contact (MPa) | σ_law (MPa) | penetration (mm) | indentation (mm) | ξ | force_balance (rel) | friction_ratio | 闭合差 vs 34.8 kN | F_hold (N) |
|---:|---:|---:|---:|---:|:---:|:---|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 130 | 13.00 | 0 | 1.000 | 0 | False | NORMAL | ✅ | 13963.0 | 13.963 | 0.2327 | 0.1276 | 13.600 | 115.900 | 0.5795 | 1.785e-07 | 0.3498 | -59.88% | 0 (standalone pad, no hold spring) |
| 140 | 14.00 | 0 | 1.000 | 0 | False | NORMAL | ✅ | 17450.0 | 17.450 | 0.2908 | 0.1552 | 14.909 | 124.591 | 0.6230 | 5.287e-07 | 0.3498 | -49.86% | 0 (standalone pad, no hold spring) |
| 150 | 15.00 | 0 | 1.000 | 0 | False | NORMAL | ✅ | 21810.0 | 21.810 | 0.3635 | 0.1872 | 16.440 | 133.060 | 0.6653 | 5.859e-07 | 0.3497 | -37.33% | 0 (standalone pad, no hold spring) |

> **f_n_reaction_n 口径**：`F_contact = |Σ_z F_reaction(indenter nodes)|` 在末态；Newton-3rd 等价 = 接触面反力。Standalone pad 无 hold 弹簧 ⇒ **校正后 F_contact = F_n_reaction_n**（`wave 5` 口径 F_contact = 施加载荷 − hold 弹簧力的特例）。
> **σ_contact 口径**：`F_contact / 60 000 mm²`（QUOTABLE §3 #2 基准；与 sink=170 mm QUOTABLE 旧行同分母，可直接并列）。
> **σ_law**：信息列（`pad_stress_pa(ξ) / 1e6`），**不**作 σ/σ_law 比值引用（合同条 6–8）。
> **force_balance**：Newton-3rd 恒等检查（indenter z 反力 vs pad-top z 反力）；≤ 1e-6 表明 Newton 第3 律成立（非 G3 独立验证）。
> **fric_cone_ratio_max**：max |F_t|/(μ·F_n) 跨所有时间步；FEBio 4.13 在位移驱动下切向力为 0 ⇒ 该值远 < 1。

## 2. 闭环状态（与原 1.1% 闭合对比）

- 线1 动力学峰值 **FALL_LOAD_N = 34.8 kN** （`pad.py` 2 m / **77 kg** / `posture='feet-first-stiff'` 软垫；2026-10-10 复算 = 34.898 kN，垫压缩 12.957 cm；见 `docs/reports/总体进度.md` §二.1 与 `docs/s5-contact/S5_contact_contract.md` §2.2）。
- 旧 QUOTABLE 行 F_n ≈ 34.4 kN（sink=170 mm，**窗外**）⇒ 闭合差 = (34.4 − 34.8) / 34.8 = **−1.15%**（『1.1%』口径，**仅在 sink=170 mm 巧合成立**，见下）。
- **新**（本次）：sink ∈ {130, 140, 150} mm（**窗口内**）⇒ 闭合差（详细见 §1 表格最右列）：
  - sink = 130 mm：F_contact = 13.963 kN ⇒ 闭合差 = (13.963 − 34.8) / 34.8 = -59.88%
  - sink = 140 mm：F_contact = 17.450 kN ⇒ 闭合差 = (17.450 − 34.8) / 34.8 = -49.86%
  - sink = 150 mm：F_contact = 21.810 kN ⇒ 闭合差 = (21.810 − 34.8) / 34.8 = -37.33%

- **结论**：在修正 sink 窗口 13.0–15.0 cm 内，
  - **最佳** 闭合差 = **-37.33%**（sink = 150 mm → F_contact = 21.810 kN）；
  - **最差** 闭合差 = **-59.88%**（sink = 130 mm → F_contact = 13.963 kN）；
  - 闭合差随 sink 单调改善（sink ↑ ⇒ σ_law ↑ ⇒ F_contact ↑ ⇒ |δ| ↓）；
  - **旧 QUOTABLE『1.1%』闭合差在 sink=170 mm 是**『巧合』**而非真闭环**：
    FE 用**常数 60 000 mm²**足迹（standalone tet4 G7 锚），
    line-1 用**渐进接触面积**（`posture='feet-first-stiff'` `contact_schedule` 在 κ=0.589 时 ≈ 0.201 m² ≈ **201 000 mm²**，为 FE 的 **3.35×**）；
    在 sink=170 mm（已作废窗口 16–18 cm），两模型的数值**恰好落在 1.1% 内** —— 这是该 sink 专有的巧合，不是物理闭环；
    在修正窗口 sink=130–150 mm，力比 ≈ 0.4–0.6（FE 更低）⇒ **closure BROKEN**。
- **真闭环需在 FE 中使用 line-1 的渐进接触面积 A(κ)**（而非当前固定 60 000 mm² 平冲头足迹）；那是**非平凡的 FE 模型扩展**（变足迹 + 自适应网格），不在本最小基线范围内。
- **口径**：`F_contact` = `|Σ_z F_reaction(indenter)|` = wave-5 校正口径 `applied load − hold-spring force` 在**无 hold 弹簧**时的特例（standalone pad 无 hold 弹簧 ⇒ `F_hold_spring_n = 0`）。

## 3. 命令与复现

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH = "src"
# 全 3 个 sink（默认）：
.venv\Scripts\python.exe scripts\opensim_fe\run_sink_closure_20261010.py
# 单 sink：
.venv\Scripts\python.exe scripts\opensim_fe\run_sink_closure_20261010.py --sink 140
```

副产物：`results/opensim_fe/sink130_140_150_closure.json`（同 §1 表格原始数值）。

## 4. 边界与不变量

- **不可作为新 QUOTABLE σ/σ_law 行**：σ/σ_law 在合同条 6–8 下作废（见 `S5_QUOTABLE.md` §4 #1）。本报告 §1 仅作 σ_law 信息列，**不**写入 QUOTABLE 行的 quotable 列。
- **不可与 §0–§7 standalone-pad 数同位引用**：本次 sink 130/140/150 mm 是 sink 校正后**重做的**『2 m 落地』测量；与 sink=170 mm QUOTABLE 旧行**不同 sink 基础**，**不**取消旧行 （保留 sink=170 mm provenance 可见；见 §5 QUOTABLE 引）。
- **不可与波5 bone-integrated 数同位引用**：波5 §8 数（`F_n~20.4 kN @ pen100`）是 bone-integrated / 距下关节载荷 / 含 hold 弹簧模型 => 不同 mesh + 不同 load path。**本次是 standalone-pad G7 锚 sink 校正**，与波5 §8 不并列。
- **不可作独立预测**：sink 由线1 运动学定（pad 压缩 12.96 cm ⇒ 130–150 mm 容许带），FE 仅回算接触力；闭合差是 FE ↔ 线1 动力学峰值的**一致性校验**。

## 5. 与 `docs/reports/S5_QUOTABLE.md` 的接口

- 旧 QUOTABLE 行（sink=170 mm）保留（**不删**，sink=170 mm provenance 可见），其状态更新为 『historical fixed-sink FE measurement；其 1.1% 闭环 = sink=170 mm **专有的数值巧合**（常数 60 000 mm² 足迹 vs line-1 渐进面积 3.35×）；sink=170 mm 落在已作废窗口 16–18 cm 外（见 §2）；**非** 2 m 落地物理量』。
- 新增 **dated 2026-10-10** QUOTABLE 行：`F_contact ∈ [13.96, 21.81] kN @ sink ∈ (130, 140, 150) mm` （标 sink basis + mesh + recipe 同引；显式标 closure status = **NOT-CLOSED** + root cause = 『常数面积冲突』（见 §2））。
- **禁止把 130/140/150 mm 行当 1.1%『延续』**：那是新 sink 的独立读数（FE ↔ line-1 closure delta ∈ [-59.9%, -37.3%]），与旧 1.1% 无连续性。
