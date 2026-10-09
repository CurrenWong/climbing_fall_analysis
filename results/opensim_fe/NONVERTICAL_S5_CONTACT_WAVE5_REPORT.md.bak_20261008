# S5-contact 波5 · 真实泡沫垫矩阵报告（G7 配方 + hold 弹簧 + rigid 加载路径）

> 生成时间 2026-10-07 00:49:38；仓库 `D:\Project\climbing_fall_analysis`；单位 mm-N-MPa-s；FEBio 4.13。
> 计划 `docs/S5_wave5_matrix.md`；契约 `docs/S5_contact_contract.md` §4/§5/§6。
> 引用限定：仅 `docs/S5_QUOTABLE.md` 的 QUOTABLE 行。F_n 一律 **INFERRED**（见 §6）。

## 0. 一句话结论

1. 真实泡沫垫矩阵（anatframe）上，**7 条先前失败 contact 行现全部 CONVERGED+CARRYING**（rc=0、end_t=1.0）。
2. **hold 弹簧必需**：同一 case 去 hold ⇒ 2–3 s 内发散（接触闭合前 3 刚体模态）。（converged+carrying contact 行数 = 8）
3. **加载路径是剩余阻塞**：波4 contact 行用 `use_rigid=False`（PressureLoad）⇒ 在接触闭合 t≈0.42 发散（阴性对照 `_press`）；本波改用 `use_rigid=True`（= proof 与 S4 主口径）后收敛。
4. penalty 敏感性：pen0.1 → penetration ≈ 40.45 mm；pen100 → 2.19 mm；两者均承载（≪ 无接触 sink ≈110 mm）。

## 1. 关键发现：加载路径 + hold 弹簧

| 配置 | 网格 | use_rigid | hold | rc | end_t | 结果 |
|---|---|---:|---:|---:|---:|---|
| rigid + hold | anatframe | True | 1.0 | 0 | 1.0 | **CONVERGED+CARRYING** |
| pressure + hold | anatframe | False | 1.0 | 1 | 0.419308 | **DIVERGED** |
| rigid, no hold | anatframe | True | None | 1 | — | **DIVERGED** |

- rigid-body 加载路径收敛且承载；PressureLoad 路径在接触闭合 (t≈0.42) 处 FEBio 收敛判据退化 ⇒ ERROR。波4 contact 行用 use_rigid=False，这是矩阵与 proof 的真实差异之一。

## 2. 逐 case 结果表（全部 17 行）

| case id | scenario | use_rigid | penalty | hold_k | febio_rc | end_t | negJac | plantar_faces | log_invalid | penetration (mm) | contact_area (mm²) | carries | verdict |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---|
| g0_fixed_axial_ref | on_pad | False | — | — | 0 | 1.0 | 0 | — | 2 | — | — | None | PASS |
| h2_fixed_pad | on_pad | False | — | — | 0 | 1.0 | 0 | — | 2 | — | — | None | PASS |
| h3_fixed_pad | on_pad | False | — | — | 0 | 1.0 | 0 | — | 2 | — | — | None | PASS |
| h3.5_fixed_pad | on_pad | False | — | — | 0 | 1.0 | 0 | — | 2 | — | — | None | PASS |
| h4.5_fixed_pad | on_pad | False | — | — | 0 | 1.0 | 0 | — | 2 | — | — | None | PASS |
| h2.0_spring_k100_pad | on_pad | False | — | — | 0 | 1.0 | 0 | — | 2 | — | — | None | PASS |
| h2.0_contact_rigid_mu0_pad | on_pad | True | 0.1 | 1.0 | 0 | 1.0 | 51 | 137q/0t | 0 | 40.454 | 2953.8 | True | CONVERGED+CARRYING |
| h2.0_contact_rigid_mu0.3_pad | on_pad | True | 0.1 | 1.0 | 0 | 1.0 | 29 | 137q/0t | 0 | 40.454 | 2953.8 | True | CONVERGED+CARRYING |
| h2.0_contact_rigid_mu0.6_pad | on_pad | True | 0.1 | 1.0 | 0 | 1.0 | 29 | 137q/0t | 0 | 40.454 | 2953.8 | True | CONVERGED+CARRYING |
| h3_contact_rigid_mu0.6_pad | on_pad | True | 0.1 | 1.0 | 0 | 1.0 | 19 | 137q/0t | 0 | 42.285 | 2953.8 | True | CONVERGED+CARRYING |
| h3.5_contact_rigid_mu0.6_pad | on_pad | True | 0.1 | 1.0 | 0 | 1.0 | 27 | 137q/0t | 0 | 43.989 | 2953.8 | True | CONVERGED+CARRYING |
| h4.5_contact_rigid_mu0.6_pad | on_pad | True | 0.1 | 1.0 | 0 | 1.0 | 44 | 137q/0t | 0 | 48.120 | 2953.8 | True | CONVERGED+CARRYING |
| h2.0_fixed_hard_surface | hard_surface | False | — | — | 0 | 1.0 | 0 | — | 2 | — | — | None | PASS |
| h2.0_contact_rigid_mu0.0_hard_surface | hard_surface | True | 100.0 | 1.0 | 0 | 1.0 | 0 | 137q/0t | 0 | 2.315 | 146.3 | True | CONVERGED+CARRYING |
| h2.0_contact_rigid_mu0.6_pad_nohold | on_pad | True | 0.1 | — | 1 | — | 21 | 137q/0t | 0 | — | — | None | DIVERGED |
| h2.0_contact_rigid_mu0.6_pad_pen100 | on_pad | True | 100.0 | 1.0 | 0 | 1.0 | 0 | 137q/0t | 0 | 2.186 | 130.5 | True | CONVERGED+CARRYING |
| h2.0_contact_rigid_mu0.6_pad_press | on_pad | False | 0.1 | 1.0 | 1 | 0.419308 | 0 | 137q/0t | 2 | — | — | None | DIVERGED |

> `plantar_faces`：emitted deck 的 `<Surface name="plantar">` 内 `<quad4>/<tri3>` 计数（**跖面专属**；任务书写的 146 是 proof 的 calcnframe 网格，矩阵 anatframe = 137）。`log_invalid` = FEBio 日志全局 `invalid facets` 计数（收敛参照行基线 = 2）。

### 2.1 emitted deck 证据（代表 case `h2.0_contact_rigid_mu0.6_pad`）

```xml
<!-- plantar Surface: 137 quad4, 0 tri3 -->
<Surface name="plantar">
<quad4 id="1">
<quad4 id="2">
    ... (137 quad4 total, 0 tri3)
</Surface>
<time_steps>2400</time_steps>
<dtmax>0.0004166666666666667</dtmax>
<cutback>0.125</cutback>
<max_retries>20</max_retries>
<opt_iter>15</opt_iter>
```

## 3. hold 弹簧必要性（同 case 对照）

- `with hold k=1.0 N/mm`：rc=0, end_t=1.0, penetration=40.45384979248047, verdict=**CONVERGED+CARRYING**
- `without hold (None)`：rc=1, end_t=None, penetration=None, verdict=**DIVERGED**
- ⇒ **hold_spring_required = True**

## 4. penalty 敏感性

| 配置 | penetration (mm) | verdict |
|---|---:|---|
| pen0.1 | 40.454 | **CONVERGED+CARRYING** |
| pen100 | 2.186 | **CONVERGED+CARRYING** |

## 5. 验收标准（承载 ≠ rc=0）

- **CONVERGED**：`febio_rc==0` 且 `end_t==1.0`。
- **CARRYING**：`plantar_drop ≪ 无接触 sink`（= 载荷/191 N·mm⁻¹ ≈ 110 mm @ 21086 N）。
- **UNLOADED**：收敛但 drop ≈ sink。**DIVERGED**：rc≠0 / end_t<1。

## 6. F_n 诚实边界

- F_n 标 INFERRED（静力平衡）；FEBio 4.13 本 deck 不报固定 DOF 反力，reaction-sum 仅记 unreliable_diagnostic。
- 量级锚：`docs/S5_QUOTABLE.md` 的 **F_n ≈ 34.4 kN**（2 m 落地，tet4 pad / 固定 sink / NEW PENALTY，网格收敛 <0.6%）。本矩阵载荷是**距下关节**载荷 （h=2.0 → 21086 N），与 34.4 kN 是**不同 mesh 与不同 load path** ⇒ 只作量级锚，不并列。

## 7. 未解决的风险与下一步

- 波5 把 contact 行切到 `use_rigid=True`（proof 已验证的加载路径）。若下游要求 `use_rigid=False` 的 PressureLoad 口径，接触闭合处仍需换数值策略（如分段 ramp / AUGLAG / 更细接触网格）——本波只证明 rigid 路径可行，**未**解决 pressure 路径。
- `penalty=0.1`（G7 配方，为泡沫标定）在 bone-vs-rigid 上穿透 40 mm（仍承载）；pen100 降到 ~2 mm。建议报告 σ/穿透时同引 penalty。
- F_n 全靠力平衡推断；未解析接触压力场。

## 8. 复现命令

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\nonvertical_s5_contact_wave5.py
```

