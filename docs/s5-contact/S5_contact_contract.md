# S5（FE 脚-垫接触）接口契约 —— 协调层文档

> 由主协调（Sisyphus）写于 2026-10-05。**所有参与 S5 的 subagent 必须先读本文件。**
> 规格来源：`docs/S5接触方案.md`（权威）、`docs/非垂直落地扩展方案.md` §4/§6。
> 状态（2026-10-06 晚间）：✅ 核心求解器已实现（93 测试）· ✅ 6 项 G7 验证完毕 ·
> ✅ `docs/S5_QUOTABLE.md` 已定为取数唯一入口 · ⚠️ **波3 能力已接通但真实网格未验证承载**
> （合成网格门禁 93/19+1/11/11 全绿；**真实 THUMS 上 `invalid facets`**，见 `docs/S5_wave3_contact.md` 顶部补注）·
> ✅ **波4 已交付**（`nonvertical_s5_contact.py` + 报告 + JSON；**5 套门禁 93/19+1/11/11/13 全绿**；
> 门禁码 G0 `PASS` / G1 `N/A_minimum_version` / G7 `RETIRED_by_contract_条6-8` 符合 §5.1）
> **但 G2–G6 全 FAIL**（7 条 contact 行收敛 0 条）→ `docs/S5_wave4_impl.md`。
> ✅ **波5 已交付（2026-10-06 深夜）—— 真实网格从「不收敛」到「收敛且承载」**：三缺陷（跖面 **tri3→quad4**、
> 内向法向、FEBio 默认 `<Control>`）+ 加载路径（`PressureLoad` → **`use_rigid=True`**）全部修复；
> **17/17 case 跑完，8/10 条 contact 行 CONVERGED+CARRYING**（`nohold`/`press` 对照 DIVERGED），门禁
> **162 passed / 1 skipped**，`F_n` 口径修正为 **`施加载荷 − hold 弹簧力`**（原式在 pen0.1 时高报 1.53×）
> → `docs/S5_wave5_matrix.md` / `results/opensim_fe/NONVERTICAL_S5_CONTACT_WAVE5_REPORT.md`。
> ⚠️ 波5 的接触 `F_n`（含 hold 弹簧建模补充，INFERRED）与 penalty 政策仍属**未决建模决策**，
> **尚未**进入 `S5_QUOTABLE.md` 的 QUOTABLE 行。
> 启动口径见 §4 **条 6–8**（原「不启动」已按实测修订为「可启动但限定可信口径」）。

---

## 0. 三个 "S5" 的消歧（防串线，务必先读）

| 代号 | 含义 | 位置 | 本次做不做 |
|---|---|---|---|
| **S5-contact** ✅ | **FE 脚-垫接触 + 摩擦** | `docs/S5接触方案.md` | **本次做这个** |
| S5-swap | 换软垫 → 抱石（只换 GRF 输入） | `OpenSim_FE复现方案.md` §五 S5.1–S5.3 | ❌ 不做 |
| S5-supination | Phase-S5 踝旋后/内翻（R1） | `results/opensim_fe/NONVERTICAL_S5_SUPINATION_REPORT.md` | ✅ 已完成，勿动 |

**禁止**把 S5.1/S5.2/S5.3 编号当成 S5-contact 的子步（那是 S5-swap 的）。

---

## 1. 目标

把当前**规定式跖面 BC**（`plantar_bc="fixed"` / `"roller"` / `"spring"`）升级为**真实接触**：
在 FEBio 里让「足底（跟骨跖面）」与「垫」发生 contact + Coulomb 摩擦，从而得到
**接触压力分布 / 真实压入量 / 摩擦切向力**，并**消除 A7 边界奇异**。

两个版本（都要，分先后）：
- **§3.1 最小版**：垫 = **刚性面**（rigid wall）+ μ。不建泡沫实体。
- **§3.2 完整版**：垫 = **可变形泡沫实体** + 接触 + 摩擦。← 这一版才能给出"压入量"。

---

## 2. 关键物理输入（已从源码核实，直接用）

### 2.1 垫的本构（必须复现的曲线）
`src/climbing/pad.py:55,104-110`：

```
σ(ξ) = s_pl·(1 − e^(−ξ/ξ_pl)) + s_d·(ξ/ξ_d)^n + E_bot·max(ξ − ξ_d, 0)
s_pl = 35_000 Pa   ξ_pl = 0.15
s_d = 365_000 Pa   n = 3.56   ξ_d = 0.85
E_bot = 5.0e6 Pa
thickness = 0.20 m
```
（`pad.py` 里 ξ 用 `xi_dense_max` 记作 0.85；公式见其 docstring。）

### 2.2 验收靶（交叉校验）

> ⚠️ **2026-10-07 重算 —— 原 16–18 cm 窗口作废，修正窗口 = 13–15 cm（ξ ≈ 0.65–0.74）。**
> 原推导用**初始**接触面积 **0.06 m²**（实为 0.045 m²）反解，但 `pad.py` 的接触面积是
> **渐进**的 `A(κ)`：κ=0→0.045、κ=0.7→0.228、κ=1→0.300 m²。峰值压缩处面积已达 **0.201 m²**（3.4×），
> 故应力实为 **173 kPa**（不是 580 kPa），ξ 实为 **0.648**（不是 0.90）。

**重算依据**（`pad.py` 动态解，2 m / 77 kg，见 `scripts/opensim_fe/make_fall_animation.py` 现用口径）：

| posture | 最大压缩 | ξ |
|---|---:|---:|
| `feet-first-stiff`（S1 默认） | **12.96 cm** | 0.648 |
| `controlled-drop` | 13.99 cm | 0.699 |
| `toe-point` | 14.40 cm | 0.720 |
| `one-leg-awkward` | 14.76 cm | 0.738 |

⇒ **修正窗口 = 13–15 cm（ξ ≈ 0.65–0.74）**；S1 默认单值 **12.96 cm（ξ=0.648）**。

> **原（已作废，保留为时间线）**：2 m 坠落峰值 34.8 kN、双足面积 0.06 m² → 应力 580 kPa
> → 反解 §2.1 → ξ ≈ 0.90（≈18 cm）；FD 身体下沉 16.2 cm ⇒ 窗口 16–18 cm。
> ⚠️ OpenSim FD 的 **16.2 cm** 是**无接触约束**的刚体下沉（比 pad.py 的垫压缩大 ~3.2 cm），
> **不**用作窗口 —— 垫压缩的诚实来源是 `pad.py`（见 `S5接触方案.md` §6 / `opensim_grf.py` docstring）。

---

## 3. 接口（函数签名 —— 必须严格照此实现，便于拼装）

### 3.1 `src/climbing/coupling/pad_foam.py`（材料拟合，T1 负责）
```python
def pad_stress_pa(xi) -> np.ndarray        # 复刻 pad.py 的 σ(ξ)，向量化
def fit_pad_material(...) -> dict          # 拟合 FEBio 可表达的材料 → 参数 + 拟合报告
```
返回 dict 至少含：`{"febio_type": str, "params": {...}, "max_rel_err": float,
"xi_range": (0.0, 0.85), "table": [(xi, sigma_fe, sigma_pad), ...]}`。
**必须**说明 FEBio 侧到底用什么材料块（`hyperfoam` / `ogden` / `neo-Hookean` / 手写 XML），
以及**为什么**它能拟合"平台 + 压实"两段（PU 泡沫的压实段对常规超弹材料是难点）。

### 3.2 `src/climbing/coupling/pad_mesh.py`（垫网格，T2 负责）
```python
def pad_block_tet(nx, ny, nz, size_mm=(300.0, 300.0, 200.0), *, out_msh=None) -> dict
```
返回与 `meshing.vtp_to_tet` **同 schema** 的 dict（`nodes`/`tets`/`surfaces`/`node_sets`/`volume_mm3`/`min_jacobian`/`jacobian_frac_lt_0_3`/`method`），
`surfaces` 至少含 `"top"`（接触面）、`"bottom"`（固定面）、`"sides"`。单位 mm。

### 3.3 `src/climbing/coupling/fe_post.py`（后处理补全，T3 负责）
```python
def read_displacement(xplt_path, *, hdf5_path=None, last=True) -> np.ndarray   # (N,3) mm
def read_reaction_forces(xplt_path, *, hdf5_path=None, last=True) -> np.ndarray # (N,3) N
def contact_force_sum(...) -> np.ndarray | float   # 接触/支承反力的合力（用于守恒检查）
```
HDF5 路径：`states/{last}/node_data/displacement/all` 与 `.../node_data/reaction forces/all`。
**必须**复用 `scripts/opensim_fe/bone_feb._xplt_to_hdf5_patched` 的 FEBio 4.13 兼容性修补
（缺了这个会因无名 shell domain 而解析失败）。

### 3.4 `docs/S5_contact_notes.md`（接触配方提取，T4 负责，只读研究）
从 `scripts/ankle_fe/build_feb.py:652-686` + `temp/pyfebio_demo/gen_load_xml.py:47`
提取**在 FEBio 4.13 上真正跑通**的 contact XML 配方。

---

## 4. 验收（硬门槛，逐条可自动判定）

| # | 判据 | 来源 |
|---|---|---|
| G0 | **默认回归逐位不变**；`pytest tests\test_opensim_fe.py tests\test_febio_run.py -q` → **19 passed / 1 skipped** | `S5接触方案.md:107` |
| G1 | §2.1 的 σ(ξ) 被 FE 复现：ξ∈[0,0.85] 内 **max_rel_err ≤ 5%**（FOAM 版） | `S5接触方案.md:67-74` |
| G2 | 极限自检：`μ=0` + 垫极刚 → 收敛到 rigid-support / `spring k→∞` | 同上 |
| G3 | 接触守恒：接触合力 = 施加外力（≤2%）；穿透量有界 | 同上 |
| G4 | 摩擦锥：任何时刻 `|F_t| ≤ μ·F_n` | 同上 |
| G5 | **A7 消失**：`contact` 的 gauge 峰值远低于 `fixed`（边缘奇异被消除） | 同上 |
| G6 | 敏感性单调/有界：μ、垫刚度、penalty 扫描 | 同上 |
| G7 | **压入量 13–15 cm 窗口**（§2.2，2026-10-07 重算；原 16–18 cm 作废） | 本契约 |
| G7' | **实测修正（2026-10-06）：FE 给 15.14 cm，G7 判 FAIL** —— 见下方「G7 修正说明」 | 实测 |

> **⚠️ 引用限定（2026-10-06，针对 §G7 + §G7' 及其提及的 1.56/1.59/1.62/1.64、2.70×、15.14 cm）**
>
> **(i) 配方标注**：1.56/1.59/1.62/1.64 是 **OLD AUGLAG pen=1** 配方在 tet4 (3,3,2) 单网格上对
> `docs\S5_contact_convergence.md` §3.1 的历史重现（数值偏差 ≤0.1%）；2.70× 与 15.14 cm 来自
> **NEW PENALTY 0.1** 配方在 sink=170 mm tet4 (6,6,4) 上的单次测量（事后已证不可复现：
> 同配置 3 次给出 negJac 27/3/15 全部失败）。
>
> **(ii) 网格未收敛**：1.56–1.64 未做 mesh-convergence 扫描（仅是同配方下重现历史值 ≤0.1%）；
> 2.70× 是单网格（6,6,4）测量，3 网格 sweep 下 σ/σ_law 漂移 **+30.62 %** 且单调递增
> （`S5_g7_verify_mesh.md` §3）；NEW PENALTY 配方在 sink ≤50 mm 同样在 3 网格上漂移
> **−0.75 %…−5.47 %**（`S5_g7_verify_anchors_mesh.md` §3）；固定 ξ ≈ 0.70 时漂移
> **+43.13 %**（`S5_g7_verify_fixedxi.md` §4）。
>
> **(iii) 不可跨配方对比**：OLD AUGLAG pen=1 与 NEW PENALTY 0.1 是两条不同的曲线，σ/σ_law 在
> 4 个 anchor sinks 上单调偏移 **−5.85 %…−9.38 %**（`S5_g7_verify_anchors.md` §4.1），
> 约为 ≤0.1% 重现噪声的 60–100× ⇒ OLD 与 NEW 数字**不可并列画在同一坐标轴**。OLD 配方
> 装上 NEW stepper 立即在 closure `negJac=55` 失败 ⇒ 两配方**不可互换**。2.70× 与 15.14 cm
> 同属 NEW PENALTY 配方、同一单网格点 ⇒ 不构成跨配方对比，但其单 mesh 抽水结果**已被
> 自身不可复现证伪**。
>
> **(iv) 权威裁决**：`docs\S5_QUOTABLE.md` —— 本契约 §G7 修正说明保留为**事件时间线记录**，
> 但所有可引用数字以 `S5_QUOTABLE.md` 的 QUOTABLE / NOT-QUOTABLE 列为准。
> 当前唯一在 2 m 落地可引用且**网格收敛**的量是 **`F_n ≈ 34.4 kN`** 与
> **`σ_contact ≈ 0.574 MPa`**（NEW PENALTY 配方，tet4，固定 sink，跨 3 网格漂移 <0.6%）。

### G7 修正说明（2026-10-06，实测；2026-10-07 窗口重算）

> ✅ **2026-10-07 G7 重跑（窗口修正后）—— 判定仍为 FAIL。**
> **澄清**：FE 施加的 **580 kPa 是压头压力**（= 34.8 kN / 0.06 m² **压头**面积），**不是**垫应力 ——
> 所以 FE 的**加载目标是对的**；作废的只是 **§2.2 的窗口**（它把 580 kPa 误当**垫应力**反解）。
> 窗口修正为 **13–15 cm**（§2.2）后重跑（`python -m climbing.coupling.pad_contact`，2026-10-07）：
> FE 在 580 kPa 下压入 **15.14 cm** → **15.14 > 15 ⇒ FAIL**（边缘，差 0.14 cm）。
> **旧窗口 16–18 与修正窗口 13–15 都判 FAIL** ⇒ 该 FAIL 结论**对窗口不敏感**。
> FE 扫描（prescribed sinkage → 压入量）：50→4.94、100→9.63、150→13.84、170→15.27、200 mm→16.28 cm。

**原窗口 16–18 cm 本身就有问题**：它是我从 `pad.py` 的 **1D 双质点模型**反解出来的，而 `pad.py` 的泡沫参数在
它自己的 docstring 里就写明「⚠️ 这些是**文献量级**，不是本项目的实测标定」。

**3D FE 接触实测**（`python -m climbing.coupling.pad_contact`；2026-10-07 重跑，窗口已修正）：
- 载荷驱动路径：施加 **580 kPa**（= 34.8 kN / 0.06 m²）→ **压力驱动 deck 本身不收敛**（`FebioRunError`）；
  回退用「规定沉入扫描 → 在 FE σ = 580 kPa 处插值」（响应单调，插值合法，**非调参凑数**）。
- **实测压入量 = 15.14 cm**（xi = 0.697）→ 相对**重算后**窗口（13–15 cm）仍在上沿之外 ⇒ G7 **FAIL**。
- 同一压入量下 **FE σ / pad 律 σ = 2.70×** —— FE 的垫子比 1D 律**硬约 2.7 倍**。

**结论（须写进 S5 报告与论文）— 2026-10-06 更新，详见 `docs/S5_contact_convergence.md`**：
1. **可靠区间**：接触解仅在名义应变 **ξ ≤ 0.25**（压头下沉 ≤ 50 mm）稳健 —— `rc=0`、**穿透 ≈ 0**（<0.25 mm）、**FE σ / 1D 律 σ = 1.56 → 1.64** 平滑单调。
2. **G7 工况（ξ≈0.85）不可解**【**旧 AUGLAG 配方下成立** —— 新 PENALTY 配方已使「不可复现」**解除**、但 σ/σ_law **仍不收敛**，现行口径见**条 6–7**】：2 m 坠落要求压掉 85% 垫厚，超出 `hex→tet4` 泡沫的稳健域 → **3 次完全相同配置给出 3 个不同结果（negJac 27 / 3 / 15，全部失败）**，即该解**不可复现**；细化网格 (6,6,4)/(8,8,6)/(10,10,8) 与提高 penalty 同样反演（`temp/opensim_fe/s5_t9/repro.py`、`mesh_refine.py`）。
3. ⚠️ **因此 15.14 cm / 2.70× 一律不得引用**（既非标定值、也非可复现值）。**此前「2.70× = 大应变平冲头物理」的结论作废** —— 它取自收敛包线**之外**的解；T8 的**材料验证**（同一 Ogden 的均质自由侧面 FE ≡ 1D 律 ≤0.5%）**仍成立**，作废的只是其对高 sink 接触比值的推断。
4. **验收口径**：同边界条件下 **FE vs FE**；可靠区间内可信的平冲头增强为 **≈1.6**。
5. **收敛配方（已固化进库）**：`laugon="AUGLAG"` + `penalty=1` + `tolerance=0.005` + `aug_controls=AUGLAG_CONVERGENT_CONTROLS`（`gaptol=0.001/minaug=3/maxaug=200/smooth_aug=1`；FEBio 4.13 的 AUGLAG 默认子参数会发散）。默认参数不传 ⇒ 行为逐字不变。
6. ~~**覆盖 G7 需换数值策略**：六面体单元（而非 hex→tet 裂解）、近不可压大应变泡沫专用公式、或自适应重划分 —— 属后续工作；在那之前**波3/波4 仍不应启动**。~~

   > ### 🔄 2026-10-06 晚间修订 —— 上述阻塞的两半已分别处置
   >
   > 原文保留为**事件时间线记录**（勿据其开工），现行状态如下：
   >
   > | 阻塞条件 | 现状 | 依据 |
   > |---|---|---|
   > | ① G7 **不可复现**（三跑三结果 `negJac 27/3/15`） | ✅ **已解除** —— 换 `laugon="PENALTY"` + `penalty=0.1`（并**保留 `search_radius=20`，大而非小**）后 `rc=0 NORMAL / end_t=1 / 2400/2400 / negJac=0 / 无 NAN`，**复现 ×4**，力平衡 ≈1e-7 | `S5_contact_udg.md` §7；`S5_g7_verify_mesh.md` §3 |
   > | ② G7 的 **σ/σ_law 不网格收敛** | ❌ **仍失败** —— 4 项独立 FAIL：自由 ξ **+30.62 %**、固定 ξ **+43.13 %**、小应变 −0.75 %…−5.47 %、跨配方 −5.85 %…−9.38 % | `S5_g7_verify_{mesh,fixedxi,anchors_mesh,anchors}.md` |

7. **⇒ 波3/波4 现在可以启动，但必须限定在可信口径**（本条**取代**原文的「不启动」）：
   - **可引用**：仅限 `docs/S5_QUOTABLE.md` 的 **QUOTABLE** 行 —— 绝对量 `F_n ≈ 34.4 kN`、`σ_contact ≈ 0.574 MPa`、`σ_1D-law(ξ)`，以及 **ξ ≤ 0.25 可靠包线**内的结果；
   - **禁止**：任何 σ/σ_law 比值（含 `2.205`）作为"收敛值"；任何跨配方对比；`15.14 cm / 2.70×`；小应变下的 `F_n` 声称网格收敛；
   - **根因归属（本条同时作废早先的"平头边缘奇异"归因）**：比值发散的真凶是 **`penalty=0.1` 软罚的压力外溢** —— 压力溢出刚性足迹约一个单元宽 ⇒ 穿透随网格变（11.6→18.2→27.7 mm）⇒「固定行程」与「固定压入量」无法同时收敛。**既非边缘奇异，也非单元公式**（同 ξ ≤ 0.65 两单元一致 ≤ ±2.5 %）。依据：`S5_g7_verify_metric.md` §5/§6。

8. **波3 的最小版本就不在 G7 阻塞范围内**：`S5接触方案.md` §3.1 明确最小版 = **刚性面 + μ、不建泡沫实体、不依赖泡沫本构标定**，且标注「**推荐先做**」。G7 阻塞针对的是 §3.2 **泡沫完整版**。⇒ **最小版可先行**。

---

## 5. 约束（红线）

- **只读勿改**：`scripts/ankle_fe/`、`temp/pyfebio_demo/`、`temp/ankle_*`、`FE_PIPELINE_TODO.md`、
  `LEARNINGS` LRN-010~019、以及所有既有 `results/opensim_fe/*` 产物（可新增文件，不可覆盖）。
- **默认行为逐位不变**：任何新参数必须 opt-in，默认路径一个字节都不许变。
- 单位 **mm–N–MPa–s**；OpenSim 模型是 m（转 mm 时 ×1000）。
- 解释器 `D:\Project\climbing_fall_analysis\.venv\Scripts\python.exe`；`$env:PYTHONPATH="src"`；
  cwd = 仓库根；PowerShell（用 `;` 不用 `&&`）。
- FEBio：`D:\Program\FEBioStudio\bin\febio4.exe`；**用 `febio_run.run_febio`**，不要自己起 subprocess。
- 新增产物一律写进 `results/opensim_fe/` 或 `temp/opensim_fe/`，不要污染其它目录。
- `.feb` 段落顺序：**Loads → Boundary → Contact**（顺序错 = `unrecognized tag`）。
- 接触面必须是**真实单元面**（tet facet），不能是 npz 里的四边形。
- 不要 `as any` / `@ts-ignore` / 空 catch；失败必须显式抛错。

---

## 6. 汇报格式（每个 subagent 返回时必须给）

1. 改了/新增了哪些文件（绝对路径）；
2. 跑了什么命令、退出码、关键输出；
3. G0/G1… 逐条判定结果（通过/未通过 + 数字）；
4. 未解决的风险与下一步；
5. 一个**可复现的最短命令**。
