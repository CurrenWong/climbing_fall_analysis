# G4 · 全量 FE 损伤材料 · 工作计划（② FE 级深水区）

> 状态：**计划稿 v2**（2026-10-11；经 **Momus 完整性评审** + **联网可行性验证**后修订；**待用户批准 → 预注册冻结**）。
> 定位：`docs/methods/G4_损伤演化设计.md`（最小路径，已跑 **A4 FAIL**）的**升级件**——把 σ(h) **后处理**损伤律升级为 **FEBio 原生损伤材料 + 逐高度 FE 重解**。
> 上游：`docs/methods/目标重构_刚性验证与垫上推断.md` §三（G4/G6 缺口）；`docs/reports/最终结论_2026-10-10.md` §八。
> 实证来源：**FEBio 4.13 官方手册**（本机 KB `C:\Users\curre\.copaw\workspaces\default\knowledge\books\engineering\`）+ **FEBio 4.13.0 源码 / 官方测试实例**（联网核实）。
>
> **v2 变更摘要（2026-10-11）**：
> - Momus 3 阻塞已修：§5 重写 D1–D6 新值 · §3/§5 加量化门 · P2/P3 补起点或降级；
> - Momus 5 次级已修：引用校正（`fracture_matrix.py:606-608`）· 映射唯一化 · FE 验证场景 · 分级 DoD · 依赖核对；
> - 联网 2 新风险已并入：**R1 判据**（改 `DC Drucker shear stress`）· **R2 正则化**（FEBio 无 → 稳定性/网格客观性受限）。
> - 逐条见 **§9 评审修订记录**。

---

## 0. 一句话目标

用 **FEBio 原生损伤本构**（`elastic damage` / `uncoupled elastic damage`，参数取自 THUMS AM50 材料卡）替换 σ(h) 后处理的**逐部位统一缩放**，取得**单元级应力重分布**——主目标 **A4 rank-2**（统计 1/4 → 2/4），附带 **A2 股骨/胸椎**、**V7 阶跃/渐进**（加分项）。

---

## 1. 价值与边界（诚实）

| 目标 | 射程 | 机制 |
|---|---|---|
| **A4 rank-2**（当前 rank-1/FAIL） | ✅ **主目标** | 单元级损伤演化 ⇒ σ(h) 形状分化（后处理拿不到） |
| **A2 股骨/胸椎** | ➕ 加分项 | 需**同时**修 FE 载荷分配（P2） |
| **V7 阶跃/渐进** | ➕ 加分项 | 需骨盆环多骨几何（P3） |
| **P-4 踝骨折率** 上界→点估计 | ➕ 潜在最大收益 | 踝 FE 并入损伤材料，破"骨折≡卡住"判据结构 |
| **A3 腓骨 OR** | ❌ **补不上** | G8/G9 跨路线 A/B 已证 ≈1.13–1.15 **结构性上限** |
| **V5 分布形状** | ❌ 非 FE | 需多体动力学自由度（上肢/膝/肘） |
| **V8b @20 m ≤5 ms** | ❌ 非 FE | 需应力波模型层 |

**天花板**：统计 **1/4 → 最多 ~3/4**（全套），单做 A4 → 2/4。**补不满。**（DoD 分级见 §6，避免加分项未达成被误判为整体 FAIL。）

---

## 2. 材料来源（实证，已确认）—— 本计划的核心

**"全量损伤材料"不需要外求**：参数已在仓库，本构 **FEBio 原生自带**。

### 2.1 参数：THUMS AM50 材料卡（仓库内，非调参）

- master deck：`model/AM50_V71_Occupant/main_THUMS_AM50_V71.k`（79 行）；`*MAT_PIECEWISE_LINEAR_PLASTICITY` 卡实际位于其 `*INCLUDE` 的 **`THUMS_model/mat_bone_AM50_V71_Occ_no_fracture.k`**（`read_material_damage` 走 include 树解析，功能正确）。
- 逐 PID：`SIGY`=col4、`ETAN`=col5。
- 既有解析：`scripts/opensim_fe/bone_feb.py::read_material_damage()`（`Material.yield_stress_mpa` / `tangent_modulus_mpa`）。
- 逐部位取值：`results/opensim_fe/G4_DAMAGE_REPORT.md` §1（腓骨 `81000900` σ_y=37 / E_t=16700；胫骨 `81000700` …）——**G4 最小实现已用同一套参数**（故全量的增量只在"单元级重分布"，不在参数）。

### 2.2 本构：FEBio 4.13 原生损伤材料（**已确认可用**）

> 来源：FEBio 4.13 官方手册（KB）+ **4.13.0 源码**（`FEBioMech/FEBioMechModule.cpp` 注册表）+ 官方测试实例 `dm01/dm05/dm09/dm14.feb`。

| 要素 | FEBio 实证 | 出处 |
|---|---|---|
| **材料类型** | `elastic damage`（约束基体）/ `uncoupled elastic damage`（近不可压，带 `k`） | User Manual §4.6.1（L13005） |
| 损伤框架 | Reactive damage：`σ=(1−D)σ₀`，D = 断键质量分数，`Dmax∈[0,1]` | §4.6（L12970） |
| **损伤判据（v2 修订）** | **主：`DC Drucker shear stress`**（压敏；α=0 时退化 von Mises）；**对照：`DC von Mises stress`**、`DC max normal stress`、`DC max normal Lagrange strain` 等 | Library §4.6.3（L13276 `DC von Mises`） |
| **CDF（软化）** | `CDF Weibull` / `Simo` / `log-normal` / `quintic` / `gamma` / `step` / `user` | §4.6.2（L13028+） |
| 弹塑 + 损伤 | 区分 intact bond 损伤 vs **yielded bond** 损伤 | Theory Manual（L15382） |

**⚠️ 硬性约束（源码实证）**：
1. **coupled/uncoupled 不匹配 = hard error**：`elastic damage` 的 `Init()` 拒绝 uncoupled 基体（报 `"Elastic material should not be of type uncoupled"`）。⇒ 近不可压骨（带有限 `k`）须用 **`uncoupled elastic damage`** + uncoupled 基体。
2. CDF 只接受**非负**判据值 ⇒ 配合 ≥0 的判据（von Mises / Drucker / SED / max normal Lagrange strain）。
3. 软化须时间步退让（`max_retries`/`dtmin`/`cutback`）。

**唯一映射表（冻结，消除"跑后调参"空间）**：

| 要素 | 冻结值 | 来源（非调参） |
|---|---|---|
| 基体本构 `<elastic>` | **neo-Hookean**（有限变形下病态单元更稳健） | 官方示例 + FEBio 作者建议 |
| 损伤判据 `<criterion>` | **`DC Drucker shear stress`** | 骨拉压不对称（见 §4 R1）；α 见下行 |
| Drucker 压敏常数 α | **由骨拉压强度比（~1.4–1.6）反推，取自文献**（非拟合） | Morgan 2018 / Bayraktar 2004 |
| CDF 族 | **`CDF Simo`**（单参数可覆盖饱和/渐进） | G4 设计 §4 选型 |
| 阈值（SIGY→） | 逐 PID THUMS `SIGY` | THUMS 卡 |
| 软化（ETAN→） | CDF 参数（**模型类别变更，非参数搬运**，见 §4 R1） | THUMS 卡 + 文献 |

示意卡（结构，非最终值）：

```xml
<material id="1" type="elastic damage">
  <elastic type="neo-Hookean">
    <density>DENS</density><E>E_THUMS</E><v>NU</v>
  </elastic>
  <criterion type="DC Drucker shear stress"> <a>ALPHA_FROM_LIT</a> </criterion>
  <damage type="CDF Simo"> <a>..</a><b>..</b><Dmax>DMAX</Dmax> </damage>
</material>
```

⇒ **P0 的"FEBio 无等价材料"风险已解除**；剩余未知 = **单解成本**（P0 计时）、**判据适用性**（R1）、**网格客观性**（R2）。

---

## 3. 分阶段工作包

### 🔴 P0 · 可行性 spike（≤1 天，先做，决定投不投）
> 目的：把"能不能做 / 值不值得做"从假设变事实。**同时验 R1/R2**。

1. `febio4.exe` 材料库枚举，确认 `elastic damage` / `uncoupled elastic damage` + `DC Drucker shear stress` 在 4.13 本构建可用（对照 `dm05.feb`）。
2. 取**一个**部件（腓骨或胫骨）+ **1 高度**，用 THUMS `(E, ν, SIGY, ETAN)` + §2.2 唯一映射表建模型，跑通。
3. 对照既有弹性件 σ_ref；记录**单解墙钟 + 收敛 + 网格依赖**（做 2 个网格密度）。

- **产物**：`scripts/opensim_fe/g4_fe_spike.py` + `results/opensim_fe/G4_FE_SPIKE_REPORT.md`
- **DoD（量化门，v2）**：
  - 单骨单高跑通，σ_ref 与弹性件**网格收敛后**相对偏差 ≤ 5%（容差，v2 明确）；
  - **单解墙钟 ≤ 5 min**；
  - 网格 2 密度下首折高度偏差 ≤ 10%（**R2 网格客观性初判**）。
- 🛑 **止损线（量化，v2）**：若单解 > 5 min **或** 9 部件 × 目标高度数总墙钟 > **12 h** **或** 网格客观性不可接受（偏差 > 10%）⇒ ② 记"已验证不可行/缩水"，回退 ③ 或仅做非软化版。

### 🟠 P1 · 全量损伤材料 FE（核心，**主目标**）
- 全 **9 部件**建 FEBio 损伤模型（替换 `bone_feb.py` 弹性近似）。
- `fracture_matrix.py` 的**线性缩放**（**L606–608**：`ratio = load_h/ref_load; sig = …*ratio`）→ **逐高度 FE 重解**。
- 重算 σ(h) → `a4_rank2_from_sigma` → **如实**判 A4。
- **产物**：`g4_fe_matrix.json` + `G4_FE_REPORT.md`（新件；冻结输入只读）

### 🟡 P2 · FE 载荷分配修复（➕ 加分项，**只治股骨/胸椎**）
- **起点（v2 补齐）**：`scripts/opensim_fe/` 既有 OpenSim 关节反力链（`run_dead_drop` / `joint_reaction`）+ `s4_alignment_bc2.py` 的载荷向量；对股骨/胸椎引入**各关节真实分配**（替代当前"缺失分配"）。
- **验收（v2 补齐）**：A2 股骨/胸椎从 `all_zero` 变为可分类状态（阶跃/渐进）。
- **不做**：腓骨 OR（永久范围外）。

### 🟢 P3 · 骨盆环多骨几何（➕ 加分项，治 V7）
- **起点（v2 补齐）**：THUMS 抽取机器 `extract_part.py` / `part_meshing.py`；源部件 = 髋骨（`R_HIPBONE`）+ 骶骨。
- **验收（v2 补齐）**：V7 逐部位 `classify_step_progressive` 与 [Y25] `PAPER_STEP`/`PAPER_GRADUAL` 一致数 **≥3/5**（v2 建议冻结值；批准时确认）。

### ⚫ 非 ②（另立预注册）
- V5 → 多体动力学（上肢/膝/肘自由度）；V8b → 应力波模型层。

---

## 4. 风险 / 成本 / 止损（含 2 新风险）

| 风险 | 概率 | 缓解 / 结论 |
|---|---|---|
| ~~FEBio 无等价材料~~ | ~~中~~ → **低（§2.2 已确认有）** | 材料已实证 |
| **R1 · 判据（von Mises 对骨不合适）** | **高** | 骨拉压不对称+各向异性；von Mises 等向/对称/压敏无关 ⇒ **模型形式错误风险**（Morgan 2018："isotropic and symmetric criteria… not capable of describing the multiaxial strength"）。**缓解**：改用 `DC Drucker shear stress`（压敏，α 取自骨拉压比 ~1.4–1.6），von Mises 仅作敏感性对照。诚实注：Keyak&Rossi 2000 发现 von Mises 对**整体股骨骨折载荷回归**表现好——仅限"载荷回归"，不含损伤定位。 |
| **R2 · 正则化（FEBio 无）** | **高** | FEBio reactive damage **纯局部**（源码：无特征长度/无空间平均/无梯度）；**无** nonlocal/gradient/crack-band/phase-field；仅 `erosion`（**非**正则化，网格依赖）+ `Dmax<1`。**完全软化 → 单元倒置 → 负 Jacobian 终止**（FEBio 作者亲述）。**缓解**：`Dmax<1` + neo-Hookean 基体 + P0 mesh-refinement 初判；**若须严格网格客观 ⇒ FEBio 不可行**（须自写 plugin 或换 ABAQUS）。 |
| 全量 FE **不比后处理强**（同参数） | 中 | P0 单骨对照先证；弱则如实报限 |
| 逐高度重解成本爆炸 | 中 | P0 计时门（§3 P0）；超预算砍高度 |
| 网格依赖（软化病态） | 中 | = R2 |

**投入估算**：P0 = 1 天；P1 = 数天–1 周（含调试）；P2/P3 各数天。**全部 = 多周**，天花板 ~3/4。**R2 若触发硬约束 ⇒ 投入翻倍或不可行。**

---

## 5. 预注册冻结决策（**本件新值，v2**）

> v2 关键修正：**不再"沿用" G4 设计 §6**——本件在 D1/D5 上**取代**该设计的值，逐条写新值。

| # | 决策 | **本件冻结值（v2）** | 相对 G4 设计 §6 |
|---|---|---|---|
| **D1** | 本构 | **FEBio 原生 `elastic damage` / `uncoupled elastic damage`**（§2.2 唯一映射表） | **取代**（原：双线性弹塑·降阶后处理） |
| **D2** | 屈服来源 | 逐 PID THUMS `SIGY`（**非调参**） | 沿用 |
| **D3** | 软化/压敏来源 | CDF 参数 + Drucker α（α 取自**文献**拉压比，**非拟合**） | 沿用+细化 |
| **D4** | rank-2 度量 | PC2 解释方差比 vs **平行分析 null / 断棍 null**（沿用 `a4_rank2_from_sigma` 口径） | 沿用（口径不变） |
| **D5** | 集成点 | **逐高度 FE 重解**（替换线性缩放 `fracture_matrix.py:606-608`）；冻结输入只读 | **取代**（原：σ(h) 后处理层·不改 FE） |
| **D6** | 载荷链 | 单一冻结 `PassiveStiffness.nominal_landing()`；**不扫不调** | 沿用 |

> **冻结即锁**：跑后**禁止**调 `SIGY`/`ETAN`/CDF 参数/α 去凑 A4/A2/V7。

---

## 6. 分级 DoD（v2）

**主目标（决定 ② 成立与否）**
1. P0 spike 落盘 + 量化门通过（§3 P0）；
2. 全 9 部件 FEBio 损伤模型跑通；
3. **A4 rank-2/LP 如实 PASS/FAIL** 落盘（新件）；
4. 字节锁（`test_pad2d.py`）+ G1 + 波杆 + G4 全绿，全回归 **≥709 passed / 8 skipped** 不减；
5. 纯增量（冻结输入只读，`git diff --stat` 确认）。

**加分项（不达成不算整体 FAIL）**
6. A2 股骨/胸椎获区分度（P2）；V7 ≥3/5（P3）；
7. FE 验证场景通过（见 §7）。

---

## 7. FE 验证场景（v2 新增）

新 FE 损伤求解须有**等价回归锁**（参照 G4 最小实现"转录规则逐值等于基线"）：
1. **低载一致性**：在 `σ_lin < SIGY`（未屈服）区，FE 损伤解 σ(h) 与既有弹性解**逐值一致**（容差 ≤1e-6 相对）；
2. **D→0 验证**：低载下 `<var type="damage"/>` 全场 D≈0；
3. **收敛性**：时间步退让收敛（无负 Jacobian，除全损终止外）；
4. **网格无关性**：2 密度下结果偏差 ≤10%（与 P0 门一致）。

---

## 8. 依赖与复现

- **FEBio**：`D:\Program\FEBioStudio\bin\febio4.exe`（4.13.0）✅
- **手册**：本机 KB `...\knowledge\books\engineering\febio-{theory,user,studio}-manual\` ✅
- **解释器**：`.venv\Scripts\python.exe`（`PYTHONPATH=src`）✅
- **材料卡**：`model/AM50_V71_Occupant/main_THUMS_AM50_V71.k`（含 include 树）✅
- **参考实例**：FEBio TestSuite `dm01/dm05/dm09/dm14.feb`

---

## 9. 评审修订记录（2026-10-11）

**Momus（完整性，判定 [REJECT]）→ 已修**
- 阻塞1 §5 冻结决策矛盾+缺新值 → **§5 重写 D1–D6**（标注"取代/沿用"）；
- 阻塞2 量化门缺失 → **§3 P0 止损线 + §6 V7 命中数**（≥3/5）；
- 阻塞3 P2/P3 无起点 → **§3 P2/P3 补起点+验收**（保留为加分项）；
- 次级A 引用不精 → 校正 `fracture_matrix.py:606-608`、卡在 include 树、"A4 rank-2"；
- 次级B 映射歧义 → **§2.2 唯一映射表**；
- 次级C 缺 FE 验证 → **§7 新增**；
- 次级D 范围 → **§1/§6 分级 DoD**；
- 次级E 依赖 → 已核对 OK。

**联网验证（可行性）→ 已并入**
- FEBio 材料 API：✅ 确认（4.13.0 源码 + `dm*` 实例）→ §2.2；
- **R1 判据**：von Mises 对骨 = 模型形式错误风险 → §2.2 改 **`DC Drucker shear stress`**（+von Mises 对照）、§4；
- **R2 正则化**：FEBio 无正则化 → §4 + §3 P0 网格门 + §7。

---

*本件为 ② 的活跃工作计划（v2）；P0 结果用于决定 P1 是否上马；R2 若触发硬约束则整体可行性重估。*
