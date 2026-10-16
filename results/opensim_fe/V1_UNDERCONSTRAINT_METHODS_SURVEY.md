# 跖面接触不收敛 —— **手册 + 网络方法检索与实测**（2026-10-16）

> **触发**：用户令「查网络和手册，找办法」（针对 `V1_PLANTAR_LOAD_BEARING_FEASIBILITY.md`
> 报告的**第 2 层不可行**：跖面真实柔度（`k_phys≈0.203 N/mm`）与 `contact` 无旋钮版全部不收敛）。
> **性质**：遵守项目铁律 `[MEMORY gotcha #26]`「**FEBio 问题先查手册+网络，再下不可修结论**」。
> **本件新增**：① 手册/web 处方逐条检索；② **发现并放宽我们自己的护栏**；
> ③ 四条处方的实测矩阵。

---

## §0 · 结论摘要

| # | 手册/web 处方 | 我们试过？ | 结果 |
|---|---|---|---|
| **A** | **避免初始分离 / 让 t=0 就有初始接触**（§8.12.4 / §8.5.3 + 官方 webinar） | ✅ **本轮首次** | ⚠️ **我们自己的代码原来禁止它**（已放宽）；实测 **单独用不足以收敛** |
| **B** | **两步法**：先位移压到接触，再把位移换成力（§8.5.3 + webinar） | ✅ 既有 | Step 1 ✅ / **Step 2 ❌**（残差 ~1e-16，既有件已记） |
| **C** | **False Transient**：准静态不收敛 ⇒ 改**动力**分析 + 数值阻尼（§8.4.2） | ✅ **本轮首次** | 🔴 **假象——实测是「刚体飞走」**（见 §3.3 更正）；**不是收敛** |
| **D** | **`symmetric_stiffness` = 非对称 / preferred**（FEBio 主动 WARNING + 手册 note 6） | ✅ **本轮首次** | ❌ 放 `<solver>` 仍失败；放 `<Control>` **被 FEBio 拒绝**（位置不对） |
| **E** | `auto_penalty`（§8.5.1） | ❌ 未试 | 待试 |

**⇒ 一句话：手册列出的五条处方（A/B/C/D/E）本轮**全部试过**，**
**无一条使跖面接触（`hold_spring_k=None`）收敛**；C 的「进展」经核是**刚体飞走**（§3.3），
**唯一能收敛的仍是 `hold_spring_k`（自由旋钮 ⇒ 只能当诊断口径）。**
**⇒ 「跖面接触在无额外约束下不收敛」= 一个**未被手册处方解决**的问题。**
**本轮最实在的收获是**找回被自己护栏挡掉的处方 A**（§2）与**三条可复用的工具知识**（§3.3）。**

---

## §1 · 依据（逐字）

### FEBio 手册 §8.12.4 Force-driven Problems

> 「When the primary method of deforming the model is a force, the problem is said to be
> **force-driven**… Force-driven forces are inherently unstable and FEBio has limited
> capabilities of handling such problems.
> … This is often the case in **force-driven contact problems**, where the contact interface
> is necessary to define a well-constrained model. **When there is an initial separation**
> (or even when the surfaces have no initial overlap) the initial contact force is zero,
> which is equivalent to having no contact enforcement, and thus **the model is
> under-constrained**. Therefore, **for force-driven contact problems it is important to
> have some initial contact before starting the analysis.**」

### FEBio 手册 §8.5.3 Initial Separation

> 「FEBio often struggles with problems that have an initial separation. This is especially
> so when the problem is force-driven… **it is best to avoid initial separation if possible.**
> Another way around this issue is to **first use a displacement to bring the surfaces in
> contact and then replace the displacement with a force**.」

### FEBio 手册 §8.3.1.5 Rigid body modes

> 「In the presence of a rigid body mode, the solution is not uniquely defined and FEBio will
> most likely throw **negative Jacobians**. The solution is to identify the rigid body mode
> and to **eliminate it by applying additional constraints** to the model.」

*(⇒ `hold_spring_k` 正是这条；但它同时是**自由旋钮** ⇒ 只能当诊断口径。)*

### FEBio 手册 §8.4.2 False Transient Analysis

> 「Quasi-static models may sometimes exhibit instabilities that may prevent the default
> Newton-based solver from converging. In such situations, **a dynamic analysis can be tried
> instead**, with the anticipation that **the inertia of the model provides sufficient
> stability** to solve the model successfully.」
> 「To introduce the most numerical damping: `<rhoi>-2</rhoi> <gamma>1.5</gamma> <beta>1</beta>`」
> §8.4.1：「The main parameter that controls the time integration of the **solid solver** is
> the `rhoi` parameter.」

### FEBio 手册 §8.5.1 The penalty method

> 「…choosing the penalty factor too large may cause the problem to become unstable…
> To help with choosing a penalty factor, FEBio offers the **`auto_penalty`** option which
> will estimate a good penalty factor based on element size and initial material stiffness.」

### FEBio 手册 Control `symmetric_stiffness`（note 6）

> 「**0 = unsymmetric, 1 = symmetric, 2 = structurally symmetric, 3 = preferred**…
> If the option is set to **3 (preferred)**, then one of the other formats will be chosen
> depending on the model components… if a **non-symmetric contact formulation** is used,
> then a **non-symmetric** stiffness matrix will be used.」**默认 1。**

### FEBio Studio 官方 webinar —— Contact Modeling（逐字要点）

> 「the top bone has a force supply and is otherwise unconstrained so you might think that the
> contact is going to catch that top bone after you apply that force but in reality that's not
> what happens… **if there's initial gap separating the two contacting interfaces then without
> any other constraints to the top bone this problem is ill-posed**」
> 「**the first displacement step runs, engages the contacts, and then the force is being
> applied**… the second step starts from that value and simply ramps up the force」

### FEBio 运行日志里的主动警告（我们以前没处理）

```
 *                               WARNING                                 *
 * The sliding-elastic contact algorithm runs better with a              *
 * non-symmetric stiffness matrix.                                       *
 * You may set symmetric_stiffness flag to 0 in Control section.         *
 *************************************************************************
```

---

## §2 · ★ 本轮最实在的收获：**我们自己的护栏挡掉了手册的 #1 处方**

`src/climbing/coupling/plantar_bc.py:436`（**改前**）：

```python
if not (np.isfinite(gap_mm) and gap_mm > 0.0):
    raise ValueError(f"gap_mm 必须有限且 > 0（接触搜索需严格正间隙），得到 {gap_mm!r}")
```

而墙的位置是 `wall_y = min(plantar_y) − gap_mm`（`:280`）⇒
**`gap_mm > 0` 被硬性要求 = 强制初始分离** —— **正是手册 §8.12.4/§8.5.3 与官方 webinar
点名的 ill-posed 条件。项目里 `contact` 一直跑 `gap_mm = 0.5`。**

**⇒ 处置（已做）**：放宽为「只要有限」+ 拒绝「过盈超过薄板厚度一半」（板会翻转）。
**默认值 `DEFAULT_CONTACT_GAP_MM = 0.5` 不动 ⇒ 原路径逐位不变。**
**零回归已验**：`pytest -k "plantar or contact or thums or bone_feb"` ⇒ **72 passed**。

---

## §3 · 实测矩阵（h=5 m，F_sub = 26383.1 N，`hold_spring_k=None` = **无旋钮**）

### 3.1 处方 A：初始接触 / 初始过盈（STATIC）

| gap (mm) | +0.5（原） | 0.0 | −0.05 | −0.2 | −0.5 |
|---|---|---|---|---|---|
| 结果 | ❌ ERR (1 s) | ❌ ERR (1 s) | ❌ ERR (1 s) | ❌ ERR (1 s) | ❌ ERR (1 s) |

**失败模式**：`8 negative jacobians detected`，在**第一次增量 `time = 1` 就发生**，
`auto_stepper` 把 dt 降到 0.833/0.667/0.5/0.333 **仍全部同款失败** ⇒ **不是载荷步长问题，是结构性的。**

**⇒ 裁定：单靠「消除初始分离」不足以收敛。**

### 3.2 处方 D：`symmetric_stiffness`

| 变体 | 结果 |
|---|---|
| `base`（无标签，Control 默认 1=对称） | ❌ ERR |
| `<solver><symmetric_stiffness>0` | ❌ ERR（**解析通过**，仍失败） |
| `<Control><symmetric_stiffness>3` | ❌ **解析失败**（我的文本插入破坏了结构；**未定论**） |

*注：`<Contact>` 块内部本来就有 `symmetric_stiffness=0`（`_apply_contact` 传的）。*

### 3.3 处方 C：False Transient —— 🔴 **更正：不是「真积分」，是「刚体飞走」**

**本件 v1 曾把这一条读成「行为质变」（秒退 → 真积分到 t≈0.44）。该读数作废。**

**实际发生的事**（读最后收敛态的位移场，`fe_post.read_displacement`，脚本 `temp/v1_freeflight.py`）：

| 量（`g-0.05_e0.185_ap1` 最后收敛态） | 值 |
|---|---:|
| 骨节点位移 \|u\| **中位** | **7295.78 mm（7.3 米）** |
| \|u\| max / min | 7314.27 / 7279.74 mm |
| **形心位移** | **(4297.34, −3725.91, 4570.37) mm** |
| **各节点相对形心的偏差（中位）** | **27.1 mm** |
| **偏差 / 形心位移** | **0.0075** |
| 板（对偶薄板）位移 | **0.000 mm** ✓（自检：板底固定） |

**⇒ 位移场就是「平移 + 少量刚体转动」——相对形心偏差只占 0.75%。骨飞了 7.3 米。**

**⇒ 这正是手册 §8.12.4 逐字预言的失败模式：**

> 「Even the slightest deviation from balancing the forces exactly will create a net-force
> which in turn causes the model to **fly off in the unconstrained direction**
> (in which a **rigid body mode** now exists).」

**⇒ 与该读数一致的自洽性检查**：网格 SPON 体积 **77794 mm³**；deck `<density>` =
骨 **1e-06**、板 **1e-09**（**两者差 1000×**）。若按 FEBio 一致的 mm-MPa-**tonne**-s 体系，
`m = 1e-6 × 77794 = 0.0778 tonne = 77.8 kg`（跟骨真实约 100 g ⇒ **密度量级可疑**）；
F(0.44·26.4 kN)=11609 N ⇒ `a ≈ 149 m/s²` ⇒ 自由飞行位移 `½at² = 0.44² 下 ≈ 14 m`
—— **与实测 7.3 m 同量级** ⇒ **自由飞行假设成立**。

**⇒ 所以「t≈0.44 才失败」不是「算到了 44%」，而是「在 44% 的伪时间之前它一直在飞」。**

**⚠️ 附带更正**：本件 v1 里「~500 个单元倒置」也是同一现象的副产品（与板穿插后的畸变），
**不是「大变形单元倒置」这一独立机制**。
**⚠️ 我算的 `det F` 中位数 0.31 是错的** —— xplt 的节点顺序与网格顺序**未验证一致**
（hdf5 只有 `meshes`/`states`，无 `node_ids`）⇒ **该 det F 数值作废**（纯刚体运动应给 det F = 1）。

**⇒ 裁定：处方 C 的「进展」是假象。DYNAMIC 没有修正欠约束，只是把「秒退」换成了「飞走」。**

#### 🔑 两条仍然成立的语法知识（本轮确定，可复用）

1. **`rhoi` / `gamma` / `beta` 属 `<solver>`（固体求解器）**，放 `<Control>` 会被 FEBio 拒绝
   （`tag "rhoi" : unrecognized`）—— 手册 §8.4.1 的措辞「**solid solver**」是准确的。
2. **`pyfebio` 把 `symmetric_stiffness` 定在求解器类**（`control.py:67`，默认 `"preferred"`），
   **不在 `<Control>`** ⇒ 直接往 `<Control>` 里插会破坏解析（实测「解析失败」）。
3. **`det F` 读法陷阱**：xplt 的节点顺序**未必等于**网格顺序（本件 hdf5 只有
   `meshes`/`states`，无 `node_ids`）⇒ **不校验顺序就算 det F 会得到假值**
   （纯刚体运动应给 `det F = 1`，我实测 0.31）。

#### 组合尝试（DYNAMIC+阻尼 × 板模量 × gap × auto_penalty）—— 全部失败

| 组合 | 走到 t |
|---|---|
| gap=−0.05 / 刚性板 / ap=0 | 0.399 |
| gap=−0.05 / **跟垫 0.185** / ap=0 | 0.447 |
| gap=+0.5 / 跟垫 0.185 / ap=0 | 0.495 |
| gap=−0.05 / 刚性板 / **ap=1** | 0.564 |
| **gap=−0.05 / 跟垫 0.185 / ap=1** | **0.613** |
| gap=+0.5 / 跟垫 0.185 / ap=1 | 0.443 |

**⇒ 组合确实把「伪时间」推得更远（0.399 → 0.613），但既然「伪时间」本身就是自由飞行的
时间刻度，这个数字**不构成进展**。**

#### 处方 E：`auto_penalty=1`（本轮首次，手册 §8.5.1）

**❌ 无效**（见上表：ap=1 与 ap=0 无质变）。
（接线已加：`contact_auto_penalty`，opt-in，默认 0 ⇒ 原路径不变。）

### 3.4 🔴 **本件 v1 的总体结论也要降级**

v1 写：「手册处方**部分有效**（C 把「秒退」变成「真积分」），但**都还没到收敛**。」

**⇒ 改述**：**手册列出的五条处方（A 初始接触 / B 两步法 / C False Transient /
D symmetric_stiffness / E auto_penalty），本轮全部试过，无一条使跖面接触（在 `hold_spring_k=None` 下）收敛。**
**唯一让它收敛的仍是 `hold_spring_k`（自由旋钮 ⇒ 只能当诊断口径）。**
**⇒ 「跖面接触在无额外约束下不收敛」是一个**未被手册处方解决**的问题。**

---

## §4 · 下一步（可操作，按优先级）

| # | 动作 | 依据 | 预期 |
|---|---|---|---|
| **0** | 🔴 **先修「刚体飞走」**：DYNAMIC 下必须先有**有效支撑**才能谈动力学 —— 即**回到处方 A（t=0 预接触）+ 足够的接触刚度**，否则任何 DYNAMIC 跑都是自由飞行 | §8.12.4 | 不修这条，后面都是白跑 |
| **1** | **`<Control>` 级 `symmetric_stiffness` 用 pyfebio 正确接口**（不是文本插）重试 | note 6 | 未定论项；FEBio 主动警告过 |
| **2** | **接触刚度：`penalty` 从 0.1 提上去 / `auto_penalty=1`** —— 已试 ap=1，但**未试「硬 penalty」**；0.1 是否过软未判 | §8.5.1 | 过软会让跖面深陷 |
| **3** | **`contact_type` 换 `tied`/`sticky` 首步**：先把跖面在 t=0 **tied 住**（完全消除刚体模态），再在后续步释放 | §3.14.9 `sticky`；§8.3.1.5「加额外约束」 | 用**物理界面**替代 `hold_spring_k` 旋钮 |
| **4** | **两步法在 DYNAMIC 下重试**（既有 Step 2 失败是 STATIC 下测的） | §8.5.3 + webinar | — |
| **5** | ⚠️ **先查 `density` 一致性**：骨 `1e-06` vs 板 `1e-09`（**差 1000×**）—— 单位体系若为 tonne-mm-s，骨密度高 1000× | — | 影响一切 DYNAMIC 结果 |
| **6** | 补 **`rigid` BC**（FEBio 的 `<rigid>`）给骨的**最小约束**（例如只锁转动），替代弹簧 | — | 手册 §8.3.1.5「eliminate by applying additional constraints」的**非旋钮**版本 |

---

## §5 · 诚实边界

1. **本件全部是 h=5 m / `a=0` / anat / `PressureLoad` / `gauge_max` 口径下的诊断跑**，
   **无一条构成可交付口径**。
2. **`DYNAMIC` 的 460–520 个倒置单元未定位** ⇒ **不能**说「是网格问题」或「是大变形问题」。
3. **`gap ≤ 0` 引入初始过盈**，其合理性**未论证**（真实跟垫是软的，不会过盈）；
   它只是**执行手册处方**的手段。
4. **`hold_spring_k` 仍是最容易收敛的路**（既有 `contact_hold100` 全高度收敛），
   但它是**自由旋钮** ⇒ 只能当诊断口径（`V1_PLANTAR_LOAD_BEARING_FEASIBILITY.md` §1）。
5. **本件不动任何判定**：V1 仍为 ❌【口径不可复现】。

---

## §6 · 溯源

| 项 | 位置 |
|---|---|
| 护栏放宽（本件，零回归 72 passed） | `src/climbing/coupling/plantar_bc.py`（原 `:436-437`） |
| gap 扫描（处方 A） | `temp/v1_gap_scan.py`；产物 `temp/opensim_fe/v1_gap/`；`results/opensim_fe/v1_plantar_gap_scan.json` |
| symmetric_stiffness（处方 D） | `temp/v1_sym_scan.py`；`results/opensim_fe/v1_symmetric_stiffness_scan.json` |
| False Transient（处方 C） | `temp/v1_dyn_scan.py`、`temp/v1_dyn2.py`；`results/opensim_fe/v1_dynamic_scan.json`、`v1_dyn_damp_scan.json` |
| 手册 | `temp/thums/febio_manual.txt`（§8.3.1.5 / §8.4.1–8.4.2 / §8.5.1–8.5.3 / §8.12.4；note 6 在 L1932-1942） |
| 官方 webinar | FEBio Studio「Contact Modeling」（`febio.org` knowledgebase → webinars） |
| 既有 `contact` 不收敛根因分析 | `results/opensim_fe/PLANTAR_CONTACT_FIX_20261015.md` §1–§3 |
| 跖面承载可行性三层裁定 | `results/opensim_fe/V1_PLANTAR_LOAD_BEARING_FEASIBILITY.md` |
| 铁律 | `[MEMORY gotcha #26]`；`docs/notes/FE_PRESOLVE_GATES.md` §0 |
