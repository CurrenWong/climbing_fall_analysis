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
| **C** | **False Transient**：准静态不收敛 ⇒ 改**动力**分析 + 数值阻尼（§8.4.2） | ✅ **本轮首次** | 🔵 **行为质变**：从「1 s 秒退、8 negJ」→「**真积分到 t≈0.44、460~520 单元倒置**」⇒ **未收敛，但瓶颈转移** |
| **D** | **`symmetric_stiffness` = 非对称 / preferred**（FEBio 主动 WARNING + 手册 note 6） | ✅ **本轮首次** | ❌ 放 `<solver>` 仍失败；放 `<Control>` **被 FEBio 拒绝**（位置不对） |
| **E** | `auto_penalty`（§8.5.1） | ❌ 未试 | 待试 |

**⇒ 一句话：手册处方**部分有效**（C 把「秒退」变成「真积分」），但**都还没到收敛**；
而本轮最实在的收获是**找回了 A 这条被我们自己护栏挡掉的路**，并**把瓶颈从「欠约束」转移到了「大变形单元倒置」**。**

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

### 3.3 处方 C：False Transient（**行为质变**）

| 变体 | 结果 |
|---|---|
| `DYNAMIC`（无阻尼） | ❌ ERR（5–6 s，比 STATIC 走得远） |
| `DYNAMIC + rhoi=-2/gamma=1.5/beta=1`（放 `<Control>`） | ❌ **解析失败**：`tag "rhoi" : unrecognized` |
| **`DYNAMIC + rhoi=-2/gamma=1.5/beta=1`（放 `<solver>`）** | ❌ ERR，**但跑了 120–125 s**，积分到 **t ≈ 0.443–0.446** 才失败，且失败模式变成 **460~520 negative jacobians** |

**⇒ 这是本轮最重要的读数**：**失败模式从「秒退 + 8 个 negJ（刚体模态）」变成
「真积分到 44% 载荷 + ~500 个单元倒置（大变形）」** ⇒
**瓶颈已从「欠约束」转移到「大变形下单元倒置」。**

**⇒ 两条有用的语法知识（本轮确定）**：
1. **`rhoi` / `gamma` / `beta` 属 `<solver>`（固体求解器），放 `<Control>` 会被 FEBio 拒绝**
   （`tag "rhoi" : unrecognized`）—— 手册 §8.4.1 的措辞「solid solver」是准确的。
2. `pyfebio` 把 `symmetric_stiffness` 定在**求解器类**（`control.py:67` 默认 `"preferred"`），
   **不是 `<Control>`** ⇒ 直接往 `<Control>` 里插会破坏解析。

---

## §4 · 下一步（可操作，按优先级）

| # | 动作 | 依据 | 预期 |
|---|---|---|---|
| **1** | **把 DYNAMIC + `<solver>` 阻尼配方正式接进 builder**（opt-in，默认不变），并按 §8.4.2 用 `dt` 分步把载荷推完 | §8.4.2 | 已证明能真积分；剩下的问题是单元倒置，不是欠约束 |
| **2** | **`auto_penalty=1`**（我们一直用 `penalty=0.1, auto_penalty=0`） | §8.5.1 | 接触刚度过软可能放大大变形 |
| **3** | **把「单元倒置」定位**：~500 个倒置单元在哪（骨/板）？本件**未查** | 手册 §8.3.1.1–§8.3.1.3（材料失稳/畸变） | 决定是网格问题还是大变形问题 |
| **4** | 两步法（处方 B）在 **DYNAMIC** 下重试 | §8.5.3 + webinar | 既有 Step 2 失败是 STATIC 下测的 |
| **5** | `<Control>` 级 `symmetric_stiffness` **用 pyfebio 正确接口**（不是文本插）重试 | note 6 | 未定论项 |

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
