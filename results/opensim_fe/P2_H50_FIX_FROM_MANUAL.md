# P2 · T6 `h=50` —— **手册 + 网页检索 → 修法（能力门 3/4 → 4/4）**

> **日期**：2026-10-15（深夜）｜ **性质**：**修法证据件**（不自行改判；采纳须用户决定）
> **触发**：用户令「从网页和手册找办法」—— 推翻本会话前件 `P2_H50_RESIDUAL.md` §6 的
> 「模型承载能力上限」**推理结论**（该结论**错**：真正缺口在**载荷施加方式**）。
> **能力门现状**：**0.20 冻结值 ⇒ 3/4**；**0.35 / 0.50 ⇒ 4/4（rc=0、零负 Jacobian）**。

---

## 0. 一句话

**`h=50` 不收敛的真因＝52.6 kN 被压成 59 个 891.8 N 的离散点力**；
**把载荷铺到更多节点（每节点力 892 → 191 / 120 N）即 4/4 收敛。**
**这不是「模型承载上限」，是「点载荷过于集中」——可由手册与网页双双指认，并被单变量实验钉死。**

---

## 1. 手册检索（`temp/thums/febio_manual.txt`，18593 行）

### 1.1 第 8 章 Troubleshooting 全文要点

| 手册节 | 原文要点（要点摘录） | 本轮是否已试 |
|---|---|---|
| **§8.1.1** | 「when elements are **too thick** (thickness **same order as the element size**)… may be particularly a problem in **areas of high curvature**. The only solution… might be to **remesh**」 | ✅ 已按「壳下实体单元高」量测（比 0.506）⇒ 不成立 |
| **§8.1.1** | 「If FEBio discovers negative jacobians **before the first time step**… the nodes of the elements are **not defined in the proper order**」 | ✅ 我们不是首步前 |
| **§8.1.3** | 「**over-constrained**… conflicting boundary conditions. **a force is applied to a node that is fixed**」 | ⚠️ 未查（见 §5 待办） |
| **§8.2** | 调试：`febio -i file.feb -g` / `-g2`（只输出导致失败的那一个非收敛状态） | ✅ 改用 `-output_negative_jacobians` 达成 |
| **§8.3.1** | 「set `output_negative_jacobians` … to print out a **detailed list of the inverted elements**」 | ✅ **本轮据此定位** |
| **§8.3.1.1** | 材料不稳定（大变形软化）⇒「no easy solution… consider a **different constitutive model**」 | —（本构为线弹性，非软化） |
| **§8.3.1.2** | 「The most common issue is that the **minimum time step is set too large**. **Cutting back the minimum** in addition to the initial time step can often help」 | ✅ 试了 `dtmin=1e-8`+50 retries ⇒ 8 步 / t=0.83，**仍不收敛** |
| **§8.3.1.3** | 「When elements get **too distorted**… adjust the **mesh**」 | ⚠️ 见 §5 |
| **§8.3.1.4** | 壳过厚 ⇒「remesh / **adjust the shell thickness** / **replace the shells with solid elements**」 | ✅ 已试（薄壳更差） |
| **§8.3.1.5** | 刚体模态 ⇒ 加约束 | ✅ 已由注册修正解决 |
| **§8.3.2.1** | 无载荷时收敛困难 ⇒ `min_residual` | —（我们有大载荷） |
| **§8.4.2** | **False Transient**：「Quasi-static models may sometimes exhibit **instabilities** that may prevent the default Newton-based solver from converging. In such situations, **a dynamic analysis can be tried** instead, with the anticipation that **the inertia of the model provides sufficient stability**」；最阻尼 `rhoi=-2, gamma=1.5, beta=1` | ✅ **已试 ⇒ 无效**（见 §3.6.2） |
| **§8.5.1** | **`auto_penalty`**：「will **estimate a good penalty factor based on element size and initial material stiffness**」；「choosing the penalty factor **too large** may cause the problem to become **unstable** since the contact force (and corresponding stiffness) will **dominate**」 | ✅ 试了 ⇒ **更差**（0 步） |
| **§8.5.2** | **Augmented Lagrangian**：「often gives **better enforcement** of the contact constraint compared to the penalty method」 | ✅ 试了 ⇒ 22 步但只到 t=0.20 |
| **§8.5.3** | 初始分离：力驱动问题最难；「first use a **displacement** to bring the surfaces in contact then **replace the displacement with a force**」 | ⚠️ 我们已是「位移驱动压板 + 力驱动骨」混合 |

### 1.2 命令行 / 配置项（§2.3、§2.5、§9.1）

- **`-output_negative_jacobians[=n]`**：打印倒置单元清单（**本轮的核心工具**）
- **`-g` / `-g2`**：调试模式，把非收敛状态写入 plot 文件
- 配置项（`set <name> <0|1>`）：`output_negative_jacobians` / `print_model_params` / `show_warnings_and_errors`
- **`§8.3.1.2` 的 `dtmin`** 与 **`§8.3.1.4` 的壳厚** 均已试

### 1.3 手册的一条**未能满足的前提**

手册 §8.1.1 的「too thick（same order as the element size）」**没有给量化阈值**。
本轮据此量到 T6 **壳 1.5 / 壳下实体高 2.965 = 0.506**（**确属「same order」**），
故**不能仅凭此排除**该条 —— 但**薄壳实测更差**（2 步），**排除依据是实验而非比值**。

---

## 2. 网页检索

| 来源 | 要点 |
|---|---|
| **audilab（FEBio 常见问题）** | 「This may mean that your **deformations have become so large that some elements have become too distorted**… **You could try reducing the pressure or force loads.**」——即**载荷量级/施加方式**是首查项 |
| **GIBBON issue #51（FEBio 实践）** | ① 「Negative Jacobians may also occur if we are **forcing FEBio to take a big leap** to the next time increment」 ② 罚因子「**I tend to start with a contact penalty parameter that is too low**… starting at 0.1, then 1, then 10」 ③ 「if the element is initially of **bad quality**… we need to look at **improving the mesh quality**」 |
| **FEBio 官方调试教程（YouTube）** | 「**set output_negative_jacobians 1**… in the log file it will tell you **which elements are inverted at every time step**」；★「if I have a **high penalty** then… it's going to push really hard and **that might invert my elements**… **simply editing the contact parameters can really have a huge impact**」；**element inspector 查 Jacobian 0~0.0001 的单元** |
| **FEBio Studio Webinar: Contact Modeling** | ★「**auto penalty flag**… takes the **stiffness** into account, takes the **mesh size** into account… with that flag on you can still set the penalty factor yourself but in that case **the penalty factor becomes a scale factor**」；★ 负 Jacobian 时推荐 **swap primary/secondary surface**；「turn on the **augmentations**… two or three」 |
| **FEBio Studio Webinar: Contact Modeling Part 2** | 「**lower the penalty factor** or… **swap the primary and secondary surface**」 |

**⇒ 网页的一致指向：负 Jacobian 优先怀疑 (a) 接触参数（罚），(b) 载荷量级/施加方式，(c) 网格质量。**

---

## 3. 本轮实测（全部单变量，`ramp0` + 注册修正）

### 3.1 ★ 决定性诊断：**倒置单元的空间位置**

用 `-output_negative_jacobians=200` 取得清单，再映射到网格域 + 沿轴位置
（`temp/fe_out/negjac_probe.py` / `negjac_map.py` / `negjac_spatial.py`）：

**第一次大失败（首步 t=0.1）归属**：

| 域 | 类型 | 唯一倒置单元 | 占比 |
|---|---|---:|---:|
| **t6_spon** | tet4（实体松质） | **247** | **70 %** |
| t6_cort | quad4（皮质壳） | 80 | 23 % |
| **plate** | hex8（E=200 GPa 压板） | 24 | 7 % |

**`t6_spon` 倒置单元沿轴分布**（f：0 = 载荷端，1 = cap 接触端；**载荷带 = f ∈ [0, 0.197]**）：

| f 区间 | 倒置数 |
|---|---:|
| 0.0–0.1 | 1 |
| 0.1–**0.2**（载荷带边缘 **0.197**） | 26 |
| **0.2–0.3** | **166** |
| **0.3–0.4** | **152** |
| 0.4–1.0 | 7 |

**⇒ 90 % 挤在载荷带边缘的**紧邻上侧**（f∈[0.2,0.4)）；且径向偏**内侧**（中位 7.22 vs 全体 9.58）。**

### 3.2 ★ 机制钉死：**每节点集中力的大小**（而非边缘突变）

| 试验 | 节点数 | 每节点力 | 步数 | 收敛 |
|---|---:|---:|---:|---|
| 对照（LF=0.20，等值） | 59 | **891.8 N** | 4 | ❌（275 negJ） |
| **平滑衰减**（内 60 % 满值、边缘线性→0，ΣF 守恒） | 6 带 | 最大 **1393 N** | **3** | ❌（**更差**） |
| **平滑衰减**（内 30 % 满值） | 6 带 | 最大 **1826 N** | **3** | ❌（**更差**） |
| **加宽 LF=0.35** | **276** | **190.7 N** | 4 | ✅ **rc=0，零 negJ** |
| **加宽 LF=0.50** | **440** | **119.6 N** | 4 | ✅ **rc=0，零 negJ** |

**⇒ 平滑化把总力归一化到更少节点 ⇒ 每节点力升到 1393/1826 N ⇒ 更糟。**
**⇒ 结论：`h=50` 的失效由「离散点力过于集中」驱动；点力 ↓ ⇒ 收敛。**

**闸门 3（patch 生效验证）**：逐次打印载荷端节点数（59 → 276 → 440）与逐带每节点力，
并校验 `ΣF = 52 618.62 N`（目标一致）✅

### 3.3 能力门全量（`ramp0`，四档高度 × 三种 LOAD_FRAC）

| LOAD_FRAC | h=10 | h=20 | h=30 | h=50 | **收敛** |
|---|---|---|---|---|---|
| **0.20（冻结）** | ✅ | ✅ | ✅ | ❌ rc=1 | **3/4** |
| **0.35** | ✅ | ✅ | ✅ | ✅ rc=0 | **4/4** |
| **0.50** | ✅ | ✅ | ✅ | ✅ rc=0 | **4/4** |

`bear_frac`：4/4 各档 **1.00000–1.00002**（冻结值时 h=50 为 **0.77778**）。

### 3.4 σ 后果（**这是本次修法的代价 —— 必须与修法并列引用**）

`σ_median`（皮质域，MPa；**σ_c = 150**）：

| h | LF=0.20 | LF=0.35 | LF=0.50 |
|---:|---:|---:|---:|
| 10 | 209.626 | 168.670 | **141.661** ⚠️ |
| 20 | 319.769 | 261.238 | 217.278 |
| 30 | 416.246 | 340.174 | 280.840 |
| 50 | 450.995 | 482.725 | 392.282 |

`σ_max`：3848.5 → 1193.5 → 1157.9（**载荷带边缘伪影被摊开**，与 `P2_V2_…` §2.2 同向）

**⇒ 影响预估**：
- **LF=0.35**：h=10 仍 >150 ⇒ **首折仍成立**，首折高度由 **6 m** 移至 **≈8–9 m** ⇒ T6 仍为**阶跃** ⇒ **A2 预计仍 5/5**（须重跑确认）
- **LF=0.50**：h=10 σ=141.661 **< 150** ⇒ 首折移至 **≈10–11 m** ⇒ **仍为阶跃**（h=50 σ=392>150）⇒ A2 预计仍 5/5
- ⚠️ **A1 的 #2 要求 T6 首折 ∈[16,18] m** ⇒ 两个方向都**更远**，对 A1 **不利**（A1 已定性【定义差异】，不因此改变判定）

### 3.5 已排除的候选（单变量实跑）

| 候选（依据） | 结果 | 判词 |
|---|---|---|
| `auto_penalty=1`（§8.5.1 / webinar） | **0 步**，negJ 1323 | ❌ 更差（auto 估的罚**过大** ⇒ 印证 §8.5.1「too large 会不稳定」） |
| `auto_penalty=1` + 删 `<penalty>` | **0 步**，negJ 1937 | ❌ 更差 |
| `laugon=AUGLAG`（§8.5.2） | 22 步、**只到 t=0.204** | ❌ 步更细但推进更浅 |
| `two_pass=1`（webinar） | 4 步、到过 t=1.0 的 attempt | ❌ 仍失败 |
| `dtmin=1e-8`+50 retries+cutback 0.25（§8.3.1.2） | 8 步、**t=0.8305**（最远） | ❌ 仍失败 |
| 壳厚 1.5→0.5 mm（§8.3.1.4） | 2 步 | ❌ 更差 |
| `LOAD_FRAC` 平滑衰减（本轮自创） | 3 步 | ❌ 更差（每节点力反升） |
| step 级 `<time_stepper>`（承前） | 与基线逐位同 | ❌ |
| **false transient**（§8.4.2，`DYNAMIC`+`solver rhoi=-2`） | 4 步、**275 negJ（与静态逐位同）** | ❌ |

### 3.6 两个新增工具坑（已验证）

1. **FEBio 4.13 拒绝 step 级 `<time_stepper type="default">`**（`invalid attribute "type"`）；
   **去掉 `type` 属性即被接受**（旧记录「拒绝裸 `<time_stepper>`」的真因）。
2. **`rhoi` / `gamma` / `beta` 属于 `<solver>` 块，不属于 `<Control>`** ——
   手册 §3.3.5（L1985-1989）把它们列为 **solver 参数**。
   写在 `<Control>` 下报 `tag "rhoi" (line 5827) : unrecognized tag`。
   正确写法：`<Control><solver><rhoi>-2</rhoi><gamma>1.5</gamma><beta>1</beta></solver></Control>`

3. **false transient 实跑 ⇒ 无效**（`temp/fe_out/h50_dyn.py`）：

   | 配置 | rc | 步数 | 覆盖 t | negJ |
   |---|---:|---:|---|---|
   | `DYNAMIC` + `step_size=0.1`（覆盖 t∈[0,1]）+ `solver rhoi=-2, gamma=1.5, beta=1` | 1 | 4 | 至 1.0 | **275**（与静态基线同） |
   | `DYNAMIC` + `step_size=0.1` + `solver rhoi=0`（最阻尼） | 1 | 4 | 至 1.0 | **275** |

   **⇒ 惯性未提供手册预期的稳定化；失败签名与静态逐位相同。**

   ⚠️ **本轮自错（假阳性）**：另有一组 `step_size=0.02` × `time_steps=10` 得 `rc=0`——
   **但那只覆盖 t∈[0,0.2]（= 仅 20 % 载荷）**，`rc=0` **无意义**，**不得**记为「false transient 有效」。
   （教训：`rc=0` 必须与「覆盖的 t 范围」一起看 —— 即闸门 3 的延伸。）

---

## 4. 判定与建议

### 4.1 事实

- **能力门可 4/4**：把载荷从 59 个 891.8 N 点力改为 276（190.7 N）或 440（119.6 N）个即收敛。
- **机制清楚**：倒置单元 90 % 在载荷带边缘紧邻区；点力↓⇒收敛（平滑化反升点力⇒更差）。
- **代价**：T6 σ 变 **−20 % ~ −33 %**（衰减到 LF=0.50 时 h=10 跌破 σ_c）⇒ **首折高度移动**。

### 4.2 建议（**须用户拍板**）

| 选项 | 内容 | 后果 |
|---|---|---|
| **(A) 不动** | 能力门维持 **3/4**，h=50 记为「冻结载荷施加方式下的上限」 | A2 5/5 不变、所有已入账数字不动、零重跑 |
| **(B) 采纳 LF=0.35** | 能力门 **4/4**；T6 σ 及首折高度重算 | A2 预计仍 5/5（**须重跑官方管线确认**）；A1 不变 |
| **(C) 改用端面分布压强** | 物理上最正当（椎体载荷经终板） | ⚠️ 实测端面仅 18 面片 / 75.80 mm² ⇒ **须先解决端面提取**（`axdiv` 阈值） |
| **(D) 先补预注册** | 任一采纳前，把「载荷带宽度是为股骨长骨而定」的**适用性缺口**写成预注册 | 防「调参凑判据」质疑（**推荐先做**） |

### 4.3 为何 (B) 不是「调参凑判据」

1. **有机制**：倒置位置（90 % 在带边缘紧邻区）＋点力因果链（平滑化反升点力⇒更差）⇒ **非阈值巧合**
2. **项目早已记录**：`P2_V2_PEN_SIGMA_ROOTCAUSE.md` §2.2 已实测 **LOAD_FRAC 0.20→0.50 ⇒ σ_max −66 %**，
   即 **0.20 这个值本身就是已知的伪影来源**
3. **0.20 是股骨长骨的约定**（`p2` 继承自 P5/股骨线），**未经椎体重新推导** ⇒ 缺口在「未重新推导」，不在「改它」
4. **手册 + 网页同向**：§8.3.1.3 / §8.1.1、audilab「try reducing the pressure or force loads」

---

## 5. 未完成 / 诚实边界

1. **false transient（§8.4.2）已跑通但无效** —— XML 形态已修正（`rhoi` 进 `<solver>`，见 §3.6.2），**实跑与静态基线逐位同（275 negJ）** ⇒ 该路在本问题上不成立。**未试**：更小的 `step_size` + 更多步（会用更长物理时间），以及 `rtol`/`dtol` 调整
2. **端面分布压强（选项 C）未打通** —— 端面面片提取口径未解决
3. **§8.1.3 过约束**（「force applied to a node that is fixed」）**未查**：
   压板 24 个 hex8 倒置，**未排除**「压板被刚体 ghost 与骨面夹住」的可能
4. **A2 官方管线未重跑** ⇒ 「LF=0.35/0.50 下 A2 仍 5/5」是**预估而非实测**
5. **本会话前件 `P2_H50_RESIDUAL.md` §6「模型承载能力上限」结论已错**，须替换式更正（§6）
6. **本节修法未采纳**：能力门**仍记 3/4**，未改任何冻结输入

---

## 6. 对前件的更正（**替换式**）

| 前件论断 | 更正 |
|---|---|
| `P2_H50_RESIDUAL.md` §6「残余＝**模型承载能力极限**（应变≫100 %、E=40 MPa、3186 tet）」 | **错**。真因在**载荷施加方式**（59×892 N 点力过于集中）；`h=50` **可收敛**（LF=0.35/0.50 ⇒ rc=0、零 negJ） |
| 同件 §0「**不是**求解器/配方缺口」 | 方向对（**不是配方**），但**也不是模型上限** ⇒ 是**载荷施加方式** |
| 同件 §4「薄壳更差 ⇒ 壳不是可修复变量」 | **保留**（结论仍成立） |
| 同件 §3「罚 ×10/×100/×1000 无效」 | **保留**（结论仍成立） |
| 本会话更早的 `P2_H50_RESIDUAL` §0 结论 ⑤「能力门维持 3/4（建议）」 | 改为：**维持 3/4 是「不动模型」下的现状**；**4/4 可达**，采纳与否待拍板 |

---

## 7. 证据索引

| 项 | 位置 |
|---|---|
| **倒置单元清单（决定性）** | `temp/fe_out/negjac_probe.py` → `negjac_D1_h50_ramp0.txt` |
| **域归属映射** | `temp/fe_out/negjac_map.py` / `negjac_attempts.py` |
| **空间分布** | `temp/fe_out/negjac_spatial.py` |
| **能力门 × LOAD_FRAC** | `temp/fe_out/h50_gate_lf.py` → `h50_gate_lf.json` |
| **平滑衰减（σ 保全试验）** | `temp/fe_out/h50_taper.py` → `h50_taper.json` |
| **手册修法扫描** | `temp/fe_out/h50_manual_fixes.py` → `h50_fixes.json` |
| **壳/实体几何** | `scripts/opensim_fe/shell_solid_ratio_audit.py`（闸门 7） |
| 手册 | `temp/thums/febio_manual.txt`（§8 / §2.3 / §2.5 / §9.1） |
| 网页 | audilab FEBio FAQ、GIBBON issue #51、FEBio 官方调试教程、FEBio Studio Contact Modeling webinar ×2 |
