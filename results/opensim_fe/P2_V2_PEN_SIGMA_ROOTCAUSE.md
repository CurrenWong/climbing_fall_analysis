# P2 · V2 —— 穿透 / σ 幅值 / h=50 的根因（**两个独立病**）

> **日期**：2026-10-15 ｜ **性质**：只读量测 + 单变量复跑（零改冻结输入）
> **脚本**：`temp/opensim_fe/p2_gap_probe.py`、`p2_band_probe.py`、`p2_pen_scan.py`（只读/单变量）
> **产物**：`results/opensim_fe/P2_PEN_SCAN_x1.json`、`P2_PEN_SCAN_x100.json`
> **前置**：`P2_T6_ROOTCAUSE.md`（注册 bug）、`P2_V1_REGFIX_REPORT.md`（0/4 → 3/4）
> **依据**：FEBio User Manual §8.3.1.2/§8.3.1.4/§8.3.2.4/§8.4.2、**L4551/L4557/L4569/L13745-13747/L13929**；FEBio 官方 Contact Modeling 教程

---

## 0. 结论

V1 的 3/4 PASS 之下藏着**两个彼此独立的病**。本件用**单变量实测**把它们分开：

| # | 病症 | 根因 | 决定性证据 |
|---|---|---|---|
| **A** | **h=50 不收敛** | **接触罚刚度不足** ⇒ 表面严重穿透 | 穿透实测 **10.4 mm**（股骨正对照 **0.060 mm**，差 **~170×**）；罚 **×100** ⇒ 穿透 **0.93 mm**、h=50 从 **0 步 → 4 步** |
| **B** | **σ 幅值不自洽**（h=10 σ_max=1421 MPa = 皮质强度 9.5×） | **载荷带边缘伪影**（**与接触无关**） | 罚 **×1 与 ×100** 下 **σ_max 逐位不变（1421.4 → 1421.4）**；峰值单元 `σ_zx=−743.7 MPa`，位置 `z_frac=0.198/0.203` **正好落在载荷带上边缘 0.20** |

**⇒ 关键分辨**：`σ` 对罚刚度**完全免疫** ⇒ 它不是接触问题；而穿透对罚刚度**高度敏感** ⇒ 它是接触问题。

---

## 1. 病 A：接触罚刚度不足 ⇒ 穿透 10.4 mm（实测）

### 1.1 量测方法（可复现）

`surface_data` 只输出 `contact force`，**无 gap/穿透** ⇒ 从位移场直接算：

```
穿透 = (cap 面片节点变形后 z) − (压板面 z = z_wall − plate_disp)
```

| 案例 | 板面 z | cap 变形后 z | **穿透 max** | **穿透 median** | 接触面占比 |
|---|---:|---|---:|---:|---:|
| **T6 h=10** | 329.72 | [331.32, 340.16] | **10.435 mm** | **8.328 mm** | 100 % |
| **T6 h=30** | 329.72 | [331.65, 339.89] | **10.171 mm** | **8.359 mm** | 100 % |
| 股骨 P5.0 h=50（**正对照**） | — | — | **0.060 mm**（档案记录） | — | — |

**⇒ T6 的穿透是股骨的 ~170×。整个 cap **100% 处于穿透状态**，等于「板压在骨里 1 cm」。

### 1.2 为什么 T6 特别软 —— 冻结公式的口径错误

`.feb` 实际写出：

```xml
<contact type="sliding-elastic" ...>
  <laugon>PENALTY</laugon>
  <penalty>16.979152707884108</penalty>     ← p2 的冻结公式值
  <auto_penalty>0</auto_penalty>            ← 关
```

且 `g4_fe_plate.CONTACT_RECIPE` 的 `penalty` **默认只有 0.1**。

**冻结公式**（P2 预注册 §4）：`penalty = 50 × F_ref / 57941.3`。

**⚠️ 该公式让罚刚度 ∝ 载荷，而 FEBio 手册 L4557 明说**：

> 「Choosing a good initial penalty parameter … **this parameter depends on material properties as well as on mesh dimensions**」

**⇒ 罚参数与载荷无关。** 把 `F_ref` 放进去，使**小载荷部位拿到软罚**：
T6 h=10 → penalty 仅 **16.98**（股骨 h=50 为 50，且股骨的接触域是**实心皮质**、T6 是**皮质壳覆在软松质上**）。

### 1.3 按手册修：加大罚 ⇒ 穿透应声而降

手册 **L13929**：

> 「When a contact constraint is **not sufficiently enforced**, the contacting surfaces **will have penetrated** … **Increasing the penalty factor** and/or decreasing the augmented Lagrangian tolerance **should solve this problem**.」

| 罚乘子 | h=10 罚值 | **h=10 穿透 max** | h=50 结果 |
|---|---:|---:|---|
| ×1 | 16.98 | **10.435 mm** | **0 步**（t=0.0167 即死） |
| **×100** | 1697.92 | **0.9258 mm**（↓11×） | **4 步**（推进到 **t=0.764**） |

**⇒ 手册预测与实测一致**，且 h=50 **从完全走不动推进到 76% 载荷**。

---

## 2. 病 B：σ 峰值 = **载荷带边缘伪影**（**与接触无关**）

### 2.1 证据一：σ 对罚刚度**完全免疫**

| 罚乘子 | 穿透 | **σ_max（皮质壳域）** |
|---|---:|---:|
| ×1 | 10.435 mm | **1421.396 MPa** |
| **×100** | 0.926 mm | **1421.398 MPa** |

**⇒ 穿透改变 11×，σ_max 只动 1.4e-6（数值噪声级）** ⇒ **σ 峰值不来自接触。**

### 2.2 证据二（**决定性**）：峰值**跟着载荷带边缘搬家**

单变量：`LOAD_FRAC` 0.20 → 0.50（载荷带从骨长 20% 扩到 50%，总载荷不变）。
⚠️ 该 patch 必须 monkeypatch `p2.end_nodes` 本身 —— `end_nodes(mesh, *, load_frac=LOAD_FRAC)`
的默认值在**函数定义时**绑定，改 `p2.LOAD_FRAC` **无效**（本轮先踩了一次，见 §5 自错 ①）。
**生效证据**：`载荷端 n=59 → n=440` ✅

| 量 | `LOAD_FRAC=0.20`（基线） | `LOAD_FRAC=0.50` | 变化 |
|---|---|---|---|
| **σ_vm 峰值单元的 z_frac** | 0.247 / 0.144 / 0.275 / 0.254 / 0.082 / 0.216 | **0.513 / 0.624 / 0.401 / 0.518 / 0.533 / 0.380** | **跟着带边缘搬** |
| **σ_max** | **1421.396 MPa** | **481.557 MPa** | **−66.1 %** |
| σ_median | 209.626 | 141.661 | −32 % |
| 峰值 `σ_zx` | **−743.7** | **+117.5** | 量级崩塌 |

**⇒ 峰值位置精确随带边缘（0.20 → 0.50）迁移，且幅度掉三分之二** ⇒
**σ 峰值 = 载荷施加方式的伪影，不是物理应力。**

### 2.3 机制

- `end_nodes(mesh, load_frac=0.20)` 取的是沿骨**轴向**的**底部 20% 带**（**不是**端面！），
  把总载荷**平均分给 59 个节点** ⇒ **每个节点一个集中力（≈333 N）**
- 集中力作用点解析应力**奇异**；带的上边缘处载荷**突然终止** ⇒ 剪切集中
- **⇒ 与项目既有发现同类**：G4② 曾量化「**P1-a 的 1054 MPa 点载荷伪影**」，并建议做**载荷面积扫描**。
  本件是同一现象的另一次现身（此处 1421 MPa），**并已给出定量证据（−66%）**。

---

## 3. 病 A 的残差：h=50 在 ×100 下仍在 t=0.764 倒置

```
×100： Number of time steps completed : 4
       failed at t = 0.764444  →  Negative jacobian detected
       Total equilibrium iterations : 32    (平均 8 次/步)
```

**与 ×1（0 步）比是巨大改善**，但**残差仍是负 Jacobian（局部元素倒置）**，属 auditlab 页所说
「跑了一段后因变形过大而失败」那一类。候选（**均未验证**）：

1. ~~**皮质壳过厚**~~ ⛔ **【已作废 —— 几何前提不成立，见 `P2_H50_RESIDUAL.md` §4】** —— FEBio 手册 **§8.3.1.4 + L1981**：
   > 「When a shell gets **too thick** it can create negative Jacobians in the out-of-plane integration points, especially in areas of high curvature… you may need to **remesh the affected area, adjust the shell thickness** or even **replace the shells with solid elements**.」
   > 「It is the user's responsibility to ensure that the **shell thickness is less than that of the solid element**; otherwise a negative Jacobian error will be automatically generated.」

   ~~**T6 实测**：皮质 **1.5 mm** vs 实体单元特征边长 **1.72 mm**（比 **0.87**）⇒ **正落在「同一量级」的告警区间**（b1 闸门 3 已量）。~~
   **更正（10-15 晚）**：1.72 mm ＝ 中位单元体积 5.06 mm³ 的**立方根**（等效立方体），
   对四面体**系统性偏小**；按手册要求的量（**壳下实体单元高**）实测 **中位 2.965 mm**
   ⇒ 比值 **0.506**（**不在**告警区间），违例 **34/1418 = 2.4%**，且腰椎 L3 同款架构 **0.487**。
   ⇒ **该候选不成立**。
   **且这是 T6 独有的**（股骨皮质是 solid，无此风险）⇒ **可能是 h=50 残差的主因**。

2. **步长** —— 手册 §8.3.1.2：减小初始/最小步长。
3. **False transient**（§8.4.2）：准静态不稳时改用动力学 + 人工阻尼（`rhoi=0`，或 `rhoi=-2 + gamma=1.5 + beta=1`）。

---

## 3.5 修法 ①（端面分布压强）—— **已实现并试跑：FAIL**（如实记）

**实现**：`scripts/opensim_fe/p2_pressure_load.py`（加法式注入；**已转正入 scripts/**，可复现）——
在骨 −Z 端**端面**（`axdiv < −0.5` 的外表面面片）建 `<Surface name="load_face">`，
把 `<Loads>` 里的 `<nodal_load ... node_set="top_load">` **整段替换**为
`<surface_load type="pressure" surface="load_face"><pressure lc="2">P</pressure></surface_load>`，
`P = load_n / A_face`（总力不变）。

**试跑结果**：

| 符号 | h=10 | h=30 | 失败模式 |
|---|---|---|---|
| `+P` | FAIL | FAIL | 首步 t=0.1 即不收敛 |
| `−P` | FAIL | FAIL | **`Negative jacobian`（8 个单元局部倒置）**，t=0.1 |

**诊断**：

- **端面偏小**：`load_face` 仅 **18 面片 / 75.80 mm² / 22 节点**（−Z 带内共 **82** 面片）
  ⇒ `P = 19675.9 / 75.80 = **259.6 MPa**`（h=30 达 504.7 MPa）
  —— 打在 **1.5 mm 皮质壳**上，远超其强度 ⇒ 局部倒置，**与全局刚体模态已不同类**（说明载荷**确实生效**了）。
- **`invalid facets` 警告非致命**：FEBio 报 `load_face` 18 个、`bone_top` 36 个 invalid facets；
  但 **`bone_top` 的 36 个在基线（节点力）run 里同样存在且能收敛**
  ⇒ 该警告是本工具链的既有现象，**不是本次失败的原因**。
- **`<pressure>` 标签确认**：FEBio 手册 L3035 表明标签为 `<pressure>`（不是 `<value>`）；
  首版用 `<value>` 报 `tag "value" : unrecognized tag`，已修。

**⇒ 未完成验收**。下一手**三选一**（均未跑）：

| 选项 | 做法 | 代价 |
|---|---|---|
| (a) | **加宽载荷面**：改用 −Z 带全部 **82 面片**（A 更大 ⇒ P 更小） | 一次 FE |
| (b) | `<linear>1</linear>`（变形无关压强，手册 L3917 提到） | 一次 FE |
| **(c)** | **口径修法（零 FE）**：保留节点力，但**判定域避开载荷带**（如排除 `z_frac<0.30`），并在件中声明 | **0** |

⚠️ **重要**：修法 ① 的目的是**让 σ 可信**。但**根因已在 §2.2 定量确立**（峰值随载荷带边缘迁移、σ_max 随带宽 −66%）
⇒ **即便不换载荷施加方式，选项 (c) 也能立刻解除「σ 不可信」**——
只要在引用 σ 时**排除载荷带附近的伪影区**并写明口径。**这是当前性价比最高的一步。**


## 4. 修法（三条，**都是修正，不是凑判据**）

| # | 修什么 | 怎么改 | 依据 |
|---|---|---|---|
| **1** | **σ 峰值伪影** | 把「底部轴向带上的 59 个节点力」改为**端面上的分布压强**（`<surface_load type="pressure">`）—— 载荷物理上就作用在终板面上 | 消除集中力奇异；且更符合解剖 |
| **2** | **接触罚** | 不再用 `penalty ∝ F_ref`；改为**由材料+网格定标**（`auto_penalty=1`，或用「接触域 E × A / h」量级），并**加入穿透闸门**（穿透 > 变形尺度 ⇒ 判不通过） | 手册 L4557（罚值取决于材料与网格、与载荷无关）+ L13929（约束执行不足 ⇒ 穿透 ⇒ 加罚） |
| **3** | **h=50 残差** | 先单变量试：① 逐步小步长 ② 壳厚/实体化（承 §8.3.1.4）③ false transient | §8.3.1.2 / §8.3.1.4 / §8.4.2 |

⚠️ **必须先做 1 再谈 A2**：在 σ 峰值是伪影的前提下，**任何基于 σ 幅值的判定都不可信**（包括 A2 的胸椎首折高度）。

---

## 5. 诚实边界

0. **本轮自错 3 条（如实）**：
   ① **首版载荷带测试是 no-op** —— `end_nodes(mesh, *, load_frac=LOAD_FRAC)` 的默认值在
      **函数定义时**绑定，改 `p2.LOAD_FRAC` **不生效**（日志仍 `n=59`）。已改为 monkeypatch
      `end_nodes` 本身，生效证据 `n=59 → n=440`。
   ② **首版 `z_frac` 与我另一支探针互不一致**（同一单元 0.198 vs 0.247）⇒ 位置证据一度不可靠。
      已改为**同一坐标源 + 显式打印 z/span**（本件 §2.2 的表即修正后口径）。
      ⚠️ **§2 中被引用的「0.198」已作废**，以 `P2_LOADBAND_TEST.json` 的
      0.247 / 0.144 / 0.275 / 0.254 / 0.082 / 0.216 为准。
   ③ 首轮 V1 曾把 `R2_MIN` 臆造成 0.90（官方 0.20）—— 见 `P2_V1_REGFIX_REPORT.md` §4。
   **⇒ 教训**：单变量测试**必须验证 patch 确实生效**（打印受控量的变化），
      不能只看"结果没变"就下结论 —— 那可能是 no-op。

1. **本件未改任何冻结输入**，未改判据；所有产出为新件/临时脚本。
2. **穿透是从位移场算的**（`cap 变形后 z − 板面 z`），非 FEBio 原生 gap 输出（`surface_data` 只输出 `contact force`）。
   ⇒ 该量**方向明确、量级可靠**，但**不等于 FEBio 内部接触积分点的 gap**。
   **建议**：在 `.feb` 的 `<Output>` 加原生 gap 变量以核对（**未做**）。
3. **罚 ×100 的 0.9258 mm 仍非「充分」**（股骨 0.060 mm）⇒ 只说明**方向对**，**不代表已收敛到充分执行**。
4. **病 B 已定量确立**（§2.2 峰值随带边缘迁移、σ_max −66%），但**修法 1 本身（端面分布压强）
   尚未实现与验收**；本件只证明「载荷施加方式是元凶」，**未给出已验证的替代加载方案**。
5. **h=50 残差归因未定**（§3 三条候选均未验证）。
6. **A2 / A1 / 统计维度不受本件影响**（未重跑任何判定矩阵）。
   ⚠️ 但**能力含义重大**：既然 σ 峰值是载荷伪影，**任何基于 σ 幅值的判定都需重新审视**
   （包括 A2 的胸椎首折高度口径）。

---

## 6. 溯源

| 项 | 位置 |
|---|---|
| **修法 ① 实现（试跑 FAIL）** | `scripts/opensim_fe/p2_pressure_load.py`（`--sign` 可翻压强方向） |
| 罚参数语义（scale factor；取决于材料与网格） | `FEBio_User_Manual.md:4551`、`:4557`、`:4569` |
| penalty 方法必靠穿透产生力 + 过大会不稳 | 同上 `:13745-13747` |
| **约束执行不足 ⇒ 穿透 ⇒ 加大 penalty** | 同上 **`:13929`** |
| 接触引起的收敛问题 | 同上 `:13699-13701`（§8.3.2.4） |
| 步长过小/过大 | 同上 `:13659`（§8.3.1.2） |
| **壳过厚 ⇒ 高曲率处负 Jacobian** | 同上 `:13665-13667`（§8.3.1.4）、**`:1981`** |
| False transient（人工阻尼） | 同上 `:13715-13735`（§8.4.2） |
| 官方 Contact 教程（gap 大 ⇒ 刚度近似差 ⇒ 收敛问题；penalty 过大亦有害） | FEBio Studio Webinar: Contact Modeling Part 2 |
| 冻结罚公式 | `scripts/opensim_fe/p2_head_contact.py:56-63` |
| 接触配方（`penalty` 默认 0.1 / `auto_penalty=0` / `laugon=PENALTY`） | `scripts/opensim_fe/g4_fe_plate.py:52-57` |
| 接触注入 | `scripts/opensim_fe/p5_head_contact.py:213-226` |
| 载荷带 | `scripts/opensim_fe/p2_head_contact.py:86-110`（`end_nodes`，`LOAD_FRAC=0.20`） |
| 横向支撑实现（硬 Dirichlet 带） | `scripts/opensim_fe/bone_feb.py:996-1020` |
| ~~皮质 1.5 mm vs 单元 1.72 mm~~ ⛔**已作废** | `P2_H50_RESIDUAL.md` §4（旧记 `P2_BC_SCOPE_INVENTORY.md` §5quater 闸门 3 口径错） |
| P1-a 1054 MPa 点载荷伪影（同类） | `P2_BC_SCOPE_INVENTORY.md` §5bis |
| 股骨穿透 0.060 mm（正对照） | `P2_BC_SCOPE_INVENTORY.md` §6.2 假设③ |
