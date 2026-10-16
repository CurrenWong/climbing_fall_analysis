---
title: FE · BC 屈曲与 σ 口径陷阱（2026-10-12 实证汇总）
date: 2026-10-12
tags: [FE, FEBio, BC, 屈曲, sigma口径, 陷阱, G4, G9, 腓骨]
---

> **用途**：凡做「单骨/细长件 + 边界固定 + 端部载荷」的 FE，**开工前先过一遍本件**。
> 每条都带**不可被约定翻转的证据**（数值/正对照），可复现。
> 证据索引见文末。

## 1. ⚠️ 最贵的一条：BC 屈曲 —— 载荷超 P_cr 时 FEBio 仍"收敛"

**症状**：FE 正常终止（rc=0）、σ 是个**看起来正常的大数**（1 000+ MPa），但它是**后屈曲折叠伪影**。
**比发散更危险**（发散至少会被发现）。

- **先估 P_cr**：本模型（腓骨单骨、自由长 283 mm、I≈682 mm⁴、E=12 GPa）
  `P_cr = π²EI/(2L)² ≈ 252 N`；而设计载荷 1 000 N = **4×P_cr**。
- **实测拐点 200–300 N**，与 Euler 估计独立吻合（`G4_FE_AREASCAN_REPORT.md` §2.3）。
- **哨兵**（比 σ 更早）：
  1. **特征点位移**：1 000 N 时 apex 轴向位移 **−294 mm**（自由长 283 mm 同量级 = 折叠）；
     有效域（100 N）只有 1.3–2.7 mm。
  2. **σ/F 是否恒定**：以载荷为扫描轴，`σ/F` 从 0.081（25 N）→ 1.099（1 000 N）= **×13.6**
     ⇒ 非线性 ⇒ 后屈曲。线性域上限 ≈200 N。
  3. **峰值位置**：若峰值落在**固定端边界**（本例 z_frac 0.32 = 固定带边缘）且与载荷区
     `in_patch=0.00`，那是 BC 边界层，不是接触/载荷效应。
- **后果提醒**：本课题因此**作废了三条旧结论**（§9.4 点载荷奇异判据、§9.10 复活条件前提、
  §9.9 首折 h=14 m）——全部是屈曲伪影。
- **修法**：见 §2；或缩短自由长（杠杆很弱：bf 0.30→0.60 只降 19%）。

## 2. 横向支撑：**覆盖率**才是关键参数（≥80% 自由长）

| 侧向支撑覆盖自由长 | FE 收敛 | σ/σ_c | 备注 |
|---|---|---|---|
| 100% | 50/50 | 0.54 | 唯一、线性 |
| **80%** | **50/50** | **0.60** | **最小可行** |
| 60% | 44/50 | 3.65 | 屈曲回归 |
| 40% | 21/50 | 2.91 | 屈曲回归 |
| 20% | 17/50 | 3.21 | 屈曲回归 |

- **≤60% ⇒ 后屈曲 ⇒ 解非唯一**（同一配置重跑会在两支间跳变：实测 18.2 / 167.1 MPa
  ——这就是 G9「单带分叉非唯一」的机制）。
- **踩坑实例（G9 D1a）**：支撑带取 `[0.45,0.55]`，而该 deck 的**固定区占底 30%、载荷带占顶 30%**
  ⇒ 自由区 = `[0.30,0.70]` ⇒ 覆盖率只有 **25%** ⇒ 解非唯一（同配置重跑两支跳变）；
  同 deck 的三段带覆盖 **50%** ⇒ 稳定线性。**教训：覆盖率必须按「自由区长度」算 —— 先减去固定区
  与载荷带**，别拿全长比例当覆盖率。（我第一版把机制写成"一半落在固定区"，是**读**出来的、
  不是量出来的 —— 量完才纠正。）
- **硬支撑 vs 软弹簧**：软化**不会**降低 σ——σ(k) 单调（k 从 1e6 降到 1：31.7→65.9 MPa），
  ⇒ **硬（导向）支撑是最优**；另得最小可用刚度 ≈1 N/mm/节点（合计 ~27 N/mm）。

## 3. σ 口径：横向约束 ⇒ 单轴应变 ⇒ **σ_vm 被压低**

- **⚠️ 统计量本身也会骗人（2026-10-12 实测）**：当域内应力**分布极偏斜**时
  （如 `fib_cort`：max≈1 100 MPa 而**带符号 median≈0**），**`median` 落在"接近零"的噪声地板上 ⇒
  刀锋敏感**：同一载荷、**逐位相同的 deck**、两次解 `median|σ_zz|` 差 **23%**、且 2× 载荷只升
  1.15×（应 2.0）。同两次解的稳健性：**mean 0.5% < p95 6.3% < max 8.7% ≪ median 23%**。
  ⇒ **判据：报"同载荷两次解"的相对差**；偏斜分布**弃用 median**，改 `mean`/`p95`。
- **⚠️ 但"换统计量"≠"换水平"**：`p95/max` 的绝对值可能远超强度（479/1092 MPa ≫ σ_c=70）
  ⇒ 直接换 ⇒ util≈7 ⇒ 判据变成 `all_one`（不可辨识）。**换统计量改变的是 util 的绝对水平，
  不改变 OR 的结构**（OR 由 load(h) 形状 + 交叉高度决定）。
- **✅ 首选：平衡锚定的名义应力** `σ_nom = F/A_eff`，其中 A_eff 由**载荷流实测**。
  正对照（必做）：`Σσ_zz·V/Δz ≈ 施加载荷`、**固定区 ≈ 0**；`A_eff` 应在不同解间**重现**
  （实测 76.01 vs 76.01 mm²，与独立实测皮质 76.8 mm² 吻合）。

- 长程横向约束下横向膨胀被约束 ⇒ 单轴应变 ⇒
  **`σ_vm/|σ_zz| = (1−2ν)/(1−ν)|_{ν=0.3} = 0.571`**；**实测 0.566–0.573**（全覆盖长导向），
  短支撑回到 **0.99**。
- ⇒ **压缩/骨折判据必须用 `|σ_zz|`（或主应力），不能用 σ_vm**；σ_vm 会低估 1.75×。
- **诊断量**：报 `vm/|σ_zz|` —— ≈1.0 = 自由/近单轴应力；≈0.57 = 完全侧向约束（单轴应变）。
- **载荷流自检**：`Σ(σ_zz·V)/段厚` 应 ≈ 施加载荷（本例 z_frac 0.5–0.6 处 −1549 N ≈ 载荷 ✓）；
  固定区内应力应 **= 0**（节点全固定 ⇒ 零应变）——若固定区有应力，说明带选错了。

## 4. 坐标系陷阱：锁"世界 x/y"可能夹住轴向（`R_net/load ≠ 1` 就是它）

- deck 用 `align_pca=False` + `load_axis=-pca` 时，**骨长轴不在 z 上**（实测 `pca·z≈0.734`）。
  直接加 x/y 零位移 ⇒ 约束里有 **~68% 的轴向分量** ⇒ 把骨**部分轴向夹死**。
- **铁证**：`R_net/load = 0.734`，且 `1298.9 N = 1768.6 N × 0.7343` **恰为载荷的 z 分量**
  （弹性/损伤两版相同 ⇒ 与材料无关，纯 BC 假象）。
- **判据**：**`R_net/load` 必须 ≈ 1.000**（弹性与损伤都查）。修法：改对齐帧（`align_pca=True`）。
- 另注：FEBio 的 `reaction forces` 在**载荷节点同时被约束**时会漏掉该方向的分量
  （分组检查可见 `bottom_fix Σ=0`、全部反力出现在载荷节点）⇒ **别只信"总和"**。

## 5. 截面口径：`A = V/L` 是**全体积/长度**，不是皮质截面

- 本例：§9.8 的 `A≈102 mm²` 实为**全体积/L**（皮质 76.8 + 松质 ~25）；
  **实测皮质 = 31 015 mm³ / 403.8 mm = 76.8 mm²**，且沿轴**强非均匀**
  （z 0.7–0.9 最薄 **48.7 mm²**，两端 107–117 mm²）。
- ⇒ 用截面做名义应力前，**按段实测** `A_eff = Σ|V_hex| / 段厚`（复用仓库
  `bone_feb._hex8_signed_volumes`，注意其局部索引是 **1-based**）。
- 名义 `F/A@1000 N`：总截面口径 9.8 MPa vs **皮质口径 13.0 MPa（最薄 20.5）**。

## 6. 其它已实证的约定坑（速查）

| 坑 | 判据/修法 |
|---|---|
| h5 `meshes/0/domains/<dom>` 布局 = **(n_elem, 1+n_nodes)**（第 0 列=单元 id，其余 **0-based** 节点位置） | 与 Python 侧单元**逐行**比对（本例 842/842 ✓）；当"8 节点+材料 id"会整表错位 |
| `nodesets/*` 是 **0-based 位置**；`meshes/0/nodes` 的 `id` 是 1..N | 坐标 = `nodes[id-1]` |
| Python 侧 `facet_geom()[...]['poly']` 是 **1-based 全局 id** | 取坐标须 `nodes[poly-1]`；写 NodeSet 可直接用 |
| `rigid_load` 平移 DOF 拼写是 **`Rz`**（不是 `z`） | 用 `dof=z` 会**读文件失败**（rc=1）；对照 `thums_feb.py` |
| **FEBio 弹簧接地面**：偏移 d 必须 ≫ 位移 | 长度 ≈ `d − ux + (uy²+uz²)/(2d)`；d=1e-3 mm 时刚弹簧变**刚性系绳**、把骨轴向夹死（vm 0.1 MPa、\|u\|≈0）。取 d=1000 mm 可忽略 |
| 6-tet 体积只是**符号代理**，强翘曲面片会骗过它 | 真判据 = 三线性 hex8 的 **detJ>0**（8 角点 + 2×2×2 Gauss） |

### 6b. pyfebio / deck 构造的约定坑（2026-10-15 新增，5 条全部实测踩到）

> 来源：S5/S7/S10「距下关节改刚性距骨+真实接触」与「细化两域模型」两轮实现。
> 每条都**实际导致求解失败或结果无效**，且**报错信息不指向真因**。

| # | 坑 | 判据 / 修法 |
|---|---|---|
| **B1** | 🔴 **双 `<Control>` 块**：pyfebio 会**同时**发射 top-level `Model.control_` 与 **step 内 `<Control>`**，而 **FEBio 取 step 级** | ⇒ ①`pyfebio.control.Control(...)` 传 `time_steps`/`step_size`/`opt_iter`/`time_stepper` 可能**被静默忽略**；②**后处理改 deck 必须改「所有」`<Control>` 块**（非贪婪 `re.search` 只改第一块 = 完全无效）。**判据**：FEBio 日志 `CONTROL DATA` 段的 `time_steps` / `step_size` 是否等于你意图的值。见 `thums_feb.py:964-969` 注释 |
| **B2** | **重名 `<Surface>`**：builder 已建过同名面后再建 ⇒ FEBio 直接拒绝 | `A surface with name "X" was already defined` ⇒ 复用已有面（`create_bone_surface=False`），不要重复 `add_surface` |
| **B3** | **hex8 节点序镜像** ⇒ 初始化即报 `Negative jacobian detected during mesh initialization` | 前 4 个节点 = 一个面、后 4 = 另一面；**与 `plantar_bc._build_pad_floor` 同款**（那里前 4 = y 较大的一侧）。**数值判据**：算 hex8 有符号体积，负则反转（`joint_contact._hex8_signed_volume`） |
| **B4** | 🔴 **`st.boundary` 绑定过早** ⇒ 之后 `bcs.append(...)` **静默无效** | `st.boundary = fbc.Boundary(all_bcs=bcs)` 之后追加的 BC **不会进 deck**。**判据**：数 deck 里 `<bc ` 的个数是否等于你追加的条数（本案例：应为 3、实际 1）。修法：把追加 BC 的调用移到绑定**之前**。**症状极具误导性**——帽 150 个 DOF 全自由 ⇒ 位移爆到 **1e13**、弹簧看起来"完全无效"、单元倒置 |
| **B5** | **接触面法向必须朝配偶面**；朝反 ⇒ 求解中单元倒置（`N negative jacobians detected`） | 用叉积复算：quad `[v1,v2,v3,v4]` 的法向 = `(v2−v1)×(v3−v1)`。**朝骨/朝板要显式算，不能靠"看起来对"**——本案例第一版写反，导致 16 个帽单元全部倒置 |

**配套两条（接触专属）**

* **C1 · tri3 facet 挂在 hex8 域上 ⇒ FEBio `invalid facets` ⇒ 接触完全不承载**（老坑，见
  `docs/S5_contact_facet_fix.md` §1.3）。修法 = **按 boundary_polys 重建 quad4 facet**
  （`_build_thums_joint_quads` / `_build_thums_plantar_quads`）。
  **判据**：接触路径下"接触力/总载荷"应 ≈1；若 ≈0 且位移爆 1e12 ⇒ 先查这条。
* **C2 · 单向接触 + 某个 DOF 无支撑 ⇒ 刚度矩阵奇异**（位移爆 1e12~1e13）。
  修法 = 给该 DOF 挂**很弱**的 hold 弹簧。**但 k 有窗口**：k **太大**会"**偷**"载荷
  —— 帽 296 节点 × 100 N/mm = **29600 N/mm**（比骨侧接触路径更硬）⇒ 帽坐在弹簧上、
  压根没压到骨（σ_max≈0）。实测窗口：**k=0.01 N/mm/节点可用；1e-3 不收敛；1.0 偷载荷**。
  **判据**：报「弹簧总反力 / 总载荷」之比，须 ≪1。

## 6c. 🔴 多步静态分析**必须**挂 load ramp curve（2026-10-16 新增，S10 三轮 BLOCKED 的真因）

**症状**：`time_steps=N (N>1)` 的多步 STATIC 分析，**第 1 步成功、第 2 步起判据退化**：

```
 Number of time steps completed .................... : 1
 ------- failed to converge at time : 0.0004...        (≈ dt×1.05)
 convergence norms :     INITIAL         CURRENT         REQUIRED
    residual            4.359361e-17    3.726064e-18    0.000000e+00   ← 要求恒零！
 ...
 * Max nr of iterations reached. *  /  * Max nr of reformations reached. *
 E R R O R   T E R M I N A T I O N
```

**机理**：无 `lc`（load curve）时载荷**在第 1 步就加满**；第 2 步载荷增量 = 0 ⇒
残差 ≈ 机器零、`REQUIRED = rtol×R0 = 0` ⇒ Newton 永远不满足 ⇒ 报"不收敛"。
**这不是物理不收敛，是判据退化。**

**修法**（`thums_feb.py:1019` 早有注释，见 `load_ramp_time_s`）：给载荷挂一条
`0→1` 的线性 ramp curve，`PressureLoad(..., pressure=Scale(lc=1, text=p))`：

```python
model.loaddata_.add_load_curve(fld.LoadCurve(
    id=1, interpolate="LINEAR", extend="CONSTANT",
    points=fld.CurvePoints(points=["0.0,0.0", "1.0,1.0"])))
```

**⚠️ 陷阱放大**：载荷**越小**每步迭代越多 ⇒ 在固定 `max_retries` 内走得越**少** ⇒
表现为「载荷小 100 倍反而更早失败」的**非物理**现象（实测 t=0.0129 vs 0.0032）。
**看到"非单调于载荷"的失败，先怀疑判据/预算，不要先归因大变形。**

**同族症状**：`nonvertical_s5_contact_wave5.py:8-11` 的「判据退化为零残差」、
`3a9c857` 的「载荷没进骨（σ≈7e-11）」。

**另外两处 FEBio 4.13 适配**（新写 deck 时必做，否则 FEBio 直接拒绝读文件）：
- `<solver type="...">` → 必须替换成 `<solver/>`（调 `febio_model._fix_febio413`）
- 非多步路径 `control.Control(..., time_stepper=None)`，否则落盘
  `<time_stepper type="default">` ⇒ `invalid attribute "type"`

---

## 6d. 🔴 两域分层**必须体积守恒**，否则网格对比被混淆（同上，2026-10-16）

**症状**：按「单元质心到表面距离 ≤ 固定阈值 t」分皮/松两域时，**网格越细皮质体积越大**：

| mesh | 单元边长中位 | 皮质体积 (mm³) |
|---|---|---|
| cl=3.0 | 3.847 mm | 1789.7 |
| cl=2.5 | 3.220 mm | 6823.0 |
| cl=2.0 | 2.585 mm | 13470.4 |

**⇒ 7.5 倍**。此时 σ 的任何"细化趋势"都混了材质分布的变化 ⇒ **不是纯网格效应**。
（朴素阈值下 `cl=4.5` 甚至只剩 **4 个**皮质单元 ⇒ 退化。）

**修法**：**逐档二分解出阈值 t**，使 `ΣV(质心距 ≤ t) = t_target × A_ref` 为**常数**
（本项目取 `t_target = 1.40 mm`、`A_ref = 13088.876 mm²` ⇒ 目标 **18324.4 mm³**）。
实测解出 2.7212 / 2.0474 / 1.8253 / 1.6197 mm。

**两种结果都要报**：朴素阈值下同序列给**下降**趋势（691→417→223 MPa），
守恒阈值下给**平坦**（183.8/189.6/179.7/170.6）——**只报后者是选择性报告。**

**⇒ 由此得到的判决**（`S10_MESH_CONVERGENCE.md`）：`gauge_max`(ℓ=4 mm) **收敛 ±5%**
（tet 数变 10 倍），裸 `max` **发散 +74%**，热点恒在皮/松界面（**S9b 坐实**）。

---

## 7. 三条纪律（本轮全部踩到并靠它们纠错）

1. **每个判据都要有正对照**：k→∞ 必须复现硬 BC 结果（实测 31.7 vs 31.6 MPa ✓）；
   `guide_ids(0.50,1.00)` 必须等于 deck 的节点集（1002/1002 ✓）；
   h5 连接表必须与 Python 逐行相等（842/842 ✓）；截面 V/L 必须与名义值对得上（**这条没过，才发现了 §5**）。
2. **先验产物再报结论**：`R_net/load≈1`、`max|u|` 是否物理、`σ/F` 是否恒定、
   σ 峰值落在哪里（载荷区 / 固定端 / 支撑边缘）。
3. **更正替换，不叠加**：旧结论被推翻时明确写"作废"，不要在文件里留下两条相反论断。

## 8. 证据索引

- `docs/methods/G4_FE全量_工作计划.md` §9.13（A–D 四类扫描）/ §9.14（P1 试点导向重跑）/
  §9.15（BC 声明 + ② 结论修订 + 覆盖要求）
- `results/opensim_fe/G4_FE_AREASCAN_REPORT.md`（面积 / 载荷 / 固定长度 / 支撑刚度）
- `results/opensim_fe/G4_FE_PILOT_GUIDED_REPORT.md`（§5e σ 口径、§5f 弱化导向）
- `docs/methods/G9_腓骨FE侧向支撑设计.md` §11（D1a′/D3′ 修订提案）
- 脚本：`scripts/opensim_fe/g4_fe_areascan.py`、`g4_fe_pilot_guided.py`
- 结果 json：`results/opensim_fe/G4_FE_AREASCAN_{full,valid100,loadsweep,lin100,bcscan,latscan,latarea,kscan}.json`、
  `G4_FE_PILOT_GUIDED_{ctrl2,full_phys,full_frozen,guidesweep,sharefull,sharefull_g60}.json`

**§6b（deck 构造 5 条 + 接触 2 条）的证据索引（2026-10-15）**

- 实现：`src/climbing/coupling/joint_contact.py`（`_build_talus_cap` / `_build_conforming_cap` /
  `_hex8_signed_volume` / `_apply_cap_hold_spring` / `apply_joint_contact`）、
  `scripts/opensim_fe/thums_feb.py`（`joint_bc="contact"` 分支、`_build_thums_joint_quads`）
- 提交：`76d3330`（S5 首版 + B1~B5 逐个修）、`e1c5cd4`（S7b 找到 C2 的 k 窗口）
- 报告：`results/opensim_fe/V1_LOAD_INTRO_SINGULARITY.md` **§6**（5 bug 表 + 结果表）、**§8**（S10）
- 相关闸门：`docs/notes/FE_PRESOLVE_GATES.md` **闸门 3**（D1~D4 同款坑）
- 老坑交叉：`docs/S5_contact_facet_fix.md` §1.3（C1：tri3 挂 hex8 ⇒ `invalid facets`）

