# 跖面接触 BC 修复探针（②）—— 2026-10-15

> **触发**：口令「按推荐落库」后的待办 ②「**修跖面 BC**（手册 §8.1.3 / §8.5）」。
> **工具**：`scripts/opensim_fe/plantar_contact_fix_probe.py`
> **产物**：`temp/fe_out/plantar_fix/h{2,5}/probe.json`（**temp/ 下，不入库**）

---

## §0 · 结论（**两条硬结论 + 一条工具修复**）

| # | 结论 | 证据 |
|---|---|---|
| **1** | **根因确认**：`plantar_bc="contact"` 失败 = **欠约束（刚体模态）**，不是网格/步长/罚参数 | 手册 §8.1.3 / §8.3.1.5 / §8.5.3 ＋ 残差签名（下 §2） |
| **2** | **`hold_spring_k` 能收敛，但*不是*干净修复** —— 所需 k **随载荷变**（h=5 需 ≥10、h=2 需 ≥100），且 **σ 随 k 单调趋向 `fixed`**（0.987→1.000）⇒ 它可在 contact↔fixed 之间插值，属**正则化** | 下 §3 两张表 |
| **3** | **发现并修复一个 pyfebio 缺陷**（与本题独立，价值更普适）：含接触的 `.xplt` 有**无名域** ⇒ `to_hdf5` 抛 `KeyError('name')` ⇒ 接触运行的**全部后处理失效** | `src/climbing/coupling/xplt_shim.py` |

**⇒ 手册给的*通用*正解（「先用位移压到接触、再把位移换成力」两步分析）尚未实现 ⇒ ② 未收口。**

---

## §1 · 手册依据（本次把 §8 整章读了一遍，§0 规矩的又一次收益）

| 手册节 | 原文要点 | 与本问题的关系 |
|---|---|---|
| **§8.1.3** | over-constrained 最常见原因＝**冲突的边界条件**，例如「**一个力施加在一个被固定的节点上**」 | **排除**：contact 分支只固定**pad**，不固定骨 ⇒ 不是此因 |
| **§8.1.3 / §8.3.1.5** | **underconstrained** ⇒ 刚体模态 ⇒ *"FEBio will most likely throw **negative Jacobians**"* | **命中**：contact 是**单边**约束，gap 未闭合前骨无切向支撑 |
| **§8.5.3** | *"FEBio often struggles with problems that have an **initial separation**. This is especially so when the problem is **force-driven** … when the contact interface is necessary to define a **well-constrained** problem it is best to **avoid initial separation** … **Another way around this issue is to first use a displacement to bring the surfaces in contact and then replace the displacement with a force.**"* | **命中**：我们是 force-driven＋有 0.5 mm 初始间隙 |
| **§8.2** | `-g` / **`-g2`**（只把**导致失败的非收敛态**写进 plot）；`-g` 还查**刚度矩阵对角线上的零** | 已接入探针 `--debug`（本次未用上，见 §5 待办） |

**仓库自己的注释早已写下同一判断**（`src/climbing/coupling/plantar_bc.py`，`hold_spring_k` 处）：

> 跖面接触是**单边**约束（只在压入时提供刚度），静止时切向摩擦也为零 ⇒ 在 gap 未闭合前，
> 骨只有接触支撑 ⇒ 切向刚度矩阵奇异 ⇒ Newton 解出 1e23 级位移而发散。

**⇒ 缺陷是*已知*的，只是那个开关默认 `None`（关），而 `s2_thums_sweep` 只传 `plantar_bc`。**

---

## §2 · 残差签名（支持「模型错，不是步长错」＝闸门 6）

h=2 / baseline（gap=0.5、无 hold）的 `.log`：

```
residual  2.211055e+07 → 2.782211e+08 → 7.671263e+08 → 7.003044e+09
        → 2.208146e+10 → 6.065578e+10 → 2.670964e+11
```

- **每次迭代涨 3–4×** ⇒ 发散，不是「刚度大但能收敛」
- 自动时间步从 dt=1 一路退到 **0.1667**，而 **`residual INITIAL` 恒为 `2.211055e+07`、残差序列逐位相同**
  ⇒ **与步长无关** ⇒ 手册 §8.3.1.2「时间步太大」**被排除**

---

## §3 · 变体扫描结果

### 3.1 h=5（subtalar 28520.6 N）

| 变体 | gap | hold_k | rc | negJ | gauge_max (MPa) | σ/σ_fixed | 热点质心 (mm) |
|---|---|---|---|---|---|---|---|
| **baseline** | 0.5 | — | **1** | **20** | — | — | — |
| gap 0.01 | 0.01 | — | 1 | 20 | — | — | — |
| gap 0.001 | 0.001 | — | 1 | 20 | — | — | — |
| hold 0.1 | 0.5 | 0.1 | 1 | 13 | — | — | — |
| hold 1 | 0.5 | 1 | 1 | **79** | — | — | — |
| **hold 10** | 0.5 | 10 | **0** | — | **447.0767** | **0.987** | (−6.6, 17.4, −11.2) |
| hold 100 | 0.5 | 100 | 0 | — | 446.8582 | 0.987 | (−6.6, 17.4, −11.2) |
| hold 1000 | 0.5 | 1000 | 0 | — | 448.4032 | 0.990 | (−6.6, 17.4, −11.2) |
| gap0.01＋hold10 | 0.01 | 10 | 0 | — | 447.0767 | 0.987 | (−6.6, 17.4, −11.2) |
| swap＋hold10 | 0.5 | 10 | 0 | — | 447.0763 | 0.987 | (−6.6, 17.4, −11.2) |
| 10 时间步＋hold10 | 0.5 | 10 | **1** | — | — | — | — |
| **ref: fixed** | — | — | 0 | — | **452.7664** | **1.000** | (−5.6, 11.4, −16.6) |
| ref: spring k=100 | — | — | 0 | — | 446.8582 | 0.987 | (−6.6, 17.4, −11.2) |

**h=5 下 k∈[10,1000] σ 只动 0.3% ⇒ 看起来「纯收敛修复」。但 h=2 推翻了这一点：**

### 3.2 h=2（subtalar 22723.6 N）

| 变体 | gap | hold_k | rc | gauge_max (MPa) | σ/σ_fixed | 热点质心 (mm) |
|---|---|---|---|---|---|---|
| **baseline** | 0.5 | — | **1** | — | — | — |
| **hold 10** | 0.5 | 10 | **1** | — | — | — |
| hold 100 | 0.5 | 100 | **0** | **357.9889** | **0.987** | (−6.6, 17.4, −11.2) |
| hold 1000 | 0.5 | 1000 | 0 | 359.3576 | 0.991 | (−6.6, 17.4, −11.2) |
| hold 10000 | 0.5 | 1e4 | 0 | 361.8341 | 0.998 | (−6.6, 17.4, −11.2) |
| **hold 100000** | 0.5 | 1e5 | 0 | **362.4641** | **1.000** | **(−5.6, 11.4, −16.6)** |
| gap0.01＋hold10 | 0.01 | 10 | 1 | — | — | — |
| swap＋hold10 | 0.5 | 10 | 1 | — | — | — |
| **ref: fixed** | — | — | 0 | **362.5611** | **1.000** | (−5.6, 11.4, −16.6) |
| ref: spring k=100 | — | — | 0 | 357.9889 | 0.987 | (−6.6, 17.4, −11.2) |

**⇒ 三条读数**

1. **所需 k 随载荷变**：h=5 需 ≥10，h=2 需 ≥100 ⇒ 没有单一的「安全 k」。
2. **σ 随 k 单调升向 `fixed`**（0.987 → 0.991 → 0.998 → 1.000），且 **k=1e5 时热点*回到*固定端位置** ⇒
   **k 足够大时 hold 弹簧就变成支撑、复现 `fixed` BC** ⇒ **它是 contact↔fixed 之间的插值旋钮**。
3. **`gap` 缩短与 `swap` 单独都无效**（都要配 hold 才收敛，且 h=2 下 hold=10 不够）。

**⇒ 裁定：`hold_spring_k` 只能当**诊断**用（回答「contact 能否收敛」＝能），**
**不能当交付口径**（k 是自由旋钮且影响结果 ⇒ 违反「不调参凑靶」）。**

---

## §4 · 副作用副产品：**pyfebio 的无名域缺陷**（已修）

**现象**：接触运行的后处理全部抛 `KeyError: 'name'`。

**链路**：

```
FEBio 为含接触的模型多写域（实测 id=4, etype=5, nelems=573，**无 name**）
  → pyfebio/xplt.py:675  parse_mesh:  mesh_dict["domains"][domain["id"][0]] = domain["name"]
       （同一函数的 **nodeset 分支有空名兜底**：if nodeset["name"] == "": nodeset["name"] = nodeset["id"][0]
         —— **domain 分支没有**）
  → KeyError: 'name'  ⇒ to_hdf5 只写出 meshes、**没有 states**
  → fe_post / s2_thums_sweep._read_domain_vm 全部失效
```

**实测域清单**（`temp/…/m.xplt`）：

| deck | 域 |
|---|---|
| `ref_fixed` | calcaneus(634) · trabecular(3323) —— **都有 name** |
| `hold10` | calcaneus(634) · trabecular(3323) · **pad_floor(16)** · **id=4, etype=5, nelems=573, 无 name** |

**修复**：`src/climbing/coupling/xplt_shim.py::patch_unnamed_domains()`（幂等，给 `_parse_domain_section`
加与 nodeset 分支同款的空名兜底），并在 **`fe_post` 模块级**调用（覆盖所有 `to_hdf5` 调用点）。
**改 shim 而不改 site-packages**：改 `.venv` 会在重装/升级时静默丢失，且遮蔽真实上游缺陷。

**影响面**：**所有** `plantar_bc="contact"`（及任何带接触的模型）的后处理；
**默认 `plantar_bc="fixed"` 的 deck 无无名域 ⇒ 行为逐位不变。**

---

## §5 · 待办（② 未收口）

1. **实现手册 §8.5.3 的通用正解：两步分析**（Step 1 用 prescribed displacement 把跖面压到接触 →
   Step 2 换成力）。这是**唯一不引入自由旋钮**的路。
2. 用 **`-g2`**（手册 §8.2）把**导致失败的非收敛态**写进 plot，直接看骨「飞」向哪 ⇒ 定位是哪个
   刚体模态（`--debug -g2` 已接好，未跑）。
3. 若两步法成功 ⇒ 重跑 **V1 全高度扫描**（`s2_thums_sweep --plantar-bc contact --contact-hold-spring-k …`
   **接口已就绪**），得**新的首折高度**并与现行 `2 m` 对照。
4. **`-g`（查刚度矩阵对角线的零）** 用于回答「那 6/17/20 个皮质 hex8 是否结构退化」。

## §6 · 溯源

| 项 | 位置 |
|---|---|
| 探针 | `scripts/opensim_fe/plantar_contact_fix_probe.py` |
| shim | `src/climbing/coupling/xplt_shim.py`（`fe_post` 模块级调用） |
| 接触实现 | `src/climbing/coupling/plantar_bc.py::_apply_contact`（`hold_spring_k` 默认 `None`） |
| sweep 接口 | `scripts/opensim_fe/s2_thums_sweep.py`（`--contact-hold-spring-k/--contact-gap-mm/--contact-swap-pair`） |
| 手册 | `temp/thums/febio_manual.txt` §8.1.3(L15765) · §8.2(L15783) · §8.3.1.4/.5(L15832) · §8.5.3(L15980) |
| 前序 | `results/opensim_fe/V1_RERUN_20261015.md`（`plantar_bc="contact"` 不收敛、6 个皮质 hex8） |





---



## §7 · 续做（2026-10-15 下半场）：**手册 §8.5.3 的两步法已实现，Step 1 通过、Step 2 卡住**



> **探针**：`scripts/opensim_fe/two_step_pretouch_probe.py`（对 deck 做最小 XML 手术）



### 7.1 手册依据（§6.1 / §6.1.1 / §6.1.2）



> **§6.1.2 Relative Boundary Conditions**：*"…**This makes it possible, for example, to switch

> between load and displacement controlled boundary conditions in multi-step analyses.**"*

> **§6.1.1**：放在**第一个 Step 之前**的 `Boundary`/`Loads`/`Contact` 是 **global**（各步都生效）；

> 放在 Step 内的只对那一步生效。



**⇒ 两步法 =** `Contact`＋pad 固定 BC 保持 **global**；**Step 1** 用 prescribed displacement 把

跖面压到接触；**Step 2** 去掉 prescribe、换上力载荷（＝"replace the displacement with a force"）。



### 7.2 **三个 FEBio 4.13 语法事实**（每一个都花了一次失败运行换来）



| # | 事实 | 证据 |

|---|---|---|

| 1 | **`<Steps>` 不被识别** | `tag "Steps" : unrecognized tag` ⇒ 容器只能用 **`<Step>`** |

| 2 | **`<Step>` 容器内必须以 `<step>` 开头**；全局段必须在**容器之外**（顶层） | `tag "Boundary" : missing tag "step"` |

| 3 | **prescribed displacement 的标签是 `<bc type="prescribed displacement">`**（手册 §3.11.6），**不是 `<prescribe>`** | `tag "prescribe" : unrecognized tag`；语法二次确认自本仓库 `scripts/opensim_fe/g4_fe_stab.py:221` |

| 4 | **`<min_residual>` 不被接受**（顶层与 step 级 Control 都报 `unrecognized tag`）—— 手册 §8.3.2.1 的「零残差」处方**在本版本不可用** | 两次失败运行；探针里已降级为 `None`（不写） |



### 7.3 **关键诊断：只控 `y` 不够**



**第一版只 prescribe `y`** ⇒ 14–44 个负 Jacobian全失败 —— 因为 **x/z 的刚体模态仍在**

（这正是 §8.5.3 要消除的东西）。§8.5.3 的目标是让问题**良约束**。



**⇒ 改成三向全控**（`x`/`z` = 0、`y` = −pretouch）后：**负 Jacobian 消失**，

失败模式从「单元倒置」变成「**不收敛**」 ⇒ **这是实质进展**。



### 7.4 **当前状态（可复现）**



| 步 | 结果 |

|---|---|

| **Step 1（pretouch, −0.6 mm, 20 步）** | ✅ **通过** —— 日志 `Number of time steps completed : 1` |

| **Step 2（force）** | ❌ **不收敛**（time 1.16667 → 2 各退一步皆失败），**无负 Jacobian**，`xplt` 0 B |



**残差签名**（step 2）：`residual ≈ 1e-16`（**实质为零**）、`REQUIRED 0.000000e+00`

⇒ 疑为手册 §8.3.2.1 的「无载荷」情形，但该节给的处方（`<min_residual>`）**在 4.13 不生效**。



**⇒ ② 仍**未收口**，但缺口已**收窄到一个具体问题**：**Step 2 的载荷/收敛口径**，

而非模型或接触本身（Step 1 已证接触能闭合、且预压步稳定）。



### 7.5 下一步（按价值）



1. **查 Step 2 的载荷是否真的进了求解**（`-g2` 导出非收敛态；或临时给 Step 2 加一条已知的小载荷做**正对照**）

2. 查 FEBio 4.13 里「零残差判定阈值」的**正确标签**（SDK 头文件 / `febio4.exe` 字符串）

3. 若 Step 2 通了 ⇒ 用 `s2_thums_sweep`（接口已就绪）跑 **V1 全高度**，得新首折并与 `2 m` 对照





---



## §8 · 续做（同日第三段）：**最小复现器证明「两步法机制有效」；跟骨 deck 的 Step 2 仍不加载**



### 8.1 为什么需要最小复现器



跟骨 2 步 deck：Step 1 通过、**Step 2 残差恒 ~1e-16**。在跟骨模型上换了 **4 种**放置都无效：



| 变体 | 结果 |

|---|---|

| `<Loads>` 放 step 内（原） | Step 2 残差 ~1e-16 ❌ |

| `<Loads>` 放**全局** | 同上 ❌ |

| **删掉顶层 `<Control>`**（手册 §6.1：多步分析里不该有它） | 同上 ❌ |

| Step 2 的 `<Control>` 换**最简形式**（只 `time_steps` + `solver`） | 同上 ❌ |



⇒ 必须**离开跟骨模型**做隔离。



### 8.2 ✅ **决定性正对照**：`min_two_step_load_test.py`（手册 §6.2 的单单元模型扩成两步）



```

[deck] min2step.feb   steps=2  prescribed=1  Loads=1  surfaces=1

[solve] rc = 0        Number of time steps completed : 4      （2 步 × 2 增量）

   Step 1（prescribe）  residual 4.160503e+02  →  5.607475e-17   ✅

   Step 2（force）      residual 3.885663e+01  →  2.726570e-11   ✅★ 非零初残差！

```



**★★ Step 2 的「初始残差」是 38.9（载荷量级），不是 ~1e-16** ⇒

**「`<Loads>` 放进 `<step>`」在 FEBio 4.13 里是生效的**，两步法机制成立。



### 8.3 本段又查出两条 FEBio 4.13 语法事实



| # | 事实 | 报错 |

|---|---|---|

| 5 | **多步模式下，每个 step 的 `<Control>` 必须带 `<solver>`** | `Component "" needs to have property "solver" defined` |

| 6 | `<NodeSet>` 必须是**逗号列表**（`<NodeSet name="x">1,2,3</NodeSet>`），**不是** `<node>` 子元素 | `tag "NodeSet" : invalid value` |

| 7 | `<zero displacement>` 的子标签是 `<x_dof>1</x_dof>`…，**不是 `<dofs>x,y,z</dofs>`**（手册 §6.2 示例用的是后者，4.13 不认） | `tag "dofs" : unrecognized tag` |



（6/7 与仓库 `scripts/ankle_fe/ankle_lig_feb.py:1372-1374` 既有记录一致。）



### 8.4 `min_residual` 的正确位置（**修正手册**）



手册 §8.3.2.1 说它是 **Control section** 的参数 —— **在 FEBio 4.13 里不对**：



- `<Control><min_residual>` ⇒ `unrecognized tag`（step 级、顶层都报）

- **正确位置是 `<solver>` 块** —— 依据：FEBio 运行日志把它印在 **`SOLVER PARAMETERS`**

  标题下；SDK 侧对应 **`FENewtonSolver::m_Rmin`**（`sdk/include/FECore/FENewtonSolver.h:229`）

- 改到 `<solver>` 后**日志确认读到**：`min_residual ... : 1e-15`



**但它没能解决 Step 2**（残差 1.6e-16 本已低于阈值）⇒ 说明 §8.3.2.1 的「零残差」判据

**不是**这里卡住的原因。



### 8.5 结论与下一步（**bisection，现在很便宜**）



**已确认**：① 两步法机制在 4.13 里**有效**（§8.2 正对照）；② 跟骨 Step 1 稳定；

③ 跟骨 Step 2 **载荷未进求解**（残差 ~1e-16，4 种放置皆然）。



**⇒ 下一步 = 用最小复现器做 bisection**（它 0.0 s 跑完，成本可忽略）：

逐个把跟骨 deck 的特征加进 `min2step`，**加到哪一步开始坏，那里就是根因**。候选增量（按怀疑度）：



1. **`<Contact>` 块**（单边约束 + `pad_floor` 域）

2. **`<LoadData>` + `lc="1"`**（跟骨 deck 有，最小复现器没有）

3. **prescribe 三向 + 释放**（最小复现器只控 `x` 一个 dof）

4. **gauge 级网格规模**（1544 节点 vs 8 节点）



**替代路径**（若 bisection 指向 Contact）：回到 §3 的 `hold_spring_k` 路线 ——

它**能收敛**且 σ 与 `fixed` 差 1.3%，只是 k 是自由旋钮；把它限定为**诊断口径**

（而非交付口径）也许是更实际的进展。





---



## §9 · bisection 两轮（**6 项特征逐项加，无一项复现失败**）



> **方法**：最小复现器（0.0 s/次）逐项吸收跟骨 deck 才有特征，**单变量**加进去，

> 看「力步初残差」是否仍为**载荷量级**（~1e1）还是塌成 ~1e-16。

> **脚本**：`min_two_step_bisect.py`（第一轮）、`min_two_step_contact_bisect.py`（第二轮）。



### 9.1 第一轮：`min_two_step_bisect.py`



| 变体 | 加什么 | rc | **力步初残差** | 结论 |

|---|---|---:|---:|---|

| `v0_base` | —（已知 ✅） | 0 | **3.885663e+01** | ✅ 载荷有效 |

| `v1_ldata` | 顶层 `<LoadData>`（lc=1，曲线 0→1） | 0 | 3.885663e+01 | ✅ |

| `v2_3dof` | Step 1 prescribe **三向** | 0 | 6.572883e+02 | ✅ |

| `v3_lc_load` | 压力载荷挂 `lc="1"` | 0 | 3.885663e+01 | ✅ |

| `v4_ldata_3dof` | v1＋v2 | 0 | 6.572883e+02 | ✅ |



### 9.2 第二轮：`min_two_step_contact_bisect.py`



两块 hex8（A 在 x∈[−2,2]、B 在 x∈[2.1,3.1]，**初始间隙 0.1**），A 右面 ⇄ B 左面 组成

`<SurfacePair>`，`<Contact>` 放**全局**（§6.1.1）；Step 1 用 prescribed displacement 把 A

压过去（+0.15 > 间隙）；Step 2 换力载荷。



| 变体 | primary/secondary | rc | **力步初残差** | 结论 |

|---|---|---:|---:|---|

| `c0_A_prim` | Aright / Bleft | 0 | **3.728475e+02** | ✅ 载荷有效 |

| `c1_B_prim` | Bleft / Aright | 0 | 3.728475e+02 | ✅ |

| `c2_rev` | Aright / Bleft（节点序反） | 1 | 3.728475e+02 | ❌ 失败（法向/配对面方向错） |



### 9.3 **结论：6 项特征全部排除 ⇒ 根因在「模型的几何/网格本体」，不在配方/结构**



**已排除**（逐项单变量，皆不破坏两步法机制）：



1. `<Loads>` 放 step 内 ↔ 放全局

2. 顶层 `<Control>` 有无

3. Step 2 的 `<Control>` 简繁

4. `<LoadData>` + `lc="1"`

5. Step 1 prescribe 的 DOF 数（1 ↔ 3）

6. **`<Contact>` + 第二块 + `<SurfacePair>`**（含 primary/secondary 两种配对）



**⇒ 剩下的是跟骨 deck 自身的几何/网格**：1544 节点、**cort hex8 634 + spon tet4 3323**（两域）

＋ pad 16、真实三维形状、以及 Step 1 后留下的**已接触状态**。



**⚠️ 本结果也削弱了「Step 2 载荷未进求解」这一解释的确定性** —— 6 项结构特征都不复现，

也可能是**载荷进了但残差本征很小**（例如 Step 1 后骨已压在垫上、Step 2 的压力被垫反力部分抵消）。

**⇒ 下一轮应先做「真·正对照」**：在跟骨 Step 2 里放一个**独立的小集中力**

（与接触/压力面无关），看残差是否响应 —— 这能一刀切开「没进求解」与「进了但小」。



### 9.4 若上述仍不收敛 ⇒ **转用 `hold_spring_k` 路线**（§3）



`hold_spring_k` **能收敛**、且 h=5/h=2 下 σ/σ_fixed = **0.987**（k=100）；

它的缺点只是 **k 是自由旋钮**。把它**限定为诊断口径**（明确声明「k 用于消刚体模态，

不承载，且已在 k∈[100,1e5] 上验 σ 对 k 的敏感性」）——

比继续追一个只在跟骨几何上出现、且已排除 6 项结构因素的现象**更可能是正路**。





---



## §10 · **收口**：改用「跖面 BC 敏感性」路线（结果件另立）



§9.4 的第二条路已执行 ⇒ **结果件：`results/opensim_fe/V1_PLANTAR_BC_SENSITIVITY.md`**



**一句话**：同一条 `s2_thums_sweep` 流水线内，只把跖面 BC 从 `fixed` 换成

`contact + hold k=100`，**解剖配准网格的首折从 3 m 变成 4 m**（σ −8%~−16%）；

`calcnframe` 网格上**比值符号相反**（+10.6% @h=1）⇒ **跖面 BC 是 V1 的第三个口径自由度**

（前两个：统计量口径、网格/配准）。**V1 仍 ❌【口径不可复现】**。



**本件（① 两步法）的最终状态：**



| 项 | 状态 |

|---|---|

| 根因定位 | ✅ **欠约束/刚体模态**（手册 §8.1.3/§8.3.1.5/§8.5.3） |

| `hold_spring_k` | ✅ **能**收敛；**非干净修复**（k 随载荷变、σ 随 k 趋 `fixed`）⇒ **限为诊断口径** |

| 手册 §8.5.3 两步法 | ⚠️ **Step 1（预压）通过、Step 2 不加载**；bisection **7 项结构特征全排除** ⇒ 根因在跟骨几何/网格本体 |

| pyfebio 无名域缺陷 | ✅ **已修**（`xplt_shim.py`） |

| 跖面 BC 敏感性 | ✅ **已量化**（首折 3 m ↔ 4 m） |

| **② 是否收口** | **是（以「诊断口径 + 敏感性」收口）** —— 不再追「让 contact 完全替代 fixed」 |



**⇒ 未做的（如需继续）**：① 跟骨 Step 2 不加载的根因（需在**跟骨网格**上做 bisection，

而非最小体系）；② 两步法的完整落地（若将来需要「无自由旋钮」的接触口径）。

