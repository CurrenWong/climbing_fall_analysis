# 跖面接触 BC 修复探针（②）—— 2026-10-16

> **触发**：口令「按推荐落库」后的待办 ②「**修跖面 BC**（手册 §8.1.3 / §8.5）」。
> **工具**：`scripts/opensim_fe/plantar_contact_fix_probe.py`
> **产物**：`temp/opensim_fe/plantar_fix/h{2,5}/probe.json`（**temp/ 下，不入库**）

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
| 前序 | `results/opensim_fe/V1_RERUN_20261016.md`（`plantar_bc="contact"` 不收敛、6 个皮质 hex8） |
