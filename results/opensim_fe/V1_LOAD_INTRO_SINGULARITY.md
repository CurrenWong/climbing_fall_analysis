# V1 求解：**首折高度由「载荷引入奇异」决定** —— 实测证明（2026-10-15）

> **触发**：用户要求「想办法解决 V1」。本件是一次**零调参**的归因实验。
> **产物**：`temp/opensim_fe/v1solve/ab_frame.json`（A/B 原始数据）。

---

## §0 · 一句话结论

**V1 的首折高度不是骨的材料/几何响应，而是 `PressureLoad` 施加在距下关节面**子集**上、
在该子集边界产生的**载荷引入奇异**。⇒ 实测：名义压力降 11.24×，峰值应力只降 2.18×。**

---

## §1 · 实验设计（预注册，跑前写死）

**同一网格**（`calcaneus` 1544 节点 / `elements_cort` 634 hex8）、**同一载荷**（h=2 m，
`F_sub = 22723.586 N`，`achilles = 0`）、**同一材料、同一阈值 150 MPa、同一后处理**
（`s2_thums_sweep._read_domain_vm(xplt,'calcaneus')` 的 `vm.max()`）。

**唯一变量**：网格所在的**配准帧** ⇒ 因而改变的是 `subtalar_joint` 具名面的**面积/面数**
⇒ 因而改变 `PressureLoad p = F/A` 的分布。

| 帧 | npz | `subtalar_joint` |
|---|---|---|
| **calcnframe**（`thums_feb.DEFAULT_NPZ`） | `calcaneus_r_calcnframe.npz` | **222.08 mm² / 14 面** |
| **anatframe**（`registration_anatomical.py` 产物） | `calcaneus_r_anatframe.npz` | **2495.33 mm² / 119 面** |

> `registration_anatomical.py:11-14` 逐字记载旧帧：`subtalar_joint` 规则
> （`dist < 25 mm` AND `n·(+Y) > 0.2`）退化成「**a tiny posterior-superior patch
> (14 quads / 222 mm2)**」，在那块上压 `p = F/A` 会给出「**a non-physical "fracture at
> <= 1 m"**」。本件的 222.08 / 14 与之**逐字吻合** ⇒ 口径确认。

---

## §2 · 结果（`ab_frame.json`）

| 帧 | `joint_area` (mm²) | `p_avg` (MPa) | **`σ_max` (MPa)** | **SC = σ_max/p_avg** | σ/150 |
|---|---|---|---|---|---|
| calcnframe | 222.08 | **102.32** | **405.35** | 3.96 | 2.70 |
| anatframe | 2495.33 | **9.11** | **185.64** | **20.39** | 1.24 |

**比值**：`σ_max` **×0.458** ／ `p_avg` **×0.089**。

---

## §3 · 判读（这是本件的核心）

### 3.1 若 σ_max 是「体力响应」，它必须 ∝ p_avg

线性弹性、同一几何族、同一 BC ⇒ 若峰值由**整体载荷水平**决定，则
`σ ∝ p` ⇒ σ_max 应**随 p_avg 一起降 11.24×**（405.35 → 36.0 MPa）。

**实测只降到 185.64（×0.458）。** ⇒ **假设被证伪** ⇒ 峰值**不由**名义压力决定。

### 3.2 SC 的绝对量级本身就是伪影的标记

| 帧 | SC |
|---|---|
| calcnframe | **3.96** |
| anatframe | **20.39** |

骨的体力应力集中系数在关节载荷下是 **2–4 量级**。**SC = 20.39 意味着峰值是平均值的 20 倍**
⇒ 物理上不可能 ⇒ 是**补片边界**的数值奇异（面内牵引从 `p` 突降到 0）。

### 3.3 大面积补片**反而** SC 更高 ⇒ 与「边界」而非「面积/力」相关

补片从 14 面扩到 119 面：平均压力降 11.24×，但**峰/均比上升 5.1×**。
⇒ 峰值由**边界的几何/网格密度**决定。这是**奇异**的定义性质（网格越细/边界越长，峰值越依赖离散）。

### 3.4 与既有结果的交叉一致

* `s2_thums_height_sweep_anat_refixed_now.json` @ h=2 m 的 `max` = **185.64** ⇒ 与本件
  `anatframe` 结果**逐位吻合** ⇒ 确认 `anat` 命名**确实指解剖帧**，
  本实验的 A/B 配对成立。
* `anat_contact_hold100`（仅跖面改接触）SC 由 20.64 降到 14.62 ⇒ 消除**一侧**边界使 SC 下降，
  与「边界奇异」机制**方向一致**。

---

## §4 · 修复方案（预注册 · 物理性 · 非调参）

**把距下关节的载荷施加方式从「面上压 `p = F/A`」改为「刚性距骨（talus cap）+ 真实接触」。**

**依据（三条，逐条可查）**：
1. **[Y25] 的做法**：原文 §1.1「**面面接触**，摩擦系数 **μ = 0.5**」⇒ 载荷经接触传递，
   **没有「压在一个面的子集上」这回事** ⇒ 无补片边界 ⇒ 无该奇异。
2. **项目已有同构先例**：`plantar_bc="contact"`（`apply_plantar_bc` + `_build_thums_plantar_quads`）
   已把**跖面**从 `fixed` 换成「固定数值刚体板 + 接触」⇒ 同一套机制、同一套代码模式。
3. **奇异消除的数学充分条件**：接触问题的牵引场由 **KKT 条件自洽解出**，
   在接触区内外**连续过渡**（`t_n → 0` 于接触边界），不产生**人为的牵引突降**。

**实现要点**（复用现有工具，不新写求解器）：
* `thums_feb.build_thums_feb(..., joint_bc="contact")` —— 新增分支，参照 `plantar_bc="contact"`：
  在 `subtalar_joint` 面**上方**造一块**数值刚体板（talus 替身）**，法向朝下，
  与补片构成接触对；把 `RigidForceLoad` / 总力施加在该板上（不是压在补片上）。
* 接触配方沿用已验收参数：`type='sliding-elastic'`、`laugon='PENALTY'`、
  `fric_coeff`（跖面实测 0.6；[Y25] 原文 μ=0.5 —— **两者都报**）、`contact_hold_spring_k` 作诊断。
* **判据、阈值（150 MPa）、材料、网格、后处理全部不动**；两帧都跑。

**预注册判据（跑前写死）**：
* **H1（主）**：接触路径下 `SC` 落到 **2–5** 区间 ⇒ 奇异消除成立。
* **H2（主）**：同 h 下 `σ_max` 相对 `anat_contact_hold100` 再降 ⇒ 首折抬高。
* **H3（否证）**：若 `SC` 仍在 10 以上或 `σ_max` 不降 ⇒ 奇异另有来源
  （如皮质/松质模量比 15000/73.4 = 204× 的界面），本方案需退回重做归因。

---

## §5 · 边界与未决

1. **未测**：`Achilles`（h=0 时 `achilles_kN = 0`）在更高高度/激活下的贡献未纳入本实验。
2. **未测**：网格收敛性（FEBio 无纯实体重网格器，`R2 风险`）。奇异判据的**充分性**因此
   依赖「SC 落到 2–5」而非「网格收敛」。
3. **未解**：`σ(h) ∝ h^0.281` 的**亚线性**来自载荷链（`F(h)` 本身亚线性），
   与 [Y25] 的能量法 `F ∝ h` 不同。**这是层①的另一条独立缺口**，本件不处理。
4. **不采纳**：任何以「换统计量」（max → p95/mean）进带的路线 —— 已由
   `V1V2_CALIBER_RULING.md` 裁定为口径 shopping。

---

## §6 · S5/S7 实施进展 —— **接触跑通了；但结果证伪了本件的修复假设（H3）**

**代码**：新增 `src/climbing/coupling/joint_contact.py`（`_nodal_normals` /
`_build_conforming_cap` / `_build_talus_cap`(平面) / `_apply_cap_hold_spring` /
`apply_joint_contact`）；`thums_feb.build_thums_feb` 新增 opt-in kwarg
`joint_bc="contact"`（**默认 `"pressure"` ⇒ deck 逐位不变**）、`joint_contact_mu` /
`joint_contact_gap_mm` / `joint_cap_ndiv` / `joint_cap_hold_k`；
新增 `_build_thums_joint_quads`（距下关节 quad4 facet 重建，**119 quad / 2496.3 mm²**，
与 `registration_anat.json` 的 2495.33 mm² 吻合）。

**修掉的 5 个真 bug**（每个都使求解失败或结果无效）

| # | 现象 | 根因 | 修法 |
|---|---|---|---|
| 1 | `A surface with name "subtalar_joint" was already defined` | builder 已建过同名 Surface | `create_bone_surface=False` 复用 |
| 2 | 初始化 `Negative jacobian` | 帽 hex8 节点序镜像 | 前 4 = y 大的一侧（同 `_build_pad_floor`） |
| 3 | 求解中单元负 Jacobian | 帽底 quad 法向算成 **+Y（朝内）** | 反序 ⇒ −Y（朝骨） |
| 4 | **BC 静默丢失**（帽 150 dof 全自由 ⇒ 位移 1e13） | `st.boundary = fbc.Boundary(all_bcs=bcs)` 在 `apply_joint_contact(..., bcs)` **之前** | joint 分支移到绑定**之前** |
| 5 | **接触"传力 ~0"**（σ_max 0.46–9.0 MPa） | **`hold_k` 太大**：平面帽 50×100=5000、共形帽 296×100=**29600 N/mm**，比骨侧接触路径更硬 ⇒ 帽坐在弹簧上、没压到骨 | 降到 **k=0.01 N/mm/节点**（总 ~3 N/mm，可忽略） |

> 第 5 条是**本轮最强的一条**：它同时解释了我先前的「平面/共形」误判 —— 几何换对了
> （`A_cap` 2493.9 vs 真值 2496.3 mm²）σ 仍 ≈0，因为**载荷被弹簧吃掉了**。

**结果（h=2 m，anatframe，F_sub = 22723.6 N，同一网格/材料/后处理）**

| 载荷路径 | `σ_max` (MPa) | `p_avg` (MPa) | **SC** | vs pressure |
|---|---|---|---|---|
| `joint_bc="pressure"`（现行 V1 口径） | **185.64** | 9.106 | 20.39 | 1.00 |
| `joint_bc="contact"`（平坦帽 + hold k=100） | 4.7–9.0 | 9.106 | 0.5–1.0 | ~0.03（**弹簧偷载荷，无效**） |
| `joint_bc="contact"`（共形帽 + hold k=100） | 0.46–1.30 | 9.106 | 0.05–0.14 | ~0.007（同上，无效） |
| **`joint_bc="contact"`（共形帽 + hold k=0.01）** | **454.43** | 9.106 | **49.90** | **2.45** |

### ⇒ **H3 触发：本件的修复假设被证伪**

`§4` 预注册的 **H2**（接触路径下 `σ_max` 再降）**不成立** —— 真实接触把峰值**抬高 2.45×**。
`§4` 的 **H1**（SC 落到 2–5）同样**不成立**（SC 升到 **49.90**）。

**物理上说得通，而且方向可预期**：均匀 `PressureLoad` 把力**摊平**在整个关节面上；
而**刚性距骨压在凹关节面上会把力集中在最高点** ⇒ **真实接触比均匀压力更狠**。

**⇒ 因此**：「补片边界奇异主导 V1 首折」这一诊断**被本实验证伪**。
`§3` 的观测（`σ_max` 不 ∝ `p_avg`）仍然成立，但它**不能**推出「换成接触就会降」。
**高 SC 不是载荷施加方式的产物，而是「骨网格分辨率 + 刚性/曲面配偶面」的固有属性。**

**⇒ 下一个归因对象（`§4` H3 已预注册）**：**皮质/松质模量比 15000/73.4 = 204× 的界面**，
以及**接触罚刚度的网格相关放大**（`penalty=0.1`、`auto_penalty=0`、`node_reloc=0`）。

**零回归已验**：`pytest -k "thums or plantar or contact or bone_feb or joint"`
⇒ **76 passed / 2 skipped**；`j0_pressure` 复现 **185.64 MPa** 逐位不变。

**⚠️ 交付边界**：上面 454.43 MPa 是在 **`hold_k=0.01` 这一窄窗口**内取得的
（`1e-3` 不收敛、`1.0` 被弹簧偷载荷）。**窗口窄本身是未解风险**，交付时必须一并报告。

---

## §7 · 溯源

| 项 | 位置 |
|---|---|
| A/B 原始数据 | `temp/opensim_fe/v1solve/ab_frame.json` |
| A/B 脚本 | `temp/e3.py`（草稿区） |
| 旧帧面积/面数记载 | `scripts/opensim_fe/registration_anatomical.py:11-14` |
| 解剖帧产物 | `temp/opensim_fe/thums_calcaneus/registration_anat.json`（`joint_face_area_mm2 = 2495.33`） |
| 默认帧 | `scripts/opensim_fe/thums_feb.py`（`DEFAULT_NPZ = calcaneus_r_calcnframe.npz`） |
| 载荷日志 | deck 内 `[thums_feb] 载荷：PressureLoad p=… （F=… / A=… mm²）` |
| 跖面接触先例 | `scripts/opensim_fe/thums_feb.py`（`plantar_bc="contact"` 分支） |
| [Y25] 接触口径 | `paper/MinerU_markdown_…_2106751303256973312.md` §1.1（面面接触 μ=0.5） |
