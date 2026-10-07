# Anderson 双论文关键信息（真值 / 已验证做法）

> 来源（MinerU 重解析）：
> - `paper/MinerU_markdown_anderson2007_2106695115320401920.md`（32,701 B）
>   Anderson DD, Goldsworthy JK, Li W, Rudert MJ, Tochigi Y, Brown TD.
>   **Physical validation of a patient-specific contact FE model of the ankle.**
>   *J Biomech* 2007;40(8):1662-1669. ｜ PMC1945165
> - `paper/MinerU_markdown_Intra-articular_Contact_Stress_Distributions_at_th_2106695331008294912.md`（38,362 B）
>   Anderson DD, Goldsworthy JK, Shivanna K, Grosland NM, Pedersen DR, Thomas TP,
>   Tochigi Y, Marsh JL, Brown TD.
>   **Intra-articular contact stress distributions at the ankle throughout stance phase
>   — patient-specific FE analysis as a metric of degeneration propensity.**
>   *Biomech Model Mechanobiol* 2006;5(2-3):82-89.
>
> 解析日：2026-10-04

---

## §1 两篇的关系

| | 2006 BMMB | 2007 J Biomech |
|---|---|---|
| 定位 | **方法学**（怎么建、怎么算整步态） | **物理验证**（模型 vs 尸体实测） |
| 标本 | 7 具正常踝 + 1 例胫骨 pilon 骨折 | **2 具尸体踝**（Tekscan 实测） |
| 载荷 | 13 步，**10 → 2800 N**，5° 跖屈~9° 背屈 | **600 N** 单次 |
| 软骨厚 | **1.5 mm** | **1.7 mm** |
| 求解器 | ABAQUS 6.4-1 | ABAQUS 6.6 |

---

## §2 🔴 与我模型直接冲突的四条（按影响排序）

### 2.1 软骨厚度：Anderson **1.5~1.7 mm**，我 **0.703 mm** ⚠️

> 2006: 「extruded along these point normals to define a zoned cartilage volume with
> user-specified thickness of **1.5 mm**」
> 2007: 「extruded along local bone surface normals at the tibio-talar articulation to a
> **uniform thickness of 1.7 mm**（Anderson et al., 2006）」

- 厚度依据：**El-Khoury et al. 2004**（Radiology 233(3):768-773，尸体踝软骨厚度测量）
- **✅ 已查清（2026-10-04，读 `temp/thums/cartilage.json` + `joint_gaps.json`）**：
  **我的软骨厚度中位数就是 1.5 mm**，与 Anderson 2006 完全一致！
  `cartilage.json` 实测：tibiotalar_tibia **median 1.5** / talus **median 1.5**；
  subtalar 两侧 median 1.49 / 1.5。
- ❌ **更正**：此前据「最小距 1.406 mm」推断「软骨做薄 2 倍」是**误读**。
  1.406 mm 只出现在**边缘 1~4 个节点**（`joint_gaps.json`：
  `n_within_1.5mm = 1`、`n_within_2mm = 4`）；**0.703 mm 是最小厚度不是中位厚度**
- **⇒ 唯一要改的**：`T_TARGET` **1.5 → 1.7**（对齐 Anderson 2007），一行改动
- ⚠️ **不需要重新摆位来"加厚软骨"**（此前的推断作废）

### 2.2 ⭐⭐⭐ 边界条件：Anderson 用**刚性骨 + reference node**，不是"固定距骨"

> **★★ 2006（`paper/MinerU_markdown_Intra-articular_...md`，§"Modeling constraints in
> whole duty cycle FE solutions"）逐字：**
> 「**Reference nodes, used to control translations and rotations of the assumed rigid
> tibial and talar subchondral surfaces**, were defined so as to coincide (in a loaded
> configuration) at the midpoint of the two talar condyles, along the provisional ankle
> flexion/extension axis.」
> 「the talus was **weakly restrained translationally by linear springs
> (stiffness = 500N/mm) connected to the ground**. It was otherwise allowed to
> **freely seat** relative to the tibia」
> 「to achieve stability during the **seating/pre-loading phase**, the talus was
> **fixed against rotation about the inversion-eversion and internal-external rotation
> axes**. These restraints were **removed**... so that the **tibio-talar articulation
> itself controlled talar rotation** thereafter.」
> 「Talar **displacement normal to the floor was constrained**, and talar
> **flexion/extension was prescribed** throughout the simulation.」
> **★ 落位（seating）协议，逐字：**
> 「In the **first provisional loading step**, the tibia was brought into **minimal
> contact** with the weakly constrained talus. To allow the ankle to **settle into its
> precise anatomic apposition**, guided by the tibio-talar articulation, **a 200 N
> compressive load was then applied** to the tibia in a **second provisional loading
> step** (after **releasing** restraints against talar rotation about the
> inversion-eversion and internal-external rotation axes).」
>
> **2007（`MinerU_markdown_anderson2007_...md`）逐字：**
> 「**a single-linear spring (k = 50 N/mm) placed to resist anterior–posterior
> translation**」（**代替软组织约束**，因纯轴向加载下相对平动很小）
> 实验夹具：「Flexion/extension rotations were **fixed**..., but **all other
> translational and rotational degrees of freedom were free**」

#### ⚠️ 口径澄清（**两个数字是两回事，不能混用**）

| | Anderson **2006** | Anderson **2007** |
|---|---|---|
| 弹簧 | **linear springs，500 N/mm，接地**，约束距骨**平动** | **单根** spring，**k=50 N/mm**，抗**前后(A–P)平动** |
| 用途 | 替代**缺失的韧带/后足**（距骨可自由落位） | 替代**软组织约束**（仅轴向加载场景） |
| 挂点 | 距骨**刚性**软骨下表面的 **reference node** | 同 |
| 转动 | 落位期固定**内外翻 + 内旋**两轴 ⇒ **duty cycle 前解除** | 屈伸固定 |
| 全程保留 | ① 距骨**垂直地面**位移约束 ② 距骨**屈伸** prescribed | 屈伸固定 |

#### ⭐⭐ 机制要点（**这决定了怎么实现**）

**Anderson 能"固定/解除内外翻与内旋"，是因为他用了 rigid 刚体 + reference node**
（6 个自由度，含 3 个**转动**）。
**⇒ 纯可变形实体网格没有转动自由度（`[LRN-20261005-065]`），因此这个约束在
可变形模型里根本无法表达。**

**⇒ 要忠实复现 Anderson 的边界条件，必须把距骨（至少其软骨下表面）设为
FEBio 的 `<rigid>` 刚体 + reference node**（FEBio 有 `<Rigid>` 段：
`rigid_bc` / `rigid_constraint` / `rigid_load`）。
**「只加平动弹簧」是不够的。**

**⇒ 对当前实现的判定**：`TALUS_SUPPORT=spring`（挂在 323 个距骨节点上的横向弹簧）
**语法是对的、机制是通的，但不是 Anderson 的口径** —— Anderson 是
「刚体参考点 + 500 N/mm 平动弹簧 + 转动约束按需解除」。
另：`FIX_TALUS=transverse`（距骨横向**硬**约束）**不对应** Anderson 的任何一条
（全程保留的只有"垂直地面位移"和"屈伸 prescribed"）。

**⇒ 物理依据（2007 实验）**：夹具只固定**屈伸**，**其余平动+转动全自由**
⇒ 内外翻/内旋自由是**物理事实**，不是数值技巧。


**⇒ 我的 `FIX_TALUS=all`（把距骨完全固定）是对已验证做法的重大偏离。**

### 2.3 ⭐⭐ 韧带：Anderson **完全不建韧带**，且实验证明「韧带不是主导」

> 2006: 「As **neither ligamentous structures nor distal structures of the hindfoot were
> explicitly included** in the computational model...」
> 2006: 「Absence of ligamentous restraint in the FE model is consistent with the results
> of **Tochigi et al. (2005a)** that showed that the **kinetics of normal ankle motion
> are dominated by the topography of articular surface contact**.」
> 2006: 「We recognize... that for many situations of ankle joint mechano-pathology,
> it **will be necessary to include ligament representation**.」
> 2007: 韧带不建，仅用 **单根 A-P 弹簧 k=50 N/mm** 代替软组织约束

**⇒ 直接结论：**
1. 用于**接触应力/面积验证**时，**不需要韧带**（Anderson 的模型没有韧带，仍然 1.5%~15% 精度）
2. **Tochigi et al. 2005a** 用实验证明：**正常踝运动的动力学由关节面接触地形主导，
   不是韧带** ⇒ 我在纯压缩下算出的韧带 ΔL≈0 是**符合这个结论的**
3. 吴恺的韧带位移排序 **没有实验验证**（Anderson/Tochigi 没测韧带位移排序）
4. ⚠️ 但 Anderson 也说：**病理/失稳情形必须建韧带**（我的攀岩跌落场景属于此类）

### 2.4 ⭐ 骨视为**刚体**，只让软骨变形

> 2006: 「**Reference nodes**, used to control translations and rotations of the assumed
> **rigid tibial and talar subchondral surfaces**」
> 2006: 「layers of continuum hexahedral elements zoned **outwardly (cartilage)** or
> **inwardly (subchondral bone)** from quadrilateral surface meshes」

**⇒ 整个关节的柔度 = 软骨柔度**（骨不变形）。

**粗算对照**（E=12 MPa、A≈300 mm²、两层共 3.4 mm）：
`k ≈ E·A/t = 12×300/3.4 ≈ 1059 N/mm`

**我的实测刚度**：
- `tet4` @DISP 0.1472 → 259.3/0.1472 ≈ **1762 N/mm**
- `penta6` @DISP 0.1472 → 116.9/0.1472 ≈ **794 N/mm**

**⇒ Anderson 估算值 1059 N/mm 落在我的两个值之间（偏 penta6 侧）。**
（注意：我的软骨更薄 1.406 mm，且骨可变形，不能直接等价换算 —— 仅作量级参照。）

---

## §3 🔧 三条可直接改进模型的做法

### 3.1 ⭐⭐ 全**六面体**网格（Anderson 的核心做法）

> 2006: 「the intervening volume is then meshed with continuum **hexahedral elements**」；
> 「layers of continuum hexahedral elements zoned outwardly (cartilage) or inwardly
> (subchondral bone) from **quadrilateral surface meshes**」
> 2007: 「generate an **all-hexahedral** FE mesh」（TrueGrid）

**生成流程**：三角面 → 投影平面上撒种子网格 → 沿射线投到三角面 → **四边形面网格**
→ 算每节点法向（相邻四边形外法向平均）→ **沿法向挤出**成软骨砖层。

**⇒ 这直接指向我的 2.2 倍问题：软骨应该是「沿法向挤出的六面体层」，
不是 tet4/penta6 任一种。** 六面体既无 tet4 的体积锁死，也无 penta6 的沙漏。

### 3.2 接触：主从面 + **单元尺寸比 1:3**

> 「the **tibia was assigned as the master surface**, and an element length ratio of
> roughly **1:3** was maintained between the tibial and talar meshes」
> （目的是「ensure robust performance of the master/slave contact paradigm」）

**⇒ 两侧单元尺寸应刻意拉开约 3 倍。** 我的模型两侧可能尺寸相近 ⇒ 接触抖动。

### 3.3 中和负重位姿需**主动摆位**（不是一个姿态直接用）

> 2006: 用 Data Manager 交互摆位；「A local coordinate reference frame, centered within
> the talus... A provisional ankle flexion/extension axis was defined along a line
> connecting the **centers of circles fitted to the condylar arcs of the talar dome**
> (Bottlang et al. 1999)... aligning the provisional flexion/extension axis along the
> global x direction, and rotating about this axis to align the **primary axis of the
> talus with the first metatarsal, 15° below horizontal**. With the talus thus fixed,
> the tibia was rotated about the provisional ankle flexion/extension axis so that the
> **angle between the tibial shaft and the floor was 85°**.」

**⇒ 我直接用 THUMS 原始姿态 ⇒ 关节间隙/相对位姿可能不是负重位。**
这条与 §2.1 互为因果（间隙 1.406 mm 可能就是因为没摆位）。

---

## §4 数值与设置速查（可直接抄）

| 项 | 值 |
|---|---|
| 软骨 E / ν | **12 MPa / 0.42**（两篇一致） |
| 摩擦 | Coulomb **μ = 0.01**（Linn 1967） |
| 软骨厚度 | 1.5 mm（2006）→ **1.7 mm**（2007） |
| 网格 | **全六面体**；收敛需 **约 50,000 单元**（扫了 2,997~293,402，13 次） |
| 主从比 | 胫骨=master，单元尺寸比 **1:3** |
| 平滑 | 圆形平均滤波，半径覆盖 **6 mm**；2007 说「宁少勿多」 |
| 距骨约束 | 弱弹簧 **500 N/mm**（2006）/ A-P 弹簧 **50 N/mm**（2007） |
| 预载流程 | 轻接触 → **200 N** 落位（松开内翻/内旋约束） → 正式加载 |
| 载荷 | 2006: 13 步 **10→2800 N** + 5°跖~9°背屈；2007: **600 N** |
| 计算成本 | 50,000 单元 ≈ **90 min**（HP 4-CPU Itanium） |

---

## §5 🎯 新的验证靶子 / 可下载文献

| 靶子 | 内容 | 价值 |
|---|---|---|
| **Tochigi et al. 2005a** | *Tensile engagement of the peri-ankle ligaments in stance phase.* Foot Ankle Int 26(12):1067-1073 | ⭐⭐ **韧带实验真值** —— 直击我的韧带排序问题 |
| El-Khoury et al. 2004 | *Cartilage thickness in cadaveric ankles.* Radiology 233(3):768-773 | 软骨 1.5~1.7 mm 的依据 |
| Bottlang et al. 1999 | J Biomech 32(1):63-70 | **踝关节屈伸轴 / 螺旋位移轴（SDA）** |
| 2006 峰值应力 | 正常踝 **9~14 MPa**（载荷到 2800 N）；骨折例 18 MPa | 高载荷下的量级参照 |
| **SDA 验证** | 模型算出的**螺旋位移轴**应与 Bottlang 1999 的实验一致（穿过距骨体、角行程相近） | ⭐ 运动学验证指标（我还没用过） |
| 精度声明 | 「conservatively, **within 10-15%** of those measured experimentally」 | 我的 −14% 正好在这个区间内 ✓ |

---

## §6 对我「关键问题」的影响（结论）

| 原判断 | 更新 |
|---|---|
| 「2 倍刚度差是主瓶颈，先交换接触面」 | **方向改为：软骨改六面体层（§3.1）+ 骨可刚性化（§2.4）** |
| 「松开顶面夹持做实验 1」 | **改为：距骨改弱弹簧约束 + 解除内翻/内旋约束（§2.2）**，这正是 Anderson 的做法 |
| 「韧带排序要复现吴恺」 | **吴恺的韧带排序无实验支撑；Anderson/Tochigi 说正常踝不由韧带主导** ⇒ 降级 |
| 「软骨 0.703 mm」 | **疑似偏薄 2 倍，需先量关节间隙分布（§2.1）** |
| 「我的 −14% 是缺陷」 | **Anderson 自己声明精度就是 10~15%** ⇒ 我的偏差在正常范围 |

---

# §7 ⭐⭐ 补：Tochigi et al. 2006 FAI —— 韧带**实验真值**（直击我的韧带问题）

> 来源：`paper/MinerU_markdown_2006FAIankleligaments_2106701085136199680.md`
> Tochigi Y, Rudert MJ, Amendola A, Brown TD, Saltzman CL.
> **Tensile Engagement of the Peri-Ankle Ligaments in Stance Phase.**
> *Foot Ankle Int* 2006;26(12):1067-1073.（= Anderson 引的「Tochigi et al. 2005a」）

## 7.1 方法

- **11 具新鲜冰冻尸体踝**（平均 71 岁，57~85）
- **DVRT 微型位移传感器**直接缝在韧带中段，连续测**应变**
- **6 条韧带**：ATFL（距腓前）、CFL（跟腓）、PTFL（距腓后）、
  以及浅层三角韧带三束 **ADL / MDL / PDL**
- 动态：**600 N 轴向载荷保持**，踝在 **15° 跖屈 ↔ 10° 背屈**间摆动，1 Hz，4 周期
- 扩展：其中 8 具另做准静态 **30° 跖屈 / 30° 背屈**（载荷降到 300 N）

## 7.2 🎯 核心结论（原文）

> 「In the dynamic loading tests, **none of the ligaments monitored showed a reproducible
> strain pattern** indicating a role in ankle stabilization.」
>
> 「**none of the ligaments studied were reproducibly recruited to be a primary
> stabilizing structure.** The peri-ankle ligaments are likely to be **secondary
> restraining structures** that serve to resist motion to avoid extreme positions.」
>
> 「**Stance phase ankle motion appears to be primarily controlled by articular
> congruity, NOT by peri-ankle ligament tension.**」

**「持续松弛（slack）贯穿整个运动弧」的标本数**：
ATFL **7/11**、CFL 4/11、PTFL 5/11、ADL 5/9、**PDL 9/9**
（只有 MDL 在部分标本持续紧张）

**站姿运动弧内的应变范围（表 1）**：
| 韧带 | 范围 (%) |
|---|---|
| ATFL | 0.2 ~ 6.2 |
| **CFL** | **0.2 ~ 2.1** |
| PTFL | 0.3 ~ 5.8 |
| ADL | 5.0 ~ 15.2 |
| MDL | 0.6 ~ 9.0 |
| PDL | 0.1 ~ 16.4 |

「近等长（范围 ≤2%）」很常见：CFL **9/11**、PTFL 7/11、MDL 7/11、ATFL 5/11。

## 7.3 🔥🔥 这直接**反转**我对自己韧带结果的判定

| 事实 | 出处 |
|---|---|
| 我的韧带 ΔL = **0.001~0.08 mm**（≈0.004~0.3% 应变 @25 mm 韧带） | 我的模型 |
| CFL 站姿应变范围 **0.2~2.1%**；多条韧带**整个弧内持续松弛** | **Tochigi 实测** |
| 吴恺韧带位移 **1.92~5.24 mm** ⇒ 换算 **≈8~21% 应变** | 吴恺论文 |

**⇒ 两条结论：**
1. **我的「ΔL≈0」不是负结果，而是与实验一致** —— 站姿相下大部分踝周韧带
   本来就是**松弛/近等长**的（<1% 应变的韧带**几乎不产生力**，Toe 区）
2. **吴恺的 1.92~5.24 mm 与实验严重不符** —— 实验测到的**整个站姿弧最大**才
   2.1%(CFL) / 6.2%(ATFL)，吴恺的数值隐含 8~21%
   ⇒ **它的韧带位移不能当靶子，这次有了实验背书**（此前只是"物理上不可能"的推断）

## 7.4 ⭐⭐ 边界条件**第三方确证**：内翻/外翻、内外旋、前后位移**都不约束**

> 原文：「Simultaneously, **inversion and eversion, internal and external rotation,
> and anterior and posterior translation were unconstrained** to reproduce natural
> ankle motion.」
> 图 1 注：「**Inversion and eversion of the foot are unconstrained**... the axial load
> plate incorporates a linear or rotational bearing that allows **free anterior and
> posterior translation and internal and external rotation of the tibia**.」

**⇒ 与 Anderson 2006（500 N/mm 弱弹簧 + 撤销内翻/内旋约束）、
Anderson 2007（A-P 弹簧 50 N/mm）**完全一致**：
三次独立来源都证明「**距骨/胫骨的内外翻、内外旋、前后位移必须放开**」。

**⇒ 我的 `FIX_TALUS=all`（完全固定距骨）被三方证据判定为错误边界条件。**

## 7.5 ⭐ 韧带只在**极端角度**才真正受力

扩展测试（30° 跖屈 / 30° 背屈）：
- 30° **跖屈**：ATFL 3/8 紧张、CFL/PTFL 2/8、**ADL 8/8 全部紧张**
- 30° **背屈**：**CFL/PTFL 6/8 紧张**、MDL 2/8、PDL 4/8

> 「Those ligaments therefore seem likely to function primarily in
> **restraining extreme flexion of the ankle**.」

**⇒ 韧带是「极端位的次级约束」，正常站姿弧（15°跖~10°背）内不主导。**

**⇒ 对我的意义（重要）**：
- 复现吴恺那种「正常站姿 600 N」的场景 ⇒ **韧带本来就不该是主导**，
  我的 ΔL≈0 合理
- **但我的真实目标是攀岩跌落的损伤场景（极端内外翻）** ⇒ 那正是韧带**真正受力**的区间
- ⇒ **韧带线不该废掉，而是换场景：从「复现站姿」改为「极端位/损伤位」**

## 7.6 其他可用的点

| 项 | 内容 |
|---|---|
| 1% 阈值 | 应变 >1% 记为 taut、<1% 记 slack；**<1% 应变的韧带刚度极低（Toe 区），几乎不产生稳定力** |
| 装载 | 胫腓骨中段 / 跟骨 / 末节趾骨 **三块 PMMA 分别固定**；测试轴对齐踝-距下复合体以减小扭矩 |
| 600 N | 」one body weight (600N)」 |
| 局限（作者自己说） | **没包含剪切与扭转力**（体内站姿存在，但远小于轴力）⇒ 吴恺引的「0~0.5° 内旋」这类小效应**他们没测** |
| Leardini 四杆模型对照 | Leardini 认为胫跟(MDL)+CFL 是主控；但那是**无关节面接触加载**的实验。**有接触加载时，两条虽常近等长，却都不被持续招募为主控** ⇒ **关节面负载后接管** |

---

# §8 🔴 新首要嫌疑：**接触区严重欠分辨**（P0 实测，2026-10-04）

## 8.1 实测数据（`temp/thums/joint_gaps.json` + `cartilage.json`）

| 项 | 我的模型 | Anderson |
|---|---|---|
| 胫距关节面节点数（tibia 侧） | **73** | — |
| 胫距关节面节点数（talus 侧） | **89** | — |
| 胫骨软骨单元 | **93** | — |
| 距骨软骨单元 | **129** | — |
| 全模型单元 | **26,794** | 收敛需 **约 50,000** |
| 接触区节点间距（估） | **≈2 mm** | — |
| 网格收敛研究 | **没做过** | **13 次，2,997→293,402 单元** |

`joint_gaps.json` 里 tibia__talus 的分布（分母是 tibia **全部 4836 个表面节点**）：
`min 1.406 / p1 3.66 / p5 9.61 / median 223.8 mm`
**⇒ 中位数毫无意义（大部分节点在胫骨干上）；有意义的只有：**
`n_within_1.5mm = 1`、`n_within_2mm = 4`、`n_within_3mm = 16`、`n_within_5mm = 54~114`

## 8.2 为什么这可能是 2.2 倍刚度差的根因

1. **接触斑只有 ~220 个软骨单元**（93+129）覆盖约 300 mm² 的接触面
   ⇒ 单元尺度 ≈ 接触斑尺度的 1/10，**接触压力分布根本分辨不出来**
2. 用 **master/slave 接触**时，从面节点少 ⇒ 主面被"点接触"离散化
   ⇒ 有效接触刚度强烈依赖离散 ⇒ **tet4 与 penta6 的接触面片数不同
   （148/191 vs 138/171）可能正是差异来源**
3. **Anderson 专门做了 13 点网格收敛研究**才敢用 —— 我一次都没做
4. 解释了为什么 `CART_QUALITY=1` 无效（`LRN-020`）：它改善的是单元**形状**，
   **没有增加单元数量**

## 8.2b 🔥 P1 诊断（2026-10-04 实测）：**hex8 曾失败的真因已定位**

**基准复现**（`temp/pyfebio_demo/pconv.py`，逐字一致）：

| 软骨单元 | rc | F(N) | 保留/剔除 | 合计 | 相对 |
|---|---|---|---|---|---|
| **tet4** | 0 | **259.3** | 456+596 / 6+16 | **1074** | 基准 |
| **penta6** | 0 | **116.9** | 140+175 / 14+29 | **358** | ×0.45 |
| penta15 | 0 | 49.3 | 104+128 / 50+76 | 232 | ×0.19 |

**关键结构发现**：`1074 / 358 = 3.00` 整
⇒ **tet4 路线 = 每个楔形拆成 3 个四面体**，与 penta6 **同网格、同节点、同体积**。
⇒ 两者的差别**不是网格密度**，而是**同一楔形的两种插值**。

**但 patch test 早已证明二者等价**（单楔形 vs 解析解，ν=0 下均 9.702 N；
薄歪斜扫描 h 0.2~5 mm / skew 0~1.5 mm 5 种情形**比值全 1.000**）。
⇒ **差别来自装配层，不是单元公式。**

### ⭐ 从代码注释里定位到的真因（`ankle_lig_feb.py:163-167` 与 `:253-256` 拼起来）

> `:163-167`「…**实测直接用 hex8 会因偏移方向与面不垂直而歪斜（det(J)≈0.5），
> 被 DET_MIN 闸门剔除 42/61 个 ⇒ 软骨破洞、接触面从 150/197 掉到 51/64**。
> 拆楔形后 0 剔除。」
> `:253-256`「★ 偏移方向用「指向对侧骨表面」的单位向量（dirA），
> **不用面积加权节点法向**。实测：**粗网格上相邻节点的面积加权法向严重不一致
> （相加几乎抵消）**，会让 offset 面几乎与 base 面重合 → 单元退化」

**⇒ hex8 失败的两个真因：**
1. 底网格是**三角形**（THUMS STL）⇒ 顶点法向互相抵消，无法得到一致法向
2. 偏移方向取「指向对侧最近点」⇒**不垂直于面** ⇒ hex 歪斜 ⇒ 被闸门剔除

### ✅ 而 Anderson 的配方恰好逐条对症（2006「Finite element mesh generation」）

> ① 生成**四边形**面网格（mesh projection plane 撒种子 → 投影到三角面）
> ② 算每节点法向：**(i) 找共享该节点的四边形 → (ii) 算各自外法向 → (iii) 平均**
> ③ **沿这些点法向挤出**成软骨层 ⇒ **all-hex**（TrueGrid）

**⇒ 结论：hex8 不是「不行」，是「三角形底 + 非垂直偏移」不行。**
只要按 Anderson 做「**四边形底 + 顶点法向挤出**」，就能拿到干净的 hex 层，
从根上消掉 tet4 vs penta6 的二选一。

## 8.3 ⇒ 下一步（替代原「交换接触面」）

**A. 软骨改「四边形底 + 顶点法向挤出的 hex8 层」**（= Anderson 配方，§8.2b）
   —— 消掉 tet4/penta6 二选一，并顺带解决接触面定向问题
**B. 做网格收敛研究**：逐级加密胫距接触区，看
   ①应力/面积/支反力是否收敛 ②各单元类型是否收敛到同一值
**C.（兜底）交换接触面 / 逐单元贡献对比**

---

# §9 P1 诊断实测（2026-10-04）——「关节面太粗」是根问题

探针：`temp/pyfebio_demo/norm_probe.py`（只读几何，不改模型）

## 9.1 三条硬结论（其中两条**推翻旧注释**）

| # | 结论 | 证据 |
|---|---|---|
| 1 | ✅ 底网格**已是四边形为主** | tibia 外表面 19,022 面 = **tri 4,832 / quad 14,190（74.6%）**；talus 75.0% |
| 2 | ✅ **顶点法向完全可用，没有"相加抵消"** | 合成模长中位 **31.28 ≈ 12 面 × 2.61 mm²**（= Σ面积）⇒ 法向几乎完美对齐。**旧注释「面积加权法向严重不一致、相加几乎抵消」经实测为误判** |
| 3 | 🔴 `gap<6mm` 是**纯距离**判据 ⇒ 关节面选择**严重掺假** | 73 个"关节面"节点中：夹角>60° 占 **40/73（55%）**；**>90°（法向反向）占 26/73（36%）**，p10 cos = **−0.91** |

**偏移方向实测差**（同节点、同 1.5 mm 厚度）：
中位 **1.587 mm** / p90 **2.932 mm** ⇒ p90 ≈ 2×1.5 ⇒ **法向与"指向对侧"在 1/3 节点上近乎相反**
⇒ 旧代码用"指向对侧"是**为了躲开反向节点**，属于**对的妥协**，但代价是偏移不垂直于面 ⇒ hex8 歪斜

## 9.2 但「按法向对齐过滤」**行不通**（实测）

新增开关 `CART_OFFSET=normal` + `ART_ALIGN=0.5` 后重跑 `cartilage.py`：

| 关节 | 关节面节点 | 软骨单元 | 体积 |
|---|---|---|---|
| tibiotalar_tibia | 73 → **33** | 93 → **18** | 1.057 → 0.263 cm³ |
| tibiotalar_talus | 89 → **36** | 129 → **27** | 1.109 → 0.486 cm³ |
| subtalar_talus | 212 → 92 | 397 → 95 | 3.059 → 1.132 cm³ |
| subtalar_calcaneus | 253 → 119 | 500 → 149 | 2.706 → 1.285 cm³ |

**⇒ 过滤把关节面砍掉一半，软骨只剩 18/27 个单元**（原 93/129）；
且退化丢弃率仍约 33%（18 保留 / 9 丢弃）。**这不是可用状态，已回滚。**

## 9.3 🎯 结论：**根问题是「关节面网格太粗」，不是「偏移方向选错」**

`GAP_MAX` 敏感性（`cartilage.py` 自带输出）：

| bone | <3mm | <4mm | **<6mm** | <8mm | <10mm |
|---|---|---|---|---|---|
| tibia | 16 | 31 | **72** | 103 | 135 |
| talus | 18 | 38 | **87** | 126 | 165 |

- 胫骨下关节面约 300~500 mm²，却只由 **72 个节点**代表 ⇒ 节点间距 **≈2.0~2.6 mm**
- Anderson 收敛需 **~50,000 单元**、并专门做了 **13 点收敛研究**
- ⇒ **必须先把关节面加密**，才谈得上"用法向挤 hex8"

## 9.4 ⇒ P1b 正确顺序（修正）

1. **先加密胫距（及距下）关节面的骨表面网格**（THUMS hex8 面级细分 or 局部重划）
   —— 目标：关节面节点数从 ~72 提到 **≥300**（间距 ~1 mm 级）
2. **再**按 Anderson：四边形底 + 顶点法向挤出 ⇒ **hex8 层**（此时反向节点会自然减少）
3. 用 `CART_OFFSET=normal` / `ART_ALIGN` 已就位的开关做 A/B
4. 最后做网格收敛研究（逐级加密，看 tet4/penta6/hex8 是否收敛到同一值）

**已就位的工具**（本轮新增，默认关闭不影响现有流程）：
`CART_OFFSET=toward|normal`、`ART_ALIGN=<cos 阈值>`（`scripts/ankle_fe/cartilage.py`）

## 9.5 ⚠️ 操作提醒

本轮曾用 `CART_OFFSET=normal ART_ALIGN=0.5` 跑 `cartilage.py`，
**覆盖了默认的 `cartilage.json` / `cartilage_mesh.npz`**。
**已立即用默认参数重跑恢复**（`tibiotalar_tibia 93 / 1.057cm³`、
`talus 129 / 1.109cm³`），并跑通回归 `FEBio rc=0 / 17.8 s`、接触面 `152/200`。
⇒ **教训：`cartilage.py` 直接写死输出路径，试参数前应先备份或改用独立 OUT。**

---

# §10 P1c 探查（2026-10-04）——骨网格构成 + 法向偏移**反证**

## 10.1 🔴 脚本混淆坑（浪费了一次实验）

**项目里有**两个**组装脚本，功能不同：**

| 脚本 | 用途 | 是否有 `CART_ELEM` 等开关 |
|---|---|---|
| `scripts/ankle_fe/build_feb.py` | A1–A4 冒烟组装（**定死** tibiotalar 单关节） | ❌ **没有** |
| `scripts/ankle_fe/ankle_lig_feb.py` | 全功能（韧带/截断/单元类型/驱动方式） | ✅ 有 |

`pconv.py` 跑的是 `ankle_lig_feb.py`；而 `FE_PIPELINE_TODO.md` 的「一键复现」写的是
`build_feb.py`。**我改了前者却跑了后者**，导致 `CART_ELEM=hex8` 完全没生效
（日志里软骨仍是 `tet4`），白跑一轮。

**⇒ 教训：改脚本前先确认**运行路径**用的是哪个文件；两个同名功能的脚本必须互相标注。**

## 10.2 骨网格构成（`npz_probe.py` 实测）

| part | 单元数 | 唯一节点数 | 实际类型 |
|---|---|---|---|
| `R_TIBIA_CORT` (81000700) | 2,417 | 8（2 个为 6） | **hex8** |
| `right_tibia_end_spon` (81000600) | 6,651 | **4** | **tet4** |
| `right_tibia_center_spon` (81000601) | 14,482 | **4** | **tet4** |
| `R_FIBULA_CORT` (81000900) | 842 | 8 | hex8 |
| `R_TALUS_CORT` (81001100) | 394 | 8 | hex8 |
| `R_TALUS_SPON` (81001000) | 1,766 | **4** | tet4 |
| `R_CALCANEUS_CORT` (81001300) | 634 | 8 | hex8 |
| `R_CALCANEUS_SPON` (81001200) | 3,323 | **4** | tet4 |

**⇒ 皮质骨 = hex8，松质骨 = tet4**（THUMS 的 `*ELEMENT_SOLID` 变长行被零填充成 8 列）。
**踝部合计 33,781 单元，其中约 87% 是 tet4。**

**这解释了 §9.1 里 25.4% 三角面的来源**（松质骨的四面体面）。

**节点共享（`share_probe.py`）：**
- **皮质 ↔ 松质：重合 50%**（tibia 2,418；fibula 844；talus 396；calc 636）
  ⇒ **共形界面，不能单独加密一侧**
- **各骨之间：0 重合** ⇒ 各骨是独立体（靠接触/韧带交互）✓

## 10.3 🔴 「沿顶点法向挤出」**反证：比默认更差**

受控对比（`ctrl_cmp.py`，用**正确的** `ankle_lig_feb.py`，
同一份 `CART_OFFSET=normal` 软骨，env 与基准一致）：

| `CART_ELEM` | rc | 用时 | 软骨 保留/剔除 | 接触面 |
|---|---|---|---|---|
| `tet4` | **1** ❌ | 9.0 s | 406/8、610/20 | 129/205 |
| `penta6` | 0 | 4.9 s | 88/**50**、146/**64** | 85/146 |
| **`hex8`** | 0 | 4.4 s | hex8 **仅 6/53**；penta6 25/7 | 31/60 |

**⇒ 三条结论：**
1. **`CART_OFFSET=normal` 的剔除率远高于默认 `toward`**
   （penta6 剔 50/88 = **57%**，默认约 14%；hex8 剔 **89%**）
2. **tet4 在法向偏移下直接失败**（负 Jacobian）—— 默认偏移下 tet4 反而能跑
3. **原代码注释是对的**：「直接用 hex8 会因偏移方向与面不垂直而歪斜」
   —— 我 §9.1 的探针只证明了**法向彼此对齐**（`mag ≈ Σ面积`），
   **不等于法向适合当挤出方向**。边缘节点的法向与「指向对侧」夹角可达 90°+（§9.1 的 36% 反向）

**⇒ 机制**：楔形（penta6）对倾斜**宽容**（偏移三角面只是平移，仍成楔形）；
而 hex8 要求偏移面**平行于**底面，一旦倾斜就退化成歪斜砖 ⇒ 被闸门剔除。

**⇒ 所以「先加密关节面」不是可选项，而是 hex8 路线的**前提**：**
只有把关节面做成**筛选干净、法向一致、足够细**的曲面片，
沿法向挤出才可能得到干净的 hex 层。

## 10.4 本轮收尾

- 已用备份恢复默认软骨（`cartilage.json` 93/129/397/500，与原值一致）
- 回归：`build_feb.py` → **`FEBio rc=0 / 17.2 s / 接触面 152/200`** ✓
- 新增只读探针：`npz_probe.py`、`share_probe.py`、`ctrl_cmp.py`、
  `norm_probe.py`（均在 `temp/pyfebio_demo/`）

---

# §11 P3 受控结果（2026-10-04）——**2.2 倍与材料参数无关**

脚本：`temp/pyfebio_demo/cmp_lit2.py`（4 组材料 × 2 种单元，各 `rc=0`）
固定：`DISP_MM=0.1472 / TRUNC_MM=60 / FIX_TALUS=all / TIME_STEPS=20 / MATSET=WK`

| 软骨参数 | tet4 F(N) | penta6 F(N) | **比值** |
|---|---|---|---|
| 基准（WK 10 / 0.40） | 259.3 | 116.9 | **2.22** |
| **Anderson（12 / 0.42）** | **309.5** | **131.8** | **2.35** |
| 仅 ν=0.42 | 267.2 | 118.2 | 2.26 |
| 仅 E=12 | 300.2 | 130.3 | 2.30 |

## 11.1 三条结论

1. ✅ **2.2 倍比值与材料参数基本无关**（2.22 → 2.35，波动仅 6%）
   ⇒ **不是材料假象**，是**离散/几何**层面的问题
   ⇒ **改用 Anderson 的 E=12/ν=0.42 不会解决 2.2 倍**（两种单元等比例抬高）
2. ✅ **ν 0.40 → 0.42 影响很小**（tet4 +3.0%、penta6 +1.1%）
   ⇒ 此前对「体积锁死」的担忧在**这个 ν 区间**不显著（与 `LRN-032` 的 ν 扫描一致：
   ν=0 时比值就已 ~2.0）
3. ⭐ **F 随 E 次线性增长** ⇒ 反解出各部件柔度占比

## 11.2 ⭐ 由次线性反解「骨柔度占比」

把模型看成「软骨串骨」两弹簧串联：`1/F = 1/(a·E) + 1/b`

由 E=10 → F=259.3 与 E=12 → F=300.2 解得：

| 量 | 值（同一比例常数下） |
|---|---|
| 软骨串刚 `a·E`（E=10） | **317** |
| **骨刚 `b`** | **1425** |
| ⇒ **骨贡献的柔度占比** | **317/1425 ≈ 22%** |

**⇒ 两个推论：**
- **细观到 `Anderson` 把骨视为刚体** ⇒ 对齐他后模型会**再硬约 28%**
  （`1/(1/317 + 0) = 317` vs 现在 259.3 ⇒ 317/259.3 = **1.22×**）
- 这也说明**我的模型里骨不是可以忽略的**：软骨与骨的柔度是同一量级

## 11.3 对 2.2 倍问题的意义

**排除了材料维度**，剩余的离散嫌疑集中在（按可能性）：
1. **接触离散**（接触斑仅 ~220 个软骨单元，主从面片数 148/191 vs 138/171）
2. **单元剔除**（tet4 剔 2%、penta6 剔 14%）
3. **3-tet 分解 vs 楔形**的插值差异（patch test 在 ν=0 单元素上等价，但整模不等价）

⇒ **仍指向 P1d：先把关节面加密，再谈哪个单元类型对。**

---

# §12 P1d 加密（2026-10-04）——**机制成功，但求解失败**

## 12.1 ✅ 加密机制работает：关节面节点 73 → 692（9.5×）

脚本 `temp/pyfebio_demo/refine3.py`（1 级红色加密全 body）：
- `hex8 → 8 hex8`（棱中点 + 面心 + 体心）
- `tet4 → 8 tet4`（棱中点）
- `quad4 壳 → 4 quad4`
- 节点 22,960 → **219,717**（新增 棱 137,730 / 面 50,640 / 体 8,387）

**在加密网格上重跑 `cartilage.py`（`THUMS_NPZ` / `CART_OUT` 指向 r1）：**

| | 原始 | 加密 R1 | 倍数 |
|---|---|---|---|
| tibiotalar 胫侧关节面节点 | **73** | **692** | **9.5×** |
| talus 侧 | 89 | 887 | 10× |
| 软骨单元（胫/距） | 93 / 129 | **894 / 1,256** | 9.6× |
| 软骨体积 | 1.057 / 1.109 cm³ | **2.414 / 2.573 cm³** | **2.3×** |

**⇒ ⭐ 关键发现：细网格下软骨体积涨 2.3 倍** ——
说明**粗网格严重低估了关节面的真实面积**（细网格才能捕捉曲率）。
这与「接触区欠分辨」的判断一致，**是本轮最有价值的产出**。

## 12.2 ⚠️ 但加密网格**求解失败**（`rc=1`，0.1 s）

```
tag "Elements" (line 49104) : missing attribute "id"
  → <Elements type="hex8" name="talus_cort__hex8"></Elements>  ← 空块
```
**骨与软骨单元全部被 `DET_MIN` 闸门剔除**（`软骨 []`），
说明加密后的子单元 **detJ 为负**（朝向不对）。

**已尝试的修法**：在 `refine3.py` 里加 `fix_hex`（子单元与母单元都做
「交换底面/顶面 4,5,6,7,0,1,2,3」的朝向修正）——**未能解决**。
⇒ 说明问题**不只是底/顶交换**能修的朝向，更可能是**母单元节点顺序不符合 FEBio 约定**。

## 12.3 📐 支撑上面判断的旁证

用两种体积度量测**同一**网格：

| 度量 | 原始 | 加密后 | 相对差 |
|---|---|---|---|
| 6-tet 扇形分解 | 890,100 | 880,044 | −1.13% |
| **三线性 Jacobian Gauss 积分** | **980,326** | 887,871 | −9.43% |

**⇒ 两种度量在【原始】网格上就差 10%**（980k vs 890k）。
扇形分解对**扭曲**六面体不精确；两者不一致 ⇒
**母六面体存在明显扭曲**，我的「参数化 → 棱/面中点」假设不一定成立。

**⇒ P1d 的正确收尾**：
先**逐单元检查 npz 里 hex 的节点顺序约定**（对照 THUMS/LS-DYNA 与 FEBio 的差异），
再决定是「规范化顺序后重加密」还是「改用 tied-elastic 把细软骨绑到粗骨上」（§3.14.5）。

## 12.4 🔴 操作坑（本轮代价最大）

我在**同一条 cmd 命令**里先 `set "THUMS_NPZ=...r1.npz"` 再跑 `cartilage.py`，
**随后**又在同一条命令里跑 `cmp_ref.py` —— **`set` 污染了整个 shell**，
而 `cmp_ref.py` 里是 `if npz: env[...] = npz` 的写法，
**"原始"组继承到了 r1 ⇒ 四个 run 全在跑加密网格**，连原始组也失败，
一度误判为「我的改动弄坏了模型」。

**⇒ 修法**：对比脚本里必须**无条件** `env.pop("THUMS_NPZ")` 再按需设置。
**⇒ 通用教训**：cmd 的 `set` 作用域是**整个命令行**，不要在同一条命令里
既改 env 又跑对比/回归。
