"""Tochigi 2006 FAI 韧带论文：接进 Anderson 文档 + 更新 CONV 降级项。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")

# ── ① 追加到 docs/Anderson双论文_关键信息.md ──
A = ROOT / "docs/Anderson双论文_关键信息.md"
s = A.read_text(encoding="utf-8")
ADD = """

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
"""
A.write_text(s.rstrip() + ADD, encoding="utf-8")
print(f"Anderson双论文_关键信息.md -> {A.stat().st_size:,} bytes")

# ── ② CONV_TODO.md 降级项替换为更硬的结论 ──
T = ROOT / "CONV_TODO.md"
t = T.read_text(encoding="utf-8")
old = """**降级项**：韧带排序复现吴恺 —— Anderson/Tochigi 2005a 用实验证明
「正常踝运动动力学由关节面接触地形主导，**不是韧带**」，
且 Anderson 的验证模型**完全不建韧带**仍达 1.5~15% 精度（§2.3）。"""
new = """**降级项（已由实验证据结案）**：韧带排序复现吴恺 ——
- Anderson 的验证模型**完全不建韧带**仍达 1.5~15% 精度
- **Tochigi 2006 FAI 实测**（11 具尸体踝，DVRT）：站姿相下**大部分踝周韧带持续松弛**，
  CFL 应变范围仅 **0.2~2.1%**、ATFL 0.2~6.2%，**无任何韧带可被稳定招募为主控**；
  「站姿踝运动主要由**关节面契合度**控制，不是韧带」
- ⇒ **我的 ΔL≈0.001~0.08 mm 与实验一致，不是失败**；
  **吴恺的 1.92~5.24 mm 隐含 8~21% 应变，与实验严重不符** ⇒ 不可作靶子
- ⇒ **韧带只在极端角度（30°跖/背屈）才真正受力**（CFL/PTFL 6/8 紧张）
  ⇒ **韧带线不废，改场景**：从「复现正常站姿」改为「**极端位/损伤位**」
  —— 这正是攀岩跌落的目标场景
- 详见 `docs/Anderson双论文_关键信息.md` §7"""
assert old in t, "CONV downgrade anchor not found"
T.write_text(t.replace(old, new, 1), encoding="utf-8")
print(f"CONV_TODO.md -> {T.stat().st_size:,}")

# ── ③ 抽查 ──
nd = A.read_text(encoding="utf-8")
print("\n落盘抽查：")
for name, k in [("11具尸体", "11 具"), ("无韧带可被招募", "primary\n> stabilizing"),
                ("CFL 0.2~2.1%", "0.2 ~ 2.1"), ("不约束内外翻", "unconstrained"),
                ("极端位才受力", "extreme flexion"), ("toe区", "Toe")]:
    k2 = k.replace("\n> ", "\n")
    print(f"  [{'OK' if (k in nd or k2 in nd) else 'XX'}] {name}")
