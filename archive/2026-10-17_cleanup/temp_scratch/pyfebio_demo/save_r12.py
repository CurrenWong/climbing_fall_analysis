"""落盘 P2② 完成：完整采样诊断上线（检测与处置分离）+ 零回归。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
TOD = os.path.join(R, "ANKLE_FIX_TODO.md")
LEARN = os.path.join(R, ".learnings", "LEARNINGS.md")
DAY = os.path.join(W, "memory", "2026-10-06.md")

ADD = """

### ✅ 2026-10-06 第十二轮：**P2② 完成 —— 完整采样诊断上线（检测与处置分离）**

#### 做了什么
1. **`ankle_lig_feb.py` 的 `detJ_penta6` 采样集漏顶点**：`RS` 只到 **2/3**
   ⇒ 三角的两个顶点 `(r,s)=(1,0)` / `(0,1)` **从未被采样**
   （只有 `(0,0)` 被采到）⇒ 漏掉「顶点附近翻转」的楔形单元。
2. **新增 `detJ_penta6_full`**（`RS` 含 `1.0`）+ 报警行：**只诊断、不驱动剔除**。
3. **`detJ_penta6` 恢复原 RS**（≈31 个剔除的基线行为不变）。
4. **`build_feb.py` 里的同名函数也补齐**（那份是另一条模型线用的）。

#### 实测结果（`p2a_verify.py`）
```
[detJ自检·完整采样] 🚨 有 43 个楔形在顶点附近翻转 ⇒ 局部应力不可信（σ_max 尖峰来源）；
                    **不可剔除**（在关节接触面上）⇒ 走 cart_patch.py 局部重铺
!!  tibiotalar_tibia  penta6 保留 145 剔除  9  det(J) min=1.9085e-01   ← 基线行为 ✓
!!  tibiotalar_talus  penta6 保留 182 剔除 22  det(J) min=2.3914e-03   ← 基线行为 ✓
翻转修正 1 个；仍无法修正而剔除 31 个                                    ← 原 31 ✓
FEBio rc=0  用时 8.6s                                                   ← ✅
```

#### ⭐ 实验证据：**这批单元不能剔除**（不是推论，是实测）
第一次改法（用完整采样**驱动剔除**）的结果：
```
tibiotalar_tibia  保留 142（原 145）剔除 12（原 9）    ← 多剔 3 个
tibiotalar_talus  保留 174（原 182）剔除 30（原 22）   ← 多剔 8 个
共剔除 42（原 31）
FEBio rc=1   ← 在关节接触面/tie 面上打洞 ⇒ 模型坏掉
```
**⇒ 多剔的 11 个全在接触面上 ⇒ 只能局部重铺，不能删除。**
**⇒ 由此确立设计原则：检测与处置分离。** 报警负责"让人知道"，
处置必须走 `cart_patch.py` 局部重铺。

#### 完整采样为什么报 43 个（而我先前只报 12 个）
密度不同：`p2_locate.py` 用 12 个采样点，`detJ_penta6_full` 用约 105 个
（`RS` 6 值满足 `r+s≤1` 的约 21 组合 × `TS` 5 值）⇒ **采样越密，抓到的越多**。
**⇒ 再次印证「采样点就是判据的一部分」**（`[LRN-20261006-085]`）。

#### ✅ 零回归
```
两顶面合力 = [0, 0, 592.602] N  |F| = 592.6020    ← 逐位一致
|u|max = 1.99272 mm                                ← 一致
FEBio rc=0
```
**⇒ `area_proj`（P2'a）与 `detJ` 诊断（P2②）两处改动都不动位移控制结果。**

#### 剩 P2①
`cart_patch.py` 局部重铺那一小片（修复扭角/顶点翻转），
验收：`worst-det(J) ≤ 0` 计数 **= 0** ｜ 接触面 **178/204 不变** ｜ 基线 **592.6020 N 不变**。
"""
t = io.open(TOD, encoding="utf-8").read()
if "P2② 完成" in t:
    print("· TODO 已有")
else:
    io.open(TOD, "a", encoding="utf-8").write(ADD)
    print(f"ANKLE_FIX_TODO +{len(ADD)} -> "
          f"{len(io.open(TOD, encoding='utf-8').read().encode('utf-8'))} B")

lt = io.open(LEARN, encoding="utf-8").read()
LADD = """

---

## [LRN-20261006-087] 「检测」与「处置」必须分离 —— 否则检查器自己会弄坏模型

**日期**: 2026-10-06　**类别**: best_practice　**项目**: climbing_fall_analysis / 踝 FE

### 现象
把 `detJ_penta6` 的采样补全（`RS` 增加 `1.0`，覆盖三角远端顶点）后，
被判定"坏"的楔形从 31 个增到 42 个 ⇒ 自动剔除逻辑多剔了 11 个 ⇒
**`FEBio rc=1`** —— 因为**多剔的 11 个全在关节接触面/tie 面上**，
删除等于**在接触面上打洞**，破坏 tie/contact 配对。

### 规则
> **检查器只负责"报告"，处置逻辑必须独立、且默认不改动基线。**
> 新增/加严一个判据时，**先只报警**，观察它会对既有流程做什么，
> 再决定处置策略；**绝不能**让"判据变严"自动变成"删除变多"。

### 本项目的正解
- `detJ_penta6_full`：完整采样，**只打印报警**（"不可剔除，走局部重铺"）
- `detJ_penta6`：保持原采样 ⇒ 剔除数 31 不变 ⇒ 基线 `592.6020 N` 逐位复现
- 真正的修复走 `cart_patch.py` 的**局部重铺**（不改几何拓扑、不在接触面打洞）

### 附带
采样密度直接决定抓到的数量：12 点 → 12 个坏单元；约 105 点 → 43 个。
**⇒ 采样点就是判据的一部分**（`[LRN-20261006-085]`）。
"""
if "LRN-20261006-087" in lt:
    print("· LRN 已有")
else:
    io.open(LEARN, "a", encoding="utf-8").write(LADD)
    print(f"LEARNINGS +{len(LADD)} -> "
          f"{len(io.open(LEARN, encoding='utf-8').read().encode('utf-8'))} B")

D = """
## ✅ P2② 完成（第十二轮）：完整采样诊断上线，检测与处置分离

**`detJ_penta6` 的 `RS` 只到 2/3 ⇒ 漏采三角远端顶点 (1,0)/(0,1)**
⇒ 漏掉顶点翻转的楔形。新增 `detJ_penta6_full`（`RS` 含 1.0）+ 报警行。

**⭐ 实测证据：这批单元不能剔除** —— 用完整采样**驱动剔除**时多剔 11 个，
`FEBio rc=1`（在关节接触面/tie 面打洞）。⇒ **检测与处置必须分离**（LRN-087）。

**现状**：报警 43 个（105 采样点）｜剔除 31（原基线）｜**rc=0**｜
**零回归 592.6020 N + |u|max 1.99272 mm 逐位一致** ✓

**剩 P2①**：`cart_patch.py` 局部重铺（验收：worst-det(J)≤0 计数 0、接触面 178/204、592.6020 N）
"""
d = io.open(DAY, encoding="utf-8").read()
if "P2② 完成（第十二轮）" in d:
    print("· 日记已有")
else:
    io.open(DAY, "a", encoding="utf-8").write(D)
    print(f"日记 +{len(D)} 字符")
