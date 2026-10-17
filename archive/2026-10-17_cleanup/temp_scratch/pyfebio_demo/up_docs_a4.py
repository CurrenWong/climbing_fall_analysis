"""更新文档：把「突破 +1.2mm 上界」立为进行中的 A4，并写清做法与判据。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"

# ---------- ROADMAP ----------
RP = os.path.join(R, "ANKLE_ROADMAP_TODO.md")
r = io.open(RP, encoding="utf-8").read()

# 状态表：加 A4 行
OLD_ROW = "| A2 解除横平面硬约束 | ⏳ **阻塞**（依赖 A1） | | |"
NEW_ROW = ("| **A1+A2** 合并 | ✅ **已打通**（见文末大节） | 可用配置 + 分解 + 边界 | "
           "**DISP ≤ +1.2mm（519 N）**；+1.3 起负 J |\n"
           "| **A4** 突破 +1.2mm 上界 | 🚧 **进行中**（用户 2026-10-06 选定） | | "
           "目标：让 600 N 可达 |\n"
           "| A2 解除横平面硬约束 | ✅ **已做**（`TALUS_SUPPORT=spring`，−11.74%） | | |")
if OLD_ROW in r and "**A4** 突破 +1.2mm 上界" not in r:
    r = r.replace(OLD_ROW, NEW_ROW)
    print("✓ 状态表已加 A1+A2 / A4 行")
else:
    print("· 状态表已更新过或未命中")

A4 = """

---

# 🚧 A4（进行中）：突破 +1.2 mm 量程上界

**用户 2026-10-06 选定。目标：让 600 N 可达**（当前 +1.2mm ⇒ 519 N；基线可到 +1.6mm/592.6N）。

## 已知事实（A1+A2 节）
- 上界处（+1.6mm）的卡点 = **单个单元 8241**，域 `tibiotalar_talus__penta6`，gauss point 5
- **不是距下关节的问题** —— 是**原有的胫距软骨**被新增载荷路径压过头
- ⛔ `cart_m1_fix`（修掉 11 个翻转单元的那版）**反而更差** ⇒ 单纯「消除翻转单元」不是解

## 诊断顺序（先测不改）
1. **表征单元 8241**：厚度、静息 det(J)、长细比、所在位置（是否在接触区边缘/高应力区）
   → 判断它是**离群单元**还是**代表性单元**
2. **域内分布**：`tibiotalar_talus` 全体单元的厚度/质量分布，看 8241 排在什么分位
3. **对照**：能收敛的基线（E0）在同一处的变形量 vs A1+A2

## 候选修法（按「不破坏贴合度」优先级）
| # | 方法 | 理由 | 风险 |
|---|---|---|---|
| 1 | **局部加密**（只细化 `tibiotalar_talus` 接触区） | 保持关节面贴合，只提高局部分辨率 | 单元数上升；需查 `LAUGON` |
| 2 | **厚度下限提升**（把过薄单元抬到 ~0.5mm 以上） | 薄单元在大压应变下先翻负 | 改变力学（需声明） |
| 3 | **`cart_patch.py` 提高 m**（m=2/3） | 已有工具、已验证 rc=0 | 单元数 × 5~10 |
| 4 | 接触参数针对该对调优（penalty/search_tol） | 便宜 | 治标；且已证 AUGLAG 更差 |

## 硬判据
- **不得回退**：+1.0mm 处 F 必须仍在 **449.63 N ± 小量**（A1+A2 基准）
- **不得破坏贴合度**：+1.0mm 与 +1.2mm 处应均收敛
- **目标**：`DISP ≥ +1.4mm` 可用（对应 >600 N 余量）
- 改完必须与基线（418.10 N @+1.0mm）和 A1+A2（449.63 N）**双向对照**
"""
if "A4（进行中）：突破 +1.2 mm 量程上界" not in io.open(RP, encoding="utf-8").read():
    io.open(RP, "a", encoding="utf-8").write(A4)
    print(f"✓ ROADMAP 追加 A4 节 -> "
          f"{len(io.open(RP, encoding='utf-8').read().encode('utf-8')):,} B")

# ---------- 日记 ----------
DP = os.path.join(W, "memory", "2026-10-06.md")
DADD = """
## 🚧 用户选定下一步：突破 +1.2mm 量程上界（A4）

**目标**：让 600 N 可达（现 +1.2mm ⇒ 519 N）。
**已知卡点**：单元 8241（`tibiotalar_talus__penta6`, gauss pt 5）翻负 —— **不是距下关节的问题**。
**做法**：先表征 8241（厚度/detJ/位置/分位）→ 再选修法（首选局部加密，保贴合度）。
**硬判据**：+1.0mm 处 F 须仍 ≈449.63 N（A1+A2 基准）；目标 DISP ≥+1.4mm 可用。
"""
if "用户选定下一步：突破 +1.2mm 量程上界（A4）" not in io.open(DP, encoding="utf-8").read():
    io.open(DP, "a", encoding="utf-8").write(DADD)
    print("✓ 日记已追加")
