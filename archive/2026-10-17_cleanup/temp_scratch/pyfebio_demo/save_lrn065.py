"""追加 LRN-065（带撞号保护）+ 追加 daily 小节。"""
import io
import os

P = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"

# ---- LRN-065 ----
LF = os.path.join(P, ".learnings", "LEARNINGS.md")
ID = "[LRN-20261005-065]"
src = os.path.join(P, "temp", "pyfebio_demo", "lrn065_snippet.md")
t = io.open(LF, encoding="utf-8").read()
if ID in t:
    print(f"LEARNINGS: {ID} 已存在，跳过")
else:
    b = os.path.getsize(LF)
    io.open(LF, "a", encoding="utf-8").write(io.open(src, encoding="utf-8").read())
    t2 = io.open(LF, encoding="utf-8").read()
    print(f"LEARNINGS.md {b} -> {os.path.getsize(LF)} B；{ID} 次数={t2.count(ID)}")

# ---- daily ----
DF = os.path.join(W, "memory", "2026-10-05.md")
BLOCK = """
## 🔧 辨析：位移控制 vs Anderson 弱弹簧（`[LRN-20261005-065]`）

**回答用户"弱弹簧和现在的机制区别是什么"时，先核对代码纠正了自己的用词错误。**

### ⚠️ 我此前的措辞是错的
反复说"`FIX_TALUS=all` 锁死旋转自由度" —— **实体单元节点只有 3 个平动 DOF，没有转动 DOF 可锁**。
真实机理是**运动学**上的：`talus_low_ns`(x,y 置零) + 胫腓骨顶面横向约束 ⇒
**两端横向硬约束 ⇒ 两轴强制平行 ⇒ 相对转动所需的横向位移被消 ⇒ 内外翻/内旋做不到**。

### 区别的本质
同一套方程/材料，**差别只在「哪些 DOF 硬/软约束」+「哪个量是输入」**。
- 位移控制：输入位移、输出支反力、姿势是**输入**（`POSE_TILT`）
- 弱弹簧：输入**力（=OpenSim 输出）**、输出**位移 + 关节姿态**、姿势是**输出**
- 弱弹簧**一次干两件事**：消刚体模态奇异（`[LRN-058]`）+ **保留姿势自由度**

### 定量点（容易说错）
"弱"是相对**骨**而言。**实测关节轴向刚度 ≈ 370~580 N/mm**，与 500 N/mm **同量级**。
串联下：**力不受影响**（关节内部状态仍由力定 ✅），但**总位移被稀释**
（1/450+1/500 ⇒ 表观 ~240 N/mm）⇒ **必须读关节相对变形，不是试验机总行程。**

### 🎯 修正我此前的路线说法
"找根法可代替弱弹簧"**只对一半**：找根只换"力的接口"（包装），
**横向硬约束原封不动 ⇒ 踝扭伤（26~28%，最大一块）仍算不出**。
**弱弹簧才是"打开姿势自由度"的能力**。目标需要后者。

### 动手前必须核实（不靠试错，`[LRN-20261003-015]`）
1. FEBio 弹簧边界条件**语法**（向 schema/手册要）
2. Anderson §2.2 的**弹簧布置**（挂点/总刚度/释放哪些方向）
"""
dt = io.open(DF, encoding="utf-8").read()
if "[LRN-20261005-065]" in dt:
    print("daily: 本节已存在，跳过")
else:
    b = os.path.getsize(DF)
    io.open(DF, "a", encoding="utf-8").write(BLOCK)
    print(f"daily {b} -> {os.path.getsize(DF)} B")
