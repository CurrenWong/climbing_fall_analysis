"""追加 LRN-067 + daily 小节（幂等）。"""
import io
import os

P = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"

LF = os.path.join(P, ".learnings", "LEARNINGS.md")
ID = "[LRN-20261005-067]"
src = os.path.join(P, "temp", "pyfebio_demo", "lrn067_snippet.md")
t = io.open(LF, encoding="utf-8").read()
if ID in t:
    print(f"LEARNINGS: {ID} 已存在，跳过")
else:
    b = os.path.getsize(LF)
    io.open(LF, "a", encoding="utf-8").write(io.open(src, encoding="utf-8").read())
    t2 = io.open(LF, encoding="utf-8").read()
    print(f"LEARNINGS.md {b} -> {os.path.getsize(LF)} B；{ID} = {t2.count(ID)} 次")

DF = os.path.join(W, "memory", "2026-10-05.md")
BLOCK = """
## 📖 核实 Anderson §2.2 弹簧口径（`[LRN-20261005-067]`）

**来源**：`paper/MinerU_markdown_Intra-articular_...md`（= **Anderson 2006 原文**）
+ `MinerU_markdown_anderson2007_...md`。
⚠️ **`paper/2006FAIankleligaments.pdf` 不是 Anderson 2006**（7 页、无 "spring"，
是 Tochigi FAI）—— **文件名不可信，引用前先验内容**。

### 核实后的口径
| 项 | 2006 | 2007 |
|---|---|---|
| 弹簧 | **linear springs 500 N/mm 接地**，约束距骨**平动** | **单根 50 N/mm**，抗**前后 A–P 平动** |
| 用途 | 替代缺失韧带/后足 ⇒ 距骨**自由落位** | 替代软组织（纯轴向场景） |
| **挂点** | 距骨**刚性**软骨下表面的 **reference node** | 同 |
| 转动 | 落位期固定**内外翻+内旋**两轴 ⇒ **duty cycle 前解除** | 屈伸固定 |
| 全程保留 | ① 距骨**垂直地面**位移 ② **屈伸 prescribed** | 屈伸固定 |

**★ 落位协议（逐字）**：step1 把胫骨带到与弱约束距骨的**最小接触** →
step2 **加 200 N 预载**（**解除内外翻/内旋转动约束之后**）⇒ 关节自行沉降到解剖对位 → duty cycle。

**★ 实验证据（2007）**：夹具**只固定屈伸**，其余平动+转动**全自由**
⇒ 内外翻/内旋自由是**物理事实**。

### ⚠️ 两个数字是两回事
**500 N/mm**（2006 接地平动弹簧）≠ **50 N/mm**（2007 单根 A–P 弹簧）。此前并列易误读，已澄清。

### ⭐⭐ 机制要点（决定实现）
**Anderson 能"固定/解除内外翻与内旋"，是因为他用了 rigid 刚体 + reference node**
（刚体 6 DOF，含 3 转动）。
**⇒ 纯可变形实体网格没有转动自由度（`[LRN-065]`），该约束在可变形模型里无法表达。**
**⇒ 忠实复现必须把距骨（至少软骨下表面）设为 FEBio `<rigid>` + reference node。**

### 对当前实现的判定
- ✅ `TALUS_SUPPORT=spring` 语法对、机制通
- ❌ **但不是 Anderson 口径**（他是刚体参考点，我是 323 节点分布式）
- ❌ `FIX_TALUS=transverse` **不对应 Anderson 任何一条**（全程保留的只有垂直地面位移 + 屈伸）
- ⇒ **正确下一步 = 引入 FEBio `<rigid>` 刚体距骨 + reference node**，
  再按 500 N/mm 平动弹簧 + 按需解除转动 + 200 N 预载落位

### 方法论教训
**摘录会丢"机制前提"**：我此前的 §2.2 抄了结论（弱弹簧/解除内旋），
却漏了"**因为骨是刚的、有 reference node**"这个**使结论可行的前提** ⇒ 照着结论实现就偏了。
**⇒ 文献的"前提"与"结论"必须一起抄。**
"""
dt = io.open(DF, encoding="utf-8").read()
if "[LRN-20261005-067]" in dt:
    print("daily: 已存在，跳过")
else:
    b = os.path.getsize(DF)
    io.open(DF, "a", encoding="utf-8").write(BLOCK)
    print(f"daily {b} -> {os.path.getsize(DF)} B")
