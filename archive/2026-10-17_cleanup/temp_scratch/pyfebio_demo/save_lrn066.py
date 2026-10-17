"""追加 LRN-066 + daily 小节（均带幂等保护）。"""
import io
import os

P = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"

LF = os.path.join(P, ".learnings", "LEARNINGS.md")
ID = "[LRN-20261005-066]"
src = os.path.join(P, "temp", "pyfebio_demo", "lrn066_snippet.md")
t = io.open(LF, encoding="utf-8").read()
if ID in t:
    print(f"LEARNINGS: {ID} 已存在，跳过")
else:
    b = os.path.getsize(LF)
    io.open(LF, "a", encoding="utf-8").write(io.open(src, encoding="utf-8").read())
    t2 = io.open(LF, encoding="utf-8").read()
    print(f"LEARNINGS.md {b} -> {os.path.getsize(LF)} B；{ID} 次数={t2.count(ID)}")

DF = os.path.join(W, "memory", "2026-10-05.md")
BLOCK = """
## 🔩 弱弹簧实现 + FEBio `<Discrete>` 语法（`[LRN-20261005-066]`）

**用户要求：先把经验保存，然后尝试弱弹簧实现。两件都做了。**

### 语法（从 `FEBio_User_Manual.pdf` §3.6.7/§3.16 取，且最小算例实测）
```xml
<DiscreteSet name="s"><delem>3186,4030</delem></DiscreteSet>     <!-- Mesh 内 -->
<Discrete>                                                       <!-- MeshDomains 之后 -->
  <discrete_material id="1" type="linear spring"><E>1.548</E></discrete_material>
  <discrete dmat="1" discrete_set="s"/>
</Discrete>
```
**校验**：Δk=99、δ=0.1 ⇒ ΔF 实测 **104.70 N** = 理论 105 N ✓

### 3 个实测语法坑
1. `<NodeSet>` 值是**逗号列表**（不是 `<node>` 子元素）
2. `zero displacement` 用 `<x_dof>1</x_dof>`；`prescribed displacement` 用
   `<dof>x</dof><value lc="1">…</value>`（**需 LoadData 段**）
3. **`<Discrete>` 必须在 `<MeshDomains>` 之后**（否则 `tag "discrete" : invalid value`）

**⇒ 教训：小算例只证明"该顺序可行"，不证明"顺序无关"。**

### 实现（`ankle_lig_feb.py`，默认关 ⇒ 零回归）
`TALUS_SUPPORT=spring` + `TALUS_K`(500) + `TALUS_SPRING_L`(5mm)
⇒ 323 距骨节点 × 2 横向方向 = **646 根接地弹簧**，单根 k=1.548 N/mm，
**`converged at time : 1`** ✓；零回归（209.5589 / 592.6068 逐字一致）。

### 硬结果：弹簧 vs 硬约束（轴向）
| DISP | 硬约束 | 弱弹簧 | 差 |
|---|---|---|---|
| 0.25 | 124.04 | 111.26 | −10.3% |
| 1.00 | 418.10 | 369.26 | −11.7% |

**⇒ 差 10~12% 是「边界条件真不同」，不是误差**（横向释放后关节能横向让开）。

### ⛔ 关键负面结果
`DRIVER=force` 扫 `TALUS_K=50/500/5000` ⇒ **全部 11 负 Jacobian**。
**弹簧解掉了刚体模态奇异，但没解决力控制** —— 真因是 `[LRN-058]` 那条：
**关节初始张开（~1 mm 间隙）** ⇒ 需要 **Anderson 的「先落位、后加载」**，
而手册 §6.1.2 明确支持**跨步切换 位移控制 ↔ 载荷控制**。

### 下一步
**多步：step1 位移落位 → step2 力控制到 600 N**；之后才是 OpenSim 力驱动 + 姿态输出。
仍待核实：Anderson §2.2 的弹簧布置口径（本实现是合理猜测，未对标原文）。
"""
dt = io.open(DF, encoding="utf-8").read()
if "[LRN-20261005-066]" in dt:
    print("daily: 本节已存在，跳过")
else:
    b = os.path.getsize(DF)
    io.open(DF, "a", encoding="utf-8").write(BLOCK)
    print(f"daily {b} -> {os.path.getsize(DF)} B")
