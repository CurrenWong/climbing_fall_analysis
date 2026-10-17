"""落盘 LRN-069 + daily + TODO 更新。"""
import io
import os

P = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"

LF = os.path.join(P, ".learnings", "LEARNINGS.md")
ID = "[LRN-20261005-069]"
src = os.path.join(P, "temp", "pyfebio_demo", "lrn069_snippet.md")
t = io.open(LF, encoding="utf-8").read()
if ID in t:
    print(f"LEARNINGS: {ID} 已存在，跳过")
else:
    b = os.path.getsize(LF)
    io.open(LF, "a", encoding="utf-8").write(io.open(src, encoding="utf-8").read())
    import re
    t2 = io.open(LF, encoding="utf-8").read()
    ids = re.findall(r"^##\s*\[LRN-(\d{8}-\d{3})\]", t2, re.M)
    print(f"LEARNINGS.md {b} -> {os.path.getsize(LF)} B；重复={sorted({i for i in ids if ids.count(i)>1})}")

DF = os.path.join(W, "memory", "2026-10-05.md")
BLOCK = """
## 🦴 路线 A 落地：距骨刚体 + 接地弹簧（`[LRN-20261005-069]`）

**用户选定路线 A（忠实复现 Anderson）并要求专注解决。**

### 已实现（env 开关，默认关）
`TALUS_RIGID=1` + `TALUS_SPRING_K=500` + `TAL_ROTLOCK` ⇒
距骨 916 节点 → `rigid-solid` 刚体；地面 2 哑节点；`<rigid_connector type="rigid spring">`。

### 位移控制跑通 ✅
| 配置 | 状态 | Fz |
|---|---|---|
| `TAL_ROTLOCK=1` | NORMAL | **202.29 N** |
| `TAL_ROTLOCK=0`（解除内外翻/内旋） | NORMAL | **187.90 N** |

⇒ **转轴解除使结果变 −7.1% ⇒ 机制真的在起作用。**

### 两个真 bug
1. **`body_a`/`body_b` 是材料 ID 不是名字** ⇒ 写死 1/2 会引用到 `mat_tibia_cort`
   ⇒ `Nonlinear constraint 1 failed to initialize`。正解 `len(MATERIALS)+1`。
   **最小算例"恰好成立"的常量是陷阱。**
2. **`<Rigid>` 必须在 `</Boundary>` 之后**（手册推荐序反而失败）⇒ 用 `_RIGID_LINES` 缓冲。

### ⛔ 力控制仍失败（本轮到此）
负 Jacobian 已被刚体化消除，但变成**不收敛**（t=0.0~0.235）；5 N 载荷却读出
胫骨顶 177.5 N、距骨反力 5004 N ⇒ 远未平衡。
**根因 = 关节初始张开（`[LRN-058]`）⇒ 加载时关节未闭合，近似机构。**
**⇒ 正解 = Anderson 落位协议：step1 最小接触 → step2 200 N 预载 → 再 duty cycle。**

### 下一步（明确）
**多步分析**：`<Step>` 里 step1 位移落位 / step2 力控制到 600 N（手册 §6.1.2 支持跨步切换）。

### 教训
- **"最小算例通过"只证明语法，不证明搬到大模型也成立** —— 本轮两个 bug 都属此类。
- **失败态会给出看似合理的数值** ⇒ 读数前必须确认 `t=1.0` + `NORMAL`。
"""
dt = io.open(DF, encoding="utf-8").read()
if "[LRN-20261005-069]" in dt:
    print("daily: 已存在，跳过")
else:
    b = os.path.getsize(DF)
    io.open(DF, "a", encoding="utf-8").write(BLOCK)
    print(f"daily {b} -> {os.path.getsize(DF)} B")
