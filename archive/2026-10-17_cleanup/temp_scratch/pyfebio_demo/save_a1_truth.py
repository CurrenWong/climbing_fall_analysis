"""落盘 A1 真实结果：H1 成立 + 「缺距下关节 44.8%」被否证 + A1/A2 耦合。

涉改：MEMORY.md（两处 44.8%）、ANKLE_ROADMAP_TODO.md、日记。
原则：更正要替换不要叠加。
"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
n = 0

# ---------- 1) MEMORY.md ----------
MP = os.path.join(W, "MEMORY.md")
m = io.open(MP, encoding="utf-8").read()
m0 = len(m.encode("utf-8"))
MREPL = [
    ("- 🔍 **2~4× 的来源 = 三条建模伪影（取代旧「刚度 2× 未归因」）**：① **缺距下关节接触 44.8%**\n"
     "  （`subtalar_*` 软骨**已生成、未接进 `.feb`**）② **无摩擦切向滑移 39.5%**"
     "（关节间只 0.029mm 真压扁）③ **tie 罚滑移 13.9%**；骨 <6%",
     "- 🔍 **2~4× 的来源（取代旧「刚度 2× 未归因」）**：① ~~缺距下关节接触 44.8%~~ "
     "**已否证**（2026-10-06 直接测量：接入距下关节 ⇒ F 变化 **+0.0000%**，因 "
     "`FIX_TALUS=transverse` 冻结距骨 ⇒ 接触不承力）② 无摩擦切向滑移 **39.5%**"
     "（关节间只 0.029mm 真压扁）③ tie 罚滑移 **13.9%**；骨 <6%"),
    ("- **⏭️ 下一步**：① **接入距下关节**（最大伪影 44.8%，软骨已生成⇒纯接线）② 力控找根接 OpenSim",
     "- **⏭️ 下一步**：① **A1+A2 一起做**（接入距下关节 **且** 解除 `FIX_TALUS=transverse`）"
     "—— 单接关节 **零效应**（实测 +0.0000%），因横向硬约束冻结距骨、接触不承力；"
     "而**接入后收敛需 `MAX_UPS=0`**（根因=拟牛顿更新病态，非几何）② 力控找根接 OpenSim"),
]
for a, b in MREPL:
    if a in m:
        m = m.replace(a, b); n += 1
    else:
        print("  ⚠ MEMORY 未命中：", a[:40].replace("\n", "\\n"))
io.open(MP, "w", encoding="utf-8").write(m)
print(f"MEMORY {m0} -> {len(m.encode('utf-8'))} B "
      f"({'✅' if len(m.encode('utf-8')) <= 10240 else '❌超'})")

# ---------- 2) ANKLE_ROADMAP_TODO.md ----------
RP = os.path.join(R, "ANKLE_ROADMAP_TODO.md")
r = io.open(RP, encoding="utf-8").read()
ADD = """

---

## 🎯 A1 真实结果（2026-10-06 晚，主 agent 亲测）

### ① H1 成立：收敛根因 = **拟牛顿更新病态**（不是几何/网格/tie/接触）
`temp/pyfebio_demo/a1_matrix2.py`（执行者写好、被中断未跑）已由主 agent 跑完：
```
[cal_noupd]    rc=0 NORMAL=True F=592.5697   只接跟骨侧 + MAX_UPS=0
[both_noupd]   rc=0 NORMAL=True F=592.5697   ★ 双侧全接 + MAX_UPS=0（= A1 目标配置）
[cal_noupd_r2] rc=0 NORMAL=True F=592.5697   重跑一致
```
⇒ **`MAX_UPS=0`（关 BFGS 更新）后，距下关节双侧接入完全收敛。**
⇒ 原 `rc=1`（11 负 Jacobian、元件级 `tibia_cort__tet4` jacobian −1e12）是
**BFGS 首次更新把 K 污染成近奇异 ⇒ 解沿近零刚度方向跑飞**，与软骨几何无关。

### ② ⭐「缺距下关节接触 44.8%」**被直接测量否证**
`temp/pyfebio_demo/a1_isolate.py`（正对照式隔离）：
| 配置 | 两顶面合力 |
|---|---|
| 不接，BFGS 开 | 592.6020 N（基线） |
| **不接，BFGS 关** | **592.5697 N**（`MAX_UPS=0` 自身效应 −0.0054%） |
| **接入，BFGS 关** | **592.5697 N** |
| **⭐ 接入的真实效应** | **+0.0000 N = +0.0000%** |

⇒ **接入距下关节对总刚度影响 = 0.0000%**（记录里写的 44.8% 差 4 个量级）。

### ③ ⭐ 机制与**架构级结论**：A1 与 A2 是同一件事的两半
`FIX_TALUS=transverse` 把距骨底面横向分量**硬置零** ⇒ 距骨不能平动/转动
⇒ **距下接触面根本不承力** ⇒ 载荷绕开它走（所以接入零效应）。
⇒ **不解除横向硬约束，接上去的关节等于没接。**
⇒ 因此 A1 与 A2 **必须一起做**；单独做 A1 在物理上是空操作。

### ④ 对路线的影响
- **A2 的定位升级**：从「唯一不确定的可选步骤」→ **A1 生效的必要条件**
- **A3 后「绝对刚度 ±20%」的提法暂缓**：三条伪影里最大那条已倒，
  剩下 39.5%（无摩擦切向滑移）+ 13.9%（tie 罚滑移）才是真目标
- **下一步 = A1+A2 合并**：接关节 + 解除横向硬约束 + `MAX_UPS=0`（或找更干净的求解器设置），
  再看刚度到底变多少、距下接触是否真的承力（`LOG_CONTACT` 读数通道已由执行者建好）

### ⑤ 尚未回答
- `MAX_UPS=0` 是「诊断开关」还是「可接受的求解设置」？
  ⚠️ 它改变了基线 **−0.0054%**（592.6020 → 592.5697）
  ⇒ 若要用它，**必须显式声明并更新零回归基准**（否则零回归测试会失败）
- 距下接触在解除约束后承载多少？（需 `LOG_CONTACT=1` 跑，通道已建）
"""
if "A1 真实结果（2026-10-06 晚，主 agent 亲测）" not in r:
    io.open(RP, "a", encoding="utf-8").write(ADD)
    print(f"ROADMAP +{len(ADD)} 字符")
else:
    print("· ROADMAP 已有")

# ---------- 3) 日记 ----------
DP = os.path.join(W, "memory", "2026-10-06.md")
d = io.open(DP, encoding="utf-8").read()
DADD = """
## 🎯 A1 真相（2026-10-06 晚）+ 我的一次误判更正

**我误判了 subagent**：我看到「工作区文件没提交、无报告」就写「它跳过了几何预检、
未交付预检报告」。**实际它第一步就做了预检（含正对照），且落盘在
`docs/A1_距下关节_交接报告.md`（13357B）+ `ANKLE_ROADMAP_TODO.md`「执行者补充」+
`SUBTALAR_A1_TODO.md`**。用户提醒后我才读到。
**教训：判 subagent「没做」之前，必须先找它的落盘产物**（对话里没回报 ≠ 没交付）。
已把错误陈述在 ROADMAP 里**就地替换**（不叠加），并加顶部更正横幅。

**我跑出了 A1 缺的那一步（H1 判据）**：
- `a1_matrix2.py`（执行者写好被中断）⇒ **`MAX_UPS=0` 后双侧接入完全收敛**
  ⇒ **根因 = 拟牛顿更新病态**（BFGS 首更新把 K 污染成近奇异 ⇒ 解跑飞）
- `a1_isolate.py`（我写的正对照隔离）⇒ **接入距下关节的总刚度效应 = +0.0000%**
  ⇒ **「缺距下关节 44.8%」被直接测量否证**（差 4 个量级）

**⭐ 架构级结论**：`FIX_TALUS=transverse` 冻结距骨 ⇒ 距下接触不承力
⇒ **A1 单独做是空操作；A1 与 A2 必须一起做。**
"""
if "A1 真相（2026-10-06 晚）" not in d:
    io.open(DP, "a", encoding="utf-8").write(DADD)
    print(f"日记 +{len(DADD)} 字符")
else:
    print("· 日记已有")
