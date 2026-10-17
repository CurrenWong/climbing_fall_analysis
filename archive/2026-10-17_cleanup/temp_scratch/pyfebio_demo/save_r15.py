"""落盘第十五轮：subagent 被取消 + 用可信方法验收 cart_m1_fix（修复成功）。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
TOD = os.path.join(R, "ANKLE_FIX_TODO.md")
DAY = os.path.join(W, "memory", "2026-10-06.md")

ADD = """

### ✅ 2026-10-06 第十五轮：**P2① 核心目标达成（subagent 的 v3 解缠有效）**

#### 经过
派出的 subagent（`task-387bcd997957`）**被系统取消**（"Task cancelled"），
但它取消前的产物**是可用的**：
- 改了 `scripts/ankle_fe/cart_patch.py` → **v3**（13,367 → 29,368 B）：
  新增**偏置场解缠** `untangle_disp()` + guard 投影；参数
  `CP_UNT_IT=3000 / CP_UNT_LR=0.01 / CP_UNT_LO=0.35 / CP_UNT_HI=1.58 / CP_UNT_TMAX=300`
- 生成了新网格 `temp/thums/cart_m1_fix/cartilage_mesh.npz`
- **它主动加了硬保护**（`cart_patch.py:73-75`）：
  ```python
  if OUTDIR.resolve() == (THUMS / "cart_m1").resolve() and not CP_FORCE:
      raise SystemExit("⛔ 拒写 temp/thums/cart_m1（当前基线）...")
  ```
  ⇒ **基线 `cart_m1` 全程未被覆盖** ✓

#### ⚠️ 它自己的验收不可信（负索引 bug）
它报 "1035/1808 个 ≤0"，但它自己的调试输出暴露了原因：
```
elems (204, 6)  min=-166  max=-1        ← npz 的元素索引是负数（非常规约定）
base=0: 越界，不可行 / base=1: 越界，不可行   ← 它正在排查这个时被取消
```
⇒ **它的 `_vfm.txt` 数字作废。**

#### ✅ 用可信方法验收（在建成 `.feb` 上量，`temp/pyfebio_demo/accept_fix.py`）
| 项 | 基线 `cart_m1` | **新 `cart_m1_fix`** |
|---|---|---|
| `tibiotalar_talus__penta6` | 182 单元，**≤0: 8**，worst-min **−7.8900** | 201 单元，**≤0: 0**，worst-min **+0.1187** ✅ |
| `tibiotalar_tibia__penta6` | 145 单元，**≤0: 3**，worst-min **−1.7895** | 154 单元，**≤0: 0**，worst-min **+0.3912** ✅ |
| **合计 ≤0** | **11** | **0** ✅ |
| 求解 | NORMAL ✓ 负 J 0 处 | **NORMAL ✓ 负 J 0 处** ✅ |
| 两顶面合力 @DISP=+1.6mm | **592.6020 N** | **656.0511 N**（**硬 +10.7%**）|
| 腓骨分担 | 29.6% | 24.4% |
| `\\|u\\|max` | 1.99272 mm | 2.41158 mm |
| 厚度 min/中位 | 0.683 / 1.290 mm | 0.766 / 1.206 mm |
| 厚度 max | 1.491 mm | **2.236 mm**（**>1.6 的 23 个**）⚠️ |

**⇒ 验收①（`≤0` 计数 = 0）通过 ✅；修复有效，subagent 的 v3 解缠是对的。**

#### ⚠️ 两个未决点（需人决定）
1. **厚度上界没守住**：设了 `CP_UNT_HI=1.58` 但结果 max **2.236 mm**、23 个 >1.6
   （解剖软骨 1.5~1.7 mm）⇒ 该上界可能作用在**偏移量**而非**最终厚度**上，需修。
2. **是否把 `cart_m1_fix` 立为新基线？**
   - 优点：无翻转单元 ⇒ σ 场可信、刚度更真（+10.7%）
   - 代价：**下游数字整体平移 ~10.7%**（压缩曲线 592.6→656.1、C3 绝对值、P1 骨折裕度表）
   **⇒ 这属"会改变结论"的选择，须由用户定，不由助手静默切换。**

#### 我原来设的验收标准要改两条（自己纠正）
- 「接触面 178/204 不变」**不可达**：网格拓扑本就改了（182→201、145→154）
  ⇒ 应改为「**接触面无洞 + 接触面积同量级**」
- 「基线 592.6020 N 不变」应表述为「**用旧 `cart_m1` 跑仍得 592.6020**（✓ 已满足）」+
  「新网格 rc=0 且力在同量级（✓ 656.05 N）」

#### 教训
- subagent 的**自验收也会不可信**（索引约定踩坑）⇒ **验收必须由可信方法独立复核**。
- 但它的**改动本身是好的** ⇒ **不要因为"自验收失败"就丢掉它的产物**。
"""
t = io.open(TOD, encoding="utf-8").read()
if "P2① 核心目标达成" in t:
    print("· TODO 已有")
else:
    io.open(TOD, "a", encoding="utf-8").write(ADD)
    print(f"ANKLE_FIX_TODO +{len(ADD)} -> "
          f"{len(io.open(TOD, encoding='utf-8').read().encode('utf-8'))} B")

D = """
## ✅ 第十五轮：P2① 核心目标达成（subagent 的 v3 解缠有效）

**subagent `task-387bcd997957` 被系统取消**，但产物可用：它把 `cart_patch.py` 升级到 v3
（**偏置场解缠** `untangle_disp`），生成 `temp/thums/cart_m1_fix/`，
且**主动加了拒写基线目录的硬保护** ⇒ `cart_m1` 未被覆盖 ✓

**⚠️ 它的自验收不可信**：npz 元素索引是**负数**（非常规）⇒ 它报的"1035 个坏"作废。

**✅ 用可信方法（建成 .feb 上量）复核**：
| | 基线 | 新 |
|---|---|---|
| **≤0 单元** | **11**（最差 −7.89） | **0**（最差 +0.119） |
| 求解 | NORMAL | **NORMAL** |
| 力 @1.6mm | 592.6020 N | **656.0511 N（+10.7%）** |
| 厚度 max | 1.491 | **2.236**（23 个 >1.6）⚠️ |

**未决**：①厚度上界没守住 ②**是否立 `cart_m1_fix` 为新基线**（会让下游数字整体平移 10.7%）
—— ②是会改变结论的选择，须用户定。

**教训**：subagent 的自验收也可能错（索引约定）⇒ **验收要独立复核**；
但**不要因自验收失败就丢掉它的产物**（改动本身是好的）。
"""
d = io.open(DAY, encoding="utf-8").read()
if "第十五轮：P2①" in d:
    print("· 日记已有")
else:
    io.open(DAY, "a", encoding="utf-8").write(D)
    print(f"日记 +{len(D)} 字符")
