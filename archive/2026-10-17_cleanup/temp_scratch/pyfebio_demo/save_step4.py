"""落盘 [LRN-061]：600 N 不可达 + 完整诊断链（含负面结果）。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
LT = ROOT / ".learnings/LEARNINGS.md"

t = LT.read_text(encoding="utf-8")
have = set(re.findall(r"\[LRN-(\d{8}-\d{3})\]", t))
NEW = next(f"20261005-{k:03d}" for k in range(61, 99) if f"20261005-{k:03d}" not in have)
print(f"⇒ LRN-{NEW}")

ENTRY = r"""

---

## [LRN-NEW] 踝子线：**600 N 到不了** —— 稳定窗只有轴向 0.1~0.3 mm，硬天花板 ~340 N

- **日期**：2026-10-05
- **类别**：`knowledge_gap`（**结构性天花板**，非配置问题）
- **承接**：`[LRN-20261005-060]`（`ALIGN_AXIS` 把加载方向做成无歧义的真轴向）

### 🔴 主结论：目标载荷 600 N **不可达**

固定 `ALIGN_AXIS=1 + AXIS_SRC=full + dispx + FIX_TALUS=transverse + TRUNC_MM=60 + MATSET=WK`，
只扫轴向位移：

| 轴向 `DISP_MM` | 力（顶面合力，已核对记账） | 结果 |
|---|---|---|
| 0.02 | — | **ERROR**（t=0.42，18 负 J） |
| 0.20 | **183.2 N** | NORMAL ✓ |
| 0.30 | **235.5 N** | NORMAL ✓ |
| 0.50 | — | **ERROR**（t=0.90） |
| 1.00 | — | **ERROR**（t=0.50） |

**⇒ 稳定窗 ≈ 轴向 `[0.1, 0.3] mm`，对应 **≤ 235 N**。**600 N 需要 ~0.7~0.9 mm ⇒ 远超窗口。**

### 📉 天花板是真的（已排除"力没记账全"）

一度怀疑「`fix_talus` 约束 y,z ⇒ 力从距骨旁路 ⇒ 顶面反力低估」。**实测否证**：

```
胫骨顶面反力合力 = [0, 0, 305.1271] N
距骨约束反力合力 = [0, 0, -0.0]   N   ← tal_dofs 实际是 "xy"，z 本来就自由 ⇒ 无旁路
腓骨顶面反力合力 = [0, 0,  41.3133] N
⇒ 合计 346.4404 N = 施加的关节载荷（记账完整）
```

### 📊 完整"杠杆 × 效果"矩阵（全部实测）

| 杠杆 | 在 DISP=0.50（原失败点）的结果 | 结论 |
|---|---|---|
| （基线，`DET_MIN=0.5`） | **ERROR** | 基线 |
| **`DET_MIN=0.95`**（剔近退化单元） | **NORMAL，334.4 N** | ✅ **最强杠杆** |
| `CART_E=1000`（软骨硬 100×） | **NORMAL，561.1 N** | ✅ 支持"软骨单元是限流" |
| `CART_E=100` | 494.6 N | ✅ 同上（单调） |
| `DET_MIN=0.95` + 厚软骨(`T_TARGET=3.0`) | **NORMAL，346.4 N** | 🟡 略抬 |
| `T_TARGET=3.0`（厚软骨，**穿透过 0**） | **ERROR**（0.25 mm / 94.8 N） | ❌ 单独更差 |
| `NODE_RELOC=1`（接触落位） | **ERROR**（t=0.0，11 负 J） | ❌ 更差（瞬时投影压塌薄层） |
| `FIX_TALUS=all` | ERROR（0.50 / 1.00 都失败） | ❌ 无效 |
| `LIG_MODE=none` | 发散（`|u|max=3.4e9`） | ❌ 无载荷通路 |
| `DET_MIN=0.95` @ DISP=0.80/1.00/1.20 | 全部 ERROR，**F 恒 ≈ 340 N** | ⛔ **天花板 ~340 N** |
| `DET_MIN=0.95` + 厚软骨 @ DISP=1.00 | ERROR，378.9 N | ⛔ 天花板抬到 ~380 N |

**⇒ 两条要点：**
1. **`DET_MIN` 是最强杠杆**（0.3 mm → 0.6 mm）⇒ **限流的确是**近退化单元**（软骨 `det(J) min=5.08e-01`，
   `DET_MIN=0.5` 刚好放过它们）。
2. **但所有杠杆叠加后天花板只到 ~380 N**，且 `DISP=0.8/1.0/1.2` 给出**同一个 F** ⇒
   **力不再随位移增长 ⇒ 有东西"触底"了。⇒ 单靠调参/调材料不可达 600 N。**

### ❌ 两个被否证的假设（负面结果，同等重要）

- **"两软骨面不共形（间隙 0.97 mm）是主因"** —— **不被实验支持**：
  把 `T_TARGET=1.5` 的封顶放开、重生成软骨（`T_TARGET=3.0` ⇒ 两层恰好填满，
  **四个关节`穿透过`全 0**、中位厚度 1.879/1.909 mm），单独使用**反而更早失败**（0.25 mm）。
  ⇒ 间隙不是主因。（该新网格作为可复用产物保留在 `temp/thums/cart_t3/`。）
- **"距骨轴向自由（只靠软韧带）是主因"** —— 否证：`FIX_TALUS=all` 在 0.5/1.0 仍失败。

### 🧾 本轮新增/改动的资产

| 项 | 说明 |
|---|---|
| `AXIS_SRC` / `ALIGN_AXIS` | 见 `[LRN-060]`（真轴向加载，默认关）|
| `NODE_RELOC` / `LAUGON` / `PENALTY` | 接触参数改为可测（**默认值 = 既有行为**）|
| `T_TARGET`（`cartilage.py`） | 厚度上限改为 env 门控（默认 1.5 = 既有行为）；`CART_OUT` 可写到新目录 |
| `temp/thums/cart_t3/` | `T_TARGET=3.0` 的软骨网格（**0 穿透**，可复用）|
| `acceptance.py` | 增加**距骨反力记账** + 全局平衡核对 |
| `temp/pyfebio_demo/baseline_align.py` | 位移扫描矩阵（**自动删读数防"残留冒充"**）|

### 🧠 方法论教训

- 🔴 **失败运行会留下"上一次的读数文件" ⇒ 会冒充本次结果。**
  本轮第一次跑矩阵就踩了：第一行三个数字全是上一次 `disp FIX_TALUS=all` 的残留（`|F|=326.728`
  与之前逐字相同才被识破）。**⇒ 脚本必须在每次运行前 `unlink` 读数文件。**
- **"账要平"是硬检查**：`顶面 + 距骨约束` 是否 ≈ 0（或等于施加载荷），能一次性排除
  "力从别的约束旁路"这类隐蔽错误。**建议任何新模型先做一次全约束反力求和核对。**
- **天花板要测"是否随驱动量增长"**：`DISP=0.8/1.0/1.2` 给出同一个 F，才敢说"封顶"；
  否则只是"没跑到"。
- **负面结果与正面结果同价**：本轮"间隙说"和"距骨自由说"都是靠反例否掉的，
  省下了后续几轮弯路。

### 📌 下一步（唯一剩下的强候选）

`DET_MIN=0.95` 有效 ⇒ **要一个"质量好的软骨网格"，而不是继续调参**。
现成的 `cart_patch.py`（2D 参数平面重心网格，底面翘曲恒 0，m=1..4 全 rc=0，见 `[LRN-051]`）
正是为此，**但它输出 prism6 + `tie_faces` 的"新格式"，需要 TIE 约束支持，
而 `ankle_lig_feb.py` 没有 TIE**（它靠"软骨与骨共享节点"）。
⇒ 两条路：① 给 `ankle_lig_feb.py` 移植 TIE（`build_feb.py` 已有，可参照）；
② 在共享节点格式下重做软骨的 2D 参数化网格。
**在能跑到 600 N 之前，C2~C5（姿势判据）都无法验收。**
""" .replace("LRN-NEW", f"LRN-{NEW}")

LT.write_text(t.rstrip() + ENTRY, encoding="utf-8")
print(f"LEARNINGS.md -> {LT.stat().st_size:,} B")

# ── MEMORY ──
M = WS / "MEMORY.md"
s = M.read_text(encoding="utf-8")
OLD = """- ⚠️ **加载方向不可靠**（`[LRN-060]`）：`dispx` 只 prescribe `argmax(|轴|)` 的全局单轴，且 `TRUNC_MM=60`
  把轴转 **12.18°** ⇒ 同 0.2mm 名义位移的 k 在 **290/416/916 N/mm** 跳（**3.2×**）。**正解=`ALIGN_AXIS=1`**"""
NEW2 = """- ⚠️ **加载方向不可靠**（`[LRN-060]`）：`dispx` 只按 `argmax(|轴|)` 全局单轴 + `TRUNC_MM=60` 转轴 **12.18°**
  ⇒ k 跳 **3.2×**（290/416/916）。**正解=`ALIGN_AXIS=1`**
- ⛔ **600 N 到不了**（`[LRN-061]`）：稳定窗仅轴向 **0.1~0.3mm**（≤235 N），硬天花板 **~340 N**
  （`DET_MIN=0.95` 抬到 0.6mm；叠加厚软骨 ~380 N）；`CART_E=1000` 可达 561 N ⇒ **限流=软骨单元质量**"""
TRIMS = [
    ("① 只靠接触约束的部件必有**刚体模态**⇒奇异⇒单元成批翻负 ⇒ 把加载端**整个截面全约束**（夹持端）",
     "① 只靠接触约束的部件必有**刚体模态**⇒奇异⇒单元成批翻负 ⇒ 加载端**整个截面全约束**"),
    ("""### 其他
- MCP 死代码 → 自接管：模板 `references/mcp-bypass.md`（Tavily / Flint 同款）
- 定时任务 → `copaw cron`（skill）；Server 酱推送 → `references/operational-gotchas.md`""",
     """### 其他
- MCP 死代码 → 自接管：`references/mcp-bypass.md`（Tavily / Flint 同款）；定时任务 → `copaw cron`；
  Server 酱 → `references/operational-gotchas.md`"""),
    ("""1. **MEMORY.md ≤ 10 KB 硬上限**
2. 新内容默认进 `references/`，主 MEMORY 只留指针；项目细节归对应 references
3. **搬，不是删**（迁走前确认新文件有全部信息）；整理完必验：**尺寸 + 无死链 + 信息未丢**
4. 更新同步改 `last_updated`；重大重构前备份到 `memory_archive/`""",
     """1. **MEMORY.md ≤ 10 KB 硬上限**；新内容默认进 `references/`，主 MEMORY 只留指针
2. **搬，不是删**（迁走前确认新文件有全部信息）；整理完必验：**尺寸 + 无死链 + 信息未丢**
3. 更新同步改 `last_updated`；重大重构前备份到 `memory_archive/`"""),
]
if OLD in s:
    s = s.replace(OLD, NEW2, 1)
    print("✓ MEMORY 已更新")
for a, b in TRIMS:
    if a in s:
        s = s.replace(a, b, 1)
    else:
        print(f"⚠ 未匹配: {a[:30]}...")
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'OK' if sz <= 10240 else 'XX 超'}")

# ── TODO ──
TD = ROOT / "ANKLE_GOAL_TODO.md"
s = TD.read_text(encoding="utf-8")
s = s.replace("- [ ] **Step 3** `DRIVER=force` 跑通 → 验判据 1~2",
              f"- [x] **Step 2-1** 扫出稳定窗：**轴向 0.1~0.3mm ⇒ ≤235 N**；硬天花板 **~340 N**"
              f"（`DET_MIN=0.95`），叠加厚软骨 ~380 N ⇒ **600 N 到不了**。`[LRN-{NEW}]`\n"
              "- [x] **Step 2-2** 否证两条假设：**间隙说**（`T_TARGET=3.0` 单独更差）与"
              "**距骨自由说**（`FIX_TALUS=all` 无效）；**最强杠杆=`DET_MIN`**（⇒ 限流=软骨单元质量）\n"
              "- [ ] **Step 3a** ⭐ **要一个质量好的软骨网格**（`DET_MIN=0.95` 已证明有效）：\n"
              "  `cart_patch.py` 可用但输出 prism6+`tie_faces` 新格式，**需给 `ankle_lig_feb.py` 移植 TIE**\n"
              "  （`build_feb.py` 已有可参照）；或在共享节点格式下重做 2D 参数化网格\n"
              "- [ ] **Step 3** `DRIVER=force` 跑通 → 验判据 1~2")
TD.write_text(s, encoding="utf-8")
print(f"ANKLE_GOAL_TODO.md -> {TD.stat().st_size:,} B")

# ── session + daily ──
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
SEC = f"""

---

## 十九(10-05)、**600 N 到不了** —— 稳定窗与硬天花板（2026-10-05）→ `[LRN-{NEW}]`

### 19.1 稳定窗（`ALIGN_AXIS` 下，只扫轴向位移）

| `DISP_MM` | 力 | 结果 |
|---|---|---|
| 0.02 | — | **ERROR**（t=0.42，18 负J） |
| 0.20 | **183.2 N** | NORMAL ✓ |
| 0.30 | **235.5 N** | NORMAL ✓ |
| 0.50 | — | **ERROR** |
| 1.00 | — | **ERROR** |

**⇒ 稳定窗 ≈ 轴向 `[0.1, 0.3] mm` ⇒ ≤ 235 N。** 600 N 需 ~0.7~0.9 mm ⇒ 远超窗口。

### 19.2 天花板已核实（排除"记账不全"）

`胫骨顶 305.13 + 腓骨顶 41.31 = 346.44 N`；**距骨约束反力 = 0**（`tal_dofs` 实际为 **`xy`**，z 本就自由）
⇒ **346.44 N 就是关节载荷，无旁路。**

### 19.3 杠杆矩阵（全部实测）

| 杠杆 | DISP=0.50 结果 |
|---|---|
| 基线 `DET_MIN=0.5` | **ERROR** |
| **`DET_MIN=0.95`** | **NORMAL / 334.4 N** ← 最强 |
| `CART_E=1000` | **NORMAL / 561.1 N** |
| `CART_E=100` | 494.6 N |
| `DET_MIN=0.95` + 厚软骨 | NORMAL / 346.4 N |
| `T_TARGET=3.0`（厚软骨，**0 穿透**） | **ERROR**（0.25mm / 94.8 N）← 单独更差 |
| `NODE_RELOC=1` | **ERROR**（t=0.0）← 更差 |
| `FIX_TALUS=all` | ERROR |
| `LIG_MODE=none` | 发散 |
| `DET_MIN=0.95` @ 0.8/1.0/1.2 | 全 ERROR，**F 恒 ≈340 N** |

**⇒ `DET_MIN` 最强（0.3→0.6mm）⇒ 限流=**近退化单元**（软骨 `det(J) min=0.51`，`DET_MIN=0.5` 刚好放过）。
但所有杠杆叠加只到 ~380 N，且 `F` 不随位移增长 ⇒ **有东西触底 ⇒ 调参/调材料不可达 600 N。**」

### 19.4 两条被否证的假设（负面结果）

- **"两软骨面有 0.97mm 间隙是主因"** —— 否证。放开 `T_TARGET` 封顶重生成软骨
  （`T_TARGET=3.0` ⇒ 两层填满、**四关节穿透过全 0**、中位厚 1.879/1.909mm），**单独用反而更早失败**。
- **"距骨轴向自由（靠软韧带）是主因"** —— 否证（`FIX_TALUS=all` 无效）。

### 19.5 新增资产

`AXIS_SRC`/`ALIGN_AXIS` · `NODE_RELOC`/`LAUGON`/`PENALTY`（默认=既有）· `T_TARGET` env 门控 ·
`temp/thums/cart_t3/`（0 穿透软骨，可复用）· `acceptance.py` 加**距骨反力记账** ·
`baseline_align.py`（位移扫描，**自动删读数防残留冒充**）

### 19.6 教训

- 🔴 **失败运行留下的旧读数会"冒充"本次结果** —— 本轮第一次跑矩阵就中招（三个数字全是上一次的残留）。
  **⇒ 每次运行前必须 `unlink` 读数文件。**
- **"账要平"是硬检查**：`顶面+距骨约束 ≈ 施加载荷`，一次排除"力从别的约束旁路"。
- **天花板要测"是否随驱动量增长"**：`0.8/1.0/1.2` 给同一个 F，才敢说封顶。

### 19.7 下一步

`DET_MIN=0.95` 有效 ⇒ **要质量好的软骨网格**。`cart_patch.py` 现成，但输出
prism6+`tie_faces` 新格式，**需给 `ankle_lig_feb.py` 移植 TIE**（`build_feb.py` 已有可参照）。
**在跑到 600 N 之前，C2~C5 姿势判据无法验收。**
"""
if SESS.exists():
    tt = SESS.read_text(encoding="utf-8")
    if "十九(10-05)" not in tt:
        SESS.write_text(tt.rstrip() + SEC, encoding="utf-8")
        print(f"session note -> {SESS.stat().st_size:,} B")

D = WS / "memory/2026-10-04.md"
if D.exists():
    tt = D.read_text(encoding="utf-8")
    if f"LRN-{NEW}" not in tt:
        L = ("\n- `2026-10-04/ankle-fe-solve-and-handover.md` → **§十九(10-05)** "
             f"**600 N 到不了**（`[LRN-{NEW}]`）：稳定窗仅轴向 **0.1~0.3mm**（≤235 N）；硬天花板 "
             "**~340 N**（`DET_MIN=0.95` 抬到 0.6mm、叠加厚软骨 ~380 N），且 DISP=0.8/1.0/1.2 给出**同一个 F** "
             "⇒ 确实封顶。**最强杠杆=`DET_MIN`**（⇒ 限流=**近退化单元**，软骨 det(J)min=0.51）；"
             "`CART_E=1000` 可达 561 N。**否证两条假设**：间隙说（`T_TARGET=3.0` 厚软骨单独更差，"
             "0 穿透）与距骨自由说（`FIX_TALUS=all` 无效）；`NODE_RELOC=1` 更差；`LIG_MODE=none` 发散。"
             "记账已核实（距骨约束反力=0，`tal_dofs` 实为 `xy`）。教训：**失败运行的旧读数会冒充结果 ⇒ "
             "每次运行前 unlink**\n")
        D.write_text(tt.rstrip() + L, encoding="utf-8")
        print(f"daily note -> {D.stat().st_size:,} B")
print(f"⇒ LRN-{NEW}")
