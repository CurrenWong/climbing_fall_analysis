"""保存「踝子线下一个动作」计划：LRN + ANKLE_FIX_TODO + 总体进度 + session note。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
LT = ROOT / ".learnings/LEARNINGS.md"

# ── 扫 ID（教训：分配前必须扫，否则守卫静默跳过整段）──
t = LT.read_text(encoding="utf-8")
have = set(re.findall(r"\[LRN-(\d{8}-\d{3})\]", t))
NEWID = next(f"20261005-{k:03d}" for k in range(57, 90)
             if f"20261005-{k:03d}" not in have)
print(f"已用 20261005-* : {sorted(x for x in have if x.startswith('20261005'))}")
print(f"⇒ 本次用 LRN-{NEWID}")

ENTRY = r"""

---

## [LRN-NEWID] 踝子线的真卡口 = **力控制病态**（距骨全固定）；并纠正韧带线结案的**循环论证**

- **日期**：2026-10-05
- **类别**：`knowledge_gap` / `correction`
- **触发**：用户把焦点收回到踝子线，并问「这一点下你应该做什么」

### 🎯 先对齐目标（**已定**，不是新决定）

`ANKLE_FIX_TODO.md`「⭐ 目标重定义（2026-10-04，用户明确）」原话：

> 「我希望投入到有一个可以**跟 OpenSim 结果配合**，判断**不同姿势骨折风险**的有限元模型程度。」

⇒ 子线定位 = 从**验证模型**升格为**预测模型**。（我曾一度建议"冻结子线"，**与这个已定目标冲突，撤回**。）

### 🔴 唯一硬卡口：模型**接不了 OpenSim 的载荷**

| 事实 | 出处 |
|---|---|
| `DRIVER=force`（600 N 面力）**会失败** | `FE_PIPELINE_TODO.md` 接手入口 |
| 原因：**胫骨只靠接触约束 ⇒ 刚体模态 ⇒ 奇异** | 同上 |
| 所以改用 `DRIVER=disp` + 顶面全约束（夹持端） | 同上 |

**OpenSim 输出的是力/力矩，不是位移。** 模型只能吃位移 ⇒
**「跟 OpenSim 配合」这一步在结构上就断了** —— 这才是子线的第一阻塞，不是接触罚因子。

### 🔴 「判不同姿势」被同一个选择锁死

现在**唯一可用**的配置：
```
LIG_MODE=shell DRIVER=disp TRUNC_MM=60 DISP_MM=0.20 FIX_TALUS=all
```
`DRIVER=disp` ⇒ 胫骨顶面 **xyz 全约束**；`FIX_TALUS=all` ⇒ 距骨低端 **xyz 全约束**
⇒ **踝成为"双端夹持"** ⇒ **内外翻 / 内旋在结构上不可能发生** ⇒ **"不同姿势"无从谈起。**

**⇒ 为求解稳定而选的那个配置，恰好把目标所需的两个能力全关掉了。**

### ⚠️ 纠正：韧带线结案是**循环论证**

原结案链（`FE_PIPELINE_TODO.md` / `docs/总体进度.md`）：

> 「本模型 ΔL = 0.001~0.08 mm，且 胫跟 缩短最多（最后一名）」
> → 「Tochigi 实测站姿韧带本就松弛（CFL 0.2~2.1%）⇒ **我 ΔL≈0 与实验一致**」
> → 「韧带线已降级 / 结案」

**问题**：`ΔL ≈ 0` 是被 **`FIX_TALUS=all` + 顶面 xyz 全约束**所**强制**的 ——
双端夹持 ⇒ 踝不能转 ⇒ **任何韧带都不会被拉伸**，**与模型好坏无关**。
⇒ 「与实验一致」**不是验证，是同义反复**。

（`docs/Anderson双论文_关键信息.md` §7.5 已半猜到「韧带线不该废掉，而是换场景」；
本条的贡献是给出**机制性**理由。⇒ **韧带线应从「已结案」改回「被错误配置挡住」。**）

### ✅ 好消息：能力**已经写好了**，只是不是"可用配置"

`scripts/ankle_fe/ankle_lig_feb.py` 已有全套（且**默认值就是按目标设的**）：

| 开关 | 默认 | 作用 |
|---|---|---|
| `DRIVER` | **`force`** | 力控制（`force` / `disp` / `dispx`） |
| `FORCE_N` | 600 | 载荷大小 |
| `TOPFIX=none`（配 `DRIVER=force`） | — | 顶面**完全不约束**；注释原文「靠接触 + 韧带定运动学（**允许外翻/内旋**）」 |
| `DRIVER=dispx` | — | 只约束轴向，**横向自由**（注释「交给韧带 + 接触约束」） |
| `FIX_TALUS` | **`transverse`** | 只约束横向、**放开轴向** |
| `TRUNC_MM` | 0 | 截断（实测 **160× 提速**，900s→5.4s） |
| `LIG_MODE` | `shell` | 8 条韧带 |

⇒ **作者默认值本身是按目标设的**；为了跑通，文档把配置改成 `disp` + `all`，
**把目标能力当代价付掉了**。所以任务不是"造能力"，是**让目标配置收敛**。

### ✅ 该做的（= Anderson 2006 §2.2 原文，**一件事解三/四个问题**）

| 步 | 内容 |
|---|---|
| ① | **距骨平移弱弹簧** `k = 500 N/mm`（2007 再加 A-P `k = 50 N/mm`） |
| ② | **落位阶段**先约束 about-内外翻 / about-内旋，**正式加载前移除** |
| ③ | 轻接触 → **200 N 落位** → 再上正式载荷（Anderson 三段流程） |
| ④ | 加载用 **`DRIVER=force`**（或 `dispx`），顶面**不全约束** |

**收益（一次改动解多问题）**：
1. **接上 OpenSim** —— 弱弹簧提供非接触约束 ⇒ 刚体模态消失 ⇒ **力控制可行**
2. **能判不同姿势** —— 距骨不再锁死 + 内外翻/内旋解除 ⇒ **踝能转动 ⇒ 姿势成为输入**
3. **关节真落位** —— 200 N 落位让软骨面真压上（31mm² → ~300mm²）⇒ 接触斑成形，罚因子影响下降
4. **重开韧带线** —— 踝能转 ⇒ 韧带应变非零且有意义（替换掉循环论证）

### 📏 验收标准（5 条，**全部相对量**）

子线绝对量有 **1.5~2× 不确定度**（`PENALTY≥2` 全负 Jacobian ⇒ 无法验收敛，
见 `[LRN-20261005-056]`）⇒ **不追绝对量**，逻辑同线 2 对「绝对高度 3~4× 漂移」的处理：
**只对标相对模式**。

| # | 判据 | 来源 |
|---|---|---|
| 1 | **`DRIVER=force` 能收敛**（`rc=0`），且 `\|F − 目标\| < 2%` | 自洽 |
| 2 | **末态出现非零的内外翻角** | "姿势"的前提 |
| 3 | **V9：内翻损伤 < 外翻损伤** | `[T22]` |
| 4 | **外侧(跟腓) < 内侧、胫跟 > 胫舟** | Anderson / 吴恺机制 |
| 5 | **⭐ 姿势排序对罚因子鲁棒**：`PENALTY ∈ {0.1, 1}` 两档下**排序不变** | 绕开无法收敛的绝对量 |

> **第 5 条是关键**：与其追 `k → 1059 N/mm`（`p≥2` 根本跑不了），
> 不如验证「**结论是否与罚因子无关**」。若两档罚因子给出同样排序 ⇒ 该结论可用。

### 🚫 明确不做（避免再挖进兔子洞）

- ❌ 不追 **2.2× 刚度差**的闭合（是量级问题，不在相对判据路径上）
- ❌ 不追**罚因子收敛**（`PENALTY≥2` 结构上跑不了）⇒ 记为**已知不确定度带**（同线 2 的「不可移植」）
- ❌ 不再**加密网格**（m2 不收敛，且已证明 tie 不是瓶颈）
- ❌ 不在 `build_feb.py`（软骨验证线）继续投入 —— **要动就动 `ankle_lig_feb.py`**，它才是目标的载体

### 📦 落地顺序（下一步开工时照做）

1. **先写验收脚本**（上表 5 条），对**当前配置**跑出"失败" ⇒ 自证判据有效（防止又一轮无标准的取证）
2. 在 `ankle_lig_feb.py` 加 `TAL_SPRING`（弱弹簧 k）+ `REL_IE_IR`（落位后解除转动）两个开关
3. 用 `DRIVER=force` 跑通 ⇒ 验第 1~2 条
4. 再做姿势扫描 ⇒ 验第 3~5 条
""" .replace("LRN-NEWID", f"LRN-{NEWID}")

LT.write_text(t.rstrip() + ENTRY, encoding="utf-8")
s2 = LT.read_text(encoding="utf-8")
ok = f"LRN-{NEWID}" in s2 and "循环论证" in s2
print(f"LEARNINGS.md -> {LT.stat().st_size:,} B  [{'OK' if ok else 'XX'}]")

# ── ANKLE_FIX_TODO.md：在「⭐ 目标重定义」后插入「下一步（已定）」 ──
AF = ROOT / "ANKLE_FIX_TODO.md"
s = AF.read_text(encoding="utf-8")
BLOCK = f"""
## ⭐ 下一步（已定，2026-10-05 → `[LRN-{NEWID}]`）

### 唯一硬卡口：**模型接不了 OpenSim 的载荷**
`DRIVER=force` 会失败（胫骨只靠接触约束 ⇒ **刚体模态 ⇒ 奇异**）。
⇒ OpenSim 给的是**力**，模型只吃**位移** ⇒「跟 OpenSim 配合」结构上断了。

### 且「判不同姿势」被同一个选择锁死
可用配置 `DRIVER=disp ... FIX_TALUS=all` ⇒ 顶面 xyz + 距骨 xyz 全约束
⇒ **双端夹持 ⇒ 内外翻/内旋结构上不可能**。

### ⇒ 要做（= Anderson 2006 §2.2 原文，一件事解多问题）
1. **距骨平移弱弹簧 500 N/mm**（2007 再加 A-P 50 N/mm）
2. **落位阶段**约束 about-内外翻/about-内旋，**正式加载前移除**
3. 轻接触 → **200 N 落位** → 正式载荷
4. **`DRIVER=force`**（或 `dispx`），顶面**不全约束**

**收益**：① 消刚体模态 ⇒ **力控制可行 ⇒ 接上 OpenSim** ② **踝能转 ⇒ 姿势成为输入**
③ 200N 落位 ⇒ 关节真压上（31→~300mm²）④ **重开韧带线**（ΔL≈0 是双端夹持强制的，属**循环论证**）

### 验收（5 条，全相对量 —— 绝对量有 1.5~2× 不确定度，不追）
1. `DRIVER=force` rc=0 且 |F−目标|<2%
2. 末态有**非零内外翻角**
3. **V9：内翻损伤 < 外翻**
4. **外侧(跟腓)<内侧、胫跟>胫舟**
5. **⭐ 姿势排序对罚因子鲁棒**（`PENALTY∈{{0.1,1}}` 排序不变）

### 不做
❌ 2.2× 闭合 ❌ 罚因子收敛（记为不确定度带）❌ 加密网格 ❌ 在 `build_feb.py` 上继续投入
（**要动就动 `ankle_lig_feb.py`**，它已内置 `DRIVER=force`/`TOPFIX=none`/`FIX_TALUS=transverse`/`TRUNC_MM`）

### 开工顺序
1. 先写验收脚本（对当前配置跑出"失败"以自证判据有效）
2. `ankle_lig_feb.py` 加 `TAL_SPRING` / `REL_IE_IR` 开关
3. `DRIVER=force` 跑通 → 验 1~2 条；再做姿势扫描 → 验 3~5 条

---
"""
if "下一步（已定，2026-10-05" not in s:
    i = s.find("## ⭐ 目标重定义")
    if i < 0:
        AF.write_text(s + "\n" + BLOCK, encoding="utf-8")
    else:
        j = s.find("\n---", i)
        AF.write_text(s[:j] + "\n" + BLOCK + s[j:], encoding="utf-8")
print(f"ANKLE_FIX_TODO.md -> {AF.stat().st_size:,} B  "
      f"[{'OK' if '下一步（已定，2026-10-05' in AF.read_text(encoding='utf-8') else 'XX'}]")

# ── docs/总体进度.md：子线行 + 待决策 ──
OP = ROOT / "docs/总体进度.md"
s = OP.read_text(encoding="utf-8")
OLD = "| **⚠️ F 未收敛（非伪影）** | 🟠 |"
NEW = (f"| **⭐ 下一步（已定）** | 🎯 | 目标=**接 OpenSim + 判不同姿势**。唯一卡口=**力控制病态**"
       f"（距骨全固定 ⇒ 刚体模态）；且 `FIX_TALUS=all`+顶面全约束 ⇒ **双端夹持 ⇒ 踝不能内外翻**。"
       f"解=**Anderson §2.2 距骨弱弹簧 500/50 N/mm + 解除内外翻内旋**（一件事解多问题）。"
       f"验收=**姿势排序对罚因子鲁棒**。`[LRN-{NEWID}]` |\n| **⚠️ F 未收敛（非伪影）** | 🟠 |")
if OLD in s and f"LRN-{NEWID}" not in s:
    OP.write_text(s.replace(OLD, NEW, 1), encoding="utf-8")
s = OP.read_text(encoding="utf-8")
OLD2 = "| 3 | 踝关节子线**继续投入到什么程度**？ |"
if OLD2 in s:
    pass  # 已有决策，不改
print(f"docs/总体进度.md -> {OP.stat().st_size:,} B  "
      f"[{'OK' if f'LRN-{NEWID}' in OP.read_text(encoding='utf-8') else 'XX'}]")

# ── session note §十六 ──
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
SEC = f"""

---

## 十六、踝子线的下一个动作（2026-10-05）→ `[LRN-{NEWID}]`

### 16.1 目标（已定，非新决定）
用户 2026-10-04 原话：「我希望投入到有一个可以**跟 OpenSim 结果配合**，
判断**不同姿势骨折风险**的有限元模型程度。」
⇒ 子线 = 从验证模型升格为**预测模型**。（我曾建议"冻结子线"，**与已定目标冲突，撤回**。）

### 16.2 唯一硬卡口 = **力控制病态**
- `DRIVER=force`（600N 面力）**失败**：胫骨只靠接触约束 ⇒ **刚体模态 ⇒ 奇异**
- **OpenSim 给的是力/力矩，模型只吃位移** ⇒ **「跟 OpenSim 配合」结构上断了**
- ⇒ **这才是子线第一阻塞，不是接触罚因子**

### 16.3 「判不同姿势」被同一个选择锁死
唯一可用配置 `DRIVER=disp ... FIX_TALUS=all`
⇒ 顶面 xyz + 距骨 xyz **全约束** ⇒ **双端夹持** ⇒ **内外翻/内旋结构上不可能**。
⇒ **为求解稳定而选的配置，恰好把目标所需的两个能力全关掉了。**

### 16.4 ⚠️ 纠正：韧带线结案是**循环论证**
「ΔL=0.001~0.08mm ⇒ 与 Tochigi 实验一致 ⇒ 结案」
—— 但 **ΔL≈0 是 `FIX_TALUS=all`+顶面全约束强制的**（双端夹持 ⇒ 踝不能转 ⇒
韧带必然不拉伸），**与模型好坏无关**。⇒「与实验一致」是同义反复。
⇒ **韧带线应从「已结案」改回「被错误配置挡住」。**

### 16.5 ✅ 能力已写好，只是不是"可用配置"
`ankle_lig_feb.py` 已内置（**默认值本就按目标设**）：
`DRIVER=force`（默认！）/`dispx`、`FORCE_N=600`、`TOPFIX=none`（注释原文「顶面完全不约束…
**允许外翻/内旋**」）、`FIX_TALUS=transverse`（默认，放开轴向）、`TRUNC_MM`（**160× 提速**，5.4s）、
`LIG_MODE=shell`。⇒ 任务不是造能力，是**让目标配置收敛**。

### 16.6 要做的（Anderson 2006 §2.2 原文）
① 距骨平移弱弹簧 **500 N/mm**（2007 加 A-P **50 N/mm**）
② 落位阶段约束 about-内外翻/about-内旋，**正式加载前移除**
③ 轻接触 → **200 N 落位** → 正式载荷
④ **`DRIVER=force`**（或 `dispx`），顶面**不全约束**

**收益**：①消刚体模态⇒**力控制可行⇒接上 OpenSim** ②**踝能转⇒姿势成为输入**
③200N 落位⇒关节真压上（31→~300mm²）④**重开韧带线**

### 16.7 验收（5 条，全相对量）
1. `DRIVER=force` rc=0 且 |F−目标|<2%
2. 末态有**非零内外翻角**
3. **V9：内翻损伤 < 外翻**
4. **外侧(跟腓)<内侧、胫跟>胫舟**
5. **⭐ 姿势排序对罚因子鲁棒**（`PENALTY∈{{0.1,1}}` 排序不变）

> 第 5 条关键：与其追 `k→1059`（`p≥2` 跑不了），不如验证「**结论与罚因子无关**」。

### 16.8 不做
❌ 2.2× 闭合　❌ 罚因子收敛（记为**已知不确定度带**，同线 2「不可移植」）❌ 加密网格
❌ **在 `build_feb.py` 继续投入**（要动就动 `ankle_lig_feb.py`）

### 16.9 开工顺序
1. **先写验收脚本**（5 条），对当前配置跑出"失败" ⇒ 自证判据有效
2. `ankle_lig_feb.py` 加 `TAL_SPRING`（弱弹簧 k）+ `REL_IE_IR`（落位后解除转动）
3. `DRIVER=force` 跑通 → 验 1~2 条；再做姿势扫描 → 验 3~5 条
"""
if SESS.exists():
    t = SESS.read_text(encoding="utf-8")
    if "十六、踝子线的下一个动作" not in t:
        SESS.write_text(t.rstrip() + SEC, encoding="utf-8")
    print(f"session note -> {SESS.stat().st_size:,} B")

DAILY = WS / "memory/2026-10-04.md"
if DAILY.exists():
    t = DAILY.read_text(encoding="utf-8")
    if "踝子线的下一个动作" not in t:
        DAILY.write_text(
            t.rstrip() + f"\n- `2026-10-04/ankle-fe-solve-and-handover.md` → §十六 "
            f"踝子线下一步（`[LRN-{NEWID}]`）：卡口=**力控制病态**（距骨全固定⇒刚体模态），"
            f"解=**Anderson §2.2 弱弹簧 500/50 + 解除内外翻内旋**；"
            f"并纠正韧带线结案的**循环论证**（ΔL≈0 是双端夹持强制的）\n", encoding="utf-8")
    print(f"daily note -> {DAILY.stat().st_size:,} B")
print(f"\n⇒ 本次 LRN = [LRN-{NEWID}]")
