"""落盘 [LRN-060]：加载方向不可靠 + ALIGN_AXIS 正解；并更新 TODO / MEMORY / session / daily。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
LT = ROOT / ".learnings/LEARNINGS.md"

t = LT.read_text(encoding="utf-8")
have = set(re.findall(r"\[LRN-(\d{8}-\d{3})\]", t))
NEW = next(f"20261005-{k:03d}" for k in range(60, 99) if f"20261005-{k:03d}" not in have)
print(f"⇒ LRN-{NEW}")

ENTRY = r"""

---

## [LRN-NEW] 踝子线：**加载方向不可靠** —— `dispx` 只按全局单轴 + `TRUNC_MM` 把轴转 12.18°；`ALIGN_AXIS` 是正解

- **日期**：2026-10-05
- **类别**：`correction` / `best_practice`
- **承接**：`[LRN-20261005-059]`（"未解释的反力方向异常"）

### 🔴 先把上一条的"异常"结掉：**不是模型问题，是我自己的运行器不一致**

`[LRN-059]` 记录了一个"未能解释"的现象：`FIX_TALUS=transverse` 的反力在 **x**、
`all` 的反力在 **z**，而顶面 BC 同为 `dispx`。

**真因（已定论）**：`acceptance.py` 当时**没有给 `TRUNC_MM` 兜默认值**
⇒ 一次跑 `TRUNC_MM=60`、一次跑 `TRUNC_MM=0`（脚本默认）⇒ **两次的加载轴不同**：

| 运行 | TRUNC_MM | 加载轴 | `argmax(\|轴\|)` | `dispx` prescribe |
|---|---|---|---|---|
| 第 1 次 | **60** | `[0.8153, 0.135, 0.5631]` | 0 = **x** | x ⇒ 反力在 x ✓ |
| 第 2 次 | **0** | `[0.6765, 0.1363, 0.7237]` | 2 = **z** | z ⇒ 反力在 z ✓ |

⇒ **代码是对的**（`_kax = argmax(|LAX|)`，prescribe 哪个轴反力就在哪个轴）；
**是我的对照实验没控制变量**。`[LRN-059]` 里那条"未解释异常"**由此撤销**。

### 🔴 但顺带暴露出两个**真**问题

#### 问题 A：`TRUNC_MM` 会把"加载轴"**转 12.18°**

截断逻辑本身用的是**全网格**轴（正确）；但加载轴在第 ~370 行**又从截断后的 `dom_nodes`
质心重算了一次**：

```python
tib_pts = XYZ[[i-1 for i in sorted(dom_nodes["tibia_cort"])]]   # 截断后
tal_pts = XYZ[[i-1 for i in sorted(dom_nodes["talus_cort"])]]
axis = tib_pts.mean(0) - tal_pts.mean(0)
```

**实测**（同一网格，只差 `TRUNC_MM`）：

| | 轴 |
|---|---|
| 截断后（`TRUNC_MM=60`） | `[0.8153, 0.135, 0.5631]` |
| **全网格（解剖学正确）** | `[0.6765, 0.1363, 0.7237]` |
| **夹角** | **12.18°** |

⇒ **同一具模型，加不加截断会改变"沿哪个方向加载"。**
**已加开关**：`AXIS_SRC=trunc`（**默认，保持既有基线**）/ `full`（用全网格轴）；
两种模式下都会**打印两轴及夹角**，让偏差可见。

#### 问题 B（更要紧）：`dispx` **只按全局单轴** 近似轴向 ⇒ **k 差 3.2 倍**

```python
_kax = int(np.argmax(np.abs(LAX)))
_dofs = "xyz" if DRIVER == "disp" else "xyz"[_kax]
```

`dispx` 只 prescribe `argmax(|轴|)` 对应的那**一个全局自由度**。对
`axis=[0.8153, 0.135, 0.5631]`（x 分量最大）它 prescribe **全局 x**——
这既不是真轴向，**其"轴"还随 `TRUNC_MM` 切换**（问题 A）。

**实测同一"0.2 mm 名义位移"下的轴向刚度**（`TRUNC_MM=60`, `FIX_TALUS=transverse`, t=1.0）：

| 配置 | 真轴向位移 | 轴向力 | **k = F/δ** |
|---|---|---|---|
| `AXIS_SRC=trunc`（prescribe 全局 x） | 0.2453 mm | 71.12 N | **290 N/mm** |
| `AXIS_SRC=full`（prescribe 全局 z） | 0.2764 mm | 114.95 N | **416 N/mm** |
| **`AXIS_SRC=full` + `ALIGN_AXIS=1`（真轴向 z）** | **0.2000 mm** | **183.20 N** | **916 N/mm** |

⇒ **3.2 倍**。文档可用的 `DRIVER=disp` 同样是"按全局 x 平移顶面"（实测
`两顶面合力=[217.95,-45.09,239.20]N`，也**不是**沿轴）。

**⇒ 结论：这条线的"轴向位移/轴力"一直带着一个**未被记录的 ~2~3 倍方向误差****。

### ✅ 正解：`ALIGN_AXIS=1`（已实现）

把**整个模型旋转到「加载轴 = +z」**（Rodrigues 公式构造旋转阵，作用到全部节点坐标），
于是 `dispx` 的"按全局单轴"**恰好等于真轴向**：

```
[align] ALIGN_AXIS=1：模型已旋转，加载轴 -> +z（旋转角 47.34deg）
加载轴（胫骨->转骨）= [0. 0. 1.]  [已对齐 +z]
▪ 两顶面合力 = [0. 0. 183.2044] N   |F|=183.2044
▪ 投影到加载轴 = 183.2044 N   (轴=[0. 0. 1.])   ← 投影==模，无歧义
```

- **默认关**（`ALIGN_AXIS=0`）⇒ 不改变既有行为
- 开启后**所有输出坐标/位移都在旋转后的坐标系**（轴=+z，便于比较与扫姿势）
- 与 `AXIS_SRC=full` 搭配使用才自洽（先取全网格轴，再对齐）

> 附带一致性观察：`k=916 N/mm` 与 Anderson 独立估计 `E·A/t≈1059 N/mm` **同量级**（但材料集
> `MATSET=WK` 骨 7300 与 Anderson「骨刚性」不同、BC 也不同）⇒ **仅记为弱佐证，不作结论**。

### 🧠 方法论教训

- 🔴 **对照实验必须回显"权威运行参数（含默认值）"**：本轮"异常"纯粹来自"一次用了默认
  `TRUNC_MM=0`、另一次显式 60"。**已把脚本自身 stdout 的 `[mode]/[axis]/[trunc]/[data]`
  行在验收报告里原样回显** ⇒ 以后任何两次运行是否同参，一眼可查。
- **"近似"必须有名字并打印**：`dispx` 的 argmax 近似此前没有任何提示，导致 3.2 倍的方向
  误差长期不可见。**⇒ 凡近似，必须在 stdout 里显式声明它选了哪条路径。**
- **开关默认值不要悄悄改基线**：`AXIS_SRC`/`ALIGN_AXIS` 一律默认"保持旧行为"，
  并在打印里说明"置 X=Y 可换"。
""" .replace("LRN-NEW", f"LRN-{NEW}")

LT.write_text(t.rstrip() + ENTRY, encoding="utf-8")
print(f"LEARNINGS.md -> {LT.stat().st_size:,} B")

# ── TODO ──
TD = ROOT / "ANKLE_GOAL_TODO.md"
s = TD.read_text(encoding="utf-8")
s = s.replace(
    "- [ ] **Step 2a2** 修 `dispx` 语义：当前只 prescribe `argmax(|LAX|)` 的**全局单轴**，非真轴向",
    f"- [x] **Step 2a2** ~~修 `dispx` 语义~~ **已做**：新增 `ALIGN_AXIS=1`（旋转模型到轴=+z，"
    f"`dispx` 即精确轴向）+ `AXIS_SRC=trunc|full`；实测方向误差致 k 差 **3.2 倍**（290/416/916 N/mm）。`[LRN-{NEW}]`")
s = s.replace(
    "- [ ] **Step 2a3** 查清异常：`transverse` 反力在 x、`all` 反力在 z（顶面 BC 同为 prescribe x）",
    f"- [x] **Step 2a3** ~~查清反力方向异常~~ **已结案**：**不是模型问题** —— `acceptance.py` 当时"
    f"没兜 `TRUNC_MM`，一次 60 一次 0 ⇒ 加载轴不同 ⇒ `dispx` prescribe 的全局单轴不同（x vs z）。"
    f"**同时暴露真问题**：`TRUNC_MM=60` 把轴转 **12.18°**。`[LRN-{NEW}]`")
TD.write_text(s, encoding="utf-8")
print(f"ANKLE_GOAL_TODO.md -> {TD.stat().st_size:,} B")

# ── MEMORY ──
M = WS / "MEMORY.md"
s = M.read_text(encoding="utf-8")
OLD = """  ⚠️ 韧带线结案是**循环论证**（ΔL≈0 由双端夹持强制）⇒ 待重开"""
NEWLINE = """  ⚠️ 韧带线结案是**循环论证**（ΔL≈0 由双端夹持强制）⇒ 待重开
- ⚠️ **加载方向不可靠**（`[LRN-060]`）：`dispx` 只 prescribe `argmax(|轴|)` 的全局单轴，且 `TRUNC_MM=60`
  把轴转 **12.18°** ⇒ 同 0.2mm 名义位移的 k 在 **290/416/916 N/mm** 跳（**3.2×**）。**正解=`ALIGN_AXIS=1`**"""
TRIMS = [
    ("""- **中文优先**（专有名词留英文）；**先行动后报告**；**宁缺毋滥**（搜索硬 3-5 条）；
  **保留原文**（搬运不「创造性总结」，早报除外）
- **emoji 用在合适处**（Discord/Slack 表情回应）；正文中不用
- **重要数字全列**（人名/航班号/ID）；**负面结果也如实报**（标"假设错误 + 实测数据"）""",
     """- **中文优先**（专有名词留英文）；**先行动后报告**；**宁缺毋滥**（搜索硬 3-5 条）；**保留原文**
- **emoji 用在合适处**（Discord/Slack 表情回应）；正文中不用
- **重要数字全列**；**负面结果也如实报**（标"假设错误 + 实测数据"）"""),
    ("""- **状态**：v1 模型（1D 纯压缩，97 测试全绿）+ **P0、P1 已完成**；
  下一步 **P3 二维运动学**（升格为必要项）或 **P5 抱石专用标定**（待用户选）""",
     """- **状态**：v1（1D 纯压缩，97 测试全绿）+ **P0、P1 已完成**；下一步 **P3 二维运动学**（必要项）或 **P5 抱石标定**（待选）"""),
    ("""- **P1 硬结果**：真实抱石场景**一例骨折都报不出来**（1000 次 MC，最高利用率 0.725）⇒ 损伤部位分布**空**（覆盖 100%→无旋转 32.6%→可信 13.8%→可产类型 19.3%）""",
     """- **P1 硬结果**：抱石场景**一例骨折都报不出来**（1000 MC，最高利用率 0.725）⇒ 损伤分布**空**（无旋转 32.6%→可产类型 19.3%）"""),
]
if OLD in s:
    s = s.replace(OLD, NEWLINE, 1)
    print("✓ MEMORY 踝段已加 加载方向 条目")
else:
    print("⚠ 未匹配 MEMORY 踝段锚点")
for a, b in TRIMS:
    if a in s:
        s = s.replace(a, b, 1)
    else:
        print(f"⚠ 未匹配削减: {a[:32]}...")
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'OK' if sz <= 10240 else 'XX 超'}")

# ── session note ──
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
SEC = f"""

---

## 十八(10-05)、加载方向不可靠 → `ALIGN_AXIS` 正解（2026-10-05）→ `[LRN-{NEW}]`

### 18.1 先结掉 §十七 留下的"未解释异常"：**是我运行器的锅**

§十七 记录 `transverse` 反力在 x、`all` 反力在 z（顶面 BC 同为 `dispx`）。**定论**：
`acceptance.py` 当时**没给 `TRUNC_MM` 兜默认值** ⇒ 一次 `TRUNC_MM=60`、一次 `0`（脚本默认）
⇒ **加载轴不同** ⇒ `_kax=argmax(|轴|)` 选中不同的全局自由度：

| 运行 | TRUNC_MM | 加载轴 | `argmax` | prescribe | 反力 |
|---|---|---|---|---|---|
| 1 | **60** | `[0.8153, 0.135, 0.5631]` | 0 = x | x | x ✓ |
| 2 | **0** | `[0.6765, 0.1363, 0.7237]` | 2 = z | z | z ✓ |

⇒ **代码正确；对照实验没控制变量。§十七 的"未解释异常"撤销。**

### 18.2 但暴露两个真问题

**A. `TRUNC_MM` 把加载轴转 12.18°**
截断用的是全网格轴（对），但加载轴在截断后**又重算了一次**（第 ~370 行从截断后 `dom_nodes`）。
实测：截断后 `[0.8153,.135,.5631]` vs 全网格 `[0.6765,.1363,.7237]`，**夹角 12.18°**。
⇒ **加不加截断会改变"沿哪个方向加载"**。已加 `AXIS_SRC=trunc`（默认，保基线）/ `full`，
两种模式都打印两轴与夹角。

**B. `dispx` 只按全局单轴 ⇒ k 差 3.2 倍**

| 配置 | 真轴向位移 | 轴向力 | k=F/δ |
|---|---|---|---|
| `trunc`（prescribe 全局 x） | 0.2453 mm | 71.12 N | **290 N/mm** |
| `full`（prescribe 全局 z） | 0.2764 mm | 114.95 N | **416 N/mm** |
| **`full` + `ALIGN_AXIS`（真轴向）** | **0.2000 mm** | **183.20 N** | **916 N/mm** |

文档"可用"的 `disp` 也是按全局 x 平移顶面（实测合力 `[217.95,-45.09,239.20]N`，非沿轴）。
⇒ **这条线的"轴向位移/轴力"一直带着未被记录的 ~2~3 倍方向误差。**

### 18.3 正解：`ALIGN_AXIS=1`（已实现，默认关）

Rodrigues 旋转阵把**全模型**转到「加载轴 = +z」，`dispx` 的全局单轴**恰好等于真轴向**：
```
[align] ALIGN_AXIS=1：模型已旋转，加载轴 -> +z（旋转角 47.34deg）
▪ 两顶面合力 = [0. 0. 183.2044] N    ▪ 投影到加载轴 = 183.2044 N (轴=[0,0,1])  ← 投影==模
```
默认关 ⇒ 不改基线；开启后**输出坐标/位移在旋转后坐标系**。

> `k=916 N/mm` 与 Anderson 独立估计 **1059 N/mm** 同量级（材料集/BC 不同）⇒ **弱佐证，不作结论**。

### 18.4 方法论（新增两条硬规则）

1. **对照实验必须回显"权威运行参数（含默认值）"** —— 验收报告已改为**原样回显脚本 stdout
   的 `[mode]/[axis]/[trunc]/[data]` 行** ⇒ 两次运行是否同参一眼可查。
2. **"近似"必须有名字并打印** —— `dispx` 的 argmax 近似此前零提示，导致 3.2 倍方向误差
   长期不可见。**凡近似，必须在 stdout 显式声明它选了哪条路径。**
3. **开关默认值不悄悄改基线** —— `AXIS_SRC`/`ALIGN_AXIS` 默认保持旧行为，并打印"置 X=Y 可换"。
"""
if SESS.exists():
    tt = SESS.read_text(encoding="utf-8")
    if "十八(10-05)" not in tt:
        SESS.write_text(tt.rstrip() + SEC, encoding="utf-8")
        print(f"session note -> {SESS.stat().st_size:,} B")

D = WS / "memory/2026-10-04.md"
if D.exists():
    tt = D.read_text(encoding="utf-8")
    if f"LRN-{NEW}" not in tt:
        L = ("\n- `2026-10-04/ankle-fe-solve-and-handover.md` → **§十八(10-05)** 加载方向不可靠 "
             f"（`[LRN-{NEW}]`）：**§十七 的反力方向异常结案=我运行器没兜 `TRUNC_MM`**"
             "（一次 60 一次 0 ⇒ 加载轴不同 ⇒ `dispx` 选的全局单轴不同，代码是对的）；"
             "但暴露真问题：**`TRUNC_MM=60` 把加载轴转 12.18°**（截断后又重算了一次轴），"
             "且 **`dispx` 只按 `argmax(|轴|)` 的全局单轴** ⇒ 同 0.2mm 名义位移的 k 在 "
             "**290/416/916 N/mm** 跳（**3.2×**）。**正解 `ALIGN_AXIS=1`**（旋转模型到轴=+z，"
             "投影==模无歧义，默认关）；新增规则：**验收报告回显权威运行参数**、**近似必须打印**\n")
        D.write_text(tt.rstrip() + L, encoding="utf-8")
        print(f"daily note -> {D.stat().st_size:,} B")
print(f"⇒ LRN-{NEW}")
