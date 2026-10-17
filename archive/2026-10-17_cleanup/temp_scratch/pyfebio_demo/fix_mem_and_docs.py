"""① 把 C2/应变细节搬进 references/ankle-fe-criteria.md
② MEMORY 里换成紧凑指针（回到 ≤10240）
③ 修 ANKLE_FIX_TODO 的残留错误说法
"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
REF = os.path.join(W, "references", "ankle-fe-criteria.md")
MEM = os.path.join(W, "MEMORY.md")
TOD = os.path.join(R, "ANKLE_FIX_TODO.md")
USE = os.path.join(R, "docs", "踝FE可用模型使用说明.md")

REF_TXT = """# 踝 FE 判据（C1-C5）口径与实测

> 2026-10-06 建立。主 `MEMORY.md` 只留指针；细节在此。
> 关联：`[LRN-20261006-074/075/076]`、`[ERR-20261006-001/002]`

## C2 重定义（已完成，PASS）

**原 C2「末态出现非零内外翻角」为什么废**：姿势是把足部组**烘焙进参考构型**的
⇒ t=0 已在位 ⇒ 测"t=0→末态转动"结构上得 0（实测 0.0014~0.0033°）。
**架构错位**（Anderson 的姿势是**输出**，我的是**输入**），不是模型错。

**新 C2 = 验"姿势作为输入被忠实实现且语义方向正确"**：

| 子判据 | 内容 | 实测 |
|---|---|---|
| C2a-1 幅值 | \\|实测相对转角 − 请求角\\| ≤ 0.5° | **0.000°** ✓ |
| C2a-2 轴向 | 实测轴与 AP **平行**（`\\|轴·AP\\|→1`） | **1.0000** ✓ |
| C2a-3 语义 | 三锚定判据一致 | **全一致、单调线性** ✓ |

**三个锚定判据**（关键：至少一个不依赖模型自身符号约定）
1. 足底法向外侧分量减小 ⇒ 内翻
2. 足外侧缘下沉 ⇒ 内翻
3. 外侧韧带(ATFL/CFL/PTFL)被拉长 ⇒ 内翻

**"外侧"方向由腓骨定义**：`lat = 腓骨顶质心 − 胫骨顶质心`（⊥加载轴，实测 `|lat·axis|=0.0000`）。
腓骨解剖上必在胫骨外侧 ⇒ 与模型左右/镜像无关 ⇒ 不会被模型自身符号约定污染。

**⇒ `+IE = 内翻`（代码注释正确）。** 先前"标签反了"是错的，根因是循环论证（见 `[ERR-20261006-002]`）。

**线性度**（ΔL vs 角度）：CFL ≈ 0.50 mm/度、ATFL ≈ 0.35、PTFL ≈ 0.13
（IE=+1/+3/+5/+10 ⇒ CFL +0.34/+1.45/+2.47/+5.07 mm）。

## 韧带应变的三个口径（`[LRN-20261006-075]`）

| 口径 | 定义 | 用途 |
|---|---|---|
| **(a) 预应变** | (L_姿势 − L_中性)/L_中性 | 姿势本身拉了多少 ← **内翻扭伤机制** |
| **(b) 增量** | (L_载荷 − L_姿势)/L_姿势 | **= FEBio t=1 元素应变**（参考构型=姿势网格）|
| **(c) 总应变** | (L_载荷 − L_中性)/L_中性 | **✅ 阈值设在它**（生理参考态=中性姿势）|

### ⚠️ (b) 与损伤风险**反相关**，不能当指标
松弛的韧带被压缩"拉直"时增量大；**已绷紧**的韧带压缩时增量≈0
⇒ 增量最大的恰是松弛（无风险）那根。

### 实测（`POSE_IE=+5` 内翻，`NORMAL` 7.6 s，0 负 J，`DISP_MM=-1.0`）

| 韧带 | 侧 | 预应变%(a) | 增量%(b) | **总应变%(c)** |
|---|---|---|---|---|
| PTFL | 外 | +4.58 | −0.17 | **+4.42** |
| **CFL** | 外 | +8.98 | +2.58 | **+11.80** |
| ATFL | 外 | +5.12 | −0.04 | **+5.08** |
| PTCL | – | 0.00 | +2.32 | +2.32 |
| SYND | – | 0.00 | −0.13 | −0.13 |
| SPRING | – | 0.00 | −0.04 | −0.04 |
| MD1 tibiotalar_post | 内 | −4.11 | +2.33 | **−1.88** |
| MD2 tibiocalcaneal | 内 | −4.30 | +3.73 | **−0.73** |
| MD3 tibionavicular | 内 | −3.44 | +0.38 | **−3.07** |
| MD4 tibiotalar_ant | 内 | −3.59 | +1.44 | **−2.20** |

- 外侧汇总：预应变 **+6.22%**、增量 +0.79%｜内侧：预应变 **−3.86%**、增量 +1.97%
- **⇒ 按 (c)：内翻下外侧 +4.42~+11.80%、内侧 −0.73~−3.07% ⇒ 解剖学正确、区分度好**

**旧结论"ATFL 几乎不动 (0.17~0.20%)"的解**：那是 (b)。ATFL 在 +5 已预应变 +5.12%
⇒ 压缩不再增加。**不是"不受力"，是"已经绷住了"。**

## 已知限制

1. **内外翻不对称**（`[LRN-20261006-076]`）：`+5` 内翻 `NORMAL`；**`−5` 外翻不收敛（1 负 J，t=0）**
   ⇒ **目前只有内翻侧可用**，"判姿势"覆盖范围受限。
2. **韧带建模长度偏差**：ATFL **36.99 mm** vs 解剖 15~20（**≈2× 过长**）、CFL 34.33 vs ~28.5
   ⇒ (c) 的**绝对值被稀释**（这是 ② 要修的第一项）。
3. **姿势预应变是模型化的**：韧带零应变态定义在**中性姿势** ⇒ 任何姿势位移都会产生预应变。
   真实韧带会适应。**⇒ (c) 的绝对值需在修完插入点几何后再谈定量。**

## 工具

| 脚本 | 作用 |
|---|---|
| `scripts/ankle_fe/check_c2.py` | C2a 全检（`NO_RUN=1` 只建不解，零求解成本）|
| `scripts/ankle_fe/lig_strain.py` | **权威应变读数**（三口径，从几何直算）|
| `scripts/ankle_fe/acceptance.py` | 验收/读数（`.feb` 解析 + Kabsch + 反力记账）|
"""

io.open(REF, "w", encoding="utf-8").write(REF_TXT)
print(f"✓ {REF}  {len(REF_TXT.encode('utf-8'))} B")

# ── MEMORY：紧凑替换 ──
m = io.open(MEM, encoding="utf-8").read()
OLD = """- ✅ **C2 已重定义并 PASS**（`[LRN-074]`）：姿势是**输入** ⇒ 不能测"末态转动"；改为
  幅值(误差 **0.000°**)/轴向(`|轴·AP|=1.0000`)/**语义三锚定**（足底法向·足缘高度·韧带长度，
  "外侧"由**腓骨**定义）。**⇒ `+IE=内翻`、代码注释本来就对**；我先前"标签反了"是错的（`[ERR-002]`，
  根因=用数据反推标签的**循环论证**）
- 🎯 **韧带应变三口径**（`[LRN-075]`）：**(b)载荷增量与损伤风险反相关**（松弛韧带被拉直时增量最大）
  ⇒ **阈值只能用 (c) 总应变（相对中性姿势）**。内翻(+5, NORMAL)：外侧 CFL **+11.80%**/ATFL +5.08%/PTFL +4.42%、
  内侧 −0.73~−3.07% ⇒ 解剖学正确。旧"ATFL 几乎不动"= (b) 小（它本就已绷住 +5.12% 预应变），不是不受力。
  **⚠️ 内外翻不对称：`−5`(外翻) 不收敛（1 负 J）⇒ 目前只有内翻侧可用**（`[LRN-076]`）"""
NEW = """- ✅ **C2 已重定义并 PASS**：`+IE=内翻`（代码注释本来就对，我先前"标签反了"是**错的**）
- 🎯 **韧带应变三口径：阈值只能用 (c) 总应变**；(b) 载荷增量**与损伤风险反相关**。
  **⚠️ 外翻 `−5` 不收敛 ⇒ 目前只有内翻侧可用**。
  **判据细节全部 → `references/ankle-fe-criteria.md`**（C2 三锚定 / 三口径实测表 / 限制 / 工具）"""
if OLD in m:
    m = m.replace(OLD, NEW)
    print("· MEMORY 紧凑替换完成")
else:
    print("⚠ MEMORY 目标串未命中（可能已被上一步写坏），改为直接截断检测")
io.open(MEM, "w", encoding="utf-8").write(m)
print(f"MEMORY.md -> {len(m.encode('utf-8'))} B "
      f"{'✅' if len(m.encode('utf-8')) <= 10240 else '❌ 仍超限'}")

# references 表加一行
OLD_R = "| `references/fe-modeling-gotchas.md` | **FE 陷阱**"
NEW_R = ("| `references/ankle-fe-criteria.md` | **踝 FE 判据（C1-C5）**：C2 三锚定 / "
         "韧带应变三口径 / 实测表 / 限制 |\n| `references/fe-modeling-gotchas.md` | **FE 陷阱**")
m2 = io.open(MEM, encoding="utf-8").read()
if "ankle-fe-criteria" in m2:
    print("· references 表已含该行")
elif OLD_R in m2:
    m2 = m2.replace(OLD_R, NEW_R)
    io.open(MEM, "w", encoding="utf-8").write(m2)
    print(f"· references 表已加行 -> {len(m2.encode('utf-8'))} B")

# ── ANKLE_FIX_TODO 残留更正 ──
t = io.open(TOD, encoding="utf-8").read()
FIXES = [
    ("- **⚠️ 符号更正**：`POSE_IE` 打印写\"内翻+/−外翻\"，实测 **+IE = 外翻**，**−IE 才是内翻**",
     "- **⚠️ 符号更正（2026-10-06 已再更正）**：先前写\"+IE = 外翻\"是**错的** —— "
     "`check_c2.py` 的两个**解剖锚定**判据（足底法向 / 足缘高度，\"外侧\"由**腓骨**定义）"
     "证明 **+IE = 内翻**，**代码注释正确**。错因=用韧带数据反推标签的**循环论证**"
     "（`[ERR-20261006-002]`）。真正的问题是引用了**错的口径**（旧表用的是 (b) 载荷增量，"
     "见 `[LRN-20261006-075]`）"),
]
n = 0
for old, new in FIXES:
    if old in t:
        t = t.replace(old, new)
        n += 1
    else:
        print(f"  ⚠ 未命中：{old[:44]}")
io.open(TOD, "w", encoding="utf-8").write(t)
print(f"· ANKLE_FIX_TODO 改写 {n} 处 -> {len(t.encode('utf-8'))} B")

for f, nm in ((MEM, "MEMORY"), (TOD, "ANKLE_FIX_TODO"), (USE, "使用说明")):
    x = io.open(f, encoding="utf-8").read()
    for pat in ("实际是外翻", "−IE 才是内翻", "-IE 才是内翻"):
        if pat in x:
            print(f"  ⚠ {nm} 仍残留『{pat}』×{x.count(pat)}")
