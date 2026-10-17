"""① 把 FEBio 三坑搬进 references/fe-modeling-gotchas.md；② 主 MEMORY 改短指针；③ 校验尺寸。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
REF = r"C:\Users\curre\.copaw\workspaces\default\references\fe-modeling-gotchas.md"
M = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"

APPEND = """

## 🧨 FEBio 静默坑（2026-10-05 踝线实测，三条都会**静默**出错）

### 1. 载荷曲线必须挂 `<scale>`，不能挂 `<traction>`（`[ERR-20261005-008]`）
- 错：`<traction lc="2">-0,-0,-9.93</traction>` ⇒ **`lc` 被静默忽略 ⇒ 满载荷在第一增量瞬加**
- 对（手册 §3.13.2.2 / p121）：
  ```xml
  <surface_load type="traction" surface="tibia_top">
    <scale lc="2">1.0</scale>
    <traction>0,0,-9.93</traction>
  </surface_load>
  ```
- **诊断信号**：失败点**对总载荷量级不敏感**（100/600/5000 N 都在第一增量失败）
  ⇒ 必是"满值瞬加"而非"斜坡"。

### 2. `symmetric_stiffness` 必须在 `<solver>` 内，不能在 `<Control>` 直下
- 放 Control 直下 ⇒ `tag "symmetric_stiffness" : unrecognized tag`
- **刚体连接件（`<Rigid>`/rigid spring）的硬要求**，FEBio 自己会警告
  *"Rigid connectors require non-symmetric stiffness matrix."*
- 缺失时：弹簧约束力算错 ⇒ 实测 step2 立刻 47 个负 Jacobian。修好后刚度矩阵 nnz 258222→507067。
- 顺带也是 sliding-elastic 接触推荐的设置。

### 3. runner 必须"跑前删旧 `.feb` + 构建后查 rc"（`[ERR-20261005-010]`）
- **构建失败后 FEBio 会读上一次留下的 `.feb`** ⇒ 多个不同配置给出**逐字相同**的结果
  （曾在姿势矩阵里被骗：8 个 `POSE_TILT` 全给 `Fz=−608.82`）。
- **"结果逐字相同"永远是第一嫌疑。**
- 另：**取反力必须看符号**（压缩为负）；只看幅值会把**拉伸**当成成功（`[ERR-20261005-009]`）。

### 4. 多步 `<Step>`（手册 §3.17 / §6.1）
- **主干 `<Boundary>` 贯穿所有 step；step 内的只在该 step 生效**
- **多步时 `<Control>` 不再放文件顶部**，改放**文件底部** `<Step>` 内
- `<SEAT_HOLD>`（step2 同时保持位移约束）是**过约束**，会给出垃圾反力（实测 28,599 N）
"""

txt = io.open(REF, encoding="utf-8").read() if os.path.exists(REF) else ""
if "FEBio 静默坑" in txt:
    print("references 已有该节，跳过")
else:
    io.open(REF, "a", encoding="utf-8").write(APPEND)
    print(f"references/fe-modeling-gotchas.md: {len(txt)} -> {len(txt)+len(APPEND)}")

# 主 MEMORY：长块 → 短指针
s = io.open(M, encoding="utf-8").read()
OLD = """- 🧨 **FEBio 三个坑（都会静默出错）**：① `lc` 必须挂 `<scale>`（挂 `<traction>` 被忽略 ⇒ 满载荷瞬加）
  ② `symmetric_stiffness` 必须在 `<solver>` 内（刚体连接件必需）③ 跑前删旧 `.feb`、构建后查 rc（否则 FEBio 读旧文件 ⇒ 结果"逐字相同"假象）；**取反力必须看符号**。`[ERR-20261005-008/009/010]`"""
NEW = "- 🧨 **FEBio 静默坑**（`lc`→`<scale>`｜`symmetric_stiffness` 在 `<solver>`｜跑前删旧 `.feb`+查 rc｜反力看符号）\n  → `references/fe-modeling-gotchas.md`；`[ERR-20261005-008/009/010]`"
assert s.count(OLD) == 1, s.count(OLD)
s = s.replace(OLD, NEW)
io.open(M, "w", encoding="utf-8").write(s)
n = len(s.encode("utf-8"))
print(f"MEMORY.md -> {n} / 10240  余量 {10240 - n}")
