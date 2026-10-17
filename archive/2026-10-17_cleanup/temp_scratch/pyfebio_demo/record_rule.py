"""把反思写成 durable 规则：LRN + references + MEMORY 指针 + SOUL 行为准则。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
LEARN = os.path.join(R, ".learnings", "LEARNINGS.md")
GOT = os.path.join(W, "references", "operational-gotchas.md")
MEM = os.path.join(W, "MEMORY.md")

# 先找 SOUL.md / AGENTS.md 在哪
cands = []
for root in (W, R, os.path.dirname(W)):
    for f in ("SOUL.md", "AGENTS.md"):
        p = os.path.join(root, f)
        if os.path.exists(p):
            cands.append(p)
print("找到：")
for c in cands:
    print(f"  {c}  {os.path.getsize(c)} B")

# ───────── 1) LRN：流程规则（含事故史） ─────────
L = """
## LRN-20261006-083 — 🧭 流程铁律：**方向/符号类结论必须来自"不可被约定翻转的证据"**

**类别**: best_practice / process（**跨项目通用**）  **关联**: `[ERR-20261005-009]`、`[ERR-20261006-002/003]`

### 事故史（用户指出"这条路之前走过"）
| 时间 | 我下的结论 | 结果 |
|---|---|---|
| 10-04 | 「力控制病态」（胫骨只靠接触约束） | 后被推翻 |
| 10-05 | `[ERR-20261005-009]`「**`DISP_MM<0` 才是压实**」+ 拉伸标定表 + 「600 N 到不了」+ 整套姿势/判据 | **错** |
| 10-06 | `[ERR-20261006-003]`「**`DISP_MM>0` 才是压实**」 | **回到最初约定** |

**净结果：我离开了本来正确的约定，在一个"反过来的世界"里做了一到两轮完整工作
（标定表、姿势矩阵、判据、"越阈点 +3~+4°"、7 个文件的更正横幅），然后回到原点。
最初的约定一直是对的。**

### 根因（两层）
1. **用"约定相关量"推断"方向"**：10-05 我给的"决定性证据"是
   「`DISP_MM=-0.25 → −132.56 N（压缩）` vs `+0.25 → +187.77 N（拉伸）`」——
   **压缩/拉伸这两个标签是从反力正负读出来的**，与圆规论证同类。
   **而当天本可以做三条与约定无关的检查**（各只花 1~2 次运行）：
   - **结构**：`.feb` 里到底写了什么（实测 `z = −DISP_MM`）
   - **几何**：材料点往哪边动（顶面 dz 是否朝距骨）
   - **物理**：距骨固定时哪个方向有接触力（`+1.0→149.59 N` vs `−1.0→0.0000 N`）
2. **不验证自己的验证器**：今天我的检查脚本自己把 `DISP=+1.0` 判成"拉开"
   （脚本内符号写反）；`compress_curve.py` 的 NORMAL/ERROR 标签是**假阴性**。
   **⇒ 没有已知正对照的判据，不算证据。**

### 三级放大（比错本身更糟）
- 该结论写进了**7 个文件顶部的"务必先读这里"横幅** + 使用说明 §1 环境变量表的「**必开**」
  ⇒ **主动误导**；
- 本次"修"它时，我只**前置**了新横幅，**没删旧的** ⇒ 3 个文件出现**两条相反的横幅叠在一起**。

### ⭐ 规则（三条，写下来防止再翻）
1. **方向/符号类结论，必须先给出一条"不可被约定翻转"的证据**（几何 / 结构 / 解剖 / 独立参照），
   才准写进文档。**只凭"反力正负"或"哪个数值更像"一律不算。**
2. **每个方向判据都要配已知正对照**：跑一个"我知道它是什么方向"的算例，
   断言检查器报得对。**没跑对照 ⇒ 不是证据。**
3. **翻转过的结论必须带"翻转史"**，引用前用同一套对照复核；**更正必须"替换"而不是"叠加"**
   （改文件时要删掉被取代的段落，不能只前置新说明）。

### 审后仍安全的部分（本次审计结论）
- **C2/C3 的应变表是"几何量"**（由节点位移算得）⇒ **与约定无关，安全** ✓
- `.feb` 的 z 符号已验证 ✓；压缩曲线由三证锚定 ✓
- 其余 7 个文件的方向论断已清干净（仅存 1 处故意的删除线更正）✓
"""
t = io.open(LEARN, encoding="utf-8").read()
if "LRN-20261006-083" in t:
    print("· LRN 已存在")
else:
    io.open(LEARN, "a", encoding="utf-8").write(L)
    print(f"LEARNINGS.md {len(t)} -> {len(t)+len(L)}")

# ───────── 2) references/operational-gotchas.md ─────────
G = """

## 🧭 方向/符号结论的铁律（`[LRN-20261006-083]`）

**症状**：同一件事（加载方向）在两轮里给出**相反**结论，各配一套完整文档 ⇒ 反复返工。

**根因**：用"约定相关量"（反力正负、数值大小）推断"方向"，而不是用不可被约定翻转的证据。

**规则**
1. 方向/符号类结论 ⇒ **必须**先有下面四类证据之一，才准写进文档：
   - **结构**：产物文件里到底写了什么（如 `.feb` 的 `<value>`）
   - **几何**：材料点往哪边动（位移记录）
   - **物理**：构造一个只有单一通路能解释的算例（如固定骨 + 测接触力）
   - **解剖/独立参照**：如"腓骨必在胫骨外侧"
2. **每个判据配已知正对照**（跑一个方向已知的算例，断言检查器报对）。没跑 ⇒ 不算证据。
3. **更正要"替换"不要"叠加"**：改文档时必须删掉被取代的段落；只前置新说明会留下两条相反论断。
4. 翻转过的结论**带翻转史**，引用前复核。

**实例（2026-10-05 → 10-06）**：`[ERR-20261005-009]` 判"DISP<0 是压实" ⇒
10-06 `[ERR-20261006-003]` 用结构+几何+物理三证改回"**DISP>0 是压实**"（即最初约定）。
"""
g = io.open(GOT, encoding="utf-8").read()
if "方向/符号结论的铁律" in g:
    print("· gotchas 已有")
else:
    io.open(GOT, "a", encoding="utf-8").write(G)
    print(f"operational-gotchas.md +{len(G)} 字符")

# ───────── 3) MEMORY 加一行指针（需腾空间） ─────────
m = io.open(MEM, encoding="utf-8").read()
n0 = len(m.encode("utf-8"))
TRIM = [
    ('- **`taskkill /IM` 会杀 agent 自己** —— 用 PID | `references/operational-gotchas.md` |',
     '- **`taskkill /IM` 会杀自己** —— 用 PID | `references/operational-gotchas.md` |'),
    ('| 20 | **Git 历史敏感值清除** → `digest/procedure/git-history-secret-purge.md` |',
     '| 20 | **Git 历史敏感值清除** → `digest/procedure/git-history-secret-purge.md` |\n'
     '| 23 | 🧭 **方向/符号类结论必须先有「不可被约定翻转」的证据**（结构/几何/物理/解剖）'
     '＋**判据要有已知正对照**；**更正要替换不要叠加** | '
     '`references/operational-gotchas.md` ＋ `[LRN-20261006-083]` |'),
]
for old, new in TRIM:
    if old in m:
        m = m.replace(old, new, 1)
    else:
        print(f"  ⚠ MEMORY 未命中：{old[:46]}")
io.open(MEM, "w", encoding="utf-8").write(m)
n1 = len(m.encode("utf-8"))
print(f"MEMORY.md {n0} -> {n1} B  {'✅' if n1 <= 10240 else f'❌ 超 {n1-10240}'}")
