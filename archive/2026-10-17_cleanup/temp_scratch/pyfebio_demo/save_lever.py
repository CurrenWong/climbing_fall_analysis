"""落盘：有效力臂验证 + 失效角预测（C3 绝对化的收口）。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
LEARN = os.path.join(R, ".learnings", "LEARNINGS.md")
REF = os.path.join(W, "references", "ankle-fe-criteria.md")
USE = os.path.join(R, "docs", "踝FE可用模型使用说明.md")
TOD = os.path.join(R, "ANKLE_FIX_TODO.md")
DAY = os.path.join(W, "memory", "2026-10-06.md")

L = """
## LRN-20261006-082 — 姿势预应变是**运动学必然**（力臂解剖吻合）⇒ **C3 可绝对化** + 失效角预测

**类别**: method / 定量收口  **关联**: `[LRN-20261006-081]`、`[ERR-20261006-003]`

### 待决问题的答案
`[LRN-081]` 遗留：+5° 时 CFL 预应变 **+8.98%** 是物理还是建模产物？
**判定法**：刚体转动下 `ΔL = θ · r_eff` ⇒ 由**实测 ΔL/θ** 反推力臂，与解剖比对。

### 实测有效力臂（`check_c2.py`，`NO_RUN=1`，纯几何，±1°~15° 全一致）

| 韧带 | 1° | 3° | 5° | 8° | 10° | **收敛值** | 解剖 | 判定 |
|---|---|---|---|---|---|---|---|---|
| **CFL** | 20.6 | 27.9 | 28.5 | 29.1 | 29.4 | **≈28.5 mm** | ~25–30 mm | ✅ 吻合 |
| **ATFL** | 17.8 | 20.8 | 21.0 | 20.9 | 21.3 | **≈21 mm** | ~15–22 mm | ✅ 吻合 |
| **PTFL** | — | 7.8 | 8.0 | 8.1 | 8.2 | **≈8 mm** | 后方近轴 | ✅ 合理 |
| MD1/MD2 | — | ~12/~13 | ~11/~13 | — | ~12/~13 | **≈−12 mm**（反号，内侧）| — | ✅ |

- **⇒ 力臂与解剖量级一致** ⇒ 预应变是「足绕踝轴刚性旋转」的**正确运动学后果**，
  **不是建模伪影** ⇒ **`[LRN-081]` 的待决问题结案：C3 可以绝对化**（对长度正确的 CFL 而言）
- 1° 处的偏离（20.6 vs 28.5）是**落位偏移**在小角度时占比变大所致（大角度时被淹没）

### ⭐ 失效角预测（零载，`ε(θ) = θ·r_eff/L`）

| 韧带 | 用的 L | 失效带 | **推得失效角** |
|---|---|---|---|
| **CFL** | **27.9 mm（正确 ✓）** | 9% / 13% | **5.05° / 7.29°** |
| ATFL（现状）| 34.2 mm（**≈2× 过长**）| 11% / 15% | 10.3° / 14.0° |
| **ATFL（若修到 20 mm）** | 20 mm | 11% / 15% | **6.0° / 8.2°** |

**⇒ ⇒ 修完 ATFL 长度后，模型预测 **ATFL 先失效、CFL 随后**：
5° 时 ATFL 9.15%（1.83/20）> CFL 8.98%（2.49/27.9）—— 与**临床内翻扭伤的顺序一致** ✓
**⇒ 这是一个真正的验证靶子：修 ATFL 长度 ⇒ 看能否复现「ATFL 先断」的临床排序。**

### ⚠️ 载荷会让失效角**变大**（压缩关上关节 ⇒ 卸载韧带）
实测：+5° 时 CFL 预应变 **+8.98%**（零载）→ 固定压缩 +1.0mm **+6.31%** → **500 N `+5.03%`**
⇒ **⇒ 内翻扭伤是「角度主导、轴压卸载」的机制** ⇒ 判据应是 (载荷, 角度) 平面里的**边界曲线**，
不能说"某个角度就断"。在低载时边界角最小（CFL ≈5°）。

### 局限（必须声明）
- 本节的失效角是**零载**下的；带载要查表（`[LRN-081]` 的恒定载荷表）
- ATFL 的绝对值仍受 **2× 长度误差**污染 ⇒ 其失效角（10.3°）**偏高**，修长度后会降到 ~6°
- 姿势上限 ±5°（模型约束）⇒ 失效角是**外推**（在 1~10° 的线性区内插到 5~7° ⇒ 内插，非远外推 ✓）
"""
t = io.open(LEARN, encoding="utf-8").read()
if "LRN-20261006-082" in t:
    print("· LRN 已存在")
else:
    io.open(LEARN, "a", encoding="utf-8").write(L)
    print(f"LEARNINGS.md {len(t)} -> {len(t)+len(L)}")

# ── references：C3 段追加"失效角"小节 ──
ADD = """
## ⭐ 失效角预测（`[LRN-20261006-082]`）—— C3 绝对化的收口

**有效力臂**（实测 `ΔL/θ`，与解剖吻合 ⇒ 预应变是**运动学必然**、非伪影）：
**CFL ≈28.5 mm｜ATFL ≈21 mm｜PTFL ≈8 mm｜内侧 MD ≈−12 mm**

**零载失效角** `ε(θ) = θ·r_eff/L`：

| 韧带 | 用的 L | 失效带 | **失效角** |
|---|---|---|---|
| **CFL** | 27.9 mm（正确 ✓）| 9% / 13% | **5.05° / 7.29°** |
| ATFL（现状）| 34.2 mm（**2× 过长**）| 11% / 15% | 10.3° / 14.0° |
| **ATFL（修到 20 mm 后）** | 20 mm | 11% / 15% | **6.0° / 8.2°** |

**⇒ 修完 ATFL 长度后应能复现「ATFL 先失效」的临床排序**（5° 时 9.15% > 8.98%）
—— **这是模型的验证靶子。**

**⚠️ 载荷会让失效角变大**：+5° 时 CFL 零载 +8.98% → 固定压缩 +1.0mm **+6.31%** → 500 N **+5.03%**
⇒ **判据应是 (载荷, 角度) 平面的边界曲线**；低载时边界角最小（CFL ≈5°）。

**⚠️ 姿势上限 ±5°（收敛约束）⇒ 失效角在 1~10° 线性区内插，不是远外推。**
"""
r = io.open(REF, encoding="utf-8").read()
io.open(REF, "a", encoding="utf-8").write(ADD)
print(f"ankle-fe-criteria.md +{len(ADD)} 字符")

# ── 使用说明 §5 下一步更新 ──
u = io.open(USE, encoding="utf-8").read()
OLD = "3. **修韧带插入点几何**（ATFL 优先，36.99 → 15~20 mm 量级）"
NEW = ("**✅ C3 已收口**：力臂实测与解剖吻合（CFL 28.5 / ATFL 21 / PTFL 8 mm）⇒ 预应变是运动学必然\n"
       "   ⇒ **C3 可绝对化**；零载失效角 **CFL 5.05°~7.29°**、**ATFL（修长度后）6.0°~8.2°**\n"
       "   ⇒ 应复现「ATFL 先失效」的临床排序（验证靶子）。详见 `references/ankle-fe-criteria.md`\n\n"
       "3. **修 ATFL 长度**（34.2 → 20 mm 量级）—— 修完才能验证上面的排序预测")
if OLD in u:
    u = u.replace(OLD, NEW)
    io.open(USE, "w", encoding="utf-8").write(u)
    print(f"使用说明下一步已更新 -> {len(u.encode('utf-8'))} B")
else:
    print("⚠ 使用说明下一步未命中")

TADD = """

### ✅ 2026-10-06 第五轮：C3 收口（力臂验证 + 失效角）

- **力臂实测**（`ΔL/θ`，纯几何）：**CFL ≈28.5 mm｜ATFL ≈21 mm｜PTFL ≈8 mm** ⇒ **与解剖吻合**
  ⇒ **预应变是「足绕踝轴刚性旋转」的运动学必然，不是建模伪影** ⇒ **C3 可绝对化**
- **零载失效角**（`ε=θ·r/L`）：**CFL 5.05°(保守)/7.29°(中值)**；
  ATFL 现状 10.3°/14.0°，**修到 20mm 后 → 6.0°/8.2°**
- **⇒ 验证靶子**：修完 ATFL 长度后应复现「**ATFL 先失效**」的临床排序（5° 时 9.15% > 8.98%）
- **⚠️ 载荷使失效角变大**（压缩关上关节）：+5° 时 CFL 零载 8.98% → 压缩 6.31% → 500N 5.03%
  ⇒ **判据是 (载荷, 角度) 平面的边界曲线**，低载时边界角最小
"""
t2 = io.open(TOD, encoding="utf-8").read()
if "第五轮：C3 收口" in t2:
    print("· TODO 已有")
else:
    io.open(TOD, "a", encoding="utf-8").write(TADD)
    print(f"ANKLE_FIX_TODO +{len(TADD)} 字符")

DADD = """
## ✅ 第五轮：C3 收口 —— 力臂验证 + 失效角

**待决问题结案**：+5° 时 CFL 预应变 +8.98% 是**物理的**。

**判定**：刚体转动 `ΔL = θ·r_eff` ⇒ 由实测反推力臂：

| 韧带 | 实测 r_eff | 解剖 | 判定 |
|---|---|---|---|
| CFL | **28.5 mm** | ~25–30 | ✅ |
| ATFL | **21 mm** | ~15–22 | ✅ |
| PTFL | **8 mm** | 后方近轴 | ✅ |

**⇒ 力臂解剖吻合 ⇒ 预应变是运动学必然、非伪影 ⇒ C3 可绝对化。**

**零载失效角**（`ε=θ·r/L`）：**CFL 5.05°/7.29°**；ATFL 现状 10.3°/14.0° ⇒ **修到 20mm 后 6.0°/8.2°**。

**⇒ 验证靶子**：修完 ATFL 长度应复现「ATFL 先失效」（5° 时 9.15% > 8.98%）—— 与临床一致。

**⚠️ 载荷让失效角变大**：+5° 时 CFL 零载 8.98% → 压缩 6.31% → 500 N 5.03%
⇒ **判据是 (载荷, 角度) 平面的边界曲线**；低载时边界角最小（CFL ≈5°）。

**落盘**：`[LRN-20261006-082]`｜`references/ankle-fe-criteria.md` +失效角小节｜使用说明 §5 更新。

**下一步**
1. **修 ATFL 长度**（34.2 → 20 mm）⇒ 验证「ATFL 先失效」预测
2. 查接触刚度为何仅 ~150 N/mm
3. 接 OpenSim 姿势输入
"""
d = io.open(DAY, encoding="utf-8").read()
if "第五轮：C3 收口" in d:
    print("· 日记已有")
else:
    io.open(DAY, "a", encoding="utf-8").write(DADD)
    print(f"日记 +{len(DADD)} 字符")

# ── 校验 ──
print("\n" + "=" * 68)
import ast
import re
ok = True
for p, cap in ((os.path.join(W, "MEMORY.md"), 10240), (REF, None), (USE, None),
               (TOD, None), (LEARN, None),
               (os.path.join(R, ".learnings", "ERRORS.md"), None),
               (os.path.join(R, "scripts", "ankle_fe", "ligaments.py"), None)):
    n = len(io.open(p, encoding="utf-8").read().encode("utf-8"))
    f = ""
    if cap:
        f = " ✅" if n <= cap else f" ❌ 超 {n-cap}"
        ok &= n <= cap
    print(f"  {os.path.basename(p):<32} {n:>8} B{f}")
for s in ("ligaments.py", "ankle_lig_feb.py", "check_c2.py", "lig_strain.py"):
    ast.parse(io.open(os.path.join(R, "scripts", "ankle_fe", s),
                      encoding="utf-8").read())
print("  AST ✓ ×4")
for f in (".learnings/LEARNINGS.md", ".learnings/ERRORS.md"):
    x = io.open(os.path.join(R, f), encoding="utf-8").read()
    ids = re.findall(r"^#+\s*\[?([A-Z]{3}-\d{8}-\d{3})", x, re.M)
    dup = sorted({i for i in ids if ids.count(i) > 1})
    print(f"  {os.path.basename(f)}: {len(ids)} 条，重号 {dup or '无 ✓'}")
    ok &= not dup
r2 = io.open(REF, encoding="utf-8").read()
miss = [k for k in ("28.5 mm", "5.05°", "6.0° / 8.2°", "运动学必然") if k not in r2]
print(f"  references 缺失：{miss or '无 ✓'}")
ok &= not miss
u2 = io.open(USE, encoding="utf-8").read()
miss2 = [k for k in ("C3 已收口", "5.05°", "ATFL 先失效") if k not in u2]
print(f"  使用说明缺失：{miss2 or '无 ✓'}")
ok &= not miss2
print("\n结论：", "PASS ✅" if ok else "有问题 ⚠")
