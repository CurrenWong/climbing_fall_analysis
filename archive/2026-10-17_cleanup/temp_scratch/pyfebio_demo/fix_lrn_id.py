"""修正：① 用新 ID LRN-056 写入更正 ② 修 TIE_DIVERGE_TODO ③ 改交叉引用。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
LT = ROOT / ".learnings/LEARNINGS.md"

# ── ① 先确认 ID 未被占用 ──
t = LT.read_text(encoding="utf-8")
have = set(re.findall(r"\[LRN-(\d{8}-\d{3})\]", t))
NEWID = None
for k in range(56, 70):
    if f"20261005-{k:03d}" not in have:
        NEWID = f"20261005-{k:03d}"
        break
print(f"现有 ID {len(have)} 个；选用新 ID = LRN-{NEWID}")
assert NEWID, "找不到空 ID"

CORR = r"""

---

## [LRN-NEWID] ⚠️ **更正 `[LRN-20261005-054]`**：关节几何没问题、F 不是伪影 —— 真堵点是**罚因子上限**

- **日期**：2026-10-05
- **类别**：`correction`（自我更正；`[LRN-054]` 有两处**过度结论**）
- **为什么记**：先记下这条，防止未来再犯「**测量伪影 → 惊天结论**」的错误。
- ⚠️ **本条自身也是教训**：第一版用了 `[LRN-20261005-055]`，
  但那个 ID **早已被攀岩跌落线占用**（"骨折判据按实际部位取强度"）⇒ 写入守卫
  `if ID not in 文本` 提前满足 ⇒ **整段被静默跳过，没写进去**。
  ⇒ **分配 LRN ID 前必须先扫一遍现有 ID**（本文件还另有 6 个重复 ID 待清理）。

### ❌ 撤回 `[LRN-054]` 第 ① 处：「关节两软骨面不共形 ⇒ 关节从未闭合」

**错因**：用**节点到节点**距离（tibia 69 个节点 → talus 200 面）得到「中位 1.62mm、
只有 1~3 点接触」。**节点稀疏 ⇒ 接触斑落在节点之间，被完全漏检。**

**正确测量**（`gap_area.py`：talus 外面稠密 100 子三角/面 = 20,000 采样点，
**无符号**距离，按面积加权）：

| gap | 面积 | 占比 |
|---|---|---|
| < 0.02 mm（FEBio `tolerance`） | **31.0 mm²** | **1.5%** |
| < 0.10 mm | 134.3 mm² | 6.7% |
| < 0.20 mm | 230.6 mm² | 11.4% |
| < 0.50 mm | 509.2 mm² | 25.2% |
| 中位 | **0.9732 mm** | |

⇒ **t=0 只是「相切于一点」（31mm²），压 0.2mm 即长到 231mm²** ——
这是**未预载关节的正常表现**（Anderson 也是先「轻接触 → **200 N 落位**」才进正式加载）。

**独立裁决（无假设）**：`point_in_tet.py` 直接做**点在四面体内**测试
⇒ **互穿节点仅 6 个**（4 talus 外侧节点在 tibia 软骨内 + 2 反向）
⇒ 与 `cartilage.json` 的 `n_penetrating = 0` **一致** ✓
⇒ **`cartilage.py` 的防穿透逻辑（`t = min(1.5, FILL_FRAC×g)`，`FILL_FRAC=0.5`）是正确的。**

**⚠️ 我自己的中间错误（一并记下）**：`contact_area.py` 用**有符号**距离算出
「312 mm² 干涉、最深 2.9 mm」，并据此提出「软骨厚度 1.5mm 与骨最小间隙 1.406mm 冲突」的
新根因 ⇒ **那是法向定向不一致造成的伪影**，换成无符号距离 + 点在-tet 双重验证后消失。

### ❌ 撤回 `[LRN-054]` 第 ② 处：「F 100% 是罚因子伪影，AUGLAG 也不收敛」

**错因**：只看 F 从 0.83 涨到 419.56 就断定「无收敛值」，**既没看增量倍率，也没算物理量**。

**正确做法：算关节刚度 `k = F / δ`**（δ = 关节实际压缩量，`k_joint.py`）：

| PENALTY | F (N) | δmax (mm) | **k = F/δ (N/mm)** | 加载后接触面积 |
|---|---|---|---|---|
| 1e-3 | 6.99 | 0.14229 | **49.1** | 544 mm² |
| 1e-2 | 50.65 | 0.33557 | **150.9** | 476 mm² |
| 1e-1 | 220.08 | 0.48767 | **451.3** | 374 mm² |
| **1** | **419.56** | 0.47680 | **880.0** | **344 mm²** |

- **k 单调趋近物理值**：Anderson §2.4 的**独立**估计
  `k ≈ E·A/t = 12×300/3.4 ≈ **1059 N/mm**` ⇒ p=1 时 880，**仅差 −17%**
- **F 增量倍率递减**：8.4 → 7.2 → 4.3 → **1.9** ⇒ **在收敛，不是发散**
- **加载后接触面积 344 mm² ≈ Anderson 生理值 300 mm²** ✓ ← **一条真正的验证**

⇒ `F` 不是伪影，是「罚函数法在逼近硬接触」的正常数列。
（`AUGLAG` 给 2.00/15.31/110.63/315.06 而不趋同 = **AUGLAG 自身增广没收敛**，属独立问题，
不能用来否定 PENALTY 法的收敛趋势。）

### ✅ 更正后仍然成立的真堵点：**罚因子推不上去 ⇒ F 无法验收敛**

`pen_converge.py`（baseline，其余默认）：

| PENALTY | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| rc | **0** | 1 | 1 | 1 | 1 |
| 正常终止 | Y | n | n | n | n |

⇒ **p≥2 全部 `Negative jacobian`** ⇒ 只能停在 p=1。而 p=1 时 k 仍以 **~1.9×/十倍** 爬升
⇒ **真收敛值很可能比 419.56 N 高 1.5~2×** ⇒ **`436.8 N` 是偏低值**（方向与 `[LRN-054]` 相反：
不是"失效"，是"未收敛/偏低"）。

### 📌 更正后的完整图景

| 项 | 结论 |
|---|---|
| `\|u\|max` 超指定 10× | ✅ 是 **400mm 胫骨悬臂弯曲**（刚体拟合 0.17° 转动占 79%），非故障 |
| tie 界面 | ✅ **正常**（界面滑移 = 自身运动的 **3%**） |
| 关节几何 | ✅ **正常**（互穿 6 节点；`cartilage.py` 防穿透正确） |
| 接触面积 | ✅ t=0 31mm² → 加载后 **344mm²**，对齐 Anderson **300mm²** |
| `F` 数值 | ⚠️ **未收敛**，p=1 偏低，真值可能高 **1.5~2×** |
| `m2 + PENALTY` rc=1 | 🔴 **真故障**：由关节罚因子引起，与 tie 无关 |
| **首要待解决** | 🔴 **让接触在更高罚因子 / 细网格下收敛** |

### ➡️ 更正后的下一步（替代 `[LRN-054]` 的"下一步"）

1. **提高接触求解器的罚因子上限**（唯一真堵点）：
   `symmetric_stiffness` / `node_reloc` / 更小步长 / `max_refs`，
   或把 AUGLAG 的增广调到真收敛（`maxaug`/`augtol`）再验证 `k → 1059`
2. **收敛后重跑**：`k` 是否落到 1059 附近 ⇒ 这才是验收判据
3. 之后再看 m2 vs baseline、`tet4 vs penta6` 的 2.2×（很可能是**同一族**问题）
4. **路线3（骨刚性/截断/Anderson 边界）仍值得做，但优先级降到第 2** ——
   几何既然没问题，它主要收益是「去掉长骨弯曲」+「对齐 Anderson 设定」，不是救火

### 📦 产物（全部带自检）
`temp/pyfebio_demo/`：`closest_vec.py`（**向量化 pt-tri，自检 ALL PASS**）、
`gap_area.py`、`point_in_tet.py`、`k_joint.py`、`pen_converge.py`、`align_diag.py`、
`bone_gap.py`、`check_npz.py`、`contact_area.py`（**已作废：法向定向有 bug**）
""".replace("LRN-NEWID", f"LRN-{NEWID}")

t2 = t.rstrip() + CORR
LT.write_text(t2, encoding="utf-8")
print(f"LEARNINGS.md -> {LT.stat().st_size:,} B  "
      f"[{'OK' if f'LRN-{NEWID}' in LT.read_text(encoding='utf-8') else 'XX'}]")

# ── ② TIE_DIVERGE_TODO.md 重写 ──
(ROOT / "TIE_DIVERGE_TODO.md").write_text(f"""# TIE_DIVERGE — 查 tie 界面发散的真因

## ✅ 已结案（2026-10-05）→ `[LRN-20261005-054]` + `[LRN-{NEWID}]`（更正）

**结论：「tie 发散」是误诊。** tie 没有问题。

### 三条实测结论
1. **`|u|max` 超指定 10× 不是 tie 滑移** —— 是 **400mm 胫骨悬臂弯曲**
   （最小二乘刚体拟合：转动 **0.17°** 解释 **79%** 位移，残差/原始 = 20.9%；
   `F=3EIδ/L³` 算约 20N 侧向力即可弯 1.2mm）。
   ⇒ **`|u|max ≈ 指定位移` 是无效判据，已作废。**
2. **tie 本身正常** —— 界面相对位移中位 **0.0476mm = 软骨自身运动(1.5857mm) 的 3.0%**。
3. **`m2 + PENALTY` rc=1 的真因 = 关节接触罚因子，与 tie 无关**
   （关 tie `auto_penalty` ❌、tie `penalty=1000` ❌ 都无效；
   关节 `PENALTY=0.001` → rc=0 且自动步长主动放大）。

### ⚠️ 但「真根因」要按 `[LRN-{NEWID}]` 的更正理解
- 关节几何**没问题**：点在-四面体内测试互穿仅 **6 节点**；
  t=0 接触 **31mm²** → 加载后 **344mm²** ≈ Anderson 生理 **300mm²**
- `F` **不是伪影**：`k = F/δ` 随罚因子趋近物理值 **49/151/451/880** N/mm（罚 1e-3→1），
  对照 Anderson 独立估计 `E·A/t ≈ 1059 N/mm`（p=1 差 **−17%**）
- **真堵点 = 接触求解器的罚因子上限**：`PENALTY ≥ 2` 全部负 Jacobian
  ⇒ 罚推不上去 ⇒ **`F` 无法验收敛**，真值可能高 **1.5~2×** ⇒ `436.8N` **偏低**

### 关键实验（全部实跑）
| 实验 | 结果 |
|---|---|
| 4 因子矩阵（网格×TIE×LAUGON） | `LAUGON` 单独改变 F 25%；m2+PENALTY rc=1 |
| D 失败归因 | `Negative jacobian` @第一步；关 tie `auto_penalty`/加硬 tie penalty **无效** |
| 关节 `PENALTY=0.001` | ✅ rc=0 且自动步长主动放大 |
| 关节罚扫描 | F：0.83→6.99→50.65→220.08→**419.56 N**（1e-4→1） |
| **关节刚度** | **k = 49.1→150.9→451.3→880.0 N/mm**（罚 1e-3→1）vs 物理 **1059** |
| 罚上限 | `PENALTY=1/2/4/8/16` → rc **0/1/1/1/1** ⇒ **推不上去** |
| AUGLAG 初始罚扫描 | 2.00→15.31→110.63→315.06 N ⇒ **AUGLAG 增广自身没收敛**（独立问题） |
| 位移扫描（P=0.001） | 2.0 mm 后饱和在 ~80 N（低罚下软） |
| 载荷截面 H=60/30 mm | 移到关节上方也救不了（说明不是载荷引入问题） |

### 下一步
1. **提高接触求解器罚因子上限**（唯一真堵点）
2. 收敛后验证 `k → 1059` 才算通过
3. 之后再重做 m2 vs baseline、tet4 vs penta6 的 2.2× 对比
""", encoding="utf-8")
print(f"TIE_DIVERGE_TODO.md -> {(ROOT/'TIE_DIVERGE_TODO.md').stat().st_size:,} B")

# ── ③ 交叉引用 055 → NEWID ──
FIX = [
    (WS / "MEMORY.md", ["`[LRN-054/055]`", "LRN-054/055"]),
    (ROOT / "ANKLE_FIX_TODO.md", ["LRN-20261005-055", "LRN-055"]),
    (ROOT / "docs/总体进度.md", ["LRN-055", "LRN-20261005-055"]),
    (ROOT / "CART_HEX_TODO.md", ["LRN-20261005-055", "LRN-055"]),
]
for p, olds in FIX:
    if not p.exists():
        continue
    s = p.read_text(encoding="utf-8")
    n = 0
    for o in olds:
        if o in s:
            s = s.replace(o, o.replace("055", NEWID[-3:]))
            n += 1
    p.write_text(s, encoding="utf-8")
    print(f"  {p.name:24s} 改了 {n} 处  -> {p.stat().st_size:,} B")

M = WS / "MEMORY.md"
sz = M.stat().st_size
print(f"\nMEMORY.md {sz:,}/10240 余量 {10240-sz:+,} {'✓' if sz<=10240 else '✗'}")
