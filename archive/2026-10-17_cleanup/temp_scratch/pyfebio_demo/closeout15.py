"""把 LRN-055 的更正传播到各文档 + 全量校验。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"

# ── ① ANKLE_FIX_TODO banner 更正 ──
AF = ROOT / "ANKLE_FIX_TODO.md"
s = AF.read_text(encoding="utf-8")
OLD_BANNER_START = "> ## 🚨 2026-10-05 重大更正"
i = s.find(OLD_BANNER_START)
if i >= 0:
    j = s.find("\n\n", s.find("修关节几何", i))
    if j < 0:
        j = s.find("\n\n#", i)
    NB = """> ## 🚨 2026-10-05 更正（`[LRN-20261005-054]` + `[LRN-20261005-055]`）
> **`F` 未收敛，但不是伪影；关节几何没问题。** 实测 `k = F/δ` 随罚因子单调趋近物理值：
> **49 / 151 / 451 / 880** N/mm（PENALTY 1e-3→1），对照 Anderson §2.4 的独立估计
> `E·A/t ≈ 1059 N/mm`（p=1 时仅差 −17%）。**但 `PENALTY ≥ 2` 全部负 Jacobian ⇒ 无法验收敛**，
> 趋势外推真值可能**高 1.5~2×** ⇒ `436.8 N` 是**偏低值**。
> 另：`|u|max` 超指定 10× 是**400mm 胫骨悬臂弯曲**（刚体拟合 0.17° 转动占 79%），非故障；
> 关节几何**正常**（点在-tet 互穿仅 6 节点；加载后接触 **344 mm²** ≈ Anderson 生理 **300 mm²**）。
> **唯一真堵点 = 接触求解器的罚因子上限。**
"""
    AF.write_text(s[:i] + NB + s[j:], encoding="utf-8")
print(f"ANKLE_FIX_TODO.md -> {AF.stat().st_size:,} B  "
      f"[{'OK' if 'LRN-20261005-055' in AF.read_text(encoding='utf-8') else 'XX'}]")

# ── ② docs/总体进度.md 行更正 ──
OP = ROOT / "docs/总体进度.md"
s = OP.read_text(encoding="utf-8")
OLD = ("| **🚨 关节几何不共形** | 🔴 | **F 全是罚因子伪影**（PENALTY 1e-4→1 ⇒ 0.83→419.56 N；"
       "AUGLAG 也不收敛 2.00→315.06 N）⇒ 关节从未闭合（中位间隙 1.6 mm、1~3 点接触）"
       "⇒ `436.8N`/「3.4≈吴恺3.74」**失效**。`[LRN-054]` |")
NEW = ("| **⚠️ F 未收敛（非伪影）** | 🟠 | `k=F/δ` 随罚趋近物理值 **49/151/451/880** N/mm"
       "（罚 1e-3→1）vs Anderson 独立估计 **1059**；但 **`PENALTY≥2` 全负 Jacobian ⇒ 无法验收敛**，"
       "真值可能高 **1.5~2×** ⇒ `436.8N` **偏低**。几何**正常**（互穿 6 节点、接触 344mm²≈生理 300mm²）。"
       "`[LRN-055]` |")
if OLD in s:
    OP.write_text(s.replace(OLD, NEW, 1), encoding="utf-8")
print(f"docs/总体进度.md -> {OP.stat().st_size:,} B  "
      f"[{'OK' if 'LRN-055' in OP.read_text(encoding='utf-8') else 'XX'}]")

# ── ③ CART_HEX_TODO.md 更正段 ──
CH = ROOT / "CART_HEX_TODO.md"
ch = CH.read_text(encoding="utf-8")
i = ch.find("## 🚨 状态修正")
if i >= 0:
    j = ch.find("## ✅ 已完成", i)
    NB = """## ⚠️ 状态修正（2026-10-05，`[LRN-20261005-055]` 更正）

**先前的「关节两软骨面不共形 ⇒ 关节从未闭合 ⇒ F 全是罚因子伪影」判断有误，已更正：**

- ❌ 「中位间隙 1.6 mm、仅 1~3 点接触」是**节点稀疏导致的漏检**。稠密面积加权测量：
  t=0 接触 **31 mm²**（0.02mm 容差内），压 0.2mm → **231 mm²**；加载后 **344 mm²**
  ≈ **Anderson 生理值 300 mm²** ✓
- ❌ 「F 全是伪影」也错：`k = F/δ` 随罚因子**单调趋近物理值**（49/151/451/**880** N/mm，
  罚 1e-3→1），对照 Anderson 独立估计 `E·A/t ≈ 1059 N/mm`（p=1 差 −17%），
  且 F 的增量倍率在**递减**（8.4→7.2→4.3→1.9）= **在收敛**
- ✅ **几何/`cartilage.py` 是对的**：点在-四面体内测试互穿仅 **6 节点**，与
  `cartilage.json` 的 `n_penetrating = 0` 一致

**⇒ 真堵点：`PENALTY ≥ 2` 全部负 Jacobian ⇒ 罚因子推不上去 ⇒ `F` 无法验收敛。**
p=1 时 k 仍在以 ~1.9×/十倍 爬升，真收敛值**可能比 419.56 N 高 1.5~2×** ⇒ `436.8 N` **偏低**。

"""
    CH.write_text(ch[:i] + NB + ch[j:], encoding="utf-8")
print(f"CART_HEX_TODO.md -> {CH.stat().st_size:,} B  "
      f"[{'OK' if 'LRN-20261005-055' in CH.read_text(encoding='utf-8') else 'XX'}]")

# ── ④ TIE_DIVERGE_TODO.md 结论更正 ──
TD = ROOT / "TIE_DIVERGE_TODO.md"
td = TD.read_text(encoding="utf-8")
i = td.find("**结论：「tie 发散」是误诊。**")
j = td.find("### 关键实验", i)
if i >= 0 and j > i:
    NEWC = """**结论：「tie 发散」是误诊。** tie 没有问题；`|u|max` 超指定 10× 是
**400mm 胫骨悬臂弯曲**（刚体拟合 0.17° 转动占位移 79%）。

**但"真根因"要先看 `[LRN-20261005-055]` 的更正**：关节几何**没问题**、
`F` **不是**伪影（`k=F/δ` 随罚趋近物理值 880 vs 1059 N/mm）。
真堵点是 **`PENALTY ≥ 2` 全部负 Jacobian ⇒ 罚因子推不上去 ⇒ `F` 无法验收敛**。

"""
    TD.write_text(td[:i] + NEWC + td[j:], encoding="utf-8")
print(f"TIE_DIVERGE_TODO.md -> {TD.stat().st_size:,} B  "
      f"[{'OK' if 'LRN-055' in TD.read_text(encoding='utf-8') else 'XX'}]")

# ── ⑤ session note 追加 ──
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
SEC = """

---

## 十五、自我更正：关节几何没问题、F 不是伪影（2026-10-05）→ `[LRN-20261005-055]`

`[LRN-054]` 有两处**过度结论**，本轮用更好的测量自己推翻。

### 15.1 ❌ 撤回「关节不共形 ⇒ 从未闭合」
**错因**：用**节点到节点**距离（69 个 tibia 节点）⇒ 接触斑落在节点之间，**完全漏检**。
**正解**（`gap_area.py`，talus 外面稠密 20,000 点、**无符号**面积加权）：

| gap | 面积 | 占比 |
|---|---|---|
| < 0.02 mm（FEBio 容差） | **31.0 mm²** | 1.5% |
| < 0.20 mm | 230.6 mm² | 11.4% |
| 中位 | 0.9732 mm | |

**独立裁决**：`point_in_tet.py`（点在-四面体内）⇒ **互穿仅 6 节点** ✓ 与
`cartilage.json` 的 `n_penetrating = 0` 一致 ⇒ **`cartilage.py` 的防穿透是对的**。
（中间我还错过一次：`contact_area.py` 用**有符号**距离算出「312mm² 干涉、深 2.9mm」
并据此提出「软骨厚度冲突」新根因 —— **那是法向定向 bug**，已作废。）

### 15.2 ❌ 撤回「F 100% 是罚因子伪影」
**错因**：只看 F 从 0.83 涨到 419.56 就断定无收敛值，**没算物理量**。
**正解**：算关节刚度 `k = F/δ`（δ = 关节实际压缩，`k_joint.py`）：

| PENALTY | F(N) | δmax(mm) | **k=F/δ (N/mm)** | 接触面积 |
|---|---|---|---|---|
| 1e-3 | 6.99 | 0.14229 | **49.1** | 544 |
| 1e-2 | 50.65 | 0.33557 | **150.9** | 476 |
| 1e-1 | 220.08 | 0.48767 | **451.3** | 374 |
| **1** | **419.56** | 0.47680 | **880.0** | **344** |

⇒ k **单调趋近** Anderson §2.4 的独立估计 `E·A/t ≈ **1059** N/mm`（p=1 差 −17%）；
F 的增量倍率递减（8.4→7.2→4.3→**1.9**）= **收敛**。
**加载后接触面积 344 mm² ≈ Anderson 生理 300 mm²** ✓ ← 真正的验证。
（AUGLAG 的 2.00/15.31/110.63/315.06 不趋同，是 **AUGLAG 增广自身没收敛**，属独立问题。）

### 15.3 ✅ 更正后仍成立的真堵点
`pen_converge.py`：`PENALTY = 1/2/4/8/16` → rc = **0/1/1/1/1**（p≥2 全负 Jacobian）
⇒ **罚因子推不上去 ⇒ F 无法验收敛**；p=1 时 k 仍以 ~1.9×/十倍 爬升
⇒ **真收敛值可能高 1.5~2×** ⇒ `436.8 N` 是**偏低值**（方向与 LRN-054 相反）。

### 15.4 更正后完整图景
| 项 | 结论 |
|---|---|
| `|u|max` 10× | ✅ 400mm 胫骨悬臂弯曲（0.17° 占 79%），非故障 |
| tie | ✅ 正常（滑移 = 自身运动 3%） |
| 关节几何 | ✅ 正常（互穿 6 节点） |
| 接触面积 | ✅ 31mm²@t0 → **344mm²**@加载 ≈ Anderson **300mm²** |
| `F` | ⚠️ **未收敛**，p=1 偏低，真值可能高 1.5~2× |
| m2+PENALTY rc=1 | 🔴 真故障（关节罚因子引起，与 tie 无关） |
| **首要待解决** | 🔴 让接触在**更高罚因子/细网格**下收敛 |

### 15.5 ➡️ 下一步（更正后）
1. **提高接触求解器的罚因子上限**（唯一真堵点）：`symmetric_stiffness`/`node_reloc`/
   更小步长/`max_refs`，或把 AUGLAG 的增广调到真收敛
2. 收敛后验证 `k → 1059`，这才算通过
3. 之后再谈 m2 vs baseline、tet4 vs penta6 的 2.2×（很可能是同一族问题）
4. 路线3（骨刚性/截断/Anderson 边界）仍值得做，但**优先级降到第 2**（几何没问题了，它不是救火）
"""
if SESS.exists():
    t = SESS.read_text(encoding="utf-8")
    if "十五、自我更正" not in t:
        SESS.write_text(t.rstrip() + SEC, encoding="utf-8")
    print(f"session note -> {SESS.stat().st_size:,} B")

# ── ⑥ 校验 ──
print("\n" + "=" * 76)
print("全量校验")
print("=" * 76)
lt = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
ms = M.read_text(encoding="utf-8")
CHK = [
    ("LRN-055 已入库", "LRN-20261005-055" in lt),
    ("LRN-055 撤回不共形", "节点稀疏导致的漏检" in lt),
    ("LRN-055 撤回伪影", "没看增量倍率" in lt),
    ("LRN-055 k 表", "880.0" in lt and "1059" in lt),
    ("LRN-055 罚上限", "p≥2 全部 `Negative jacobian`" in lt or "p≥2 全部" in lt),
    ("LRN-055 记中间错误", "法向定向不一致导致的伪影" in lt),
    ("MEMORY 更正", "F 未收敛但非伪影" in ms),
    ("MEMORY 几何没问题", "几何没问题" in ms),
    ("ANKLE_FIX 更正", "LRN-20261005-055" in AF.read_text(encoding="utf-8")),
    ("总体进度 更正", "LRN-055" in OP.read_text(encoding="utf-8")),
    ("CART_HEX 更正", "LRN-20261005-055" in CH.read_text(encoding="utf-8")),
    ("TIE_DIVERGE 更正", "LRN-055" in TD.read_text(encoding="utf-8")),
]
for k, v in CHK:
    print(f"  [{'OK' if v else 'XX'}] {k}")
sz = M.stat().st_size
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb))`", ms)))
print(f"\n  MEMORY.md {sz:,}/10240 余量 {10240-sz:+,} {'OK' if sz<=10240 else '超'}")
for f in ("closest_vec.py", "gap_area.py", "point_in_tet.py", "k_joint.py",
          "pen_converge.py", "align_diag.py", "bone_gap.py", "check_npz.py"):
    p = ROOT / "temp/pyfebio_demo" / f
    print(f"  [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>7,}  {f}")
