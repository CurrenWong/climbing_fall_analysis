"""收口验证 + session 记忆笔记。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
DAILY = WS / "memory/2026-10-04.md"

SEC = """

---

## 十四、「tie 发散」结案：误诊，真根因是关节几何（2026-10-05）→ `[LRN-20261005-054]`

### 14.1 任务
`[LRN-053]` 说「tie 界面是细网格的新瓶颈」（`|u|max` 超指定 10×、F 不收敛、
`TIE_PENALTY` 加硬发散）。本轮核查这个结论。

### 14.2 ✅ 撤回三条
| 原结论 | 实测 |
|---|---|
| `|u|max` 10× ⇒ tie 在滑移 | **错**。胫骨位移刚体拟合 **0.17° 转动占 79%**（残差 20.9%）⇒ 是 **400mm 悬臂梁弯曲**；20 N 侧向力即可弯出 1.2mm。`axis` 400×3.0e-3 = 1.2mm ✓ |
| `TIE_PENALTY=1` 太软 | **错**。界面相对位移中位 **0.048mm** = 软骨自身运动(1.586mm)的 **3%** ⇒ tie 抓住 |
| 细网格 ⇒ 接触求解器问题 | 归因错。判据（失败点恒 `1/(11·TS)`）对，但真因是**关节罚因子** |

### 14.3 🔬 受控矩阵（4 因子）
| | 网格 | TIE | LAUGON | rc | F(N) |
|---|---|---|---|---|---|
| B | baseline | 0 | PENALTY | 0 | 419.56 |
| C | baseline | 0 | AUGLAG | 0 | 315.06 |
| D | m2 | 1 | PENALTY | **1** | — |
| E | m2 | 1 | AUGLAG | 0 | 190.80 |
⇒ `LAUGON` 单独改变 F **25%** ⇒ 之前 baseline↔m2 对比被污染。

**D 为何 rc=1（与 tie 无关）**：`Negative jacobian` @第一步。
关 tie `auto_penalty` ❌、tie `penalty=1000` ❌ **都无效**；
**关节 `PENALTY=0.001` → rc=0 且自动步长主动放大** ✓

### 14.4 💥 真根因：关节两软骨面不共形 ⇒ 从未闭合
（用**自检通过**的 pt-tri 实现测）
| 工况 | 最小间隙 | <0.05mm | 中位间隙 |
|---|---|---|---|
| 基线前向 | 0.0148 | **1/69** | **1.6202mm** |
| 基线反向 | 0.0224 | 1/69 | 1.6565 |
| m2 前向 | 0.0146 | 3/313 | 1.2771 |
初始 t=0 也是（中位 2.3810）。两软骨厚都 1.5mm（和 3.0），但两**外面**中位相距 2.38mm ⇒ 只在**一个点**相遇。

### 14.5 💥 后果：F 100% 罚因子伪影
| 关节 PENALTY | 1e-4 | 1e-3 | 1e-2 | 1e-1 | 1 |
|---|---|---|---|---|---|
| PENALTY 模式 | 0.83 | 6.99 | 50.65 | 220.08 | **419.56** |
| AUGLAG 模式 | — | 2.00 | 15.31 | 110.63 | **315.06** |
⇒ 跨 4 个数量级变 ~500×；**AUGLAG 也不收敛**。
位移扫描：2.0mm 后**饱和在 ~80N**（中位间隙停在 0.6mm）⇒ 几何上无法共形。
载荷截面移到关节上方 60/30mm：罚敏感比仍 **67×/69×**（`|u|max` 降到 0.68 但 F 没救）。

### 14.6 ⚠️ 影响面
`436.8N`/`419.56`/`190.80`/`315.06` **都不是收敛解** ⇒
「外推 600N 胫侧 ≈3.4MPa 对吴恺 3.74MPa 吻合」**失效**。
`[LRN-046~050]` 的 -0.99% 是 **tie 的** penalty 扫，关节 penalty 从未扫过。
「tet4 259.3 vs penta6 116.9（2.2×）」有了新解释路径（点接触+罚主导）⇒ **先修几何**。

### 14.7 ✅ 新验收判据（替换作废的 `|u|max ≈ 指定位移`）
1. `F` 对关节罚因子 ×10/×0.1 不敏感（<2%）
2. 末态关节面**中位间隙 ≈ 0**、接触面积 ≥30%（现在 **1.4%**）
3. 自动步长**主动放大**（不是靠 cutback 硬撑）
4. `|u|max` 仅参考，**非判据**

### 14.8 ➡️ 下一步
1. **修关节几何让两软骨面共形**（唯一真解）：
   ① 统一用一张关节面（talus 侧投到 tibia 侧）② 或给 0.05~0.1mm 初始过盈
   ③ 或按 Anderson：骨刚性 + 芯 E=12/ν=0.42 + **只取关节局部**
2. 几何修好后**重跑罚收敛**，F 稳定才算通过
3. 之后再重做 m2 vs baseline 对比

### 14.9 产物
`temp/pyfebio_demo/`：`matrix.py`·`diag_dfail.py`·`sweep_pen.py`·`aug_test.py`·
`disp_sweep.py`·`trunc_test.py`·`diag_tiehold.py`·`diag_rigid.py`·**`pttri.py`**（自检通过）·
`diag_gap3.py`·`diag_thick.py`·`diag_share.py`·`diag_frames.py`·`diag_reac.py`·
`cmp_run.py`·`read4.py`·`diag_tie2.py`
文档：`CART_HEX_TODO.md`（重写）·`TIE_DIVERGE_TODO.md`（结案）·`ANKLE_FIX_TODO.md`（警示）·
`docs/总体进度.md`（加行）·`MEMORY.md`（更正）·`[LRN-20261005-054]`
"""
if SESS.exists():
    t = SESS.read_text(encoding="utf-8")
    if "十四、「tie 发散」结案" not in t:
        SESS.write_text(t.rstrip() + SEC, encoding="utf-8")
    print(f"session note -> {SESS.stat().st_size:,} B")
if DAILY.exists():
    t = DAILY.read_text(encoding="utf-8")
    if "十四、「tie 发散」结案" not in t:
        DAILY.write_text(
            t.rstrip() + "\n- `2026-10-04/ankle-fe-solve-and-handover.md` → §十四 "
            "tie 发散结案（误诊）：真根因=**关节两软骨面不共形 ⇒ 关节未闭合 ⇒ "
            "F 全为罚因子伪影**（PENALTY 1e-4→1 ⇒ 0.83→419.56 N，AUGLAG 也不收敛）"
            "⇒ `[LRN-054]`，`436.8N` 失效\n", encoding="utf-8")
    print(f"daily note -> {DAILY.stat().st_size:,} B")

print("\n" + "=" * 76)
print("全量校验")
print("=" * 76)
ART = [("MEMORY.md", M), ("LEARNINGS.md", ROOT / ".learnings/LEARNINGS.md"),
       ("CART_HEX_TODO.md", ROOT / "CART_HEX_TODO.md"),
       ("TIE_DIVERGE_TODO.md", ROOT / "TIE_DIVERGE_TODO.md"),
       ("ANKLE_FIX_TODO.md", ROOT / "ANKLE_FIX_TODO.md"),
       ("FE_PIPELINE_TODO.md", ROOT / "FE_PIPELINE_TODO.md"),
       ("docs/总体进度.md", ROOT / "docs/总体进度.md"),
       ("scripts/ankle_fe/cart_patch.py", ROOT / "scripts/ankle_fe/cart_patch.py"),
       ("scripts/ankle_fe/build_feb.py", ROOT / "scripts/ankle_fe/build_feb.py"),
       ("pttri.py", ROOT / "temp/pyfebio_demo/pttri.py"),
       ("matrix.py", ROOT / "temp/pyfebio_demo/matrix.py"),
       ("sweep_pen.py", ROOT / "temp/pyfebio_demo/sweep_pen.py"),
       ("aug_test.py", ROOT / "temp/pyfebio_demo/aug_test.py"),
       ("disp_sweep.py", ROOT / "temp/pyfebio_demo/disp_sweep.py"),
       ("trunc_test.py", ROOT / "temp/pyfebio_demo/trunc_test.py")]
for k, p in ART:
    print(f"  [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>9,}  {k}")

lt = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
ms = M.read_text(encoding="utf-8")
CHK = [
    ("LRN-054 已入库", "LRN-20261005-054" in lt),
    ("LRN-054 含撤回①（刚体拟合）", "刚体拟合" in lt),
    ("LRN-054 含撤回②（滑移 3%）", "3.0%" in lt or "0.0476" in lt),
    ("LRN-054 含 D1/D2 变体", "5 negative jacobians" in lt),
    ("LRN-054 含 AUGLAG 扫描", "0.83" in lt and "2.00" in lt),
    ("LRN-054 含位移扫描饱和", "80.53" in lt),
    ("LRN-054 含载荷截面 H", "552.59" in lt or "573.14" in lt),
    ("MEMORY: F 非收敛", "所有 F 都不是收敛解" in ms),
    ("MEMORY: tie 无问题", "tie **没问题**" in ms),
    ("MEMORY: 体积 2.3× 已纠", "根因=接触区欠分辨**：" not in ms),
    ("CART_HEX: 已撤回 tie", "误诊" in (ROOT / "CART_HEX_TODO.md").read_text(encoding="utf-8")),
    ("TIE_DIVERGE: 结案", "已结案" in (ROOT / "TIE_DIVERGE_TODO.md").read_text(encoding="utf-8")),
    ("ANKLE_FIX: 警示", "重大更正" in (ROOT / "ANKLE_FIX_TODO.md").read_text(encoding="utf-8")),
    ("总体进度: LRN-054", "LRN-054" in (ROOT / "docs/总体进度.md").read_text(encoding="utf-8")),
]
for k, v in CHK:
    print(f"  [{'OK' if v else 'XX'}] {k}")

sz = M.stat().st_size
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb))`", ms)))
dead = [p for p in refs if not ((WS / p).exists() or (ROOT / p).exists())
        and not p.startswith(("memory/", "references/", "digest/", "docs/"))]
print(f"\n  MEMORY.md {sz:,}/10240 余量 {10240-sz:+,} {'OK' if sz<=10240 else '超限'}；"
      f"引用 {len(refs)} 死链 {len(dead)} {dead}")
