"""记录 P1c 探查：脚本混淆坑 + 骨网格构成 + 法向偏移反证。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
ND = ROOT / "docs/Anderson双论文_关键信息.md"
T = ROOT / "ANKLE_FIX_TODO.md"

s = ND.read_text(encoding="utf-8")
s = s.rstrip() + """

---

# §10 P1c 探查（2026-10-04）——骨网格构成 + 法向偏移**反证**

## 10.1 🔴 脚本混淆坑（浪费了一次实验）

**项目里有**两个**组装脚本，功能不同：**

| 脚本 | 用途 | 是否有 `CART_ELEM` 等开关 |
|---|---|---|
| `scripts/ankle_fe/build_feb.py` | A1–A4 冒烟组装（**定死** tibiotalar 单关节） | ❌ **没有** |
| `scripts/ankle_fe/ankle_lig_feb.py` | 全功能（韧带/截断/单元类型/驱动方式） | ✅ 有 |

`pconv.py` 跑的是 `ankle_lig_feb.py`；而 `FE_PIPELINE_TODO.md` 的「一键复现」写的是
`build_feb.py`。**我改了前者却跑了后者**，导致 `CART_ELEM=hex8` 完全没生效
（日志里软骨仍是 `tet4`），白跑一轮。

**⇒ 教训：改脚本前先确认**运行路径**用的是哪个文件；两个同名功能的脚本必须互相标注。**

## 10.2 骨网格构成（`npz_probe.py` 实测）

| part | 单元数 | 唯一节点数 | 实际类型 |
|---|---|---|---|
| `R_TIBIA_CORT` (81000700) | 2,417 | 8（2 个为 6） | **hex8** |
| `right_tibia_end_spon` (81000600) | 6,651 | **4** | **tet4** |
| `right_tibia_center_spon` (81000601) | 14,482 | **4** | **tet4** |
| `R_FIBULA_CORT` (81000900) | 842 | 8 | hex8 |
| `R_TALUS_CORT` (81001100) | 394 | 8 | hex8 |
| `R_TALUS_SPON` (81001000) | 1,766 | **4** | tet4 |
| `R_CALCANEUS_CORT` (81001300) | 634 | 8 | hex8 |
| `R_CALCANEUS_SPON` (81001200) | 3,323 | **4** | tet4 |

**⇒ 皮质骨 = hex8，松质骨 = tet4**（THUMS 的 `*ELEMENT_SOLID` 变长行被零填充成 8 列）。
**踝部合计 33,781 单元，其中约 87% 是 tet4。**

**这解释了 §9.1 里 25.4% 三角面的来源**（松质骨的四面体面）。

**节点共享（`share_probe.py`）：**
- **皮质 ↔ 松质：重合 50%**（tibia 2,418；fibula 844；talus 396；calc 636）
  ⇒ **共形界面，不能单独加密一侧**
- **各骨之间：0 重合** ⇒ 各骨是独立体（靠接触/韧带交互）✓

## 10.3 🔴 「沿顶点法向挤出」**反证：比默认更差**

受控对比（`ctrl_cmp.py`，用**正确的** `ankle_lig_feb.py`，
同一份 `CART_OFFSET=normal` 软骨，env 与基准一致）：

| `CART_ELEM` | rc | 用时 | 软骨 保留/剔除 | 接触面 |
|---|---|---|---|---|
| `tet4` | **1** ❌ | 9.0 s | 406/8、610/20 | 129/205 |
| `penta6` | 0 | 4.9 s | 88/**50**、146/**64** | 85/146 |
| **`hex8`** | 0 | 4.4 s | hex8 **仅 6/53**；penta6 25/7 | 31/60 |

**⇒ 三条结论：**
1. **`CART_OFFSET=normal` 的剔除率远高于默认 `toward`**
   （penta6 剔 50/88 = **57%**，默认约 14%；hex8 剔 **89%**）
2. **tet4 在法向偏移下直接失败**（负 Jacobian）—— 默认偏移下 tet4 反而能跑
3. **原代码注释是对的**：「直接用 hex8 会因偏移方向与面不垂直而歪斜」
   —— 我 §9.1 的探针只证明了**法向彼此对齐**（`mag ≈ Σ面积`），
   **不等于法向适合当挤出方向**。边缘节点的法向与「指向对侧」夹角可达 90°+（§9.1 的 36% 反向）

**⇒ 机制**：楔形（penta6）对倾斜**宽容**（偏移三角面只是平移，仍成楔形）；
而 hex8 要求偏移面**平行于**底面，一旦倾斜就退化成歪斜砖 ⇒ 被闸门剔除。

**⇒ 所以「先加密关节面」不是可选项，而是 hex8 路线的**前提**：**
只有把关节面做成**筛选干净、法向一致、足够细**的曲面片，
沿法向挤出才可能得到干净的 hex 层。

## 10.4 本轮收尾

- 已用备份恢复默认软骨（`cartilage.json` 93/129/397/500，与原值一致）
- 回归：`build_feb.py` → **`FEBio rc=0 / 17.2 s / 接触面 152/200`** ✓
- 新增只读探针：`npz_probe.py`、`share_probe.py`、`ctrl_cmp.py`、
  `norm_probe.py`（均在 `temp/pyfebio_demo/`）
"""
ND.write_text(s, encoding="utf-8")
print(f"Anderson双论文_关键信息.md -> {ND.stat().st_size:,}")

t = T.read_text(encoding="utf-8")
t = t.replace(
    """- [ ] **P1c（新首要）加密胫距/距下关节面骨表面网格** —— 关节面节点 ~72 → **≥300**（~1 mm 间距）；
      然后再做 hex8 层（此时反向节点自然减少）
- [ ] **P1d 软骨改 hex8 层**（四边形底 + 顶点法向挤出，Anderson 配方）""",
    """- [x] **P1c 探查** ✅ **完成**（`docs/Anderson双论文_关键信息.md` §10）
  - 🔴 **脚本混淆坑**：`build_feb.py`（定死、开关无效）≠ `ankle_lig_feb.py`（全功能）；
    改脚本前必须确认运行路径用哪个文件
  - 骨网格：**皮质 hex8 / 松质 tet4**（87% 是 tet4）；皮质↔松质**共节点 50%**⇒
    共形界面，**不能单独加密一侧**；各骨之间 0 重合
  - 🔴 **「沿顶点法向挤出」反证**：hex8 被剔 **89%**（6 保留/47 剔除）、
    penta6 剔 57%、tet4 **直接失败**（负 Jacobian）⇒ **法向偏移比默认更差**，
    原注释「hex8 会歪斜」**是对的**
- [ ] **P1d（新首要）加密 + 规范筛选关节面骨表面网格** —— 关节面节点 ~72 → **≥300**；
      且需**筛选干净**（不能是 `gap<6mm` 纯距离），法向一致
      ⇒ 这是 hex8 路线的**前提**，不是可选项
- [ ] **P1e 软骨改 hex8 层**（四边形底 + 顶点法向挤出）—— 待 P1d 完成后再试""")
t = t.replace("**当前：P0 ✅ / P1a-b ✅ 诊断完成；P1c 关节面加密为下一步**",
              "**当前：P0 ✅ / P1a-c ✅；P1d 加密+规范筛选关节面为下一步（hex8 的前提）**")
T.write_text(t, encoding="utf-8")
print(f"ANKLE_FIX_TODO.md -> {T.stat().st_size:,}")

L = ROOT / ".learnings/LEARNINGS.md"
l = L.read_text(encoding="utf-8")
if "LRN-20261004-042" not in l:
    l = l.rstrip() + """

## LRN-20261004-042: 同名同功能的两个脚本 = 定时炸弹；改之前先确认"跑的是哪个"

- **类别**: `correction` / `operational`
- **场景**: `scripts/ankle_fe/` 下有 `build_feb.py`（A1–A4 冒烟组装，**把关节数、
  单元类型全写死，没有 `CART_ELEM` 等开关**）和 `ankle_lig_feb.py`（全功能，
  有韧带/截断/驱动/单元类型全部开关）。两者功能重叠、名字相似。
  `FE_PIPELINE_TODO.md` 的「一键复现」用 `build_feb.py`，
  而 `pconv.py` 用 `ankle_lig_feb.py`。
  **我给 `ankle_lig_feb.py` 加了 `CART_ELEM=hex8` 分支，却用 `build_feb.py` 跑**，
  日志里软骨仍是 `tet4`，白跑一轮还差点误判「hex8 分支没生效」。
- **教训**:
  1. **加开关前先 `grep` 这个开关名在"即将运行的那个文件"里存不存在**
     （本次若先 grep 就会发现 `build_feb.py` 里根本没有 `CART_ELEM`）
  2. 两个功能重叠的脚本，**必须在文件头互相标注**「本文件 vs 另一个的差别、谁调用谁」
  3. 改完脚本后的**自检**：跑一次并确认日志里出现了**你期望的新行为**
     （如 `按类型统计` 里出现 `hex8`），不要只看 `rc`
- **行动**: 已在 `docs/Anderson双论文_关键信息.md` §10.1 记录；
  `ANKLE_FIX_TODO.md` P1c 标注
- **关联**: `[LRN-20261004-041]`（同属"改产物前后必须验证"）

## LRN-20261004-043: 局部指标好 ≠ 全局可用——"法向彼此对齐"不代表"法向适合挤出"

- **类别**: `knowledge_gap` / `best_practice`
- **场景**: 为做 hex8 软骨层，我先用只读探针验证「顶点法向是否可用」，
  测到 `合成模长 ≈ Σ面面积`（12 面 × 2.61 mm² ≈ 31.3）⇒ 判定**法向对齐良好，
  没有"相加抵消"**，于是认定可以沿法向挤出。
- **实际**: 受控对比显示 —— `CART_OFFSET=normal`（沿法向）**hex8 被剔 89%**
  （6 保留/47 剔除）、penta6 剔 57%、**tet4 直接负 Jacobian 失败**；
  而默认的 `toward`（指向对侧骨）反而稳定。
  **原代码注释「直接用 hex8 会因偏移方向与面不垂直而歪斜」是对的。**
- **错在哪**: 我测的是「**法向之间是否一致**」（局部/内部性质），
  但挤出真正需要的是「**法向与"指向对侧"方向是否一致**」（**全局/几何关系**）。
  边缘节点的法向与后者夹角可达 90°+（36% 反向），照样能把砖挤歪。
- **教训**: ① 验证一个量"能不能用"时，**要测它在使用场景下的判据**，
  而不是测它自身的自洽性（自洽性只是必要条件）
  ② 报告"某量良好"时，**必须写明它是针对哪个判据良好**
- **行动**: 结论写入 `docs/Anderson双论文_关键信息.md` §9.1 / §10.3
- **关联**: `[LRN-20261004-040]`（同属"假设验证不充分"）
"""
    L.write_text(l, encoding="utf-8")
    print(f"LEARNINGS.md -> {L.stat().st_size:,}（+LRN-042/043）")

print("\n抽查：")
nd = ND.read_text(encoding="utf-8"); td = T.read_text(encoding="utf-8")
for nm, ok in [("§10 落盘", "§10 P1c 探查" in nd),
               ("脚本混淆", "个**组装脚本" in nd or "两个**组装脚本" in nd),
               ("hex8 剔89%", "6/53" in nd),
               ("P1c 完成", "[x] **P1c" in td),
               ("P1d 前提", "P1d（新首要）" in td),
               ("LRN-042/043", "LRN-20261004-042" in L.read_text(encoding="utf-8"))]:
    print(f"  [{'OK' if ok else 'XX'}] {nm}")
