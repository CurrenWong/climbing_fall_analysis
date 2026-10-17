"""P0 结果：更正「软骨薄 2 倍」的误判；把「接触区欠分辨」立为新首要嫌疑。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
ND = ROOT / "docs/Anderson双论文_关键信息.md"

# ── ① Anderson 文档：更正 §2.1 ──
s = ND.read_text(encoding="utf-8")
old = """- **我的现值**：`t = min(1.5, 0.5×骨间隙)`，骨间隙实测 **1.406 mm** ⇒ 每层 **0.703 mm**
- **可能根因**：我测的是「骨面最小距离」，而 Anderson 是**先摆到负重中和位**再挤出软骨
- ⇒ **待查**：把 THUMS 胫距关节面的**间隙分布**（中位/中央区）算出来，
  若中央间隙 ≫1.406 mm，说明 1.406 是边缘最靠近点，**软骨厚度需要改**"""
new = """- **✅ 已查清（2026-10-04，读 `temp/thums/cartilage.json` + `joint_gaps.json`）**：
  **我的软骨厚度中位数就是 1.5 mm**，与 Anderson 2006 完全一致！
  `cartilage.json` 实测：tibiotalar_tibia **median 1.5** / talus **median 1.5**；
  subtalar 两侧 median 1.49 / 1.5。
- ❌ **更正**：此前据「最小距 1.406 mm」推断「软骨做薄 2 倍」是**误读**。
  1.406 mm 只出现在**边缘 1~4 个节点**（`joint_gaps.json`：
  `n_within_1.5mm = 1`、`n_within_2mm = 4`）；**0.703 mm 是最小厚度不是中位厚度**
- **⇒ 唯一要改的**：`T_TARGET` **1.5 → 1.7**（对齐 Anderson 2007），一行改动
- ⚠️ **不需要重新摆位来"加厚软骨"**（此前的推断作废）"""
assert old in s, "§2.1 anchor not found"
s = s.replace(old, new, 1)

# ── ② Anderson 文档：新增 §8 接触区欠分辨 ──
s = s.rstrip() + """

---

# §8 🔴 新首要嫌疑：**接触区严重欠分辨**（P0 实测，2026-10-04）

## 8.1 实测数据（`temp/thums/joint_gaps.json` + `cartilage.json`）

| 项 | 我的模型 | Anderson |
|---|---|---|
| 胫距关节面节点数（tibia 侧） | **73** | — |
| 胫距关节面节点数（talus 侧） | **89** | — |
| 胫骨软骨单元 | **93** | — |
| 距骨软骨单元 | **129** | — |
| 全模型单元 | **26,794** | 收敛需 **约 50,000** |
| 接触区节点间距（估） | **≈2 mm** | — |
| 网格收敛研究 | **没做过** | **13 次，2,997→293,402 单元** |

`joint_gaps.json` 里 tibia__talus 的分布（分母是 tibia **全部 4836 个表面节点**）：
`min 1.406 / p1 3.66 / p5 9.61 / median 223.8 mm`
**⇒ 中位数毫无意义（大部分节点在胫骨干上）；有意义的只有：**
`n_within_1.5mm = 1`、`n_within_2mm = 4`、`n_within_3mm = 16`、`n_within_5mm = 54~114`

## 8.2 为什么这可能是 2.2 倍刚度差的根因

1. **接触斑只有 ~220 个软骨单元**（93+129）覆盖约 300 mm² 的接触面
   ⇒ 单元尺度 ≈ 接触斑尺度的 1/10，**接触压力分布根本分辨不出来**
2. 用 **master/slave 接触**时，从面节点少 ⇒ 主面被"点接触"离散化
   ⇒ 有效接触刚度强烈依赖离散 ⇒ **tet4 与 penta6 的接触面片数不同
   （148/191 vs 138/171）可能正是差异来源**
3. **Anderson 专门做了 13 点网格收敛研究**才敢用 —— 我一次都没做
4. 解释了为什么 `CART_QUALITY=1` 无效（`LRN-020`）：它改善的是单元**形状**，
   **没有增加单元数量**

## 8.3 ⇒ 下一步（替代原「交换接触面」）

**做网格收敛研究**：把胫距接触区逐级加密，看
1. 接触应力 / 接触面积 / 支反力是否**收敛**
2. **tet4 与 penta6 是否收敛到同一个值**（若收敛到同一值 ⇒ 2.2 倍就是分辨率问题，结案）
"""
ND.write_text(s, encoding="utf-8")
print(f"docs/Anderson双论文_关键信息.md -> {ND.stat().st_size:,} bytes")

# ── ③ ANKLE_FIX_TODO.md：P0 完成 + P1 改成网格收敛 ──
T = ROOT / "ANKLE_FIX_TODO.md"
t = T.read_text(encoding="utf-8")
t = t.replace(
    """- [ ] **P0 量胫距关节面「间隙分布」**（不是最小值）——
      判定软骨是否做薄 2 倍（我 0.703 mm vs Anderson 1.5~1.7 mm）；
      同时判定是否需要重新摆中和负重位。**~20 min 脚本活**""",
    """- [x] **P0 量关节间隙/厚度分布** ✅ **完成（2026-10-04）**
      **结论：我的软骨中位厚度就是 1.5 mm，与 Anderson 2006 一致**；
      1.406 mm 只是边缘 1~4 个节点 ⇒ **「做薄 2 倍」是我误读最小值，已更正**。
      **顺带暴露真问题**：胫距关节面只有 **73/89 个节点、93/129 个软骨单元**，
      而 Anderson **做了 13 点收敛研究、收敛需 ~50,000 单元**"""
)
t = t.replace(
    """- [ ] **P1 软骨改「沿骨面法向挤出的六面体层」**——Anderson 全模型 all-hex；
      从根上绕开 tet4 体积锁死 vs penta6 沙漏的二选一（== 2 倍刚度差的根因方向）""",
    """- [ ] **P1（新首要）接触区网格收敛研究** —— 逐级加密胫距接触区，看
      ①应力/面积/支反力是否收敛 ②**tet4 与 penta6 是否收敛到同一值**
      （收敛到同一值 ⇒ 2.2 倍就是分辨率问题，结案）。
      依据：Anderson 13 次收敛研究 / 收敛需 ~50,000 单元 / 主从比 1:3
- [ ] **P1b 软骨改「沿骨面法向挤出的六面体层」**——Anderson all-hex；
      绕开 tet4 体积锁死 vs penta6 沙漏（与 P1 可能同源）
- [ ] **P1c `T_TARGET` 1.5 → 1.7 mm**（对齐 Anderson 2007）——一行改动"""
)
t = t.replace("**当前：P0 进行中**", "**当前：P0 ✅ 完成；P1 网格收敛研究为下一步**")
T.write_text(t, encoding="utf-8")
print(f"ANKLE_FIX_TODO.md -> {T.stat().st_size:,} bytes")

# ── ④ CONV_TODO.md：把「交换接触面」降为兜底，网格收敛提为首选 ──
C = ROOT / "CONV_TODO.md"
c = C.read_text(encoding="utf-8")
c = c.replace(
    """1. **软骨改「沿骨面法向挤出的六面体层」** —— Anderson 全模型 **all-hex**；
   这从根上绕开 tet4 体积锁死 vs penta6 沙漏的二选一（§3.1）""",
    """1. ⭐ **接触区网格收敛研究**（**新首选**）—— 实测胫距关节面只有 **73/89 节点、
   93/129 软骨单元**（`joint_gaps.json`: `n_within_5mm` 仅 114/4836），
   而 Anderson 做了 **13 点收敛研究、收敛需 ~50,000 单元、主从比 1:3**。
   **看 tet4 与 penta6 是否收敛到同一值** ⇒ 若是，2.2 倍即分辨率问题，结案。
   （`CART_QUALITY=1` 之所以无效，是因为它改单元**形状**而非**数量**）
2. **软骨改「沿骨面法向挤出的六面体层」** —— Anderson 全模型 **all-hex**；
   这从根上绕开 tet4 体积锁死 vs penta6 沙漏的二选一（§3.1）"""
)
c = c.replace(
    """3. **先量关节间隙分布** —— 我可能把软骨做薄了 2 倍（0.703 vs Anderson 1.5~1.7 mm）；
   疑因未摆到负重中和位（§2.1 / §3.3）""",
    """3. ~~先量关节间隙分布~~ ✅ **已做**：**软骨中位厚度就是 1.5 mm，与 Anderson 一致**；
   1.406 mm 只是边缘 1~4 个节点的值 ⇒ **「薄 2 倍」是我误读，已更正**；
   只需 `T_TARGET` 1.5→1.7 对齐 2007（`docs/Anderson双论文_关键信息.md` §2.1 / §8）"""
)
c = c.replace(
    "4. 若仍差 ⇒ 交换接触面 / 逐单元贡献对比（原方案，降为兜底）",
    "4. 若仍差 ⇒ 交换接触面 / 逐单元贡献对比（原方案，**降为兜底**）"
)
C.write_text(c, encoding="utf-8")
print(f"CONV_TODO.md -> {C.stat().st_size:,} bytes")

# ── ⑤ 更正 + 新假设 记入 LEARNINGS ──
L = ROOT / ".learnings/LEARNINGS.md"
l = L.read_text(encoding="utf-8")
if "LRN-20261004-039" not in l:
    l = l.rstrip() + """

## LRN-20261004-039: 分布统计只看 min/max 会得出反向结论——必须先看中位数与支撑集

- **类别**: `correction` / `knowledge_gap`
- **场景**: 我据「踝关节骨间最小距 1.406 mm」+ `t = min(1.5, 0.5×gap)`
  推断「软骨每层只有 0.703 mm，比 Anderson 的 1.5~1.7 mm **薄 2 倍**」，
  并写进了三份文档、还列为最高优先修正项。
- **实际**（读 `temp/thums/cartilage.json`）：**厚度中位数 = 1.5 mm**，
  与 Anderson 2006 完全一致。**0.703 是最小厚度，不是典型厚度。**
  `joint_gaps.json` 证实 1.406 mm 只出现在 **边缘 1~4 个节点**
  （`n_within_1.5mm = 1`、`n_within_2mm = 4`）。
- **同类陷阱**：`joint_gaps.json` 里 tibia__talus 的 `median = 223.8 mm` ——
  因为分母是 tibia **全部 4836 个表面节点**（大部分在胫骨干上）。
  **分位数一旦跨过"支撑集"就完全无意义；有意义的只有
  `n_within_5mm = 114` 这种带阈值的计数。**
- **教训**: ① 报"厚度/间隙"必须同时给 **min / 中位 / max + 分位数**，
  只给 min 会被极端值劫持 ② 算分位数前先确认分母是**同一支撑集**
  （先按阈值筛出关节面节点，再统计） ③ **改动前先读已有的 JSON 产物**——
  这次的答案早就在 `cartilage.json` 里躺着，我却先写了两份文档去猜
- **行动**: 已在 `docs/Anderson双论文_关键信息.md` §2.1 更正 + 新增 §8；
  `ANKLE_FIX_TODO.md` P0 标记完成；`CONV_TODO.md` 表述更正
- **关联**: `[LRN-20261004-038]`（开关跑过≠生效，同属"没先读现有产物"）
"""
    L.write_text(l, encoding="utf-8")
    print(f"LEARNINGS.md -> {L.stat().st_size:,} bytes（+LRN-039）")

# ── ⑥ 抽查 ──
nd = ND.read_text(encoding="utf-8")
todo = T.read_text(encoding="utf-8")
print("\n抽查：")
for nm, ok in [("§2.1 已更正", "已查清" in nd and "误读" in nd),
               ("§8 接触区欠分辨", "接触区严重欠分辨" in nd),
               ("P0 完成", "[x] **P0" in todo),
               ("P1 改网格收敛", "接触区网格收敛研究" in todo),
               ("CONV 换首选", "新首选" in c),
               ("LRN-039", "LRN-20261004-039" in L.read_text(encoding="utf-8"))]:
    print(f"  [{'OK' if ok else 'XX'}] {nm}")
