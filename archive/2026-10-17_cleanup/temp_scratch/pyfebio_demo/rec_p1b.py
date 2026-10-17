"""P1 诊断结论固化 + 计划更新。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
ND = ROOT / "docs/Anderson双论文_关键信息.md"
T = ROOT / "ANKLE_FIX_TODO.md"

s = ND.read_text(encoding="utf-8")
s = s.rstrip() + """

---

# §9 P1 诊断实测（2026-10-04）——「关节面太粗」是根问题

探针：`temp/pyfebio_demo/norm_probe.py`（只读几何，不改模型）

## 9.1 三条硬结论（其中两条**推翻旧注释**）

| # | 结论 | 证据 |
|---|---|---|
| 1 | ✅ 底网格**已是四边形为主** | tibia 外表面 19,022 面 = **tri 4,832 / quad 14,190（74.6%）**；talus 75.0% |
| 2 | ✅ **顶点法向完全可用，没有"相加抵消"** | 合成模长中位 **31.28 ≈ 12 面 × 2.61 mm²**（= Σ面积）⇒ 法向几乎完美对齐。**旧注释「面积加权法向严重不一致、相加几乎抵消」经实测为误判** |
| 3 | 🔴 `gap<6mm` 是**纯距离**判据 ⇒ 关节面选择**严重掺假** | 73 个"关节面"节点中：夹角>60° 占 **40/73（55%）**；**>90°（法向反向）占 26/73（36%）**，p10 cos = **−0.91** |

**偏移方向实测差**（同节点、同 1.5 mm 厚度）：
中位 **1.587 mm** / p90 **2.932 mm** ⇒ p90 ≈ 2×1.5 ⇒ **法向与"指向对侧"在 1/3 节点上近乎相反**
⇒ 旧代码用"指向对侧"是**为了躲开反向节点**，属于**对的妥协**，但代价是偏移不垂直于面 ⇒ hex8 歪斜

## 9.2 但「按法向对齐过滤」**行不通**（实测）

新增开关 `CART_OFFSET=normal` + `ART_ALIGN=0.5` 后重跑 `cartilage.py`：

| 关节 | 关节面节点 | 软骨单元 | 体积 |
|---|---|---|---|
| tibiotalar_tibia | 73 → **33** | 93 → **18** | 1.057 → 0.263 cm³ |
| tibiotalar_talus | 89 → **36** | 129 → **27** | 1.109 → 0.486 cm³ |
| subtalar_talus | 212 → 92 | 397 → 95 | 3.059 → 1.132 cm³ |
| subtalar_calcaneus | 253 → 119 | 500 → 149 | 2.706 → 1.285 cm³ |

**⇒ 过滤把关节面砍掉一半，软骨只剩 18/27 个单元**（原 93/129）；
且退化丢弃率仍约 33%（18 保留 / 9 丢弃）。**这不是可用状态，已回滚。**

## 9.3 🎯 结论：**根问题是「关节面网格太粗」，不是「偏移方向选错」**

`GAP_MAX` 敏感性（`cartilage.py` 自带输出）：

| bone | <3mm | <4mm | **<6mm** | <8mm | <10mm |
|---|---|---|---|---|---|
| tibia | 16 | 31 | **72** | 103 | 135 |
| talus | 18 | 38 | **87** | 126 | 165 |

- 胫骨下关节面约 300~500 mm²，却只由 **72 个节点**代表 ⇒ 节点间距 **≈2.0~2.6 mm**
- Anderson 收敛需 **~50,000 单元**、并专门做了 **13 点收敛研究**
- ⇒ **必须先把关节面加密**，才谈得上"用法向挤 hex8"

## 9.4 ⇒ P1b 正确顺序（修正）

1. **先加密胫距（及距下）关节面的骨表面网格**（THUMS hex8 面级细分 or 局部重划）
   —— 目标：关节面节点数从 ~72 提到 **≥300**（间距 ~1 mm 级）
2. **再**按 Anderson：四边形底 + 顶点法向挤出 ⇒ **hex8 层**（此时反向节点会自然减少）
3. 用 `CART_OFFSET=normal` / `ART_ALIGN` 已就位的开关做 A/B
4. 最后做网格收敛研究（逐级加密，看 tet4/penta6/hex8 是否收敛到同一值）

**已就位的工具**（本轮新增，默认关闭不影响现有流程）：
`CART_OFFSET=toward|normal`、`ART_ALIGN=<cos 阈值>`（`scripts/ankle_fe/cartilage.py`）

## 9.5 ⚠️ 操作提醒

本轮曾用 `CART_OFFSET=normal ART_ALIGN=0.5` 跑 `cartilage.py`，
**覆盖了默认的 `cartilage.json` / `cartilage_mesh.npz`**。
**已立即用默认参数重跑恢复**（`tibiotalar_tibia 93 / 1.057cm³`、
`talus 129 / 1.109cm³`），并跑通回归 `FEBio rc=0 / 17.8 s`、接触面 `152/200`。
⇒ **教训：`cartilage.py` 直接写死输出路径，试参数前应先备份或改用独立 OUT。**
"""
ND.write_text(s, encoding="utf-8")
print(f"Anderson双论文_关键信息.md -> {ND.stat().st_size:,}")

t = T.read_text(encoding="utf-8")
t = t.replace(
    """- [ ] **P1b（首要）软骨改「四边形底 + 顶点法向挤出的 hex8 层」**（Anderson 配方）
      —— 需改 `cartilage.py`：
      ① 关节区生成**四边形**面网格（而非三角形）
      ② 顶点法向 = **相邻四边形外法向的平均**（不用「指向对侧最近点」）
      ③ 沿法向挤出 → **hex8 层**
      验收：hex8 被剔除数 ≈0（此前 42/61）、接触面完整、跑通并与 penta6 对比""",
    """- [x] **P1b 法向/偏移方向诊断** ✅ **完成**
  - ✅ 底网格已是**四边形为主**（tibia 74.6% quad）
  - ✅ **顶点法向可用**（合成模长 ≈ Σ面积）⇒ **旧注释「法向相加几乎抵消」是误判**
  - 🔴 `gap<6mm` 纯距离判据混入 **36% 法向反向**、55% 夹角>60° 的侧壁节点
  - 🔴 但「按法向对齐过滤」会把关节面**砍半**（73→33），软骨只剩 **18/27** 单元 ⇒ 不可用
  - 🎯 **结论：根问题是关节面网格太粗**（tibia 仅 72 节点 / ≈2.0~2.6 mm 间距）
  - 已加开关（默认关）：`CART_OFFSET=toward|normal`、`ART_ALIGN=<cos>`
- [ ] **P1c（新首要）加密胫距/距下关节面骨表面网格** —— 关节面节点 ~72 → **≥300**（~1 mm 间距）；
      然后再做 hex8 层（此时反向节点自然减少）
- [ ] **P1d 软骨改 hex8 层**（四边形底 + 顶点法向挤出，Anderson 配方）""")
t = t.replace("**当前：P0 ✅ 完成；P1 网格收敛研究为下一步**",
              "**当前：P0 ✅ / P1a-b ✅ 诊断完成；P1c 关节面加密为下一步**")
T.write_text(t, encoding="utf-8")
print(f"ANKLE_FIX_TODO.md -> {T.stat().st_size:,}")

L = ROOT / ".learnings/LEARNINGS.md"
l = L.read_text(encoding="utf-8")
if "LRN-20261004-041" not in l:
    l = l.rstrip() + """

## LRN-20261004-041: 只读探针要趁早；写死输出路径的脚本，试参数前必须备份

- **类别**: `correction` / `best_practice`
- **场景 A（正面）**: 「顶点法向是否可用」这个卡了很久的问题，
  写一个 **只读探针**（`norm_probe.py`，只读 npz、不改模型）15 秒就出结论：
  **合成模长 ≈ Σ面面积 ⇒ 法向对齐良好，旧注释「相加几乎抵消」是误判**。
  同时量化出真问题：`gap<6mm` 纯距离判据混入 **36% 法向反向节点**。
  **⇒ 卡住时优先写"只读、可量化"的探针，不要先改模型试。**
- **场景 B（踩坑）**: `cartilage.py` 把输出路径**写死**为
  `temp/thums/cartilage.json` + `cartilage_mesh.npz`。
  我用新开关跑了一次，**直接覆盖了默认产物**（软骨单元 93/129 → 18/27），
  会把已验证的模型搞坏。所幸立即用默认参数重跑恢复，并跑回归确认
  （`FEBio rc=0 / 17.8 s`、接触面 `152/200`）。
- **教训**: ① 诊断类探针**一律只读**，绝不写入模型产物
  ② 脚本输出路径写死时，**试参数前先备份**或加 `OUT` 环境变量
  ③ 改了产物之后**必须跑一次回归**（rc + 关键数字）再报告
- **行动**: 探针留在 `temp/pyfebio_demo/norm_probe.py`；
  结论写入 `docs/Anderson双论文_关键信息.md` §9
- **关联**: `[LRN-20261004-039]`/`[040]`（同属"先读/先验现状"）
"""
    L.write_text(l, encoding="utf-8")
    print(f"LEARNINGS.md -> {L.stat().st_size:,}（+LRN-041）")

print("\n抽查：")
nd = ND.read_text(encoding="utf-8"); td = T.read_text(encoding="utf-8")
cj = (ROOT / "temp/thums/cartilage.json").read_text(encoding="utf-8")
for nm, ok in [("§9 落盘", "§9 P1 诊断实测" in nd),
               ("旧注释推翻", "误判" in nd),
               ("P1b 完成", "[x] **P1b" in td),
               ("P1c 加密", "P1c（新首要）" in td),
               ("软骨已恢复 93", '"n_cartilage_elements": 93' in cj),
               ("LRN-041", "LRN-20261004-041" in L.read_text(encoding="utf-8"))]:
    print(f"  [{'OK' if ok else 'XX'}] {nm}")
