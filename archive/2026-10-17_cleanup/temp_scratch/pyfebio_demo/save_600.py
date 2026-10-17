"""固化 600 N 突破：更新 TODO + 论文重解析文档 + learnings。"""
import ast
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")

# ---------- 1. ANKLE_600N_TODO.md 勾掉完成项 ----------
T = ROOT / "ANKLE_600N_TODO.md"
t = T.read_text(encoding="utf-8")
t = t.replace("- [ ] 1. 诊断：从 FEBio 日志定位**具体哪些单元**翻反（域名/单元号）",
              "- [x] 1. 诊断：**FEBio 不报单元号** ⇒ 改为自己算 det(F)。"
              "★ 真凶 = **薄软骨(0.703 mm)楔形被拆成 3 个 tet4** ⇒ "
              "细长 sliver，末态 det(F) 已跌到 **0.523** ⇒ 再压就翻"
              "（tibiotalar_tibia；talus 侧 0.729；其余全部 ≥0.9994）")
t = t.replace("- [ ] 2. 试 `CART_QUALITY=1`（已有开关，此前判\"负面结果\"是针对**收敛速度**，此处目标**防翻转**，目标不同）",
              "- [x] 2. **改用 `CART_ELEM=prism`**（penta6 楔形，不再拆 sliver tet4）"
              " ⇒ **模型从「0.378 mm 就翻」变成能跑到 1.20 mm / 645.7 N / rc=0 / 7.4 s**")
t = t.replace("- [ ] 3. 若不够：软骨改 **hex 挤出**（薄壳状 hex 抗压缩远好于 tet4）或加层数",
              "- [x] 3. 已试 hex8（四边形底）：**歪斜 det(J)≈0.5，被剔 42/61 ⇒ 软骨破洞、"
              "接触面 150/198→51/64**，弃用；改「四边形也拆 2 个 penta6」⇒ 接触面回到 150/198")
t = t.replace("- [ ] 4. 若仍不够：接触调参（`two_pass` / `auto_penalty` / 减小 `tolerance`）",
              "- [x] 4. 接触调参**未做**（不需要了）：验证 `DET_MIN` 0.5→0.02 几乎无变化"
              "（117 N 不变），说明**剔除不是主因**")
t = t.replace("- [ ] 5. 跑到 600 N → 出 **3.74 MPa（应力）+ 335.5 mm²（接触面积）** 双靶对照表",
              "- [x] 5. **600 N 达到**（`DISP_MM=1.20` → 645.7 N）。**双靶对照**：\n"
              "      | 靶子 | 吴恺 | 本模型 | 偏差 |\n|---|---|---|---|\n"
              "      | 最大接触应力 | 3.74 MPa | 软骨 von Mises max **2.96** | −21% |\n"
              "      | 接触面积 | 335.5 mm² | **287.4 mm²**（间隙<0.20 mm） | −14% |")
t = t.replace("- [ ] 6. 更新 `docs/吴恺论文_重解析关键信息.md` + `ANKLE_LIG_TODO.md` + 记忆",
              "- [x] 6. 已更新文档 + 记忆\n\n"
              "## ⚠️ 遗留问题（未解决，必须知道）\n\n"
              "**楔形 vs 四面体的刚度差 2.2 倍**：同一位移 0.1472 mm 下\n"
              "- tet4 版：**259.3 N**，软骨 p99 0.765 MPa\n"
              "- prism 版：**116.9 N**，软骨 p99 0.273 MPa\n\n"
              "同几何、同材料、同边界 ⇒ 差异只能来自**离散格式**：\n"
              "线性 tet4 在 ν=0.4 下**体积锁死偏刚**，线性 penta6 **偏柔**（可能沙漏）。\n"
              "⇒「哪个更接近真值」未定。**做网格收敛研究才能定**（加密后两者应趋同）。\n"
              "**在此之前，600 N 档的数值只能当量级参考。**\n\n"
              "## 📌 本轮修好的两个**自己的** bug（重要教训）\n\n"
              "1. **hex8 的 det(F) 基向量必须从 node0 出发**：`(0,1),(0,3),(0,4)`。\n"
              "   我原先用 `(0,1),(2,3),(4,7)`，其中两对近乎平行 ⇒ det(G)≈0\n"
              "   ⇒ 纯刚体平移都算成 det=**−1.34**，**误报腓骨单元翻转**（白白追查一轮）。\n"
              "   判据：**刚体平移必须给出 det(F)=1**，不会就该怀疑公式。\n"
              "2. **楔形的 Jacobian 必须采样**（沿 r,s,t 扫），单点算会漏：\n"
              "   底面形心 det>0 但三角内部/沿 t 退化 ⇒ FEBio 在积分点报\n"
              "   `Negative jacobian detected during mesh initialization`。\n\n"
              "## 📌 FEBio 4.13 的 `penta6` 节点序（**实测确定，别再猜**）\n\n"
              "底面三角 `0,1,2` + 顶面三角 `3,4,5`（**同向对应**）\n"
              "实测 6 种排列，只有 `[1,2,3,4,5,6]` 通过；顶面反序/底三角换向全部\n"
              "报初始化负 Jacobian。\n")
T.write_text(t, encoding="utf-8")
print(f"ANKLE_600N_TODO.md -> {T.stat().st_size:,} bytes")

# ---------- 2. 论文重解析文档：补双靶结果 ----------
D = ROOT / "docs/吴恺论文_重解析关键信息.md"
s = D.read_text(encoding="utf-8")
ADD = """

---

## 八、✅ 600 N 打通 + 双靶验证结果（2026-10-04）

### 打通的关键：软骨改用 `penta6`（楔形）

原实现把薄软骨（0.703 mm）的 6 节点楔形**拆成 3 个 tet4** ⇒ 细长 sliver
⇒ 末态 det(F) 跌到 **0.523**（tibiotalar_tibia）⇒ 再压就翻负
（`DISP_MM=0.60` 在 t=0.63 / 0.378 mm 处失败）。

改用 `CART_ELEM=prism`（**penta6，不拆**）后：

```
DRIVER=disp TOPFIX=lateral TRUNC_MM=60 FIX_TALUS=all CART_ELEM=prism
TIME_STEPS=30 DISP_MM=1.20
→ rc=0 / 7.4 s / 支反力 645.7 N
```

### 🎯 双靶对照（600 N 档）

| 靶子 | 吴恺 2012 | 本模型（645.7 N） | 偏差 |
|---|---|---|---|
| **最大接触应力** | **3.74 MPa** | 软骨 von Mises max **2.96 MPa**（p99 0.97） | **−21%** |
| **接触面积** | **335.5 mm²** | **287.4 mm²**（间隙阈值 0.20 mm） | **−14%** |

接触面积对阈值的敏感性：0.10 mm → 109.8；**0.20 mm → 287.4**；0.50 mm → 722.9 mm²。

**⇒ 两个靶子都在同一量级、偏差 ~20% 以内**，比此前只比软骨应力 p99 硬得多。

### ⚠️ 遗留：楔形 vs 四面体刚度差 2.2 倍（未解决）

同一位移 `DISP_MM=0.1472`：

| 软骨单元 | 支反力 | 软骨 p99 | 接触面 |
|---|---|---|---|
| **tet4** | **259.3 N** | 0.765 MPa | 148 / 191 |
| **prism** | **116.9 N** | 0.273 MPa | 138 / 171 |

同几何/材料/边界 ⇒ 差异来自**离散格式**：线性 tet4 在 ν=0.4 下体积锁死偏刚、
线性 penta6 偏柔（可能沙漏）。**哪个更接近真值未定，需网格收敛研究**。
在此之前 600 N 档数值只能当量级参考。

（已排除：`DET_MIN` 0.5→0.02 几乎无变化 ⇒ 元素剔除不是主因。）
"""
D.write_text(s.rstrip() + ADD, encoding="utf-8")
print(f"docs/吴恺论文_重解析关键信息.md -> {D.stat().st_size:,} bytes")

# ---------- 3. learnings ----------
L = ROOT / ".learnings/LEARNINGS.md"
s = L.read_text(encoding="utf-8")
ADD = """

## [LRN-20261004-030] 刚体平移必须给出 det(F)=1 —— 这是检验变形梯度公式的判据

- **类别**: correction（自查方法）
- **场景**: 怀疑 `fibula_cort__hex8` 有 2 个单元翻转（det min=−1.34），
  追查一轮后才发现是**我自己的 det(F) 公式取错了基向量**。
- **错误做法**: hex8 用 `(0,1),(2,3),(4,7)` 三对边做变形梯度；
  这些"边"里有两对近乎平行 ⇒ `det(G)≈0` ⇒ 商爆炸。
- **正确做法**: 从 **node0** 出发的三条边 `(0,1),(0,3),(0,4)`。
- **判据（可复用）**: 拿到位移场后，先喂一个**纯刚体平移**的假位移
  （全部节点加同一向量），**det(F) 必须恰好 =1**。不是 1 就说明公式错了，
  别急着解释物理。
- **代价**: 一轮误诊 + 差点据此去"修腓骨"。
- **相关**: LRN-20261004-031

## [LRN-20261004-031] 楔形(penta6)的 Jacobian 必须**采样**，单点会漏

- **类别**: knowledge_gap
- **症状**: `CART_ELEM=prism` 后 FEBio 报
  `Negative jacobian detected during mesh initialization`，
  但我的单点检查（底面形心）**全部为正**（min 0.96/0.74）。
- **原因**: 楔形的 Jacobian 沿 r,s,t 变化。底面形心 det>0，
  但三角**内部**或沿 t 的退化会让积分点 det<0。
  另外 `dt` 列必须是**上下三角形心差** `(P3+P4+P5-P0-P1-P2)/3`，
  写成 `P3-P0` 在偏移方向沿三角变化时也会算错。
- **正解**: 在 `(r,s)∈{0,1/6,1/3,1/2,2/3}`（r+s≤1）× `t∈{0,0.25,0.5,0.75,1}`
  上扫，取 min det。
- **修好后**: tet4/prism/无韧带三种配置全部 rc=0。
- **相关**: LRN-20261004-030

## [LRN-20261004-032] FEBio 4.13 的 `penta6` 节点序 = 底三角 0,1,2 + 顶三角 3,4,5（同向）

- **类别**: best_practice（实测确定，别再猜）
- **做法**: 用最小单单元模型（底三角 (0,0,0)(1,0,0)(0,1,0)、顶面抬高 1）
  遍历 6 种排列跑 FEBio，只看是否报初始化负 Jacobian。
- **结论**: **只有 `[1,2,3,4,5,6]`（底+顶同序对应）通过**；
  顶面反序/底三角换向/整体换序全部失败。
- **顺带教训**: 最小探针 .feb 必须自带 `<solver type="solid">`，
  否则 FEBio 报 `Component "" needs to have property "solver" defined`，
  会被误当成"元素序错误"。
"""
L.write_text(s.rstrip() + ADD, encoding="utf-8")
ast.parse("pass")
print(f".learnings/LEARNINGS.md -> {L.stat().st_size:,} bytes")

# ---------- 4. .gitignore ----------
G = ROOT / ".gitignore"
g = G.read_text(encoding="utf-8") if G.exists() else ""
if "ANKLE_600N_TODO.md" not in g:
    G.write_text(g.rstrip() + "\nANKLE_600N_TODO.md\n", encoding="utf-8")
    print(".gitignore += ANKLE_600N_TODO.md")
else:
    print(".gitignore 已有")
