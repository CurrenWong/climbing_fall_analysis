"""落盘 ATFL 长度异常 + 更新使用说明的局限。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
LEARN = os.path.join(R, ".learnings", "LEARNINGS.md")
L = """
### 5. ATFL 为何在内翻下几乎不受力（`[LRN-20261005-073]` 追查）

- **现象**：内翻（`POSE_IE=−3`）时最该受力的 **ATFL 应变只有 0.197%**，
  而 CFL 3.513%、MD2 3.708%。
- **根因（实测）**：**ATFL 的建模长度 ≈ 35.1~37.2 mm，而解剖参考只有 ~15~20 mm
  ⇒ 约 2 倍过长**。同样的绝对伸长，被换算成应变时被**稀释 ~2.3×**。
- **对照（建模长度 vs 解剖参考）**：

  | 韧带 | 建模长度中位 (mm) | 解剖参考 | 判定 |
  |---|---|---|---|
  | PTFL | 13.5~15.3 | 15~20 | ✓ |
  | CFL | 24.5~34.8 | 20~30 | ~✓ |
  | **ATFL** | **35.1~37.2** | **15~20** | **❌ ~2× 过长** |
  | MD1/MD2 | 22.9~27.6 | 20~30 | ✓ |
  | MD3 | 36.4~37.5 | 30~40 | ✓ |

- **⇒ 结论**：韧带应变**不能直接**当作损伤指标 —— 先要保证各韧带的**几何/插入点**正确。
  这是"判断骨折风险"从定性走向定量之前必须先修的一项。
- **另注**：`ATFL` 的 4 个 pair 全部连到同一个节点 1497（腓骨侧），
  与 3 个不同伙伴节点配对 ⇒ 也值得核对拓扑是否合理。
"""
txt = io.open(LEARN, encoding="utf-8").read()
if "ATFL 为何在内翻下几乎不受力" in txt:
    print("★ 已存在，跳过")
else:
    io.open(LEARN, "a", encoding="utf-8").write(L)
    print(f"LEARNINGS.md {len(txt)} -> {len(txt)+len(L)}")

# 使用说明：把 ATFL 的"待查"改成"已查明"
D = os.path.join(R, "docs", "踝FE可用模型使用说明.md")
d = io.open(D, encoding="utf-8").read()
OLD = """**ATFL 内翻时几乎不受力**（0.2%）—— 最该被拉的韧带
> 却不动，怀疑几何/插入点或壳单元受力口径，值得单查；"""
NEW = """**ATFL 内翻时几乎不受力**（0.2%）—— **已查明根因**：
> 其**建模长度 ≈ 35~37 mm**，而解剖参考只有 ~15~20 mm（**约 2 倍过长**）⇒ 同样的绝对伸长
> 被稀释 ~2.3×。⇒ **韧带应变不能直接当损伤指标，须先修插入点几何**；"""
if d.count(OLD) == 1:
    d = d.replace(OLD, NEW)
    io.open(D, "w", encoding="utf-8").write(d)
    print(f"使用说明已更新 ATFL 条目（{len(d.encode('utf-8'))} B）")
else:
    print(f"⚠️ 使用说明 ATFL 条目未命中（{d.count(OLD)}）")
