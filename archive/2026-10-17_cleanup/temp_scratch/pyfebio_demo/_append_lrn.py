"""追加本次踝关节建模的经验到 .learnings/LEARNINGS.md（幂等）。"""

import re
from pathlib import Path

P = Path(__file__).resolve().parents[2] / ".learnings" / "LEARNINGS.md"

ADD = r"""
---

## [LRN-20261003-010] 「pip 装了同一个库」≠「你能用它的全部 API」

**类别**：`best_practice`（工具排查）

**背景**：要把踝关节 STL 转四面体网格。`uv pip install gmsh` 装了 4.15.2，
但所有 STL→实体的调用全失败：

| 调用 | 结果 |
|---|---|
| `occ.importShapes(stl, format="stl")` | `Unknown file type`（6 种 format 全试） |
| `occ.reclassify(1,1)` | `AttributeError`（不存在） |
| `geo.addTriangle(...)` | `AttributeError`（geo 无此函数） |
| `mesh.importStl(path)` | `takes 0 positional arguments` |
| `merge(stl)` → `healShapes` | `Unknown OpenCASCADE entity of dimension 2` |

而 `General.BuildInfo` 显示 **OCC 7.8.1 确实编译进去了**（还带 TetGen/Netgen/Mmg）。

**我当时做了错事**：在同一套残缺 API 上反复换写法试了 5、6 轮。
用户已经选了"换带完整几何内核的网格器"，我却还在原地打转。

**根因**：pip 的 gmsh wheel **只装了 Python 绑定 `gmsh.py` + 一个精简 dll**，
没有 `gmsh.exe`，部分 OCC 绑定未导出。官方 SDK 是完全不同的东西。

**教训**：
1. **"同一个库/同一个版本号"不代表"同一份二进制"**。pip wheel 和官方 SDK
   可以是两个能力集。判断依据应该是**实测 API 是否存在/可调用**，不是版本号。
2. 用户说"换一个东西"时，**要换真的不同的东西** ——
   在已确认失败的载体上换 API 写法，是伪装的重复尝试。
   **判据**：如果重试用的还是同一个文件/同一个进程，那就不是换方案。
3. 找官方示例比试 API 快得多。官方 `examples/api/stl_to_brep.py` 一读就明白
   正确做法是**逐三角形手工重建 BREP**，而不是我猜的那些函数。

**关联**：`scripts/ankle_fe/geom.py`、`ANKLE_FE_TODO.md` §关键技术记录

---

## [LRN-20261003-011] 「拓扑干净」不等于「几何可缝合」

**类别**：`knowledge_gap`

**背景**：踝关节软骨 STL 被 OCC 判为非封闭（18 条边界边），无法建实体。
我依次验证并**排除了**两个假设：

| 假设 | 检验方法 | 结果 |
|---|---|---|
| 有开口 | 统计每条边被几个三角形共享 | **13581 条全是 2 次，0 例外** → 无开口 |
| 非流形 | 统计 >2 次的边 | **0 条** |
| 绕向不一致 | 统计共享边是否被反向遍历 | **13581 反向 / 0 冲突** → 完全一致 |
| 精度不够 | 扫合并容差 1e-8…1e-3 | 18→10→5→**0**，但**体积同时被毁** |

**结论**：这是一个**拓扑完全干净**的闭壳，却无法被 OCC 缝合成实体。
加大容差确实能"封口"，但代价是体积从 2.33 cm³ 变成 **-0.12 cm³（负！）**
或 4.15 → 0.036 cm³。**封口成功 = 几何被毁**。

**教训**：
1. **"边被恰好 2 个面共享"是拓扑条件，不是几何条件。**
   几何缝合还要求相邻面在共享边处**连续**（切向/位置），
   而这需要 OCC 的实际缝合容差内可匹配。细长片（这里 19–243 μm）被合并后
   消失，就留下了 OCC 认得的洞，而拓扑统计（在合并前的点上做）看不出来。
2. **要交叉验证两种"闭壳"判据**：VTK 说闭合（`n_open_edges=0`）而 OCC 说
   非封闭时，**以 OCC 为准**（它要真缝合）。反过来不能只看 VTK。
3. **参数扫描要一直扫到"能做/不能做"的边界，并检查解的可用性**。
   这次如果只扫到 tol=1e-3 看到"✅ 封闭"就收工，
   就会拿到一个**体积为负**的实体 —— 比失败更糟。
   **看到成功要立刻质疑：这个"成功"的产物还有意义吗？**
4. 遇到缺陷源数据，**换生成方式往往比修数据更划算**：吴恺原文的软骨本来就是
   从骨面 offset 生成的，不是分割出来的 —— 绕开缺陷同时还更忠于论文。

**关联**：`temp/pyfebio_demo/{cart_diag,cart_tol,cart_wind,cart_close}.py`、
`ANKLE_FE_TODO.md` §阻塞
"""

text = P.read_text(encoding="utf-8")
ids = sorted(set(re.findall(r"LRN-\d+-\d+", text)))
print("现有:", ids)

for title, block in (
    ("## [LRN-20261003-010] 「pip 装了同一个库」", ADD.split("---\n\n## ")[1]),
    ("## [LRN-20261003-011] 「拓扑干净」不等于", ADD.split("---\n\n## ")[2]),
):
    if title in text:
        print("skip (exists):", title)
    else:
        P.write_text(P.read_text(encoding="utf-8").rstrip() + "\n\n---\n\n## " + block, encoding="utf-8")
        print("appended:", title)
print("最终:", sorted(set(re.findall(r"LRN-\d+-\d+", P.read_text(encoding="utf-8")))))
