"""再裁 MEMORY.md 到 ≤10240。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
LIMIT = 10240
s = M.read_text(encoding="utf-8")

CUTS = [
    ("- ⭐ **根因取向=接触区欠分辨**：仅 **73 节点**(~2.5mm)；加密到 692 后\n"
     "  软骨体积涨 **2.3×** ⇒ 粗网格低估面积",
     "- ⭐ **根因取向=接触区欠分辨**：仅 **73 节点**(~2.5mm)；加密到 692 后软骨体积\n"
     "  涨 **2.3×** ⇒ 粗网格低估面积"),
    ("- ⚠️ **刚度 2 倍不确定度**：tet4 **259.3** vs penta6 **116.9** N；\n"
     "  **材料已排除**（比值恒 2.22~2.35）⇒ 剩离散/接触",
     "- ⚠️ **刚度 2 倍**：tet4 **259.3** vs penta6 **116.9** N；材料已排除\n"
     "  （比值恒 2.22~2.35）⇒ 剩离散/接触"),
    ("| 5 | **cmd `set X=Y && Z` 末尾追加 space**（40013）",
     "| 5 | **cmd `set X=Y && Z` 末尾追加 space**（40013）"),
]
for a, b in CUTS:
    if a != b and a in s:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
        print(f"  -> {M.stat().st_size:,} (余量 {LIMIT - M.stat().st_size:+,})")

# 兜底：逐轮砍最短收益项
EXTRA = [
    ("**160 倍**", "160x"), ("**1.5~1.7mm**", "1.5~1.7mm"),
    ("`memory/2026-10-03/ankle-fe-modeling.md` §六–§十",
     "`memory/2026-10-03/ankle-fe-modeling.md`"),
    ("**与实验一致**；吴恺 1.92~5.24mm 隐含 8~21% **与实验不符**",
     "**与实验一致**；吴恺 1.92~5.24mm 隐含 8~21% **不符实验**"),
    ("**真值**=`docs/Anderson双论文_关键信息.md`",
     "真值=`docs/Anderson双论文_关键信息.md`"),
]
i = 0
while M.stat().st_size > LIMIT and i < len(EXTRA):
    a, b = EXTRA[i]
    if a in s:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
        print(f"  extra{i} -> {M.stat().st_size:,}")
    i += 1
sz = M.stat().st_size
print(f"\n最终 disk = {sz:,} / {LIMIT}  {'OK' if sz <= LIMIT else '超限 ' + str(sz - LIMIT)}")
