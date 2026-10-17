"""MEMORY 精修三：最后 70+ B。"""
import io

P = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
t = io.open(P, encoding="utf-8").read()
print(f"改前 {len(t.encode('utf-8'))} B")
REPL = [
    ("- ⚠️ 加载须 `ALIGN_AXIS=1`（否则 `TRUNC_MM` 转轴 12.18° ⇒ k 跳 3.2×）`[LRN-060]`",
     "- ⚠️ 加载须 `ALIGN_AXIS=1`（否则 `TRUNC_MM` 转轴 12.18°⇒k 跳 3.2×）`[LRN-060]`"),
    ("- 🔍 **2~4× 的来源（取代旧「刚度 2× 未归因」）**：① ~~缺距下关节接触 44.8%~~ "
     "**已否证**（2026-10-06 直接测量：接入距下关节 ⇒ F 变化 **+0.0000%**，因 "
     "`FIX_TALUS=transverse` 冻结距骨 ⇒ 接触不承力）② 无摩擦切向滑移 **39.5%**"
     "（关节间只 0.029mm 真压扁）③ tie 罚滑移 **13.9%**；骨 <6%",
     "- 🔍 **2~4× 的来源**：① ~~缺距下关节 44.8%~~ **已否证**（2026-10-06 实测：接入后 F "
     "**+0.0000%**，因 `FIX_TALUS=transverse` 冻结距骨⇒接触不承力）② 无摩擦切向滑移 "
     "**39.5%**（只 0.029mm 真压扁）③ tie 罚滑移 **13.9%**；骨 <6%"),
    ("- 🧨 **坑**：`lc`→`<scale>`｜`sym_stiffness` 在 `<solver>`｜删旧 `.feb`+查 rc｜"
     "**patch test 须用畸变单元**（规则单元 `J=cI` 假通过）｜"
     "**h5 nodesets ≠ `.feb` NodeSet**（差 46/35/80/59 节点）",
     "- 🧨 **坑**：`lc`→`<scale>`｜`sym_stiffness` 在 `<solver>`｜删旧 `.feb`+查 rc｜"
     "**patch test 须用畸变单元**｜**h5 nodesets ≠ `.feb` NodeSet**（差 46/35/80/59）"),
]
for a, b in REPL:
    if a in t:
        t = t.replace(a, b)
    else:
        print("  ⚠ 未命中：", a[:34].replace("\n", "\\n"))
io.open(P, "w", encoding="utf-8").write(t)
n = len(t.encode("utf-8"))
print(f"改后 {n} B  " + ("✅ 余量 %d" % (10240-n) if n <= 10240 else "❌ 仍超 %d" % (n-10240)))
