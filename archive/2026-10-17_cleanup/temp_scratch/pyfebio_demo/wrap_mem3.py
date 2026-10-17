"""收口①c：精修 MEMORY 至 ≤10240。"""
import io
import os

W = r"C:\Users\curre\.copaw\workspaces\default"
P = os.path.join(W, "MEMORY.md")
t = io.open(P, encoding="utf-8").read()

REPL = [
    ("- 入口：**`ANKLE_FIX_TODO.md` §第十三~十六轮**（最新全量状态，自包含）｜\n"
     "  运行说明 `docs/踝FE可用模型使用说明.md`｜判据 `references/ankle-fe-criteria.md`",
     "- 入口：**`ANKLE_FIX_TODO.md` §第十三~十六轮**（自包含）｜"
     "运行说明 `docs/踝FE可用模型使用说明.md`"),
    ("  （npz 里 `subtalar_*` 软骨**已生成但未接进 `.feb`**）② **无摩擦接触切向滑移 39.5%**\n"
     "  （关节间只 0.029 mm 真压扁）③ **tie 罚滑移 13.9%**。骨合计 <6%、弯曲 <0.2%",
     "  （`subtalar_*` 软骨**已生成、未接进 `.feb`**）② **无摩擦切向滑移 39.5%**"
     "（关节间只 0.029mm 真压扁）③ **tie 罚滑移 13.9%**；骨 <6%"),
    ("- 🔬 **软骨网格大修**：① `area_proj` **quad4 翘曲 bug 已修**（47.8→**555.5** mm²；力控除数错 11.6×）\n"
     "  ② **11 个顶点翻转单元、全在关节面 ⇒ 不可剔除**；`detJ_penta6_full` **只诊断**（`[LRN-087]`）\n"
     "  ③ `cart_m1_fix` 把 ≤0 降到 **0** 但**贴合度回归**（+5° 张开 0.598 mm）⇒ **保持 `cart_m1` 为基线**",
     "- 🔬 **软骨网格大修**：① `area_proj` **quad4 翘曲 bug 已修**（47.8→**555.5** mm²）"
     "② **11 个顶点翻转单元全在关节面⇒不可剔除**；`detJ_penta6_full` **只诊断**（`[LRN-087]`）"
     "③ `cart_m1_fix` ≤0 降到 **0** 但**贴合度回归**（+5° 张开 0.598mm）⇒ **保持 `cart_m1` 为基线**"),
    ("- 🧨 **坑**：`lc`→`<scale>`｜`sym_stiffness` 在 `<solver>`｜删旧 `.feb`+查 rc｜\n"
     "  **patch test 必须用畸变单元**（规则单元 `J=cI` 假通过）｜**h5 nodesets ≠ `.feb` NodeSet**\n"
     "  （差 46/35/80/59 节点）→ `references/fe-modeling-gotchas.md`",
     "- 🧨 **坑**：`lc`→`<scale>`｜`sym_stiffness` 在 `<solver>`｜删旧 `.feb`+查 rc｜"
     "**patch test 须用畸变单元**（规则单元 `J=cI` 假通过）｜"
     "**h5 nodesets ≠ `.feb` NodeSet**（差 46/35/80/59 节点）"),
    ("- **Anderson 真值**（`[LRN-067]` 核过原文）：all-hex｜距骨**弱弹簧 500 N/mm 接地**+解除内外翻内旋\n"
     "  ｜**骨刚性+reference node**｜不建韧带｜E=12 ν=0.42 厚 1.5~1.7mm｜落位 200 N｜精度 10~15%",
     "- **Anderson 真值**（`[LRN-067]`）：all-hex｜距骨**弱弹簧 500 N/mm**+解除内外翻内旋｜"
     "**骨刚性+reference node**｜不建韧带｜E=12 ν=0.42 厚 1.5~1.7mm｜落位 200 N"),
]
n = 0
for a, b in REPL:
    if a in t:
        t = t.replace(a, b)
        n += 1
    else:
        print(f"  ⚠ 未命中：{a[:50]}...")
io.open(P, "w", encoding="utf-8").write(t)
n1 = len(t.encode("utf-8"))
print(f"替换 {n}/{len(REPL)} 处 ⇒ {n1} B  "
      + ("✅ 余量 %d" % (10240-n1) if n1 <= 10240 else "❌ 仍超 %d" % (n1-10240)))
