"""修 3 处手术副作用：#21 截断 / 表中间空行 / gmsh 重复短语。"""
import io

p = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
t = io.open(p, encoding="utf-8").read()
b = len(t.encode("utf-8"))
fixes = []


def sub(old, new, tag):
    global t
    assert t.count(old) == 1, f"[{tag}] 命中 {t.count(old)}"
    t = t.replace(old, new)
    fixes.append(tag)


# ① gotcha 21 被截断的尾巴补回
sub('④ **做优化前先确认它是不是瓶颈**（我用"提高网格质量"|',
    '④ **优化前先确认它是不是瓶颈**（我曾把"网格质量"当瓶颈，结果质量提升 42 倍反而更差）|',
    "g21-tail")

# ② SecLog 表被空行截断（Markdown 里会断成两张表）
sub("只返空 URL |\n\n| 金额显示", "只返空 URL |\n| 金额显示", "table-blank")

# ③ gmsh 那条重复了「occ 绑定残缺」
sub("- **`uv pip install gmsh` ≠ 官方 gmsh SDK**：pip wheel 的 occ 绑定残缺\n"
    "  （occ 绑定残缺 ⇒ STL→实体走不通）；官方 SDK 在 `D:\\Program\\gmsh-sdk\\`，\n"
    "  DLL 在 `lib/` 不在 `bin/`。",
    "- **`uv pip install gmsh` ≠ 官方 SDK**：pip wheel 的 occ 绑定残缺 ⇒ STL→实体走不通；\n"
    "  官方 SDK 在 `D:\\Program\\gmsh-sdk\\`，DLL 在 `lib/` 不在 `bin/`。", "gmsh-dup")

# ④ 「保留原文」那条也顺手并进上一行，省一行
sub("- **中文优先**（专有名词保留英文）；**先行动后报告**（能跑就跑）；\n"
    "  **宁缺毋滥**（搜索硬 3-5 条）\n"
    "- **保留原文**：搬运/转载类任务不做「创造性总结」（每日早报除外）",
    "- **中文优先**（专有名词保留英文）；**先行动后报告**（能跑就跑）；\n"
    "  **宁缺毋滥**（搜索硬 3-5 条）；**保留原文**（搬运/转载不做「创造性总结」，每日早报除外）",
    "style-merge")

io.open(p, "w", encoding="utf-8").write(t)
a = len(t.encode("utf-8"))
print(f"修复 {len(fixes)} 处: {fixes}")
print(f"{b} -> {a} bytes  "
      f"{'OK, 余量 ' + str(10240 - a) if a <= 10240 else 'FAIL 差 ' + str(a - 10240)}")
