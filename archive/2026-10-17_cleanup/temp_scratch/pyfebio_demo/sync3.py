from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
L = 10240
s = M.read_text(encoding="utf-8")
CUTS = [
    ("- **入口**：`FE_PIPELINE_TODO.md`「🧭接手」｜`ANKLE_FIX_TODO.md`｜"
     "`docs/Anderson双论文_关键信息.md`",
     "- **入口**：`FE_PIPELINE_TODO.md`「🧭接手」｜`ANKLE_FIX_TODO.md`｜"
     "`docs/Anderson双论文_关键信息.md`"),
    ("**参考原文**（搬运/转载不「创造性总结」，早报除外）",
     "**保留原文**（搬运/转载不「创造性总结」，早报除外）"),
    ("**搜索硬 3-5 条**", "**搜索硬 3-5 条**"),
    ("| 云函数冷启动 | 云函数改动必须重新 `updateFunctionCode` **+二次 invoke 验证** |",
     "| 云函数冷启动 | 云函数改动须重 `updateFunctionCode` **+二次 invoke 验证** |"),
    ("- **emoji 用在合适处**", "- **emoji 用在合适处**"),
    ("**宁缺毋滥**（搜索硬 3-5 条）", "**宁缺毋滥**（搜索硬 3-5 条）"),
]
for a, b in CUTS:
    if M.stat().st_size <= L:
        break
    if a in s and a != b:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
        print(f"  -> {M.stat().st_size:,} 余量 {L - M.stat().st_size:+,}")
sz = M.stat().st_size
print(f"disk = {sz:,} / {L}  {'OK' if sz <= L else '仍超 ' + str(sz - L)}")
