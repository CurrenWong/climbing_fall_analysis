import io

L = io.open(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md",
            encoding="utf-8").read().split("\n")
for i, l in enumerate(L, 1):
    if (l.startswith("| 21 |") or l.startswith("| 只返空")
            or l.startswith("| 金额显示") or "occ 绑定残缺" in l
            or l.startswith("| 保留了") or "保留原文" in l):
        print(f"{i:3d}| {l[:170]}")
