import io

L = io.open(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md",
            encoding="utf-8").read().split("\n")
print("L29 尾部:", repr(L[28][-150:]))
print("L29 结尾是 | ?", L[28].rstrip().endswith("|"))
