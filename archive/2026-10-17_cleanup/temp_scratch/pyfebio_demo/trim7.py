"""MEMORY 精修：把「坑」条压成指针，腾出 200+ B。"""
import io

P = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
t = io.open(P, encoding="utf-8").read()
print(f"改前 {len(t.encode('utf-8'))} B")
A = ("- 🧨 **坑**：`lc`→`<scale>`｜`sym_stiffness` 在 `<solver>`｜删旧 `.feb`+查 rc｜"
     "**patch test 须用畸变单元**｜**h5 nodesets ≠ `.feb` NodeSet**（差 46/35/80/59）→ "
     "`references/fe-modeling-gotchas.md`")
B = ("- 🧨 **FEBio 坑**（`lc`→`<scale>`／`sym_stiffness` 在 `<solver>`／删旧 `.feb` 查 rc／"
     "patch test 须用畸变单元／**h5 nodesets ≠ `.feb` NodeSet**、差 46/35/80/59）"
     "→ `references/fe-modeling-gotchas.md`")
C = ("- 🧨 **FEBio 坑** → `references/fe-modeling-gotchas.md`"
     "（`lc`→`<scale>`／`sym_stiffness` 在 `<solver>`／删旧 `.feb` 查 rc／"
     "patch test 用畸变单元／**h5 nodesets ≠ `.feb` NodeSet**）")
for src in (A, B):
    if src in t:
        t = t.replace(src, C); print("✓ 已压缩"); break
else:
    print("⚠ 未命中，打印含「坑」的行：")
    for ln in t.split("\n"):
        if "🧨" in ln:
            print("   ", ln[:200])
io.open(P, "w", encoding="utf-8").write(t)
n = len(t.encode("utf-8"))
print(f"改后 {n} B  " + ("✅ 余量 %d" % (10240-n) if n <= 10240 else "❌ 超 %d" % (n-10240)))
