"""最后一轮精简（差 119 字节）。"""
import io

M = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
s = io.open(M, encoding="utf-8").read()
n0 = len(s.encode("utf-8"))
PAIRS = [
    ("""- ⚠️ **加载方向不可靠**（`[LRN-060]`）：`dispx` 只按 `argmax(|轴|)` 全局单轴 + `TRUNC_MM=60` 转轴 **12.18°**
  ⇒ k 跳 **3.2×**（290/416/916）。**正解=`ALIGN_AXIS=1`**""",
     """- ⚠️ 加载方向须 `ALIGN_AXIS=1`（否则 `argmax(|轴|)` 单轴 + `TRUNC_MM` 转轴 12.18° ⇒ k 跳 3.2×）`[LRN-060]`"""),
    ("""- 🧨 **FEBio 静默坑**（`lc`→`<scale>`｜`symmetric_stiffness` 在 `<solver>`｜跑前删旧 `.feb`+查 rc｜反力看符号）
  → `references/fe-modeling-gotchas.md`；`[ERR-20261005-008/009/010]`""",
     """- 🧨 **FEBio 静默坑**：`lc`→`<scale>`｜`sym_stiffness` 在 `<solver>`｜删旧 `.feb`+查 rc｜反力看符号
  → `references/fe-modeling-gotchas.md`｜`[ERR-…008/009/010]`"""),
]
for old, new in PAIRS:
    c = s.count(old)
    if c != 1:
        print(f"⚠️ 未命中({c})：{old.splitlines()[0][:50]}")
        continue
    s = s.replace(old, new)
    print(f"✓ {old.splitlines()[0][:46]}")
io.open(M, "w", encoding="utf-8").write(s)
n1 = len(s.encode("utf-8"))
print(f"MEMORY.md {n0} -> {n1} / 10240  余量 {10240 - n1}  "
      f"{'✅ 合规' if n1 <= 10240 else '❌ 仍超'}")
