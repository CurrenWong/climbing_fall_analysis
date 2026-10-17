"""MEMORY 真正压缩：把 519N 条与坑条都削短。"""
import io

P = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
t = io.open(P, encoding="utf-8").read()
print(f"改前 {len(t.encode('utf-8'))} B")
REPL = [
    # 1) 519N 条削短
    ("- 🧭 **「519 N」是收敛代价、非承载极限**（2026-10-07 同位移实测）：位移驱动 ⇒ 它只是"
     "最后一个能解出的点；**A2 后距骨转动 0.52°（基线 0.0033°，158×）⇒ 胫距软骨由均匀受压"
     "变边缘受压 ⇒ 薄楔单元(1.2mm)剪翻**。基线(锁死)可达 **615.4N@+2.5mm** ⇒ 挪向物理正确"
     "掉 ~16% 余量。⚠️ 勿写成「踝只能承 500N」",
     "- 🧭 **「519 N」是收敛代价、非承载极限**：A2 后距骨转动 **0.52°**（基线 0.0033°）"
     "⇒ 软骨边缘受压 ⇒ 薄楔(1.2mm)剪翻；**基线(锁死)@+2.5mm 可达 615.4N**。"
     "⚠️ 勿写成「踝只能承 500N」"),
    # 2) 坑条压成短指针
    ("- 🧨 **FEBio 坑** → `references/fe-modeling-gotchas.md`（`lc`→`<scale>`／`sym_stiffness` 在 `<solver>`／删旧 `.feb` 查 rc／patch test 用畸变单元／**h5 nodesets ≠ `.feb` NodeSet**）",
     "- 🧨 **FEBio 坑**（`lc`/`sym_stiffness`/删旧 `.feb`/patch test/**h5 nodesets≠.feb NodeSet**）"
     "→ `references/fe-modeling-gotchas.md`"),
]
hit = 0
for a, b in REPL:
    if a in t:
        t = t.replace(a, b); hit += 1
    else:
        print("  ⚠ 未命中：", a[:40].replace("\n", "\\n"))
io.open(P, "w", encoding="utf-8").write(t)
n = len(t.encode("utf-8"))
print(f"命中 {hit}/{len(REPL)}；改后 {n} B  "
      + ("✅ 余量 %d" % (10240-n) if n <= 10240 else "❌ 超 %d" % (n-10240)))
