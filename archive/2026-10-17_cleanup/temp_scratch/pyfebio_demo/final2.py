"""按行插入两条限制 + 复查残留 + 全量校验。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
D = os.path.join(R, "docs", "踝FE可用模型使用说明.md")

ls = io.open(D, encoding="utf-8").read().splitlines()
NEW = [
    '| **内外翻不对称**（`[LRN-20261006-076]`）| `POSE_IE=+5`(内翻) `NORMAL`；'
    '**`−5`(外翻) `rc=1`/1 个负 J / t=0 即失败** | 外翻侧落位平移需 **+0.186 mm**'
    '（内翻只需 −0.092）⇒ 该侧几何更"顶"。**⇒ 目前只有内翻侧可用**，"判姿势"覆盖范围受限 |',
    '| **韧带插入点几何** | ATFL 建模 **36.99 mm** vs 解剖 15~20（**≈2× 过长**）；'
    'CFL 34.33 vs ~28.5 | 同量绝对伸长被稀释 ⇒ **(c) 的绝对值暂不可当定量指标**；这是下一项要修的 |',
]
idx = [i for i, l in enumerate(ls) if "59,590" in l]
if not idx:
    print("⚠ 未找到锚点行")
elif any("内外翻不对称" in l for l in ls):
    print("· 已存在这两行")
else:
    i = idx[0]
    ls[i + 1:i + 1] = NEW
    io.open(D, "w", encoding="utf-8").write("\n".join(ls) + "\n")
    print(f"✓ 已插入 2 行 -> {len(io.open(D, encoding='utf-8').read().encode('utf-8'))} B")

# ── 全量校验 ──
print("\n" + "=" * 70)
print("压缩/交接校验")
print("=" * 70)
ok = True
files = [
    (os.path.join(W, "MEMORY.md"), 10240),
    (os.path.join(W, "references", "ankle-fe-criteria.md"), None),
    (os.path.join(W, "references", "fe-modeling-gotchas.md"), None),
    (D, None),
    (os.path.join(R, "ANKLE_FIX_TODO.md"), None),
    (os.path.join(R, ".learnings", "LEARNINGS.md"), None),
    (os.path.join(R, ".learnings", "ERRORS.md"), None),
    (os.path.join(R, "scripts", "ankle_fe", "check_c2.py"), None),
    (os.path.join(R, "scripts", "ankle_fe", "lig_strain.py"), None),
    (os.path.join(R, "scripts", "ankle_fe", "ankle_lig_feb.py"), None),
]
for p, cap in files:
    if not os.path.exists(p):
        print(f"  ✗ 缺失 {p}")
        ok = False
        continue
    n = os.path.getsize(p)
    flag = ""
    if cap:
        flag = " ✅" if n <= cap else f" ❌ 超 {n-cap}"
        if n > cap:
            ok = False
    print(f"  {os.path.relpath(p, R) if p.startswith(R) else p:<58} {n:>8} B{flag}")

# AST
import ast
for s in ("check_c2.py", "lig_strain.py", "ankle_lig_feb.py", "acceptance.py"):
    try:
        ast.parse(io.open(os.path.join(R, "scripts", "ankle_fe", s),
                          encoding="utf-8").read())
        print(f"  AST ✓ {s}")
    except SyntaxError as e:
        print(f"  AST ✗ {s}: {e}")
        ok = False

# LRN/ERR 重号
for f in (".learnings/LEARNINGS.md", ".learnings/ERRORS.md"):
    t = io.open(os.path.join(R, f), encoding="utf-8").read()
    ids = re.findall(r"^#+\s*\[?([A-Z]{3}-\d{8}-\d{3})", t, re.M)
    dup = sorted({i for i in ids if ids.count(i) > 1})
    print(f"  {f}: {len(ids)} 条，重号 {dup or '无 ✓'}")
    if dup:
        ok = False

# 残留错说法（用真实字符）
print("\n--- 错说法残留（用真实字符复查）---")
pats = ["实际是外翻", "才是内翻", "是反的"]
for f in (os.path.join(W, "MEMORY.md"), D, os.path.join(R, "ANKLE_FIX_TODO.md")):
    t = io.open(f, encoding="utf-8").read()
    hits = []
    for p in pats:
        for m in re.finditer(re.escape(p), t):
            seg = t[max(0, m.start() - 40):m.end() + 30].replace("\n", " ")
            hits.append(f"{p} :: ...{seg}...")
    nm = os.path.basename(f)
    print(f"  {nm}: {len(hits)} 处")
    for h in hits:
        print(f"     {h[:150]}")

print("\n结论：", "PASS ✅" if ok else "有问题 ⚠")
