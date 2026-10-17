"""压缩前最终校验：AST / 死链 / 真实重复条目 / 尺寸 / 关键内容。"""
import ast
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
os.chdir(R)
ok = True

# 1. AST + 死代码
ast.parse(io.open("scripts/ankle_fe/ankle_lig_feb.py", encoding="utf-8").read())
print("1) AST OK")
import subprocess
p = subprocess.run([".venv\\Scripts\\python.exe", "-u",
                    "temp\\pyfebio_demo\\scan_dead.py",
                    "scripts\\ankle_fe\\ankle_lig_feb.py"],
                   capture_output=True, text=True, timeout=300)
print("   " + (p.stdout or "").strip().replace("\n", " | "))

# 2. 真实重复条目（只看 ## [XXX 标题）
for f in (".learnings/LEARNINGS.md", ".learnings/ERRORS.md"):
    t = io.open(f, encoding="utf-8").read()
    ids = re.findall(r"^#+\s*\[?([A-Z]{3}-\d{8}-\d{3})", t, re.M)
    dup = sorted({i for i in ids if ids.count(i) > 1})
    print(f"2) {f}: 条目 {len(ids)} 条，真实重号 {dup or '无'}")
    if dup:
        ok = False

# 3. 死链：MEMORY.md 里 references/ 与 memory/ 的指针是否存在
M = os.path.join(W, "MEMORY.md")
ms = io.open(M, encoding="utf-8").read()
refs = set(re.findall(r"references/[A-Za-z0-9_\-]+\.md", ms))
refs |= set(re.findall(r"memory/\d{4}-\d{2}-\d{2}[^\s）\)|]*", ms))
miss = [r for r in sorted(refs) if not os.path.exists(os.path.join(W, r))]
print(f"3) MEMORY 指针 {len(refs)} 个，缺失 {miss or '无'}")

# 4. 关键内容存在性
D = os.path.join(R, "docs", "踝FE可用模型使用说明.md")
ds = io.open(D, encoding="utf-8").read()
need = ["LOAD_SIGN=+1", "−608.82", "POSE_IE", "79.71", "ERR-20261005-009",
        "ERR-20261005-010", "ATFL"]
print("4) 使用说明关键内容：" + "，".join(
    f"{k}{'✓' if k in ds else '✗'}" for k in need))
B = os.path.join(R, "ANKLE_FIX_TODO.md")
print("   ANKLE_FIX_TODO 横幅：" + ("✓" if "重大更正（本文件下方" in
                                io.open(B, encoding="utf-8").read() else "✗"))

# 5. 尺寸
print("5) 尺寸：")
for f in ("scripts/ankle_fe/ankle_lig_feb.py", "ANKLE_FIX_TODO.md",
          "docs/踝FE可用模型使用说明.md", ".learnings/LEARNINGS.md",
          ".learnings/ERRORS.md"):
    print(f"   {f:<40} {os.path.getsize(os.path.join(R, f)):>8} B")
print(f"   MEMORY.md {len(ms.encode('utf-8'))} / 10240 "
      f"{'✅' if len(ms.encode('utf-8')) <= 10240 else '❌'}")
print(f"   references/fe-modeling-gotchas.md "
      f"{os.path.getsize(os.path.join(W, 'references', 'fe-modeling-gotchas.md'))} B")
print("\n全部检查：", "PASS ✅" if ok else "有问题 ⚠")
