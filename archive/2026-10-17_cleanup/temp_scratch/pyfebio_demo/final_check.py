"""最终校验：AST / LRN 查重 / 文件尺寸。"""
import ast
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
os.chdir(R)
P = "scripts/ankle_fe/ankle_lig_feb.py"
ast.parse(io.open(P, encoding="utf-8").read())
print("AST OK")
L = io.open(".learnings/LEARNINGS.md", encoding="utf-8").read()
ids = re.findall(r"(LRN-\d{8}-\d{3})", L)
dup = sorted({i for i in ids if ids.count(i) > 1})
print("LRN 条数", len(ids), "| 重复号:", dup or "无")
print("LRN-20261005-070 出现", L.count("LRN-20261005-070"), "次")
src = io.open(P, encoding="utf-8").read()
print("MULTISTEP 出现", src.count("MULTISTEP"), "次")
for f in (P, "ANKLE_FIX_TODO.md", ".learnings/LEARNINGS.md",
          "docs/总体进度.md"):
    if os.path.exists(f):
        print(f"  {f:<34} {os.path.getsize(f):>8} B")
