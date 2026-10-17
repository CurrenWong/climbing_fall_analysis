"""审计：会话开头加的 7 个横幅里，关于"符号/拉伸"的论断现在是否已错。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
FILES = ["ANKLE_FIX_TODO.md", "ANKLE_GOAL_TODO.md", "ANKLE_600N_TODO.md",
         "TIE_PORT_TODO.md", "ANKLE_LIG_TODO.md", "FE_PIPELINE_TODO.md",
         "docs/总体进度.md"]
PAT = re.compile(r"(拉伸|拉开|DISP_MM<0|DISP_MM>0|LOAD_SIGN|已作废|符号)")
for rel in FILES:
    p = os.path.join(R, rel.replace("/", os.sep))
    if not os.path.exists(p):
        print(f"✗ {rel}")
        continue
    t = io.open(p, encoding="utf-8").read()
    ls = t.splitlines()
    hits = [(i + 1, l.strip()) for i, l in enumerate(ls)
            if PAT.search(l) and i < 30]
    print("=" * 96)
    print(f"=== {rel}   （前 30 行里的相关论断）")
    for i, l in hits:
        print(f"  {i:>3}: {l[:150]}")
