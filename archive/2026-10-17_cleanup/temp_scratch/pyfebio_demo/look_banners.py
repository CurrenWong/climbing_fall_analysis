"""看清 3 个文件里新旧横幅的分界，然后精确删除旧的。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
for rel in ("ANKLE_GOAL_TODO.md", "ANKLE_600N_TODO.md", "TIE_PORT_TODO.md"):
    p = os.path.join(R, rel)
    ls = io.open(p, encoding="utf-8").read().splitlines()
    print("=" * 90)
    print(f"=== {rel}  共 {len(ls)} 行")
    for i, l in enumerate(ls[:46], start=1):
        mark = ""
        if "2026-10-05 重大更正" in l:
            mark = "   ← 旧横幅开始"
        if "2026-10-06" in l and "翻转史" in l:
            mark = "   ← 新横幅"
        if "翻转史（写下来防止再翻）" in l:
            mark = "   ← 新横幅的'翻转史'节"
        print(f"  {i:>3}: {l[:118]}{mark}")
