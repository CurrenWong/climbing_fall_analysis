"""查基准文献里的 ATFL / 韧带长度口径。"""
import glob
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
pats = ["距腓前", "ATFL", "腓前韧带", "韧带长度", "长度约", "LIG_TALOFIB"]
files = (glob.glob(os.path.join(R, "paper", "*.md"))
         + glob.glob(os.path.join(R, "paper", "_text", "*.txt"))
         + glob.glob(os.path.join(R, "docs", "*.md")))
for p in files:
    try:
        t = io.open(p, encoding="utf-8", errors="ignore").read()
    except OSError:
        continue
    hits = []
    for i, l in enumerate(t.splitlines()):
        if any(k in l for k in pats) and re.search(r"[\d.]+\s*(mm|毫米)|长度|length", l, re.I):
            hits.append((i + 1, l.strip()[:176]))
    if hits:
        print("=" * 22, os.path.basename(p), f"{len(hits)} 处")
        for i, l in hits[:12]:
            print(f"  {i:>5}: {l}")
