"""查历史是否记录过双靶结果（3.74 MPa / 335.5 mm²）。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
PAT = re.compile(r"3\.74|335\.5|双靶")
SKIP = (".venv", ".git", "__pycache__", "node_modules", "paper", "model")
files = []
for root in (R, W):
    for dp, dirs, fs in os.walk(root):
        dirs[:] = [x for x in dirs if x not in SKIP]
        for f in fs:
            if f.endswith((".md", ".txt", ".log")):
                p = os.path.join(dp, f)
                try:
                    if os.path.getsize(p) > 4_000_000:
                        continue
                    t = io.open(p, encoding="utf-8", errors="ignore").read()
                except OSError:
                    continue
                if PAT.search(t):
                    files.append((os.path.relpath(p, root), t))
print(f"命中 {len(files)} 个文件\n")
for rel, t in files:
    print("=" * 90)
    print(f"=== {rel}")
    for i, l in enumerate(t.splitlines()):
        if PAT.search(l):
            print(f"  {i+1:>5}: {l.strip()[:150]}")
