"""找 §4 逐韧带应变数字（3.824/2.842/3.513）的来源文件。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
NEEDLES = ("3.824", "2.842", "3.513")
SKIP = (".venv", ".git", "node_modules", "__pycache__", "paper")
found = []
for root, dirs, files in os.walk(R):
    dirs[:] = [d for d in dirs if d not in SKIP]
    for f in files:
        p = os.path.join(root, f)
        if os.path.getsize(p) > 8_000_000:
            continue
        try:
            t = io.open(p, encoding="utf-8", errors="ignore").read()
        except OSError:
            continue
        hit = [n for n in NEEDLES if n in t]
        if hit:
            found.append((os.path.relpath(p, R), hit, len(t)))
for rel, hit, n in sorted(found):
    print(f"{rel:<62} {hit}  {n} 字符")
print(f"\n共 {len(found)} 个文件")
