"""找 THUMS 原始 .k 文件与任何部件名/文档，核对 8110xxxx 的权威名称。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
# 1) .k 文件
print("=== 大 .k / .key 文件 ===")
for root, dirs, files in os.walk(R):
    dirs[:] = [d for d in dirs if d not in (".venv", ".git", "node_modules")]
    for f in files:
        if f.lower().endswith((".k", ".key", ".keyword")):
            p = os.path.join(root, f)
            try:
                n = os.path.getsize(p)
            except OSError:
                continue
            if n > 500_000:
                print(f"  {n/1e6:>8.1f} MB  {os.path.relpath(p, R)}")

# 2) 任何含 8110 的元数据文件
print("\n=== 含 '8110' 的文本文件（前 200 KB）===")
cnt = 0
for root, dirs, files in os.walk(R):
    dirs[:] = [d for d in dirs if d not in (".venv", ".git", "node_modules",
                                            "paper", "model")]
    for f in files:
        p = os.path.join(root, f)
        try:
            if os.path.getsize(p) > 3_000_000:
                continue
            t = io.open(p, encoding="utf-8", errors="ignore").read(200_000)
        except OSError:
            continue
        if "81100701" in t or "81100801" in t or "81100901" in t:
            rel = os.path.relpath(p, R)
            print(f"  {rel}")
            cnt += 1
            if cnt > 25:
                break
    if cnt > 25:
        break

# 3) temp/thums 目录内容
print("\n=== temp/thums 前 40 项 ===")
d = os.path.join(R, "temp", "thums")
if os.path.isdir(d):
    for x in sorted(os.listdir(d))[:40]:
        p = os.path.join(d, x)
        sz = os.path.getsize(p) if os.path.isfile(p) else 0
        print(f"  {'DIR ' if os.path.isdir(p) else '    '}{x:<44} {sz:>12}")
