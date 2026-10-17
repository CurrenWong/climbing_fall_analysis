"""定位 LEARNINGS.md 的权威副本（项目 vs workspace）+ 下一个 LRN 编号。"""
import os
import re

P = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
for tag, base in (("项目", P), ("workspace", W)):
    f = os.path.join(base, ".learnings", "LEARNINGS.md")
    if not os.path.exists(f):
        print(f"{tag}: (无) {f}")
        continue
    t = open(f, encoding="utf-8", errors="replace").read()
    ids = sorted(set(re.findall(r"\[LRN-20261005-(\d+)\]", t)))
    print(f"{tag}: {f}\n    {os.path.getsize(f)} bytes  20261005 ids={ids}")
    print(f"    064 已用? {'[LRN-20261005-064]' in t}   062 已用? {'[LRN-20261005-062]' in t}")
for tag, base in (("项目", P), ("workspace", W)):
    f = os.path.join(base, "MEMORY.md")
    if os.path.exists(f):
        print(f"{tag} MEMORY.md: {os.path.getsize(f)} bytes")
