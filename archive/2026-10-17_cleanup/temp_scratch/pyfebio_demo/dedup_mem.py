from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
L = 10240
s = M.read_text(encoding="utf-8")
DUP = "- **详见**：`memory/2026-10-03/ankle-fe-modeling.md`"
n = s.count(DUP)
print(f"重复行出现 {n} 次")
if n > 1:
    # 保留第一次，删掉其余
    first = s.find(DUP)
    head, tail = s[:first + len(DUP)], s[first + len(DUP):]
    tail = tail.replace(DUP, "").replace("\n\n\n", "\n\n")
    s = head + tail
    M.write_text(s, encoding="utf-8")
print(f"disk = {M.stat().st_size:,} / {L}  余量 {L - M.stat().st_size:+,}")
