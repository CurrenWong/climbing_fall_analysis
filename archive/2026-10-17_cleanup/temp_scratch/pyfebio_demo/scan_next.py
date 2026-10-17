"""扫出两文档里与 C2/C3/判据/下一步 相关的行。"""
import io
import re

for f in (r"docs\踝FE可用模型使用说明.md", r"ANKLE_FIX_TODO.md"):
    t = io.open(f, encoding="utf-8").read()
    print("=" * 24, f, len(t.encode("utf-8")), "B")
    for i, l in enumerate(t.splitlines()):
        if re.search(r"C2|C3|判据|下一步|待办|重定义", l):
            print(f"{i+1:>5}: {l.strip()[:132]}")
