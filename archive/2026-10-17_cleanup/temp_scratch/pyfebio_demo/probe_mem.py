"""查 MEMORY 尾部与总行数，定位多出来的内容。"""
import io

p = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
m = io.open(p, encoding="utf-8").read()
print(f"总 {len(m.encode('utf-8'))} B / {len(m.splitlines())} 行\n")
ls = m.splitlines()
for i, l in enumerate(ls[-46:], start=len(ls) - 46 + 1):
    print(f"{i:>4}: {l}")
