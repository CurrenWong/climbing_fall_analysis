"""安全检查：被 subagent 改过的 cart_patch.py 的结构（不运行它）。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
P = os.path.join(R, "scripts", "ankle_fe", "cart_patch.py")
t = io.open(P, encoding="utf-8").read()
ls = t.splitlines()
print("=" * 96)
print(f"cart_patch.py  {len(t.encode('utf-8'))} B / {len(ls)} 行")
print("=" * 96)
print("--- 顶部 docstring（前 30 行）---")
for l in ls[:30]:
    print("  " + l[:130])

print("\n--- 所有 def / class ---")
for i, l in enumerate(ls):
    if re.match(r"^\s*(def |class )", l):
        print(f"  {i+1:>5}: {l.strip()[:120]}")

print("\n--- 输出路径 / env 开关（关键：会不会覆盖 cart_m1）---")
for i, l in enumerate(ls):
    if re.search(r"cart_m1|OUT|_DIR|environ|savez|argparse|sys\.argv", l):
        print(f"  {i+1:>5}: {l.strip()[:130]}")
