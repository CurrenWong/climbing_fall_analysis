"""校验 MEMORY.md 尺寸与关键更正是否到位。"""
import io

M = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
s = io.open(M, encoding="utf-8").read()
n = len(s.encode("utf-8"))
print(f"MEMORY.md {n} / 10240  余量 {10240 - n}")
for k in ("LOAD_SIGN", "-608.82", "592.6", "79.7", "ERR-20261005-009",
          "ERR-20261005-010", "ERR-20261005-008", "已作废", "拉伸"):
    print(f"  {k:<18} x{s.count(k)}")
