"""按 stat().st_size 精确收敛 MEMORY.md ≤ 10240。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"

raw = M.read_bytes()
print(f"disk = {len(raw):,} bytes")
print(f"CRLF = {raw.count(bytes([13, 10]))} | LF = {raw.count(bytes([10]))} | "
      f"CR = {raw.count(bytes([13]))}")
print(f"BOM = {raw[:3] == bytes([0xEF, 0xBB, 0xBF])}")
