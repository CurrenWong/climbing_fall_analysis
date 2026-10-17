"""修正：MEMORY.md 裁到磁盘≤10240；FE_PIPELINE_TODO.md 在正确锚点插入。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
L = 10240
d = lambda p: p.stat().st_size  # noqa: E731

print("--- FE_PIPELINE 入口附近 ---")
F = ROOT / "FE_PIPELINE_TODO.md"
ft = F.read_text(encoding="utf-8")
lines = ft.split("\n")
for i, ln in enumerate(lines[18:48], start=19):
    print(f"{i:3d}| {ln[:110]}")
