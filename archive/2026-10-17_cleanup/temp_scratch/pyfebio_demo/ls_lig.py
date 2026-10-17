"""列 febio_lig 目录，找 FEBio 自己的 .log（可能有单元号）。"""
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
for p in sorted(W.iterdir(), key=lambda q: -q.stat().st_mtime)[:20]:
    print(f"  {p.stat().st_size:>12,}  {p.name}")
