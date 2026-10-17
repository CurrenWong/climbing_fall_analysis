"""列出当前 .feb 的 NodeSet 名称 + 各集合规模，并给出解剖锚点可行性判断。"""
import io
import re
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
for p in sorted(W.glob("ankle_lig_*.feb")):
    t = p.read_text(encoding="latin-1", errors="replace")
    names = re.findall(r'<NodeSet name="([^"]+)">([\d,]*)</NodeSet>', t)
    print(f"\n=== {p.name}  {p.stat().st_size/1024:.0f} KB  集合 {len(names)} 个")
    for n, ids in names:
        k = len([x for x in ids.split(",") if x])
        print(f"   {n:<26} {k:>6} 节点")
    # 只要一个文件就够
    break
