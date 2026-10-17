"""把 .feb 里的 elem id 区间映射到域/单元类型（用于解读 FEBio 的负 Jacobian 元件号）。"""
import re
from pathlib import Path

feb = Path(r"D:\Project\climbing_fall_analysis\temp/thums/febio_lig/ankle_lig_WK.feb")
t = feb.read_text(encoding="latin-1")
blocks = []
for m in re.finditer(r'<Elements type="(\w+)" name="([^"]+)">(.*?)</Elements>',
                     t, re.S):
    et, name, body = m.groups()
    ids = [int(x) for x in re.findall(r'<elem id="(\d+)">', body)]
    blocks.append((name, et, min(ids), max(ids), len(ids)))
for name, et, a, b, n in blocks:
    print(f"  {a:>6}..{b:<6} n={n:<5} {et:7s} {name}")
