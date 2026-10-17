"""打印指定日志文件的 FEBio 段。用法: python show_log.py <log 路径>"""
import sys
from pathlib import Path

p = Path(sys.argv[1] if len(sys.argv) > 1 else
         r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_FIX.log")
L = p.read_text(encoding="utf-8", errors="replace").splitlines()
i = [k for k, l in enumerate(L) if "FEBio rc=" in l]
k = i[-1] if i else 0
print(f"=== {p.name} 共 {len(L)} 行，FEBio 段从第 {k+1} 行 ===")
print("\n".join(L[max(0, k - 12):k + 30]))
