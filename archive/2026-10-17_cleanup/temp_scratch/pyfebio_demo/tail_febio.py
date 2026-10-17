"""看 step2 停滞增量的逐迭代信息（残余/比率）。"""
import io
import re
from pathlib import Path

t = io.open(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_ms_febio.txt",
            encoding="utf-8", errors="replace").read()
lines = t.splitlines()
i = max(k for k, l in enumerate(lines) if "beginning time step" in l
        or "failed to converge" in l)
print("=== 最后 70 行 ===")
print("\n".join(l for l in lines[-70:] if l.strip()))
