"""看 <Rigid>/接地弹簧的实际定义（查是否有预拉伸）。"""
import io
import re
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
f = (W / "ankle_lig_WK.feb").read_text(encoding="latin-1", errors="replace")
for tag in ("Rigid", "Discrete"):
    i, j = f.find(f"<{tag}>"), f.find(f"</{tag}>")
    if i < 0:
        print(f"<{tag}> 不存在")
        continue
    blk = f[i:j + len(tag) + 3]
    print(f"===== <{tag}> ({len(blk)} 字符) =====")
    print(blk[:1800])
    print("...")
