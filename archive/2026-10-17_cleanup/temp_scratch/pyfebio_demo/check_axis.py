"""直接查 lig_read_WK.npz 的 axis + 文件时间戳。"""
import os
import time
from pathlib import Path

import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
for f in ["lig_read_WK.npz", "ankle_lig_WK.feb", "log_WK.txt", "disp_WK.txt"]:
    p = W / f
    if p.exists():
        print(f"  {f:22} {p.stat().st_size:>10,} B   mtime={time.strftime('%H:%M:%S', time.localtime(p.stat().st_mtime))}")

p = W / "lig_read_WK.npz"
d = np.load(p)
print(f"\nnpz keys: {sorted(d.files)}")
ax = np.asarray(d["axis"], float)
print(f"  axis = {ax}   |axis|={np.linalg.norm(ax):.6f}")
print(f"  argmax|axis| = {int(np.argmax(np.abs(ax)))}  → 'xyz'[{int(np.argmax(np.abs(ax)))}]"
      f" = {'xyz'[int(np.argmax(np.abs(ax)))]}")
print(f"  与 [0.8153,0.135,0.5631] 的点积 = {float(ax @ np.array([0.8153,0.135,0.5631])):.6f}")
print(f"  与 [0.6766,0.1363,0.7237] 的点积 = {float(ax @ np.array([0.6766,0.1363,0.7237])):.6f}")

# .feb 里实际写了哪些 prescribed displacement（看 dof）
t = (W / "ankle_lig_WK.feb").read_text(encoding="latin-1", errors="replace")
import re
print("\n.feb 里的 prescribed displacement BC：")
for m in re.finditer(r'<bc name="([^"]+)" type="prescribed displacement"\s*node_set="([^"]+)">\s*<dof>(\w)</dof><value lc="1">([-\d.]+)</value>', t):
    print(f"  {m.group(1):16} ns={m.group(2):18} dof={m.group(3)}  value={m.group(4)}")
print("\n.feb 里的 zero displacement BC（dof 逐行）：")
for m in re.finditer(r'<bc name="([^"]+)" type="zero displacement" node_set="([^"]+)">(.*?)</bc>', t, re.S):
    dofs = re.findall(r"<(\w)_dof>", m.group(3))
    print(f"  {m.group(1):20} ns={m.group(2):22} dofs={dofs}")
