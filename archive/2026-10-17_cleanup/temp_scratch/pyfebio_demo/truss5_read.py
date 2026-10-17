"""读 truss 模型的支反力：给 pyfebio 打容错补丁（beam 域没有 name 字段）。"""
import subprocess
from pathlib import Path
import numpy as np
import h5py
import pyfebio.xplt as X

# ── 容错补丁：parse_mesh 遇到没有 "name" 的域时用它自己的 id 兜底 ──
_orig = X.parse_mesh
import inspect
SRC = inspect.getsource(_orig)
print("原 parse_mesh 关键行:")
for l in SRC.split("\n"):
    if '"name"' in l or "domains" in l:
        print("   ", l.strip()[:110])

WORK = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\truss_probe")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
E, L, DL = 260.0, 10.0, 1.0

print("\n=== 支反力（v2.0 xplt 文本块直接扫） ===")
for tag in ("et_iso",):
    xp = WORK / f"{tag}.xplt"
    raw = xp.read_bytes()
    print(f"\n{tag}: {len(raw):,} bytes")
    txt = raw.decode("latin-1")
    # 找 REACTION / FACE 等块名
    import re
    for m in re.finditer(r"(reaction|RFT|FORCE|DISPLACEMENT|NODAL)", txt, re.I):
        i = m.start()
        seg = txt[max(0, i - 40):i + 120].replace("\n", "|")
        print(f"  @{i}: {seg[:150]}")
