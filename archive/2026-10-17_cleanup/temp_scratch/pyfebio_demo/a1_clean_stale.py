"""删掉参考工况的旧产物（temp/ 下，避免 FEBio/acceptance 读到旧文件出假象）。"""
import os
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
targets = [
    ROOT / "temp/thums/febio_lig/ankle_lig_WK.feb",
    ROOT / "temp/thums/febio_lig/log_WK.txt",
    ROOT / "temp/thums/febio_lig/disp_WK.txt",
    ROOT / "temp/thums/febio_lig/reac_tibia_WK.txt",
    ROOT / "temp/thums/febio_lig/reac_fibula_WK.txt",
    ROOT / "temp/thums/febio_lig/reac_talus_WK.txt",
    ROOT / "temp/thums/febio_lig/lig_read_WK.npz",
    ROOT / "temp/thums/febio_lig/febio_WK.log",
    ROOT / "temp/pyfebio_demo/_acc_last.log",
]
for t in targets:
    if t.exists():
        n = t.stat().st_size
        os.remove(t)
        print(f"  删 {t.name} ({n:,} B)")
    else:
        print(f"  跳过（本就不存在）{t.name}")
print("done")
