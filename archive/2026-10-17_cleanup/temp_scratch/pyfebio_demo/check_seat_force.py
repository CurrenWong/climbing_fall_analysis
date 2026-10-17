"""读 FEBio log：step1 末尾支反力到底多大（验证落位是否真生效）。"""
import os
import subprocess
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
e = dict(os.environ)
e.update(DRIVER="force", FORCE_N="600", FIX_TALUS="transverse", TRUNC_MM="60",
         AXIS_SRC="full", ALIGN_AXIS="1",
         CART_DIR=str(R / "temp" / "thums" / "cart_m1"), CART_NEW="1", TIE="1",
         DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0", MULTISTEP="1",
         SEAT_MM="0.5", TIME_STEPS="20", TIME_STEPS_SEAT="20", LOG_DATA="1")
subprocess.run([str(R / ".venv" / "Scripts" / "python.exe"), "-u",
                str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
               capture_output=True, text=True, env=e, cwd=str(R), timeout=3000)
subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "ms.xplt"],
               capture_output=True, text=True, cwd=str(W), timeout=3000)
print("log 目录文件：", sorted(p.name for p in W.glob("*.txt")))
for nm in ("reac_tibia", "reac_talus", "disp"):
    fs = sorted(W.glob(f"{nm}_*.txt"))
    if not fs:
        print(f"  {nm}: (无文件)")
        continue
    a = np.loadtxt(fs[0], delimiter=",", ndmin=2)
    print(f"\n{nm} <- {fs[0].name}  形状 {a.shape}（列: 时间 + 分量）")
    print(f"  首行 {np.round(a[0], 4)}")
    print(f"  末行 {np.round(a[-1], 4)}")
    if a.shape[1] >= 4:
        F = np.linalg.norm(a[:, 1:4], axis=1)
        print(f"  合力 |F| 前5 {np.round(F[:5], 3)}  后5 {np.round(F[-5:], 3)}")
