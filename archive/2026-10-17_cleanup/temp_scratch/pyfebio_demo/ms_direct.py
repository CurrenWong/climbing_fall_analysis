"""直接跑 FEBio（不经 acceptance），抓 step2 失败的逐迭代细节。"""
import os
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
e = dict(os.environ)
e.update(DRIVER="force", FORCE_N="600", FIX_TALUS="transverse", TRUNC_MM="60",
         AXIS_SRC="full", ALIGN_AXIS="1",
         CART_DIR=str(R / "temp" / "thums" / "cart_m1"), CART_NEW="1", TIE="1",
         DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0", MULTISTEP="1",
         SEAT_MM="0.5")
subprocess.run([str(R / ".venv" / "Scripts" / "python.exe"), "-u",
                str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
               capture_output=True, text=True, env=e, cwd=str(R), timeout=3000)
feb = W / "ankle_lig_WK.feb"
r = subprocess.run([FEB, "-i", feb.name, "-o", "ms.xplt"],
                   capture_output=True, text=True, cwd=str(W), timeout=3000)
out = (r.stdout or "") + (r.stderr or "")
Path(R / "temp" / "pyfebio_demo" / "_ms_febio.txt").write_text(out, encoding="utf-8")
lines = out.splitlines()
print(f"FEBio 输出 {len(lines)} 行")
i = [k for k, l in enumerate(lines) if "failed to converge" in l]
k = i[0] if i else len(lines) - 1
print("\n".join(lines[max(0, k - 42):k + 6]))
