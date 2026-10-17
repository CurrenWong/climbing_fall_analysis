"""dump 生成的 <Step>/<Boundary>/<Loads> 段 + FEBio step 切换处的原始输出。"""
import io
import os
import re
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
         SEAT_MM="0.5", TIME_STEPS="20")
subprocess.run([str(R / ".venv" / "Scripts" / "python.exe"), "-u",
                str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
               capture_output=True, text=True, env=e, cwd=str(R), timeout=3000)
txt = (W / "ankle_lig_WK.feb").read_text(encoding="latin-1", errors="replace")
for tag in ("Step", "Loads"):
    i, j = txt.find(f"<{tag}>"), txt.find(f"</{tag}>")
    print(f"\n===== <{tag}> ===== (从 {i})")
    print(txt[i:j + len(tag) + 3] if i >= 0 else "(无)")
print("\n===== 全局 Boundary 里的 BC 名单 =====")
b = txt[txt.find("<Boundary>"):txt.find("</Boundary>")]
for m in re.finditer(r'<bc name="([^"]+)" type="([^"]+)"', b):
    print(f"  {m.group(1):<24} {m.group(2)}")
