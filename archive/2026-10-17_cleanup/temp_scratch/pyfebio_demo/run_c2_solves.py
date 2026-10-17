"""跑姿势（真求解）+ 权威应变读数。

对每个姿势：求解 → 存 .feb / disp / lig_read → 用 lig_strain.py 分析。
中性作为基线（几何参照），故 (c) 总应变 = 相对中性姿势。
"""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
C2D = WORK / "c2"
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
STRAIN = ROOT / "scripts/ankle_fe/lig_strain.py"

BASE = {
    "DET_MIN": "0.0", "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
    "FIX_TALUS": "transverse", "DRIVER": "dispx", "DISP_MM": "-1.0",
    "LOAD_SIGN": "+1", "CART_NEW": "1", "TIE": "1",
    "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "LOG_DATA": "1",
}
vals = [float(x) for x in (sys.argv[1] if len(sys.argv) > 1
                           else "5,-5").split(",")]

for v in vals:
    tag = f"IE_{int(round(v*100))}"
    env = os.environ.copy()
    env.update(BASE)
    env["POSE_IE"] = str(v)
    env["PYTHONIOENCODING"] = "utf-8"
    print(f"\n{'#'*90}\n# 求解 POSE_IE={v}  ->  {tag}\n{'#'*90}")
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    (ROOT / f"temp/pyfebio_demo/_solve_{tag}.log").write_text(out, encoding="utf-8")
    m = re.search(r"FEBio rc=(-?\d+)\s+用时 ([\d.]+)s", out)
    nm = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in out else "ERROR"
    print(f"  rc={m.group(1) if m else '?'}  用时={m.group(2) if m else '?'}s  "
          f"{nm}  负J={len(re.findall('negative jacobian', out, re.I))}")
    for l in out.splitlines():
        if l.startswith(("[pose]", "[mode]", "[axis]", "[trunc]")) or \
           "两顶面合力" in l or "关节面最小间距" in l:
            print("   " + l.strip()[:118])
    for suf in (".feb",):
        src = WORK / f"ankle_lig_WK{suf}"
        if src.exists():
            shutil.copy2(src, C2D / f"{tag}{suf}")
    for src, dst in ((WORK / "disp_WK.txt", C2D / f"disp_{tag}.txt"),
                     (WORK / "lig_read_WK.npz", C2D / f"lig_{tag}.npz")):
        if src.exists():
            shutil.copy2(src, dst)
    # 权威应变读数
    print(f"\n--- 权威应变读数（lig_strain.py）---")
    r2 = subprocess.run(
        [str(PY), "-u", str(STRAIN), f"NEUTRAL=c2/IE_0.feb", "TAG=WK",
         f"POSED={C2D / f'{tag}.feb'}", f"DISP={C2D / f'disp_{tag}.txt'}",
         f"LIGNPZ={C2D / f'lig_{tag}.npz'}"],
        cwd=str(ROOT), capture_output=True, text=True, encoding="utf-8",
        errors="replace", timeout=1200)
    print((r2.stdout or "") + (r2.stderr or ""))
