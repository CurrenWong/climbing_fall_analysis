"""韧带几何修复后的回归 + 重测 C3（+5 内翻）。"""
import os
import re
import shutil
import subprocess
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
CASES = [("旧格式零回归(不设CART_DIR/TIE)", {"DRIVER": "disp", "DISP_MM": "1.6",
                                          "TIE": "1", "CART_NEW": "1"}),
         ("中性(POSE=0)", {"POSE_IE": "0"}),
         ("内翻 +3", {"POSE_IE": "3"}),
         ("内翻 +5", {"POSE_IE": "5"}),
         ("外翻 -3", {"POSE_IE": "-3"}),
         ("外翻 -4（原失败）", {"POSE_IE": "-4"})]

print(f"{'用例':<26} {'rc':>3} {'NORMAL':>7} {'负J':>4} {'用时':>7} {'顶面合力':>12}")
print("─" * 78)
for name, kw in CASES:
    env = os.environ.copy()
    env.update(BASE)
    env.update(kw)
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    io_ = ROOT / f"temp/pyfebio_demo/_fx_{re.sub(r'[^A-Za-z0-9]+','_',name)}.log"
    io_.write_text(out, encoding="utf-8")
    m = re.search(r"FEBio rc=(-?\d+)\s+用时 ([\d.]+)s", out)
    nm = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in out else "ERROR"
    nj = len(re.findall(r"negative jacobian", out, re.I))
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)", out)
    print(f"{name:<26} {m.group(1) if m else '?':>3} {nm:>7} {nj:>4} "
          f"{(m.group(2)+'s') if m else '?':>7} "
          f"{(float(fz.group(1)) if fz else float('nan')):>12.4f}")
    if nm == "NORMAL" and "POSE" in kw:
        tag = f"fx_IE_{kw['POSE_IE']}"
        shutil.copy2(WORK / "ankle_lig_WK.feb", C2D / f"{tag}.feb")
        shutil.copy2(WORK / "disp_WK.txt", C2D / f"disp_{tag}.txt")
        shutil.copy2(WORK / "lig_read_WK.npz", C2D / f"lig_{tag}.npz")
        r2 = subprocess.run([str(PY), "-u", str(STRAIN), "NEUTRAL=c2/IE_0.feb",
                             "TAG=WK", f"POSED={C2D / f'{tag}.feb'}",
                             f"DISP={C2D / f'disp_{tag}.txt'}",
                             f"LIGNPZ={C2D / f'lig_{tag}.npz'}"],
                            cwd=str(ROOT), capture_output=True, text=True,
                            encoding="utf-8", errors="replace", timeout=1200)
        lines = (r2.stdout or "").splitlines()
        for l in lines:
            if re.search(r"外侧（|内侧（|CFL_|ATFL_|PTFL_", l):
                print("      " + l.strip()[:118])
