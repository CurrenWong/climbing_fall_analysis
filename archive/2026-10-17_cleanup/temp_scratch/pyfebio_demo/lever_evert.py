"""外翻(−IE)不收敛的杠杆矩阵（每个 ~5-8 s，便宜）。"""
import io
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"

BASE = {
    "DET_MIN": "0.0", "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
    "FIX_TALUS": "transverse", "DRIVER": "dispx", "DISP_MM": "-1.0",
    "LOAD_SIGN": "+1", "CART_NEW": "1", "TIE": "1",
    "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "LOG_DATA": "1",
}

CASES = [
    ("基线 IE=-5",              {"POSE_IE": "-5"}),
    ("步长×10 (200 步)",         {"POSE_IE": "-5", "TIME_STEPS": "200"}),
    ("落位限 0.1mm",            {"POSE_IE": "-5", "POSE_SEAT_LIM": "0.1"}),
    ("落位限 0.0 (不落位)",      {"POSE_IE": "-5", "POSE_SEAT_LIM": "0.0"}),
    ("DET_MIN=1e-2",           {"POSE_IE": "-5", "DET_MIN": "1e-2"}),
    ("姿势小一点 IE=-3",         {"POSE_IE": "-3"}),
    ("载荷小一半 -0.5mm",        {"POSE_IE": "-5", "DISP_MM": "-0.5"}),
    ("对照 内翻 IE=+5",          {"POSE_IE": "+5"}),
]

print(f"{'用例':<24} {'rc':>3} {'NORMAL':>7} {'负J':>4} {'用时':>7}  "
      f"{'顶面合力':>12}  {'关节面间距':>10}  备注")
print("─" * 108)
for name, kw in CASES:
    env = os.environ.copy()
    env.update(BASE)
    env.update(kw)
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(-?\d+)\s+用时 ([\d.]+)s", out)
    nm = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in out else "ERROR"
    nj = len(re.findall(r"negative jacobian", out, re.I))
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)", out)
    gp = re.search(r"落位后 ([\d.]+) mm", out)
    note = ""
    if nm == "ERROR":
        mm = re.search(r"failed to converge at time : ([\d.]+)", out)
        if mm:
            note = f"失败于 t={mm.group(1)}"
    print(f"{name:<24} {m.group(1) if m else '?':>3} {nm:>7} {nj:>4} "
          f"{(m.group(2)+'s') if m else '?':>7}  "
          f"{(float(fz.group(1)) if fz else float('nan')):>12.3f}  "
          f"{(float(gp.group(1)) if gp else float('nan')):>10.4f}  {note}")
    io.open(ROOT / f"temp/pyfebio_demo/_lever_{re.sub(r'[^A-Za-z0-9]+','_',name)}.log",
            "w", encoding="utf-8").write(out)
