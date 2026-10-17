"""压缩侧姿势 × 罚因子：用 acceptance.py 取 C2（内外翻角）并检验 C5（排序鲁棒）。

两档罚因子：PENALTY=0.1 / 1，看姿势排序是否不变。
"""
import os
import re
import subprocess
import sys
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
PY = str(R / ".venv" / "Scripts" / "python.exe")
ACC = str(R / "scripts" / "ankle_fe" / "acceptance.py")
POSES = [(0, 0), (10, 0), (-10, 0), (0, 10), (0, -10)]


def run(tilt, tml, pen):
    args = [f"DRIVER=dispx", f"DISP_MM=-1.0", f"LOAD_SIGN=+1",
            "FIX_TALUS=transverse", "TRUNC_MM=60", "AXIS_SRC=full",
            "ALIGN_AXIS=1", "CART_NEW=1", "TIE=1", "DET_MIN=0.0",
            "TALUS_RIGID=1", "TAL_ROTLOCK=0", "TIME_STEPS=20",
            "CART_DIR=" + str(R / "temp" / "thums" / "cart_m1"),
            f"LAUGON=PENALTY", f"PENALTY={pen}", f"POSE_TILT={tilt}",
            f"POSE_TILT_ML={tml}"]
    p = subprocess.run([PY, "-u", ACC] + args, capture_output=True, text=True,
                       cwd=str(R), timeout=1800)
    return (p.stdout or "") + (p.stderr or "")


def grab(t, pat, cast=float):
    m = re.search(pat, t)
    return cast(m.group(1)) if m else None


print(f"{'pose':>10} {'pen':>5} {'NORMAL':>7} {'Fz':>9} {'Fshear':>8} "
      f"{'刚体转动°':>10} {'内外翻°':>9}", flush=True)
out = {}
for pen in ("1", "0.1"):
    for tilt, tml in POSES:
        t = run(tilt, tml, pen)
        ok = "NORMAL TERMINATION" in t
        fz = grab(t, r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)")
        fx = grab(t, r"两顶面合力 = \[\s*([-\d.e+]+)")
        fy = grab(t, r"两顶面合力 = \[\s*[-\d.e+]+\s+([-\d.e+]+)")
        rb = grab(t, r"全体刚体拟合转动 = ([\d.eE+-]+)")
        inv = grab(t, r"关节内外翻角[^\d]*([\d.eE+-]+)")
        sh = (fx * fx + fy * fy) ** 0.5 if (fx is not None and fy is not None) else None
        out[(tilt, tml, pen)] = (sh, inv)
        print(f"{tilt:>4}/{tml:<5} {pen:>5} {str(ok):>7} {fz if fz else 0:>9.2f} "
              f"{sh if sh else 0:>8.2f} {rb if rb else 0:>10.4f} "
              f"{inv if inv else 0:>9.4f}", flush=True)

print("\n===== C5：剪切排序对罚因子是否鲁棒 =====")
for pen in ("1", "0.1"):
    order = sorted(((k[0], k[1], v[0]) for k, v in out.items() if k[2] == pen),
                   key=lambda x: -x[2])
    print(f"  PENALTY={pen:>4}: " + " > ".join(f"({a},{b}):{s:.1f}" for a, b, s in order))
