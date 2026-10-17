"""检查「刚体+单点弹簧=近机构」的假设：锁转轴 / 加大 k 是否恢复收敛。"""
import os
import re
import subprocess

R = r"D:\Project\climbing_fall_analysis"
PY = R + r"\.venv\Scripts\python.exe"
ACC = R + r"\scripts\ankle_fe\acceptance.py"
BASE = dict(DRIVER="force", FORCE_N="5", FIX_TALUS="transverse", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1", CART_DIR=R + r"\temp\thums\cart_m1",
            CART_NEW="1", TIE="1", DET_MIN="0.0", TALUS_RIGID="1")
CASES = [
    ("rot 锁 + k=500", dict(TAL_ROTLOCK="1")),
    ("rot 开 + k=5000", dict(TAL_ROTLOCK="0", TALUS_SPRING_K="5000")),
    ("rot 开 + k=50000", dict(TAL_ROTLOCK="0", TALUS_SPRING_K="50000")),
]
for tag, ov in CASES:
    d = dict(BASE)
    d.update(ov)
    e = dict(os.environ)
    e.update(d)
    r = subprocess.run([PY, "-u", ACC] + [f"{k}={v}" for k, v in d.items()],
                       capture_output=True, text=True, env=e, cwd=R, timeout=3000)
    t = r.stdout + r.stderr
    tm = re.search(r"末态\s*:\s*step=(\d+)\s*t=([\d.]+)", t)
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)\]", t)
    tal = re.search(r"距骨约束反力合力 = \[\s*([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\]\s*N\s*\|\s*F\|=([\d.]+)", t)
    print(f"\n===== {tag} =====")
    print(f"  末态 t={tm.group(2) if tm else '?'}  "
          f"两顶面合力={float(fz.group(1)) if fz else float('nan'):.2f} N  "
          f"距骨反力|F|={float(tal.group(4)) if tal else float('nan'):.1f} N")
