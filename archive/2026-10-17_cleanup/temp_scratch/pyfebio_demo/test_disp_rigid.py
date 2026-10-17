"""对照：刚体+弹簧系统在「位移控制」下能承受多大位移（判 step1 落位是否本身有问题）。"""
import os
import re
import subprocess

R = r"D:\Project\climbing_fall_analysis"
PY = R + r"\.venv\Scripts\python.exe"
ACC = R + r"\scripts\ankle_fe\acceptance.py"
BASE = dict(FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full", ALIGN_AXIS="1",
            CART_DIR=R + r"\temp\thums\cart_m1", CART_NEW="1", TIE="1",
            DET_MIN="0.0", TALUS_RIGID="1", TAL_ROTLOCK="0")
print("位移控制 + 刚体距骨：")
for d in ("0.25", "0.5", "0.9", "1.6"):
    dd = dict(BASE)
    dd.update(DRIVER="dispx", DISP_MM=d)
    e = dict(os.environ)
    e.update(dd)
    r = subprocess.run([PY, "-u", ACC] + [f"{k}={v}" for k, v in dd.items()],
                       capture_output=True, text=True, env=e, cwd=R, timeout=3000)
    t = r.stdout + r.stderr
    nm = ("NORMAL" if "NORMAL TERMINATION" in t else
          ("ERROR" if "E R R O R" in t else "?"))
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)\]", t)
    nj = re.findall(r"(\d+) negative jacobian", t)
    print(f"  DISP={d:>4}  {nm:>6}  Fz={float(fz.group(1)) if fz else float('nan'):>8.2f} N"
          f"  negJ={nj[:1]}")
