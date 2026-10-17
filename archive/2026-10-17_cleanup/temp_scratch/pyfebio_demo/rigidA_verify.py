"""路线 A 的关键验证：
  ① TAL_ROTLOCK=1 vs 0（放开内外翻/内旋）—— 力与姿态是否改变
  ② DRIVER=force 力控制（此前因刚体模态奇异一直失败）
"""
import os
import re
import subprocess

ROOT = r"D:\Project\climbing_fall_analysis"
PY = rf"{ROOT}\.venv\Scripts\python.exe"
ACC = rf"{ROOT}\scripts\ankle_fe\acceptance.py"
BASE = dict(FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full", ALIGN_AXIS="1",
            CART_DIR=rf"{ROOT}\temp\thums\cart_m1", CART_NEW="1", TIE="1",
            DET_MIN="0.0", TALUS_RIGID="1")


def run(**ov):
    d = dict(BASE)
    d.update(ov)
    e = dict(os.environ)
    e.update(d)
    r = subprocess.run([PY, "-u", ACC] + [f"{k}={v}" for k, v in d.items()],
                       capture_output=True, text=True, env=e, cwd=ROOT, timeout=2400)
    t = r.stdout + r.stderr
    norm = ("NORMAL" if "NORMAL TERMINATION" in t else
            ("ERROR" if "E R R O R" in t or "failed to converge" in t else "?"))
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)\]", t)
    fx = re.search(r"两顶面合力 = \[\s*([-\d.e+]+)\s+([-\d.e+]+)", t)
    return norm, (float(fz.group(1)) if fz else float("nan")), t, \
        (float(fx.group(1)) if fx else float("nan"))


print("=== ① 转轴锁定 vs 解除（DRIVER=dispx, DISP=0.25）===")
for lock in ("1", "0"):
    n, fz, t, fx = run(DRIVER="dispx", DISP_MM="0.25", TAL_ROTLOCK=lock)
    print(f"  TAL_ROTLOCK={lock}  {n:>6}  Fz={fz:>9.2f} N  Fx={fx:>9.2f} N")

print("\n=== ② 力控制 DRIVER=force（此前必失败）===")
for fn in ("200", "600"):
    n, fz, t, fx = run(DRIVER="force", FORCE_N=fn, TAL_ROTLOCK="0")
    print(f"  FORCE_N={fn:>4}  {n:>6}  Fz={fz:>9.2f}  Fx={fx:>9.2f}")
    if n == "?":
        for ln in t.splitlines():
            if "onstraint" in ln or "converge" in ln or "negative" in ln:
                print(f"      {ln.strip()[:100]}")
                break
