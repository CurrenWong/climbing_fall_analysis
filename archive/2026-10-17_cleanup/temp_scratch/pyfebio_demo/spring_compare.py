"""弱弹簧支撑：① 轴向响应对比硬约束；② 力控制(DRIVER=force)可行性。"""
import os
import re
import subprocess

ROOT = r"D:\Project\climbing_fall_analysis"
PY = rf"{ROOT}\.venv\Scripts\python.exe"
ACC = rf"{ROOT}\scripts\ankle_fe\acceptance.py"
BASE = dict(FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full", ALIGN_AXIS="1",
            CART_DIR=rf"{ROOT}\temp\thums\cart_m1", CART_NEW="1", TIE="1", DET_MIN="0.0")


def run(**ov):
    e = dict(os.environ)
    d = dict(BASE)
    d.update(ov)
    e.update(d)
    r = subprocess.run([PY, "-u", ACC] + [f"{k}={v}" for k, v in d.items()],
                       capture_output=True, text=True, env=e, cwd=ROOT, timeout=2400)
    t = r.stdout + r.stderr
    norm = ("NORMAL" if "NORMAL TERMINATION" in t else
            ("ERROR" if "E R R O R" in t or "failed to converge" in t else "?"))
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)\]", t)
    um = re.search(r"\|u\|max = ([\d.]+)", t)
    return norm, (float(fz.group(1)) if fz else float("nan")), \
        (float(um.group(1)) if um else float("nan")), t


print("=== ① 轴向响应：硬约束 vs 弱弹簧 (DRIVER=dispx) ===")
print(f"{'DISP':>6} {'硬约束 Fz':>12} {'弱弹簧 Fz':>12} {'差':>9}")
print("-" * 46)
for disp in ("0.25", "0.50", "1.00", "1.60"):
    n1, f1, _, _ = run(DRIVER="dispx", DISP_MM=disp)
    n2, f2, _, _ = run(DRIVER="dispx", DISP_MM=disp, TALUS_SUPPORT="spring")
    d = (f2 - f1) / f1 * 100 if f1 else float("nan")
    print(f"{disp:>6} {n1[:4]:>4}{f1:>8.2f} {n2[:4]:>4}{f2:>8.2f} {d:>8.1f}%")

print("\n=== ② 力控制 DRIVER=force（此前因刚体模态奇异而失败）===")
for k in ("500", "50", "5000"):
    n, f, um, t = run(DRIVER="force", FORCE_N="600", TALUS_SUPPORT="spring", TALUS_K=k)
    print(f"  TALUS_K={k:>5}  {n:>6}  Fz={f:>9.2f}  |u|max={um:.4f}")
    if n == "?":  # 打印症状
        for ln in t.splitlines():
            if any(x in ln for x in ("negative", "converge", "singular", "ERROR")):
                print(f"      {ln.strip()[:100]}")
                break
