"""LRN-062 判据检验：位移控制下"力"是否可信，取决于柔度由谁提供。
物理构件（软骨/韧带）⇒ 可信 ✅    数值参数（接触罚因子）⇒ 伪影 ❌
检验：扫接触罚因子/laugon，力显著变化 ⇒ 伪影。
"""
import os
import re
import subprocess

ROOT = r"D:\Project\climbing_fall_analysis"
PY = rf"{ROOT}\.venv\Scripts\python.exe"
ACC = rf"{ROOT}\scripts\ankle_fe\acceptance.py"
BASE = dict(DRIVER="dispx", DISP_MM="1.00", FIX_TALUS="transverse", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1",
            CART_DIR=rf"{ROOT}\temp\thums\cart_m1", CART_NEW="1", TIE="1", DET_MIN="0.0")

VARIANTS = [
    ("基线 (PENALTY=1)", {}),
    ("PENALTY=0.1", {"PENALTY": "0.1"}),
    ("PENALTY=10", {"PENALTY": "10"}),
    ("PENALTY=100", {"PENALTY": "100"}),
    ("LAUGON=AUGLAG", {"LAUGON": "AUGLAG"}),
    ("TIE_PEN=10 auto=0", {"TIE_PENALTY": "10", "TIE_AUTO_PENALTY": "0"}),
    ("TIE_PEN=1000 auto=0", {"TIE_PENALTY": "1000", "TIE_AUTO_PENALTY": "0"}),
]
print(f"{'变体':<24} {'状态':>7} {'Fz(N)':>11}  相对基线")
print("-" * 62)
base = None
for name, ov in VARIANTS:
    e = dict(os.environ)
    e.update(BASE)
    e.update(ov)
    r = subprocess.run([PY, "-u", ACC] + [f"{k}={v}" for k, v in
                                          list(BASE.items()) + list(ov.items())],
                       capture_output=True, text=True, env=e, cwd=ROOT, timeout=1800)
    t = r.stdout + r.stderr
    norm = "NORMAL" if "NORMAL TERMINATION" in t else ("ERROR" if "ERROR" in t else "?")
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)\]", t)
    v = float(fz.group(1)) if fz else float("nan")
    if base is None and v == v:
        base = v
    rel = f"{(v/base-1)*100:+.1f}%" if (base and v == v) else "-"
    print(f"{name:<24} {norm:>7} {v:>11.2f}  {rel}")
