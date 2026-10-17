"""收口：零回归（旧格式 / 新格式+TIE / 弱弹簧 / 刚体）+ LRN 唯一性。"""
import io
import os
import re
import subprocess

P = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
PY = rf"{P}\.venv\Scripts\python.exe"
ACC = rf"{P}\scripts\ankle_fe\acceptance.py"
ok = True


def run(extra, disp, driver="dispx"):
    args = (f"DRIVER={driver}", f"DISP_MM={disp}", "FIX_TALUS=transverse",
            "TRUNC_MM=60", "AXIS_SRC=full", "ALIGN_AXIS=1") + tuple(extra)
    r = subprocess.run([PY, "-u", ACC] + list(args), capture_output=True,
                       text=True, cwd=P, timeout=3000)
    s = r.stdout + r.stderr
    m = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)\]", s)
    return ("NORMAL" if "NORMAL TERMINATION" in s else
            ("ERROR" if "E R R O R" in s else "?")), \
        (float(m.group(1)) if m else float("nan"))


nr, v = run([], "0.25")
ok &= nr == "NORMAL" and abs(v - 209.5589) < 0.01
print(f"[1] 旧格式 0.25 -> {nr} {v:.4f} N (期望 209.5589)  {'PASS' if ok else 'FAIL'}")

M1 = ["CART_DIR=" + rf"{P}\temp\thums\cart_m1", "CART_NEW=1", "TIE=1", "DET_MIN=0.0"]
nr, v = run(M1, "1.60")
p = nr == "NORMAL" and abs(v - 592.6068) < 0.05
ok &= p
print(f"[2] 新格式+TIE 1.60 -> {nr} {v:.4f} N (期望 592.6068)  {'PASS' if p else 'FAIL'}")

nr, v = run(M1 + ["TALUS_SUPPORT=spring"], "0.25")
ok &= nr == "NORMAL"
print(f"[3] 弱弹簧 0.25 -> {nr} {v:.2f} N (期望 NORMAL)  {'PASS' if nr=='NORMAL' else 'FAIL'}")

for lock, exp in (("1", 202.29), ("0", 187.90)):
    nr, v = run(M1 + ["TALUS_RIGID=1", f"TAL_ROTLOCK={lock}"], "0.25")
    p = nr == "NORMAL" and abs(v - exp) < 0.5
    ok &= p
    print(f"[4] 刚体 转轴锁={lock} -> {nr} {v:.2f} N (期望 {exp})  {'PASS' if p else 'FAIL'}")

lf = os.path.join(P, ".learnings", "LEARNINGS.md")
lt = io.open(lf, encoding="utf-8").read()
ids = re.findall(r"^##\s*\[LRN-(\d{8}-\d{3})\]", lt, re.M)
dups = sorted({i for i in ids if ids.count(i) > 1})
ok &= not dups
print(f"[5] LRN 共 {len(ids)} 条，重复 {dups}  {'PASS' if not dups else 'FAIL'}")

n = os.path.getsize(os.path.join(W, "MEMORY.md"))
print(f"[6] MEMORY.md {n}/10240 (仅记录，用户已说不用检查)")

print(f"\n=== 总结: {'ALL PASS ✅' if ok else 'HAS FAIL ❌'} ===")
