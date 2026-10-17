"""收口验证：MEMORY 尺寸 + 死链、LRN 唯一性、零回归（旧格式 + 新格式）。"""
import io
import os
import re
import subprocess

P = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
PY = rf"{P}\.venv\Scripts\python.exe"
ACC = rf"{P}\scripts\ankle_fe\acceptance.py"
ok = True

f = os.path.join(W, "MEMORY.md")
n = os.path.getsize(f)
ok &= n <= 10240
print(f"[1] MEMORY.md {n}/10240  余量 {10240-n}  {'PASS' if n <= 10240 else 'FAIL'}")

t = io.open(f, encoding="utf-8").read()
refs = sorted(set(re.findall(r"references/([A-Za-z0-9_\-\.]+\.md)", t)))
miss = [r for r in refs if not os.path.exists(os.path.join(W, "references", r))]
ok &= not miss
print(f"[2] references {len(refs)} 引用 / 缺 {len(miss)} {miss}  {'PASS' if not miss else 'FAIL'}")

lf = os.path.join(P, ".learnings", "LEARNINGS.md")
lt = io.open(lf, encoding="utf-8").read()
ids = re.findall(r"^##\s*\[LRN-(\d{8}-\d{3})\]", lt, re.M)
dups = sorted({i for i in ids if ids.count(i) > 1})
c66 = ids.count("20261005-066")          # ★ 只数**定义**（## 标题），不数正文引用
ok &= not dups and c66 == 1
print(f"[3] LRN 共 {len(ids)} 条；重复 {dups}；066 定义 {c66} 次  "
      f"{'PASS' if not dups and c66 == 1 else 'FAIL'}")


def run(extra, disp):
    args = (f"DRIVER=dispx", f"DISP_MM={disp}", "FIX_TALUS=transverse",
            "TRUNC_MM=60", "AXIS_SRC=full", "ALIGN_AXIS=1") + tuple(extra)
    r = subprocess.run([PY, "-u", ACC] + list(args), capture_output=True,
                       text=True, cwd=P, timeout=2400)
    s = r.stdout + r.stderr
    m = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)\]", s)
    return ("NORMAL" if "NORMAL TERMINATION" in s else
            ("ERROR" if "E R R O R" in s else "?")), (float(m.group(1)) if m else float("nan"))


nr, v = run([], "0.25")
p = nr == "NORMAL" and abs(v - 209.5589) < 0.01
ok &= p
print(f"[4] 旧格式 0.25 -> {nr} {v:.4f} N (期望 209.5589)  {'PASS' if p else 'FAIL'}")

nr, v = run(["CART_DIR=" + rf"{P}\temp\thums\cart_m1", "CART_NEW=1", "TIE=1",
             "DET_MIN=0.0"], "1.60")
p = nr == "NORMAL" and abs(v - 592.6068) < 0.05
ok &= p
print(f"[5] 新格式+TIE 1.60 -> {nr} {v:.4f} N (期望 592.6068)  {'PASS' if p else 'FAIL'}")

nr, v = run(["CART_DIR=" + rf"{P}\temp\thums\cart_m1", "CART_NEW=1", "TIE=1",
             "DET_MIN=0.0", "TALUS_SUPPORT=spring"], "0.25")
p = nr == "NORMAL"
ok &= p
print(f"[6] 弱弹簧 0.25 -> {nr} {v:.2f} N (期望 NORMAL)  {'PASS' if p else 'FAIL'}")

print(f"\n=== 总结: {'ALL PASS ✅' if ok else 'HAS FAIL ❌'} ===")
