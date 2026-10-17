"""收口验证：MEMORY.md 尺寸+死链、LRN-064 唯一性、旧/新格式双向回归。"""
import io
import os
import re
import subprocess

P = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
PY = rf"{P}\.venv\Scripts\python.exe"
ACC = rf"{P}\scripts\ankle_fe\acceptance.py"
ok = True

# --- 1. MEMORY.md 尺寸 ---
f = os.path.join(W, "MEMORY.md")
n = os.path.getsize(f)
ok &= n <= 10240
print(f"[1] MEMORY.md {n}/10240  余量 {10240 - n}  {'PASS' if n <= 10240 else 'FAIL'}")

# --- 2. MEMORY.md 引用的 references/*.md 是否都存在 ---
t = io.open(f, encoding="utf-8").read()
refs = sorted(set(re.findall(r"references/([A-Za-z0-9_\-\.]+\.md)", t)))
miss = [r for r in refs if not os.path.exists(os.path.join(W, "references", r))]
ok &= not miss
print(f"[2] references 引用 {len(refs)} 个，缺 {len(miss)} 个 {miss}  {'PASS' if not miss else 'FAIL'}")

# --- 3. LRN-064 唯一性 + 无重复号 ---
lf = os.path.join(P, ".learnings", "LEARNINGS.md")
lt = io.open(lf, encoding="utf-8").read()
c64 = lt.count("[LRN-20261005-064]")
ids = re.findall(r"^##\s*\[LRN-(\d{8}-\d{3})\]", lt, re.M)
dups = sorted({i for i in ids if ids.count(i) > 1})
ok &= (c64 == 1) and not dups
print(f"[3] LRN-064 出现 {c64} 次；重复编号 {dups}  "
      f"{'PASS' if c64 == 1 and not dups else 'FAIL'}  (共 {len(ids)} 条)")

# --- 4. 旧格式回归 ---
def run(disp, extra):
    e = dict(os.environ)
    e.update({"DRIVER": "dispx", "DISP_MM": disp, "FIX_TALUS": "transverse",
              "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1"})
    e.update(extra)
    args = [f"{k}={v}" for k, v in e.items() if k in
            ("DRIVER", "DISP_MM", "FIX_TALUS", "TRUNC_MM", "AXIS_SRC", "ALIGN_AXIS",
             "CART_DIR", "CART_NEW", "TIE", "DET_MIN")]
    r = subprocess.run([PY, "-u", ACC] + args, capture_output=True, text=True,
                       env=e, cwd=P, timeout=1800)
    s = r.stdout + r.stderr
    m = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)\]", s)
    nr = "NORMAL" if "NORMAL TERMINATION" in s else ("ERROR" if "ERROR" in s else "?")
    return nr, (float(m.group(1)) if m else float("nan"))

nr, v = run("0.25", {})
p4 = nr == "NORMAL" and abs(v - 209.5589) < 0.01
ok &= p4
print(f"[4] 旧格式 DISP=0.25 -> {nr} {v:.4f} N  (期望 NORMAL / 209.5589)  {'PASS' if p4 else 'FAIL'}")

nr2, v2 = run("1.60", {"CART_DIR": rf"{P}\temp\thums\cart_m1", "CART_NEW": "1",
                       "TIE": "1", "DET_MIN": "0.0"})
p5 = nr2 == "NORMAL" and v2 > 580
ok &= p5
print(f"[5] 新格式+TIE DISP=1.60 -> {nr2} {v2:.2f} N  (期望 NORMAL / >580)  {'PASS' if p5 else 'FAIL'}")

print(f"\n=== 总结: {'ALL PASS ✅' if ok else 'HAS FAIL ❌'} ===")
