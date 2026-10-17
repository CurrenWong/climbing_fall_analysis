"""验证「加载轴」是否可复现：同配置跑两次，比对 stdout 里的轴与网格统计。"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
SC = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
CFG = {"DRIVER": "dispx", "DISP_MM": "0.01", "FIX_TALUS": "transverse",
       "TRUNC_MM": "60", "TIME_STEPS": "1"}

for run in (1, 2):
    env = os.environ.copy()
    env.update(CFG)
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SC)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    print(f"\n===== run {run}  (rc={r.returncode}) =====")
    for pat in (r"加载轴（胫骨->距骨）= \[[^\]]+\]",
                r"\[trunc\] 截断到[^\n]+",
                r"翻转修正[^\n]+",
                r"顶面投影面积[^\n]+",
                r"单元：骨\+软骨[^\n]+",
                r"FEBio rc=[^\n]+"):
        m = re.search(pat, out)
        print("  " + (m.group(0)[:100] if m else f"(未找到) {pat[:30]}"))
