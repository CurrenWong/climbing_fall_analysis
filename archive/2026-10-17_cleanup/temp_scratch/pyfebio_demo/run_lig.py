"""踝韧带线运行器（Python 版，替代 .bat）——用 env dict 传参，绕开 cmd `set` 的坑。

用法：
    .venv\\Scripts\\python.exe temp\\pyfebio_demo\\run_lig.py DRIVER=force FIX_TALUS=all TRUNC_MM=60

功能：
  1. 以显式 env 运行 scripts/ankle_fe/ankle_lig_feb.py
  2. 打印 rc / 用时 / 关键 FEBio 结论行
  3. 顺带核对 .feb 的 <Loads> 行数（[ERR-20261005-007] 的回归闸门）
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"

args = {}
for a in sys.argv[1:]:
    if "=" in a:
        k, v = a.split("=", 1)
        args[k.strip()] = v.strip()
env = os.environ.copy()
env.update(args)
# 未被命令行覆盖时给一组可复现的默认（与 ANKLE_GOAL_TODO 一致）
env.setdefault("LIG_MODE", "shell")
env.setdefault("TIME_STEPS", "20")
env.pop("PYTHONIOENCODING", None)
env["PYTHONIOENCODING"] = "utf-8"

print(f"[run] {args or '(默认)'}")
t0 = time.time()
r = subprocess.run(
    [str(ROOT / ".venv/Scripts/python.exe"), "-u",
     str(ROOT / "scripts/ankle_fe/ankle_lig_feb.py")],
    cwd=str(ROOT), env=env, capture_output=True, text=True,
    encoding="utf-8", errors="replace", timeout=3600,
)
el = time.time() - t0
out = (r.stdout or "") + (r.stderr or "")
LOG = ROOT / "temp/pyfebio_demo/_lig_last.log"
LOG.write_text(out, encoding="utf-8")
print(f"[rc] {r.returncode}   用时 {el:.1f}s   （日志 {LOG}）")

# ── 关键结论行 ──
KEY = ("negative jacobian", "N O R M A L", "E R R O R", "failed to converge",
       "Nr of equations", "trunc", "截断", "翻转修正", "接触面", "顶面投影面积",
       "总用时", "Elapsed time")
for l in out.splitlines():
    ls = l.strip()
    if any(k.lower() in ls.lower() for k in KEY):
        print("   " + ls[:110])

# ── <Loads> 回归闸门 ──
tag = env.get("MATSET", "WK")
if env.get("BONES"):
    tag += "_" + "_".join(env["BONES"].split(","))
feb = WORK / f"ankle_lig_{tag}.feb"
if feb.exists():
    L = feb.read_text(encoding="latin-1", errors="replace").splitlines()
    i = next((k for k, l in enumerate(L) if "<Loads>" in l), None)
    j = next((k for k, l in enumerate(L) if "</Loads>" in l), None)
    n = (j - i - 1) if (i is not None and j is not None) else -1
    tr = [l.strip() for l in L if "<traction>" in l]
    ok = n > 0 if env.get("DRIVER") == "force" else True
    print(f"[Loads] {n} 行  {'OK' if ok else 'XX 空! [ERR-20261005-007]'}  "
          f"tractions={[re.sub(r'</?traction>', '', x) for x in tr]}")

# ── 读数（若有 read_lig） ──
rl = ROOT / "scripts/ankle_fe/read_lig.py"
print(f"[next] 读数脚本存在: {rl.exists()}")
