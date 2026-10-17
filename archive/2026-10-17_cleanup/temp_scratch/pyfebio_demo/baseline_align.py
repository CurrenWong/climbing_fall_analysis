"""在 ALIGN_AXIS=1（真轴向）下重扫 DISP，建立可信基线 k=F/δ。

用法：python baseline_align.py            # 默认 dispx + transverse + TRUNC=60
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
PY = ROOT / ".venv/Scripts/python.exe"
SC = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"

ARGS = {a.split("=", 1)[0]: a.split("=", 1)[1] for a in sys.argv[1:] if "=" in a}
BASE = {"DRIVER": "dispx", "FIX_TALUS": "transverse", "TRUNC_MM": "60",
        "TIME_STEPS": "20", "AXIS_SRC": "full", "ALIGN_AXIS": "1", "MATSET": "WK",
        "LOG_DATA": "1"}     # ★ 必须：否则读数文件不更新，失败时会拿上一次的值冒充
BASE.update(ARGS)
DISPS = [float(x) for x in os.environ.get("DISPS", "0.30,0.50,1.00").split(",")]


def read_records(p: Path):
    out, cur, buf = {}, None, []
    if not p.exists():
        return {}
    for l in p.read_text(encoding="utf-8", errors="replace").split("\n"):
        if l.startswith("*Step"):
            if cur is not None and buf:
                out[cur[0]] = (cur[1], np.array(buf, float))
            cur = (int(l.split("=")[1]), None)
            buf = []
        elif l.startswith("*Time"):
            if cur is not None:
                cur = (cur[0], float(l.split("=")[1]))
        elif l.strip() and not l.startswith("*"):
            try:
                buf.append([float(x) for x in l.split(",")])
            except ValueError:
                pass
    if cur is not None and buf:
        out[cur[0]] = (cur[1], np.array(buf, float))
    return out


def pick(rec):
    """取 t=1.0 的记录（[LRN-049] 纪律）。"""
    if not rec:
        return None
    k = max(rec, key=lambda x: rec[x][0])
    return rec[k] if abs(rec[k][0] - 1.0) < 1e-9 else rec[k]


print(f"{'DISP(mm)':>9} {'|F|(N)':>9} {'Fproj(N)':>9} {'k=Fproj/δ':>11} "
      f"{'|u|max':>9} {'umb':>9} {'step/t':>10} negJ  status")
rows = []
tag = "WK"
for d in DISPS:
    env = os.environ.copy()
    env.update(BASE)
    env["DISP_MM"] = f"{d}"
    env["PYTHONIOENCODING"] = "utf-8"
    # 先删读数文件：失败运行时它们不会被重写，留着会「冒充」本次结果
    for _f in (f"disp_{tag}.txt", f"reac_tibia_{tag}.txt", f"reac_fibula_{tag}.txt"):
        (WORK / _f).unlink(missing_ok=True)
    r = subprocess.run([str(PY), "-u", str(SC)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    normal = "N O R M A L   T E R M I N A T I O N" in out
    negj = len(re.findall(r"negative jacobian", out, re.I))
    m = re.search(r"FEBio rc=(-?\d+)", out)
    rc = int(m.group(1)) if m else None
    dt = pick(read_records(WORK / f"disp_{tag}.txt"))
    rt = pick(read_records(WORK / f"reac_tibia_{tag}.txt"))
    rf = pick(read_records(WORK / f"reac_fibula_{tag}.txt"))
    F = (rt[1][:, 1:4].sum(0) + (rf[1][:, 1:4].sum(0) if rf else 0)
         if rt else np.zeros(3))
    Fproj = float(F[2])                              # ALIGN_AXIS ⇒ 轴=+z
    umax = float(np.linalg.norm(dt[1][:, 1:4], axis=1).max()) if dt else float("nan")
    umb = float(np.abs(dt[1][:, 3]).max()) if dt else float("nan")   # uz 最负
    k = abs(Fproj) / d if d else float("nan")
    st = f"{dt[0] if dt else -1}/{dt[1] if False else (pick(read_records(WORK/f'disp_{tag}.txt'))[0] if dt else -1)}"
    tval = pick(read_records(WORK / f"disp_{tag}.txt"))
    print(f"{d:>9.4f} {np.linalg.norm(F):>9.3f} {Fproj:>9.3f} {k:>11.1f} "
          f"{umax:>9.5f} {umb:>9.5f} {str(tval[0]) if tval else '-':>10} "
          f"{negj:>4}  {'NORMAL' if normal else 'ERROR'} (rc={rc})")
    rows.append((d, np.linalg.norm(F), Fproj, k, umax, negj, normal))

# 线性度检查
if len(rows) >= 2:
    ks = [x[3] for x in rows]
    print(f"\n  k 范围 {min(ks):.1f} ~ {max(ks):.1f} N/mm   "
          f"非线性度 {(max(ks)/min(ks)-1)*100:.1f}%")
