"""受控矩阵：隔离 LAUGON / TIE / 网格 三个因子对 F 的影响。"""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio"
PY = ROOT / ".venv/Scripts/python.exe"
sys.path.insert(0, str(ROOT / "temp/pyfebio_demo"))

CONFIGS = [
    ("B: 基线", dict(CART_MESH="", TIE="0", LAUGON="PENALTY")),
    ("C: 基线+AUGLAG", dict(CART_MESH="", TIE="0", LAUGON="AUGLAG")),
    ("D: m2+PENALTY", dict(CART_MESH="cartilage_mesh_patch_m2.npz", TIE="1",
                           LAUGON="PENALTY")),
    ("E: m2+AUGLAG", dict(CART_MESH="cartilage_mesh_patch_m2.npz", TIE="1",
                          LAUGON="AUGLAG")),
]

base_env = {k: v for k, v in os.environ.items()}
base_env.update({"DRIVER": "disp", "DISP_MM": "0.1472", "TIME_STEPS": "20"})
res = {}
for tag, cfg in CONFIGS:
    env = dict(base_env)
    for k, v in cfg.items():
        if v == "":
            env.pop(k, None)
        else:
            env[k] = v
    print(f"\n{'='*76}\n>>> {tag}   {cfg}\n{'='*76}", flush=True)
    r = subprocess.run([str(PY), "-u", "scripts/ankle_fe/build_feb.py"],
                       cwd=str(ROOT), env=env, capture_output=True, text=True,
                       timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    rc = re.search(r"FEBio rc=(\d+)\s+用时\s+([\d.]+)s", out)
    tieline = re.findall(r"\[tie\].*", out)
    flo = re.findall(r"外层接触面.*", out)
    for l in tieline + flo:
        print("   " + l)
    print(f"   -> {rc.group(0) if rc else 'NO rc'}")
    if not rc:
        print(out[-2000:])
        continue
    nm = tag.split(":")[0].strip()
    for ext in (".feb", ".xplt"):
        src = W / f"ankle_tibiotalar{ext}"
        if src.exists():
            shutil.copy2(src, W / f"_{nm}{ext}")
    res[nm] = dict(rc=int(rc.group(1)), secs=float(rc.group(2)),
                   mesh=cfg["CART_MESH"] or "baseline",
                   tie=cfg["TIE"], laugon=cfg["LAUGON"])

print("\n\n" + "=" * 76)
print("读数（末帧）")
print("=" * 76)
import importlib

import numpy as np
import cmp_run
importlib.reload(cmp_run)
AX = cmp_run.AX
rows = []
for nm, m in res.items():
    f = W / f"_{nm}.feb"
    x = W / f"_{nm}.xplt"
    if not x.exists():
        continue
    import h5py
    import pyfebio.xplt as X
    X.to_hdf5(str(x), str(W / f"_r{nm}.h5"))
    with h5py.File(W / f"_r{nm}.h5", "r") as h:
        ns = sorted((int(k) for k in h["states"]), key=int)
        U = {s: np.asarray(h[f"states/{s}/node_data/displacement/all"][()], float)
             for s in ns}
        R = {s: np.asarray(h[f"states/{s}/node_data/reaction forces/all"][()], float)
             for s in ns}
    fm = ns[-1]
    mag = np.linalg.norm(U[fm], axis=1)
    rows.append((nm, m["mesh"], m["tie"], m["laugon"], m["rc"], m["secs"],
                 len(ns), (R[fm] @ AX).sum(), mag.max(),
                 np.linalg.norm(U[fm] @ AX, axis=1).max() if False else
                 float(np.abs(U[fm] @ AX).max())))

print(f"\n  {'':4} {'mesh':22s} {'tie':>4} {'laugon':>8} {'rc':>3} {'s':>6} "
      f"{'帧':>4} {'F(N)':>9} {'|u|max':>8} {'|u·ax|max':>10}")
for r in rows:
    print(f"  {r[0]:>4} {r[1]:22s} {r[2]:>4} {r[3]:>8} {r[4]:>3} {r[5]:>6.1f} "
          f"{r[6]:>4} {r[7]:>9.2f} {r[8]:>8.4f} {r[9]:>10.4f}")
