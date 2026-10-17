"""m2+PENALTY 为何第一步就挂？抓日志 + 罚因子变体。"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio"
PY = ROOT / ".venv/Scripts/python.exe"

VARIANTS = [
    ("D ", dict(CART_MESH="cartilage_mesh_patch_m2.npz", TIE="1", LAUGON="PENALTY")),
    ("D1", dict(CART_MESH="cartilage_mesh_patch_m2.npz", TIE="1", LAUGON="PENALTY",
                TIE_AUTO_PENALTY="0")),
    ("D2", dict(CART_MESH="cartilage_mesh_patch_m2.npz", TIE="1", LAUGON="PENALTY",
                TIE_AUTO_PENALTY="0", TIE_PENALTY="1000")),
    ("D3", dict(CART_MESH="cartilage_mesh_patch_m2.npz", TIE="1", LAUGON="PENALTY",
                PENALTY="0.001")),
]

base = dict(os.environ)
base.update({"DRIVER": "disp", "DISP_MM": "0.1472", "TIME_STEPS": "20"})
for tag, cfg in VARIANTS:
    env = dict(base)
    for k, v in cfg.items():
        env[k] = v
    print(f"\n{'='*76}\n>>> {tag}  {cfg}\n{'='*76}", flush=True)
    r = subprocess.run([str(PY), "-u", "scripts/ankle_fe/build_feb.py"],
                       cwd=str(ROOT), env=env, capture_output=True, text=True,
                       timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    log = W / "febio.log"
    print("  " + (re.search(r"FEBio rc=\d+\s+用时\s+[\d.]+s", out) or
                  re.search(r"FEBio rc=.*", out)).group(0))
    # 取 FEBio 日志里的失败信息
    if log.exists():
        L = log.read_text(encoding="utf-8", errors="replace").splitlines()
        keys = ("warning", "error", "fail", "not converge", "diverg", "singular",
                "negative", "MAX", "iterations", "reforming", "step", "TIME",
                "Reforming", "Nr of")
        hits = [l for l in L if any(k.lower() in l.lower() for k in
                                   ("error", "fail", "not converge", "diverg",
                                    "singular", "negative jacobian", "stop"))]
        print("  --- 关键行 ---")
        for l in hits[:25]:
            print("    " + l.strip()[:150])
        # 时间步信息
        st = [l for l in L if "time step" in l.lower() or "Time step" in l]
        print("  --- 时间步 ---")
        for l in st[:8]:
            print("    " + l.strip()[:150])
        print("  --- 尾部 18 行 ---")
        for l in L[-18:]:
            print("    " + l.strip()[:150])
