"""A1 失败二分矩阵：逐个配置跑参考工况，保存各自的 FEBio log + .feb。
全部只读诊断（不改模型参数，只切 SUBTALAR / NO_CONTACT 诊断开关）。
用法: python temp/pyfebio_demo/a1_matrix.py [tag ...]
"""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
OUT = ROOT / "temp/pyfebio_demo/a1_matrix"
OUT.mkdir(exist_ok=True)
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"

BASE = dict(
    DRIVER="dispx", DISP_MM="+1.6", DET_MIN="0.0", TRUNC_MM="60",
    AXIS_SRC="full", ALIGN_AXIS="1", FIX_TALUS="transverse",
    CART_NEW="1", TIE="1", CART_DIR="temp/thums/cart_m1", POSE_IE="0",
    LOG_DATA="1",
)

CASES = {
    "both": {},                       # A1 目标配置（双侧 + 接触）
    "nocontact": {"NO_CONTACT": "1"},  # 双侧软骨+tie，关掉所有关节接触
    "talus": {"SUBTALAR": "talus"},    # 只接距下距骨侧软骨
    "calcaneus": {"SUBTALAR": "calcaneus"},
}


def run(tag: str) -> dict:
    env = os.environ.copy()
    env.update(BASE)
    env.update(CASES[tag])
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    (OUT / f"{tag}.stdout.txt").write_text(out, encoding="utf-8")
    m = re.search(r"FEBio rc=(-?\d+)\s+用时 ([\d.]+)s", out)
    rc = int(m.group(1)) if m else None
    feb = WORK / "ankle_lig_WK.feb"
    log = WORK / "febio_WK.log"
    if feb.exists():
        shutil.copy(feb, OUT / f"{tag}.feb")
    if log.exists():
        shutil.copy(log, OUT / f"{tag}.febio.log")
    lg = WORK / "log_WK.txt"
    normal = False
    if lg.exists():
        shutil.copy(lg, OUT / f"{tag}.logdata.txt")
        normal = "N O R M A L   T E R M I N A T I O N" in lg.read_text(
            encoding="utf-8", errors="replace")
    # 顶面合力（末 t=1.0）
    F = None
    rt = WORK / "reac_tibia_WK.txt"
    if rt.exists() and normal:
        cur, buf, best = None, [], {}
        for l in rt.read_text(encoding="utf-8", errors="replace").split("\n"):
            if l.startswith("*Step"):
                cur = int(l.split("=")[1]); buf = []
            elif l.startswith("*Time"):
                cur_t = float(l.split("=")[1])
            elif l.strip() and not l.startswith("*"):
                try:
                    buf.append([float(x) for x in l.split(",")])
                except ValueError:
                    pass
                if cur is not None:
                    best[cur_t] = np_sum(buf)
        if best:
            F = best[max(best)]
    return dict(tag=tag, rc=rc, normal=normal, F=F)


def np_sum(buf):
    import numpy as np
    a = np.array(buf, float)
    return a[:, 1:4].sum(0)


def main():
    tags = sys.argv[1:] or list(CASES)
    for t in tags:
        d = run(t)
        print(f"[{t}] febio_rc={d['rc']} NORMAL={d['normal']} "
              f"胫骨顶面反力={d['F']}")


if __name__ == "__main__":
    main()
