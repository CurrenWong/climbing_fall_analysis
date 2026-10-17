"""回答「为什么卡在 ~519 N」：同位移 (+1.2mm) 下 基线 vs A1+A2 的对比。

要看的机制候选：
  H-a: A2（解除横向约束）让距骨倾斜/平动更大 ⇒ 胫距软骨边缘过压 ⇒ 单元翻负
  H-b: A1（距下接触）本身改变载荷分配
取 acceptance.py 报告的「距骨低端刚体转动」与 |u|max + 两顶面合力。
"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
ACC = ROOT / "scripts/ankle_fe/acceptance.py"
OUT = ROOT / "temp/pyfebio_demo/a4why"
OUT.mkdir(parents=True, exist_ok=True)

COMMON = ["DRIVER=dispx", "DISP_MM=1.20", "FIX_TALUS=transverse",
          "DET_MIN=0.0", "TRUNC_MM=60", "AXIS_SRC=full", "ALIGN_AXIS=1",
          "CART_NEW=1", "TIE=1",
          "CART_DIR=" + str(ROOT / "temp/thums/cart_m1"), "POSE_IE=0"]

CASES = {
    "base_1.20": {},
    "a1a2_1.20": {"SUBTALAR": "1", "SUBTALAR_LAUGON": "PENALTY",
                  "TALUS_SUPPORT": "spring", "TALUS_K": "500", "MAX_UPS": "0"},
    "a2only_1.20": {"SUBTALAR": "0", "TALUS_SUPPORT": "spring",
                    "TALUS_K": "500", "MAX_UPS": "0"},
}

KEYS = ["两顶面合力", "|u|max", "距骨低端刚体转动", "全体刚体拟合转动",
        "求解", "距骨约束反力"]

for tag, extra in CASES.items():
    env = os.environ.copy()
    env.update(extra); env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(ACC)] + COMMON, cwd=str(ROOT),
                       env=env, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    (OUT / f"{tag}.txt").write_text(out, encoding="utf-8")
    print(f"=== {tag} ===")
    for ln in out.split("\n"):
        if any(k in ln for k in KEYS):
            print("   ", ln.strip()[:150])
    print()
