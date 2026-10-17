"""决定性实验：同力(436.817 N)下 tet4 vs penta6 的应力场是否一致。

原理：**线性弹性下，同力 ⇒ 应力场与刚度无关**（σ 由平衡决定，ε=σ/E 才依赖 E）。
若 2.2 倍只是"位移刚度"差异 ⇒ 力控下两者应力场应一致 ⇒ 不阻塞 OpenSim 耦合。
若力控下应力场仍不同 ⇒ 是**载荷分配**错误 ⇒ 必须先修。
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio"
XPLT = WORK / "ankle_tibiotalar.xplt"
FORCE_N = "436.817"
STEPS = "30"
CASES = [("penta6", "prism"), ("tet4", "tet4")]

env = dict(os.environ)
env["PATH"] = r"D:\Program\FEBioStudio\bin;" + env["PATH"]
env.update(DRIVER="force", FORCE_N=FORCE_N, TIME_STEPS=STEPS)

results = {}
for name, elem in CASES:
    print(f"\n{'='*72}\n▶ {name}  (CART_ELEM={elem}, DRIVER=force, FORCE_N={FORCE_N})\n{'='*72}")
    sys.stdout.flush()
    e = dict(env, CART_ELEM=elem)
    r = subprocess.run([str(ROOT / ".venv/Scripts/python.exe"), "-u",
                        str(ROOT / "scripts/ankle_fe/build_feb.py")],
                       cwd=ROOT, env=e, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=1800)
    out = r.stdout + r.stderr
    for ln in out.split("\n"):
        if any(k in ln for k in ("外层接触面", "FEBio rc=", "负体积", "单元  ", "接触面")):
            print("   ", ln.strip())
    # 保存本次 xplt
    dst = ROOT / f"temp/pyfebio_demo/fc_{name}.xplt"
    if XPLT.exists():
        shutil.copy2(XPLT, dst)
        print(f"    -> xplt 存档 {dst.name} ({dst.stat().st_size:,} B)")
    results[name] = out

(ROOT / "temp/pyfebio_demo/_forcectrl.txt").write_text(
    "\n".join(f"##### {k}\n{v}" for k, v in results.items()), encoding="utf-8")
print("\n完成，日志 -> temp/pyfebio_demo/_forcectrl.txt")
