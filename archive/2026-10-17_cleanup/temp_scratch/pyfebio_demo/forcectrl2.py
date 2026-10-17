"""力控对照（用耦合就绪的 ankle_lig_feb.py：顶面横向约束 + 轴向自由）。

问题：同力下 tet4 vs penta6 的**应力场**是否一致？
  线性弹性：同力 ⇒ σ 由平衡决定，**与刚度无关**（只有 ε=σ/E 依赖 E）
  ⇒ 若一致：2.2 倍只是"位移刚度"伪影，**不阻塞 OpenSim 耦合**
  ⇒ 若不一致：载荷分配错误，**必须先修**
"""
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
TMP = ROOT / "temp/pyfebio_demo"
CASES = [("penta6", "prism"), ("tet4", "tet4")]

env = dict(os.environ)
env["PATH"] = r"D:\Program\FEBioStudio\bin;" + env["PATH"]
env.update(DRIVER="force", FORCE_N="600", TOPFIX="lateral",
           FIX_TALUS="transverse", LIG_MODE="none")

for name, elem in CASES:
    print(f"\n{'='*74}\n▶ {name}  DRIVER=force 600N  TOPFIX=lateral  FIX_TALUS=transverse"
          f"  LIG_MODE=none  CART_ELEM={elem}\n{'='*74}")
    e = dict(env, CART_ELEM=elem)
    r = subprocess.run([str(ROOT / ".venv/Scripts/python.exe"), "-u",
                        str(ROOT / "scripts/ankle_fe/ankle_lig_feb.py")],
                       cwd=ROOT, env=e, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    out = r.stdout + r.stderr
    (TMP / f"_fc2_{name}.txt").write_text(out, encoding="utf-8")
    KEEP = ("rc=", "NORMAL", "ERROR", "negative jacob", "failed to converge",
            "软骨", "接触面", "单元", "支反力", "位移", "写成", "写 ", ".feb",
            "轴", "uz", "converg")
    for ln in out.split("\n"):
        if any(k in ln for k in KEEP):
            print("   ", ln.strip()[:150])
    # 找 xplt
    for p in sorted(ROOT.glob("temp/thums/**/*.xplt"), key=lambda x: -x.stat().st_mtime)[:1]:
        if p.stat().st_mtime > 0:
            dst = TMP / f"fc2_{name}.xplt"
            shutil.copy2(p, dst)
            print(f"    -> xplt {p.name} -> fc2_{name}.xplt ({dst.stat().st_size:,} B)")
print("\n完成")
