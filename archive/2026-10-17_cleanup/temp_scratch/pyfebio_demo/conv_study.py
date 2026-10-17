"""收敛研究：CP_M = 1,2,3 扫（CP_DIAG=02 必须开，否则底面三角与骨面不匹配）"""
import os
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
THUMS = ROOT / "temp/thums"
PY = str(ROOT / ".venv/Scripts/python.exe")

env0 = dict(os.environ)
env0["PATH"] = r"D:\Program\FEBioStudio\bin;" + env0["PATH"]
env0.update(CP_DOMAIN="bary", CP_DIAG="02", TIE="1", DRIVER="disp",
            DISP_MM="0.1472", TIME_STEPS="20")

print(f"{'m':>3s} {'底面三角':>9s} {'软骨tet':>9s} {'接触面':>12s} {'剔除':>6s} "
      f"{'rc':>3s} {'FEBio':>8s} {'detJmin':>10s}")
for m in (1, 2, 3):
    out = f"cartilage_mesh_patch_m{m}.npz"
    e = dict(env0, CP_M=str(m), CP_OUT=out)
    r = subprocess.run([PY, "-u", str(ROOT / "scripts/ankle_fe/cart_patch.py")],
                       cwd=ROOT, env=e, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    ntri = re.findall(r"三角化\s+(\d+)", r.stdout)
    ncell = re.findall(r"面\s+(\d+)\(", r.stdout)
    tot_face = sum(int(x) for x in ncell) if ncell else 0

    e2 = dict(env0, CART_MESH=out)
    r2 = subprocess.run([PY, "-u", str(ROOT / "scripts/ankle_fe/build_feb.py")],
                        cwd=ROOT, env=e2, capture_output=True, text=True,
                        encoding="utf-8", errors="replace", timeout=7200)
    o = r2.stdout + r2.stderr
    cf = re.findall(r"外层接触面 (\d+)", o)
    nt = re.findall(r"保留\s+([\d,]+)\s+剔除\s+(\d+)", o)
    tet = sum(int(x.replace(",", "")) for x, _ in nt) if nt else 0
    drop = sum(int(y) for _, y in nt) if nt else 0
    rc = re.search(r"FEBio rc=(\d+)", o)
    tm = re.search(r"用时 ([\d.]+)s", o)
    dj = re.findall(r"det\(J\) min=([\d.eE+-]+)", o)
    conv = "✔" if "converged at time : 1" in o else "✘"
    print(f"{m:>3d} {tot_face:>9d} {tet:>9d} {str(cf):>12s} {drop:>6d} "
          f"{rc.group(1) if rc else '?':>3s} "
          f"{(tm.group(1)+'s') if tm else '?':>8s} {conv}")
    if rc and rc.group(1) != "0":
        for l in o.split("\n"):
            if "Negative" in l or "failed to converge" in l:
                print("       ", l.strip())
    shutil.copy2(THUMS / "febio/ankle_tibiotalar.xplt",
                 ROOT / f"temp/pyfebio_demo/conv_m{m}.xplt")
