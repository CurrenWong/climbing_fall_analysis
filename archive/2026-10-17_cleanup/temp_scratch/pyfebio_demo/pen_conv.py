"""① tie 罚因子扫描（每个 m 都要够硬，否则软骨底面滑移主导结果）
   ② 在够硬罚因子下重跑收敛研究。"""
import os
import re
import subprocess
from pathlib import Path
import sys

sys.path.insert(0, str(Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo")))
import conv_final as CF  # noqa: E402  复用其 post() 与 pyfebio 补丁

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = str(ROOT / ".venv/Scripts/python.exe")
env0 = dict(os.environ)
env0["PATH"] = r"D:\Program\FEBioStudio\bin;" + env0["PATH"]
env0.update(TIE="1", DRIVER="disp", DISP_MM="0.1472", TIME_STEPS="20",
            LAUGON="AUGLAG", CP_DOMAIN="bary", CP_DIAG="02")


def run(mesh, pen, tag):
    e = dict(env0, CART_MESH=mesh, TIE_PENALTY=str(pen))
    r = subprocess.run([PY, "-u", str(ROOT / "scripts/ankle_fe/build_feb.py")],
                       cwd=ROOT, env=e, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=7200)
    o = r.stdout + r.stderr
    rc = re.search(r"FEBio rc=(\d+)", o)
    tm = re.search(r"用时 ([\d.]+)s", o)
    if not (rc and rc.group(1) == "0"):
        return None, f"FAIL rc={rc.group(1) if rc else '?'}"
    p = CF.post()
    return p, f"{tm.group(1)}s"


print("=" * 78)
print("① tie 罚因子扫描（m=1，看是否收敛到共形解）")
print("=" * 78)
print(f"{'pen':>6s} {'耗时':>8s} {'|u|max':>9s} {'|F|(N)':>9s} {'软骨p99':>9s}")
for pen in (1, 4, 16, 64):
    p, info = run("cartilage_mesh_patch_m1.npz", pen, f"m1p{pen}")
    if p:
        print(f"{pen:>6d} {info:>8s} {p['umax']:>9.4f} {p['F']:>9.3f} "
              f"{p.get('cart_p99', 0):>9.4f}")
    else:
        print(f"{pen:>6d} {info}")

print()
print("=" * 78)
print("② 收敛研究 @ TIE_PENALTY=16")
print("=" * 78)
print(f"{'m':>2s} {'底三角':>7s} {'接触面':>14s} {'耗时':>8s} {'|u|max':>9s} "
      f"{'|F|(N)':>9s} {'软骨p99':>9s}")
for m in (1, 2, 3):
    mesh = f"cartilage_mesh_patch_m{m}.npz"
    if not (ROOT / "temp/thums" / mesh).exists():
        e = dict(env0, CP_M=str(m), CP_OUT=mesh)
        subprocess.run([PY, "-u", str(ROOT / "scripts/ankle_fe/cart_patch.py")],
                       cwd=ROOT, env=e, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    p, info = run(mesh, 16, f"m{m}p16")
    if p:
        print(f"{m:>2d} {'':>7s} {'':>14s} {info:>8s} {p['umax']:>9.4f} "
              f"{p['F']:>9.3f} {p.get('cart_p99', 0):>9.4f}")
    else:
        print(f"{m:>2d} {info}")
