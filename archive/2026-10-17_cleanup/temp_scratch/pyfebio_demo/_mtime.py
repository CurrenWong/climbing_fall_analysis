import os
import time

CANDS = [
    r"temp/thums/cart_m1/cartilage_mesh.npz",
    r"temp/thums/cart_m1_fix/cartilage_mesh.npz",
    r"temp/thums/thums_lowerlimb.npz",
    r"scripts/ankle_fe/ankle_lig_feb.py",
    r"scripts/ankle_fe/build_feb.py",
    r"scripts/ankle_fe/cart_patch.py",
    r"scripts/ankle_fe/ligaments.py",
]
print("now =", time.strftime("%Y-%m-%d %H:%M:%S"))
for p in CANDS:
    if os.path.exists(p):
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(p)))}"
              f"  {os.path.getsize(p):>12,}  {p}")
    else:
        print("缺失            ", p)
# 目录下所有 npz
d = r"temp/thums/cart_m1"
if os.path.isdir(d):
    print("\n-- cart_m1 内 --")
    for f in sorted(os.listdir(d)):
        p = os.path.join(d, f)
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(p)))}"
              f"  {os.path.getsize(p):>12,}  {f}")
