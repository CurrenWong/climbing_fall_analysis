import numpy as np, os, h5py, re

D = r"temp/thums/febio_lig"
for nm in ["reac_tibia_WK.txt", "reac_fibula_WK.txt", "reac_talus_WK.txt"]:
    p = os.path.join(D, nm)
    raw = open(p, errors="replace").read()
    lines = [l for l in raw.splitlines() if l.strip()]
    print("=== %s : %d lines, %.1f KB" % (nm, len(lines), os.path.getsize(p) / 1024))
    for l in lines[:4]:
        print("   HEAD |", l[:150])
    print("   ...")
    for l in lines[-4:]:
        print("   TAIL |", l[:150])
    print()

npz = np.load(os.path.join(D, "lig_read_WK.npz"), allow_pickle=True)
print("=== lig_read_WK.npz ===")
for k in npz.files:
    a = npz[k]
    if a.ndim == 0:
        print("  %-28s scalar = %r" % (k, a.item()))
    elif a.dtype.kind in "US" or a.ndim == 1:
        print("  %-28s shape=%s dtype=%s  %s" % (k, a.shape, a.dtype, np.array2string(a[:12], precision=6)))
    else:
        print("  %-28s shape=%s dtype=%s" % (k, a.shape, a.dtype))
        print("       last row:", np.array2string(a[-1], precision=4))
