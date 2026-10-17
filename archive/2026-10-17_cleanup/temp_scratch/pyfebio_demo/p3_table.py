import pickle, numpy as np, h5py
D = pickle.load(open(r"temp/pyfebio_demo/p3_cache.pkl", "rb"))
EL = D["EL"]
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)
z = P[:, 2]; uz = U[:, 2]
ztop = z.max()
bins = []
for e0 in np.arange(ztop, z.min() - 1, -1.0):
    m = (z <= e0 + 1e-9) & (z > e0 - 1.0 + 1e-9)
    bins.append((e0 - 0.5, e0, e0 - 1.0, int(m.sum()), uz[m].mean() if m.sum() else np.nan,
                 np.median(uz[m]) if m.sum() else np.nan))
print("TASK 2 TABLE (1 mm bins, ALL 4029 nodes) - condensed")
print("  z-range            n   mean_uz   med_uz   d(med)_prev  %of1.6")
prev = None
for (zc, hi, lo, n, mu, md) in sorted(bins, reverse=True):
    if n == 0:
        continue
    show = (md < -1.585) and (prev is None or prev < -1.585) and (zc < -660 or zc > -700)
    dr = "" if prev is None else "%+8.4f %7.2f%%" % (prev - md, 100 * (prev - md) / 1.6)
    if not show or (md < -1.585 and prev is not None and prev < -1.585):
        # compress long flat runs: print only every 10th
        if int(round(-zc)) % 10 != 0 and md < -1.585 and zc > -670:
            prev = md; continue
    print("  %7.2f..%7.2f %5d  %+.4f  %+.4f   %s" % (hi, lo, n, mu, md, dr))
    prev = md

print()
print("  5-mm windows ranked by |median drop| (ALL nodes)")
w = [b for b in sorted(bins, reverse=True) if b[3] > 0]
res = []
for i in range(len(w) - 4):
    a = w[i:i + 5]
    drop = a[0][5] - a[4][5]
    res.append((drop, a[0][1], a[-1][2], sum(x[3] for x in a)))
for d, hi, lo, n in sorted(res, key=lambda r: r[0])[:6]:
    print("   z %7.2f -> %7.2f   drop=%.4f mm  (%.1f%% of 1.6)  n=%d" % (hi, lo, d, 100 * d / 1.6, n))
print("   sum of all 1-mm median drops = %.4f mm" % (w[0][5] - w[-1][5]))

print()
print("  same, restricted to TIBIA+FIBULA domain nodes only (their own deformation)")
tf = np.unique(np.concatenate([EL[k]["conn"].ravel() for k in EL
                               if k.startswith("tibia_") or k.startswith("fibula_")]))
zt, uzt = z[tf], uz[tf]
b2 = []
for e0 in np.arange(ztop, z.min() - 1, -1.0):
    m = (zt <= e0 + 1e-9) & (zt > e0 - 1.0 + 1e-9)
    if m.sum():
        b2.append((e0, int(m.sum()), uzt[m].mean(), np.median(uzt[m]), uzt[m].std(ddof=0)))
print("   z(hi)      n   mean     med      std")
for hi, n, mu, md, sd in b2:
    if n:
        print("   %7.2f  %4d  %+.5f  %+.5f  %.2e" % (hi, n, mu, md, sd))
print("   tibia+fibula overall: mean=%.5f std=%.2e ptp=%.5f  -> absorbed %.4f mm (%.2f%% of 1.6)" % (
    uzt.mean(), uzt.std(ddof=0), np.ptp(uzt), np.ptp(uzt), 100 * np.ptp(uzt) / 1.6))
