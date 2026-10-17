import pickle, numpy as np, h5py
D = pickle.load(open(r"temp/pyfebio_demo/p3_cache.pkl", "rb"))
DATA, SIG, EL, NS = D["DATA"], D["SIG"], D["EL"], D["NS"]
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)
F_TOT = 592.60

MAT = {}
for k in EL:
    if k.startswith("lig_"): MAT[k] = (260, .4)
    elif k.startswith("tibiotalar"): MAT[k] = (10, .4)
    else: MAT[k] = (7300, .3)

print("=" * 110)
print("INDEPENDENT ENERGY CHECK:  U = 1/2 sum sigma^T C^-1 sigma * V   (uses ONLY stress)")
print("=" * 110)
totc = 0.0; per = {}
for k in EL:
    E, nu = MAT[k]; V = DATA[k][1]; s = SIG[k]
    tr = s[:, 0] + s[:, 1] + s[:, 2]
    q = ((1 + nu) * (s[:, 0] ** 2 + s[:, 1] ** 2 + s[:, 2] ** 2) - nu * tr ** 2
         + 2 * (1 + nu) * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))
    u = 0.5 * q / E * V
    per[k] = u.sum(); totc += u.sum()
print("  total (stress-only) = %.1f N.mm   |  strain-based = %.1f  |  external work = %.1f N.mm" % (
    totc, 147.90, 0.5 * F_TOT * 1.6))
print("  -> penalty (contact+tie) residual = %.1f N.mm  = %.0f%% of external work" % (
    0.5 * F_TOT * 1.6 - totc, 100 * (0.5 * F_TOT * 1.6 - totc) / (0.5 * F_TOT * 1.6)))
print()
print("  %-32s %10s %8s %12s" % ("domain", "U_Nmm", "share", "max|sigma|"))
for k in sorted(per, key=lambda x: -per[x]):
    print("  %-32s %10.3f %7.2f%% %12.4f" % (k, per[k], 100 * per[k] / totc, np.abs(SIG[k]).max()))
grp = {"tibia": [k for k in EL if k.startswith("tibia_")],
       "fibula": [k for k in EL if k.startswith("fibula_")],
       "cartilage": [k for k in EL if k.startswith("tibiotalar")],
       "foot": [k for k in EL if k.split("__")[0] in ("talus_cort", "talus_spon", "calcaneus_cort",
                                                      "calcaneus_spon", "navicular_spon")],
       "ligament": [k for k in EL if k.startswith("lig_")]}
print()
for g, ks in grp.items():
    s = sum(per[k] for k in ks)
    print("  %-12s U=%8.2f  share=%5.1f%%" % (g, s, 100 * s / totc))

print()
print("=" * 110)
print("TIE PAIR SLIP (nearest node pairing, restricted to |d_ref|<0.6 mm)")
print("=" * 110)
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    S = {k: f["meshes/0/surfaces"][k][:] for k in f["meshes/0/surfaces"].keys()}
def fn(s): return np.unique(s[:, 1:][s[:, 1:] >= 0].ravel())
for cart, bone in [("tiesurf_tibiotalar_tibia", "bonesurf_tibiotalar_tibia"),
                   ("tiesurf_tibiotalar_talus", "bonesurf_tibiotalar_talus")]:
    ia = fn(S[cart]); ib = fn(S[bone])
    d2 = ((P[ia][:, None, :] - P[ib][None, :, :]) ** 2).sum(2)
    j = d2.argmin(1); gref = np.sqrt(d2[np.arange(len(ia)), j])
    m = gref < 0.6
    du = U[ia[m]] - U[ib[j[m]]]
    print("  %-26s pairs<0.6mm: n=%d  |slip| med=%.5f mean=%.5f  duz med=%+.5f mean=%+.5f  (gap ref med=%.4f)" % (
        cart, m.sum(), np.median(np.linalg.norm(du, axis=1)), np.linalg.norm(du, axis=1).mean(),
        np.median(du[:, 2]), du[:, 2].mean(), np.median(gref[m])))

print()
print("=" * 110)
print("EXACT BODY-LEVEL BUDGET (mean uz per body)")
print("=" * 110)
def mz(pred):
    nd = np.unique(np.concatenate([EL[k]["conn"].ravel() for k in EL if pred(k)]))
    return U[nd, 2].mean(), len(nd)
steps = [("BC (tibia_top_ns, prescribed)", U[NS["tibia_top_ns"] - 1, 2].mean()),
         ("tibia bone body", mz(lambda k: k.startswith("tibia_"))[0]),
         ("tibial cartilage layer", mz(lambda k: k == "tibiotalar_tibia__penta6")[0]),
         ("talar cartilage layer", mz(lambda k: k == "tibiotalar_talus__penta6")[0]),
         ("talus bone body", mz(lambda k: k.startswith("talus_"))[0]),
         ("calcaneus body", mz(lambda k: k.startswith("calcaneus_"))[0])]
prev = None
for lab, v in steps:
    if prev is None:
        print("  %-34s uz = %+.5f" % (lab, v))
    else:
        print("  %-34s uz = %+.5f   -> absorbed %+.4f mm = %5.1f%% of 1.6" % (
            lab, v, prev - v, 100 * (prev - v) / 1.6))
    prev = v
print("  (also: fibula body uz = %+.5f, navicular uz = %+.5f)" % (
    mz(lambda k: k.startswith("fibula_"))[0], mz(lambda k: k.startswith("navicular_"))[0]))
print()
print("  implied series stiffness:")
segs = [("tibia bone (BC->body)", 0.0011), ("tie: tibia<->tib cart", None),
        ("tibiotalar joint space", None), ("tie: talus<->tal cart", None), ("talus->calcaneus", None)]
vals = [s[1] for s in steps]
for i in range(1, len(vals)):
    d = vals[i - 1] - vals[i]
    if d > 1e-5:
        print("     %-40s du=%+.4f  k=%.0f N/mm" % (
            "%s -> %s" % (steps[i - 1][0], steps[i][0]), d, F_TOT / d))
