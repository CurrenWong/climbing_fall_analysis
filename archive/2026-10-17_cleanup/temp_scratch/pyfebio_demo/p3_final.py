import pickle, numpy as np, h5py
D = pickle.load(open(r"temp/pyfebio_demo/p3_cache.pkl", "rb"))
DATA, SIG, EL, NS = D["DATA"], D["SIG"], D["EL"], D["NS"]
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)
    S = {k: f["meshes/0/surfaces"][k][:] for k in f["meshes/0/surfaces"].keys()}

print("=== stress outliers (test of the 100 MPa anomaly) ===")
for k in ["tibiotalar_talus__penta6", "tibiotalar_tibia__penta6"]:
    s = SIG[k]; m = np.abs(s).max(1)
    o = np.argsort(-m)
    print("  %-26s n=%d  |s|max: top5=%s  median=%.4f  mean=%.4f  n>|s|5MPa=%d" % (
        k, len(m), np.round(m[o[:5]], 3), np.median(m), m.mean(), int((m > 5).sum())))
    v = DATA[k][1]; print("       their volumes = %s mm^3" % np.round(v[o[:5]], 3))
    keep = m < 5.0
    E, nu = 10, .4; s2 = s[keep]; v2 = v[keep]
    tr = s2[:, 0] + s2[:, 1] + s2[:, 2]
    q = ((1 + nu) * (s2[:, 0] ** 2 + s2[:, 1] ** 2 + s2[:, 2] ** 2) - nu * tr ** 2
         + 2 * (1 + nu) * (s2[:, 3] ** 2 + s2[:, 4] ** 2 + s2[:, 5] ** 2))
    print("       energy excl |s|>5MPa: raw=%.2f  excl=%.2f N.mm (%d/%d el)  median stress=%.3f MPa" % (
        (0.5 * q / E * v2).sum() + (0.5 * ((1 + nu) * (s[~keep, 0] ** 2 + s[~keep, 1] ** 2 + s[~keep, 2] ** 2)
                                           - nu * (s[~keep, 0] + s[~keep, 1] + s[~keep, 2]) ** 2
                                           + 2 * (1 + nu) * (s[~keep, 3] ** 2 + s[~keep, 4] ** 2 + s[~keep, 5] ** 2))
                              / E * v[~keep]).sum(),
        (0.5 * q / E * v2).sum(), int((~keep).sum()), len(m), np.median(m)))

print()
print("=== paired tie budget (mean over matched pairs) ===")
def fn(s): return np.unique(s[:, 1:][s[:, 1:] >= 0].ravel())
for cart, bone, body in [("tiesurf_tibiotalar_tibia", "bonesurf_tibiotalar_tibia", "tibia bone"),
                         ("tiesurf_tibiotalar_talus", "bonesurf_tibiotalar_talus", "talar cartilage")]:
    ia = fn(S[cart]); ib = fn(S[bone])
    d2 = ((P[ia][:, None, :] - P[ib][None, :, :]) ** 2).sum(2)
    j = d2.argmin(1)
    print("  %-26s cartilage side mean uz=%+.4f | matched bone side mean uz=%+.4f | diff=%+.4f mm" % (
        cart, U[ia, 2].mean(), U[ib[j], 2].mean(), U[ia, 2].mean() - U[ib[j], 2].mean()))

print()
print("=== exact series stiffness table (mean-uz budget) ===")
def mz(pred):
    nd = np.unique(np.concatenate([EL[k]["conn"].ravel() for k in EL if pred(k)]))
    return U[nd, 2].mean()
F = 592.60
lv = [("BC prescribed (tibia_top_ns)", U[NS["tibia_top_ns"] - 1, 2].mean()),
      ("tibia bone body", mz(lambda k: k.startswith("tibia_"))),
      ("tibial cartilage layer", mz(lambda k: k == "tibiotalar_tibia__penta6")),
      ("talar cartilage layer", mz(lambda k: k == "tibiotalar_talus__penta6")),
      ("talus bone body", mz(lambda k: k.startswith("talus_"))),
      ("calcaneus body (fixed)", mz(lambda k: k.startswith("calcaneus_")))]
print("  %-32s %9s %9s %10s %10s" % ("segment", "uz", "du[mm]", "share", "k[F/du]"))
for i in range(1, len(lv)):
    du = lv[i - 1][1] - lv[i][1]
    print("  %-32s %+9.5f %9.4f %9.1f%% %10.0f" % (
        "%s -> %s" % (lv[i - 1][0], lv[i][0]), lv[i][1], du, 100 * du / 1.6,
        (F / du) if du > 1e-6 else float("inf")))
print("  total %.4f mm ; F/du_total = %.1f N/mm" % (lv[0][1] - lv[-1][1], F / (lv[0][1] - lv[-1][1])))

print()
print("=== per-domain max |stress| summary (stress field quality) ===")
for k in sorted(SIG, key=lambda x: -np.abs(SIG[x]).max())[:12]:
    print("  %-34s |s|max=%9.4f  |s|median=%.4f  MPa-per-E: E=%.0f" % (
        k, np.abs(SIG[k]).max(), np.median(np.abs(SIG[k])),
        260 if k.startswith("lig_") else (10 if k.startswith("tibiotalar") else 7300)))
