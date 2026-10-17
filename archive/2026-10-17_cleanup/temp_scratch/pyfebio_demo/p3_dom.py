import pickle, numpy as np, h5py
D = pickle.load(open(r"temp/pyfebio_demo/p3_cache.pkl", "rb"))
DATA, SIG, Udom, NS, EL = D["DATA"], D["SIG"], D["Udom"], D["NS"], D["EL"]
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)
    R = f["states/20/node_data/reaction forces/all"][:].astype(float)
F_TOT = 592.60; UTOT = 147.90

TIBIA = ["tibia_cort__hex8", "tibia_cort__tet4", "tibia_spon__tet4", "tibia_marrow__tet4"]
CART = ["tibiotalar_tibia__penta6", "tibiotalar_talus__penta6"]
FIB = ["fibula_cort__hex8", "fibula_spon__tet4"]
FOOT = ["talus_cort__hex8", "talus_spon__tet4", "calcaneus_cort__hex8", "calcaneus_spon__tet4",
        "navicular_spon__hex8"]
LIG = sorted([k for k in EL if k.startswith("lig_")])
ORDER = TIBIA + CART + FIB + FOOT + LIG

def robust_grad(k):
    Z = P[EL[k]["conn"]]; W = U[EL[k]["conn"]]
    ne, nn, _ = Z.shape; g = np.zeros(ne)
    for i in range(nn):
        for j in range(i + 1, nn):
            dx = np.linalg.norm(Z[:, i] - Z[:, j], axis=1)
            du = np.linalg.norm(W[:, i] - W[:, j], axis=1)
            g = np.maximum(g, np.where(dx > 1e-9, du / np.maximum(dx, 1e-12), 0.0))
    return g

print("=" * 132)
print("TASK 3  ---  per-domain deformation / gradient / energy")
print("=" * 132)
print(" domain                          nnd   dzmean    dzstd   dzmin   dzmax   |du|ptp |  gradmax  grad_p95 | maxEzz  |   U(N.mm)  share | u_eff(mm)")
tot_per = {}
GR = {}
for k in ORDER:
    conn = EL[k]["conn"]; nd = np.unique(conn.ravel())
    uz = U[nd, 2]; du = U[conn]
    dzptp = np.ptp(uz)
    H = DATA[k][2]; eps = DATA[k][3]
    Hn = np.linalg.norm(H, axis=(1, 2))
    g = robust_grad(k); GR[k] = g
    ezz = np.abs(eps[:, 2, 2])
    Ut = Udom[k].sum()
    print(" %-30s %4d %+8.4f %8.2e %+7.3f %+7.3f %8.4f | %8.4f %8.4f | %.5f | %9.2f %6.2f%% | %7.4f" % (
        k, len(nd), uz.mean(), uz.std(ddof=0), uz.min(), uz.max(), dzptp,
        g.max(), np.percentile(g, 95), ezz.max(), Ut, 100 * Ut / UTOT, 2 * Ut / F_TOT))
    tot_per[k] = dict(dzstd=uz.std(ddof=0), dzptp=dzptp, gmax=g.max(), Ut=Ut, nnd=len(nd))

print("\n group totals:")
for nm, ks in [("TIBIA", TIBIA), ("CARTILAGE", CART), ("FIBULA", FIB), ("FOOT BONES", FOOT), ("LIGAMENTS", LIG)]:
    U_ = sum(Udom[k].sum() for k in ks)
    nn = len(np.unique(np.concatenate([EL[k]["conn"].ravel() for k in ks])))
    print("   %-12s domains=%2d nodes=%4d  U=%8.2f N.mm  share=%5.1f%%  u_eff=%.4f mm" % (
        nm, len(ks), nn, U_, 100 * U_ / UTOT, 2 * U_ / F_TOT))
print("   TOTAL element energy %.2f  vs external work 0.5*F*u = %.2f N.mm" % (UTOT, 0.5 * F_TOT * 1.6))

print()
print("=" * 132)
print("TASK 3b  ---  where does uz change?  node groups by uz value")
print("=" * 132)
allc = np.unique(np.concatenate([EL[k]["conn"].ravel() for k in ORDER]))
def groupof(k):
    if k in TIBIA: return "tibia"
    if k in CART: return "cartilage"
    if k in FIB: return "fibula"
    if k in FOOT: return "foot"
    return "ligament"
gmap = {}
for k in ORDER:
    for i in EL[k]["conn"].ravel():
        gmap.setdefault(int(i), set()).add(groupof(k))
for lab, lo, hi in [("uz <-1.55 (BC value)", -9, -1.55), ("-1.55..-1.45", -1.55, -1.45),
                    ("-1.45..-1.0", -1.45, -1.0), ("-1.0..-0.75", -1.0, -0.75),
                    ("-0.75..-0.70 (talus level)", -0.75, -0.70), ("-0.70..-0.05", -0.70, -0.05),
                    ("-0.05..0.0 (fixed)", -0.05, 1e-9)]:
    m = (U[:, 2] > lo) & (U[:, 2] <= hi)
    ii = np.where(m)[0]
    if len(ii) == 0:
        print("  %-28s 0 nodes" % lab); continue
    from collections import Counter
    c = Counter()
    for i in ii:
        c["+".join(sorted(gmap.get(int(i), {"<none>"})))] += 1
    print("  %-28s %4d nodes  z=[%7.2f,%7.2f]  groups=%s" % (
        lab, len(ii), P[ii, 2].min(), P[ii, 2].max(), dict(c.most_common(6))))
