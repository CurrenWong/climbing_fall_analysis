import pickle, numpy as np, h5py
D = pickle.load(open(r"temp/pyfebio_demo/p3_cache.pkl", "rb"))
EL, NS, DATA = D["EL"], D["NS"], D["DATA"]
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)
    SURF = {k: f["meshes/0/surfaces"][k][:] for k in f["meshes/0/surfaces"].keys()}
Pd = P + U

print("=" * 118)
print("GEOMETRY / z-ORDER OF BODIES")
print("=" * 118)
ORDER = ["tibia_cort__hex8", "tibia_spon__tet4", "tibia_marrow__tet4", "tibiotalar_tibia__penta6",
         "tibiotalar_talus__penta6", "talus_cort__hex8", "talus_spon__tet4",
         "calcaneus_cort__hex8", "calcaneus_spon__tet4", "navicular_spon__hex8",
         "fibula_cort__hex8", "fibula_spon__tet4"]
for k in ORDER:
    nd = np.unique(EL[k]["conn"].ravel())
    print("  %-28s n=%4d  z_ref=[%7.2f,%7.2f]  uz=[%+.4f,%+.4f] mean=%+.4f" % (
        k, len(nd), P[nd, 2].min(), P[nd, 2].max(), U[nd, 2].min(), U[nd, 2].max(), U[nd, 2].mean()))
for k in sorted([x for x in EL if x.startswith("lig_")]):
    nd = np.unique(EL[k]["conn"].ravel())
    print("  %-40s n=%3d z_ref=[%7.2f,%7.2f] uz=[%+.3f,%+.3f]" % (
        k, len(nd), P[nd, 2].min(), P[nd, 2].max(), U[nd, 2].min(), U[nd, 2].max()))

print()
print("=" * 118)
print("CARTILAGE LAYER THICKNESS CHANGE (exact: centroid of tri-1 vs tri-2 of each wedge)")
print("=" * 118)
for k in ["tibiotalar_tibia__penta6", "tibiotalar_talus__penta6"]:
    conn = EL[k]["conn"]; Z = P[conn]; W = U[conn]
    # try both pairings, pick the one giving the smaller |thickness| (the true thin direction)
    res = {}
    for lab, (A, B) in [("0-2|3-5", ([0, 1, 2], [3, 4, 5])), ("0,5..|..", ([0, 2, 4], [1, 3, 5]))]:
        a = Z[:, A].mean(1); b = Z[:, B].mean(1)
        ta = (b - a); t0 = np.linalg.norm(ta, axis=1)
        ua = W[:, A].mean(1); ub = W[:, B].mean(1)
        t1 = (b + ub) - (a + ua); t1n = np.linalg.norm(t1, axis=1)
        res[lab] = (t0, t1n - t0)
    for lab in res:
        t0, dt = res[lab]
        print("  %-28s %-9s thick0: mean=%.3f med=%.3f  | d(thick): mean=%+.5f med=%+.5f min=%+.5f max=%+.5f  strain_med=%+.4f%%" % (
            k, lab, t0.mean(), np.median(t0), dt.mean(), np.median(dt), dt.min(), dt.max(),
            100 * np.median(dt) / np.median(t0)))

print()
print("=" * 118)
print("JOINT SPACE (contact) CLOSURE :  nearest-pair distance, ref vs deformed")
print("=" * 118)
def faces_nodes(s):
    return np.unique(s[:, 1:][s[:, 1:] >= 0].ravel())
for a, b in [("surf_tibiotalar_tibia", "surf_tibiotalar_talus")]:
    ia = faces_nodes(SURF[a]); ib = faces_nodes(SURF[b])
    A = P[ia]; B = P[ib]
    d2 = ((A[:, None, :] - B[None, :, :]) ** 2).sum(2)
    j = d2.argmin(1)
    gref = np.sqrt(d2[np.arange(len(ia)), j])
    gdef = np.linalg.norm(Pd[ia] - Pd[ib[j]], axis=1)
    dclos = gref - gdef
    m = gref < 3.0
    print("  %s(%d nd) vs %s(%d nd)" % (a, len(ia), b, len(ib)))
    print("    all pairs (n=%d): gap_ref mean=%.4f med=%.4f | closure mean=%+.5f med=%+.5f" % (
        len(ia), gref.mean(), np.median(gref), dclos.mean(), np.median(dclos)))
    print("    contact zone gap_ref<3mm (n=%d): gap_ref mean=%.4f med=%.4f min=%.4f | closure mean=%+.5f med=%+.5f min=%+.5f max=%+.5f" % (
        m.sum(), gref[m].mean(), np.median(gref[m]), gref[m].min(), dclos[m].mean(),
        np.median(dclos[m]), dclos[m].min(), dclos[m].max()))
    print("    -> joint-space gap change (median) = %+.5f mm  (0.873 mm of tibia->talus approach expected)" % np.median(dclos[m]))

print()
print("=" * 118)
print("TIE SLIP (cartilage inner surface vs bone articular surface)")
print("=" * 118)
for cart, bone in [("tiesurf_tibiotalar_tibia", "bonesurf_tibiotalar_tibia"),
                   ("tiesurf_tibiotalar_talus", "bonesurf_tibiotalar_talus")]:
    ia = faces_nodes(SURF[cart]); ib = faces_nodes(SURF[bone])
    A = P[ia]; B = P[ib]
    d2 = ((A[:, None, :] - B[None, :, :]) ** 2).sum(2)
    j = d2.argmin(1)
    gref = np.sqrt(d2[np.arange(len(ia)), j])
    slip = np.linalg.norm(U[ia] - U[ib[j]], axis=1)
    dzz = np.abs(U[ia, 2] - U[ib[j], 2])
    m = gref < 1.0
    print("  %-28s n=%d  |d|ref: med=%.4f | slip |du|: med=%.5f max=%.5f | |duz| med=%.5f max=%.5f" % (
        cart, len(ia), np.median(gref), np.median(slip), slip.max(), np.median(dzz), dzz.max()))
    print("      (pairs <1mm: n=%d  slip med=%.5f max=%.5f)" % (
        m.sum(), np.median(slip[m]) if m.any() else -1, slip[m].max() if m.any() else -1))

print()
print("=" * 118)
print("LIGAMENT LENGTH CHANGE (nodal, exact)")
print("=" * 118)
for k in sorted([x for x in EL if x.startswith("lig_")]):
    conn = EL[k]["conn"]; Z = P[conn]; W = U[conn]
    e_a = (Z[:, 1] - Z[:, 0] + Z[:, 2] - Z[:, 3]) / 2
    e_b = (Z[:, 2] - Z[:, 1] + Z[:, 3] - Z[:, 0]) / 2
    la = np.linalg.norm(e_a, axis=1); lb = np.linalg.norm(e_b, axis=1)
    usea = la >= lb
    e0 = np.where(usea[:, None], e_a, e_b)
    L0 = np.linalg.norm(e0, axis=1)
    # deformed direction using same node pairs
    f_a = (Z[:, 1] + W[:, 1] - Z[:, 0] - W[:, 0] + Z[:, 2] + W[:, 2] - Z[:, 3] - W[:, 3]) / 2
    f_b = (Z[:, 2] + W[:, 2] - Z[:, 1] - W[:, 1] + Z[:, 3] + W[:, 3] - Z[:, 0] - W[:, 0]) / 2
    e1 = np.where(usea[:, None], f_a, f_b)
    L1 = np.linalg.norm(e1, axis=1)
    st = (L1 - L0) / L0
    print("  %-34s n=%2d  L0 med=%.3f  dL med=%+.4f mm  strain med=%+.3f%%  max=%+.3f%%  n_tension=%d/%d" % (
        k, len(conn), np.median(L0), np.median(L1 - L0), 100 * np.median(st), 100 * st.max(),
        int((st > 0).sum()), len(st)))
