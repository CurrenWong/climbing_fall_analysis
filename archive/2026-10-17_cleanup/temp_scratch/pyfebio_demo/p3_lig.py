import pickle, numpy as np, h5py, re
D = pickle.load(open(r"temp/pyfebio_demo/p3_cache.pkl", "rb"))
SIG, NS, EL = D["SIG"], D["NS"], D["EL"]
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)
    R = f["states/20/node_data/reaction forces/all"][:].astype(float)

TIBIA = ["tibia_cort__hex8", "tibia_cort__tet4", "tibia_spon__tet4", "tibia_marrow__tet4"]
FIB = ["fibula_cort__hex8", "fibula_spon__tet4"]
FOOT = ["talus_cort__hex8", "talus_spon__tet4", "calcaneus_cort__hex8", "calcaneus_spon__tet4",
        "navicular_spon__hex8"]
LIG = sorted([k for k in EL if k.startswith("lig_")])

body = {}
for pre, g in [("tibia", "TIB"), ("fibula", "FIB")]:
    for k in TIBIA + FIB:
        for i in EL[k]["conn"].ravel():
            body.setdefault(int(i), set()).add("TIB" if k.startswith("tibia") else "FIB")
for k in FOOT:
    for i in EL[k]["conn"].ravel():
        body.setdefault(int(i), set()).add("FOOT")
fixed = np.zeros(len(P), bool); fixed[NS["calcaneus_bottom_ns"] - 1] = True
presc = np.zeros(len(P), bool); presc[NS["tibia_top_ns"] - 1] = True; presc[NS["fibula_top_ns"] - 1] = True

print("=" * 120)
print("LIGAMENT TOPOLOGY + AXIAL FORCE  (which bodies does each ligament bridge?)")
print("=" * 120)
def voigt(s):
    return np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])

tot = []
for k in LIG:
    conn = EL[k]["conn"]; Z = P[conn]; W = U[conn]; nq = conn.shape[1]
    ends = []
    for i in range(nq):
        b = set()
        for nd in conn[:, i]:
            b |= body.get(int(nd), set())
        ends.append(b)
    left = set().union(*[ends[i] for i in [0, 3]])
    right = set().union(*[ends[i] for i in [1, 2]])
    # axial direction: longer of the two opposite-edge pairs
    e_a = (Z[:, 1] - Z[:, 0] + Z[:, 2] - Z[:, 3]) / 2
    e_b = (Z[:, 2] - Z[:, 1] + Z[:, 3] - Z[:, 0]) / 2
    la = np.linalg.norm(e_a, axis=1); lb = np.linalg.norm(e_b, axis=1)
    usea = la >= lb
    nvec = np.where(usea[:, None], e_a, e_b)
    nvec = nvec / np.linalg.norm(nvec, axis=1, keepdims=True)
    wid = np.where(usea, lb, la)
    S = SIG[k]; F = np.zeros(len(conn)); sa = np.zeros(len(conn))
    for e in range(len(conn)):
        st = voigt(S[e]); nv = nvec[e]
        sa[e] = nv @ st @ nv
        F[e] = sa[e] * wid[e] * 2.0
    nfx = np.zeros(len(conn)); nfy = np.zeros(len(conn)); nfz = np.zeros(len(conn))
    # force vector along +n (tension pulls the two ends apart)
    fx = F * nvec[:, 0]; fy = F * nvec[:, 1]; fz = F * nvec[:, 2]
    tot.append((k, F, nvec, wid))
    print(" %-34s n=%2d bridge:%-22s  F_ax[N]: mean=%+7.2f min=%+7.2f max=%+7.2f  sumFz=%+8.2f  |F|sum=%7.2f" % (
        k, len(conn), "%s<->%s" % ("+".join(sorted(left)) or "?", "+".join(sorted(right)) or "?"),
        F.mean(), F.min(), F.max(), fz.sum(), np.abs(F).sum()))

mm = {}
for k, F, nvec, wid in tot:
    mm[k] = F
allF = np.concatenate([F for _, F, _, _ in tot])
print("   Σ over all 135 ligament elements: ΣF_ax = %+.2f N ; Σ|F_ax| = %.2f N" % (allF.sum(), np.abs(allF).sum()))

# net vertical force that ligaments apply to the FOOT body
net_fz_on_foot = 0.0
for k, F, nvec, wid in tot:
    conn = EL[k]["conn"]
    for i in range(conn.shape[1]):
        bg = set()
        for nd in conn[:, i]:
            bg |= body.get(int(nd), set())
        if bg == {"FOOT"} or ("FOOT" in bg and "TIB" not in bg and "FIB" not in bg):
            # this end is on the foot -> force on foot = +F*n if n points away from foot side
            pass
print("\n  (sign analysis is ambiguous element-wise; using net reaction instead)")

print()
print("=" * 120)
print("NET FORCE BALANCE  ---  who reacts 592.6 N?")
print("=" * 120)
for nm, grp in [("TIBIA domains", TIBIA), ("FIBULA domains", FIB), ("FOOT domains", FOOT)]:
    nd = np.unique(np.concatenate([EL[k]["conn"].ravel() for k in grp]))
    print("  %-16s nodes=%4d  ΣRz=%+9.4f  ΣR= (%+.4f,%+.4f,%+.4f)" % (
        nm, len(nd), R[nd, 2].sum(), R[nd, 0].sum(), R[nd, 1].sum(), R[nd, 2].sum()))
print("  fixed set   (calcaneus_bottom 194): ΣRz = %+.6f  Σ|R| = %.6f" % (
    R[NS["calcaneus_bottom_ns"] - 1, 2].sum(), np.abs(R[NS["calcaneus_bottom_ns"] - 1]).sum()))
print("  ALL nodes: ΣRz = %+.4f" % R[:, 2].sum())
print("  z of nodes with |Rz|>1N: %.2f .. %.2f  (calcaneus/talus span z=%.2f .. %.2f)" % (
    P[np.abs(R[:, 2]) > 1, 2].min(), P[np.abs(R[:, 2]) > 1, 2].max(),
    P[np.unique(np.concatenate([EL[k]['conn'].ravel() for k in FOOT])), 2].min(),
    P[np.unique(np.concatenate([EL[k]['conn'].ravel() for k in FOOT])), 2].max()))

print()
print("=" * 120)
print("TASK 4  ---  bending vs uniform compression, tibia top 10 mm")
print("=" * 120)
tn = np.unique(np.concatenate([EL[k]["conn"].ravel() for k in TIBIA]))
zt = P[tn, 2].max()
print("  tibia nodes z=[%.2f, %.2f]; top 10 mm = z in [%.2f, %.2f]" % (
    P[tn, 2].min(), zt, zt - 10, zt))
for lab, m in [("top 10 mm", (P[tn, 2] >= zt - 10)), ("middle 10 mm", (P[tn, 2] >= zt - 30) & (P[tn, 2] < zt - 20)),
               ("bottom 10 mm", (P[tn, 2] <= P[tn, 2].min() + 10))]:
    sub = tn[m]; u = U[sub]
    print("  %-14s n=%3d |ux|: max=%.5f  |uy|: max=%.5f  |uz|: mean=%+.5f std=%.2e ptp=%.5f" % (
        lab, len(sub), np.abs(u[:, 0]).max(), np.abs(u[:, 1]).max(),
        u[:, 2].mean(), u[:, 2].std(ddof=0), np.ptp(u[:, 2])))
print("  --- per-1mm: spread of uz across xy at same z (bending signature) + plane-fit tilt ---")
for zc in np.arange(zt - 0.5, zt - 10.5, -1.0):
    m = (P[tn, 2] <= zc + 0.5) & (P[tn, 2] > zc - 0.5)
    if m.sum() < 4: continue
    sub = tn[m]; u = U[sub]; X = P[sub]
    A_ = np.column_stack([np.ones(len(sub)), X[:, 0] - X[:, 0].mean(), X[:, 1] - X[:, 1].mean()])
    coef, *_ = np.linalg.lstsq(A_, u[:, 2], rcond=None)
    spanx = np.ptp(X[:, 0]); spany = np.ptp(X[:, 1])
    tilt = abs(coef[1]) * spanx + abs(coef[2]) * spany
    print("    z~%.2f n=%3d  uz mean=%+.5f std=%.2e ptp=%.5f | plane fit slope=(%.2e,%.2e) -> dz span from tilt=%.5f (%.2f%% of 1.6)" % (
        zc, len(sub), u[:, 2].mean(), u[:, 2].std(ddof=0), np.ptp(u[:, 2]),
        coef[1], coef[2], tilt, 100 * tilt / 1.6))
print("  --- per-domain transverse displacement ---")
for k in TIBIA + FIB + FOOT + ["tibiotalar_tibia__penta6", "tibiotalar_talus__penta6"]:
    nd = np.unique(EL[k]["conn"].ravel()); u = U[nd]
    print("    %-28s |ux|max=%.5f |uy|max=%.5f  |uxy|max=%.5f   |uz|max=%.5f" % (
        k, np.abs(u[:, 0]).max(), np.abs(u[:, 1]).max(),
        np.linalg.norm(u[:, :2], axis=1).max(), np.abs(u[:, 2]).max()))
