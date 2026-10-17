import pickle, numpy as np
D = pickle.load(open(r"temp/pyfebio_demo/p3_cache.pkl", "rb"))
DATA, SIG, Udom, NS, EL = D["DATA"], D["SIG"], D["Udom"], D["NS"], D["EL"]
import h5py
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)
F_TOT = 592.60; UTOT = 147.90; EXT = 0.5 * F_TOT * 1.6

TIBIA = ["tibia_cort__hex8", "tibia_cort__tet4", "tibia_spon__tet4", "tibia_marrow__tet4"]
CART = ["tibiotalar_tibia__penta6", "tibiotalar_talus__penta6"]
FIB = ["fibula_cort__hex8", "fibula_spon__tet4"]
FOOT = ["talus_cort__hex8", "talus_spon__tet4", "calcaneus_cort__hex8", "calcaneus_spon__tet4",
        "navicular_spon__hex8"]
LIG = [k for k in EL if k.startswith("lig_")]
ORDER = TIBIA + CART + FIB + FOOT + LIG

# ---------- robust per-element gradient (max over node pairs of |du|/|dx|) ----------
def robust_grad(k):
    conn = EL[k]["conn"]; Z = P[conn]; W = U[conn]
    ne, nn, _ = Z.shape
    g = np.zeros(ne)
    for i in range(nn):
        for j in range(i + 1, nn):
            dx = np.linalg.norm(Z[:, i] - Z[:, j], axis=1)
            du = np.linalg.norm(W[:, i] - W[:, j], axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                r = np.where(dx > 1e-9, du / np.maximum(dx, 1e-12), 0.0)
            g = np.maximum(g, r)
    return g

GRAD = {k: robust_grad(k) for k in ORDER}

print("=" * 108)
print("TASK 1/2  ---  BC band + axial displacement profile")
print("=" * 108)
for nm in ["tibia_top_ns", "fibula_top_ns"]:
    ix = NS[nm] - 1
    print("  %-14s n=%3d  z=[%.3f, %.3f] span=%.2f mm | ux,uy all 0 | uz: mean=%.6f std=%.2e min=%.6f max=%.6f" % (
        nm, len(ix), P[ix, 2].min(), P[ix, 2].max(), np.ptp(P[ix, 2]),
        U[ix, 2].mean(), U[ix, 2].std(ddof=0), U[ix, 2].min(), U[ix, 2].max()))

ztop = P[NS["tibia_top_ns"] - 1, 2].max()
print("\n  z_top(BC band) = %.3f ; model z=[%.2f, %.2f]" % (ztop, P[:, 2].min(), P[:, 2].max()))

tib_nodes = np.unique(np.concatenate([EL[k]["conn"].ravel() for k in TIBIA]))
print("  tibia domain nodes: %d  z=[%.2f, %.2f]" % (len(tib_nodes), P[tib_nodes, 2].min(), P[tib_nodes, 2].max()))

def zprofile(mask_nodes, z0, z1, step=1.0, tag=""):
    z = P[mask_nodes, 2]; uz = U[mask_nodes, 2]
    rows = []
    edges = np.arange(z0, z1 - 1e-9, -step)
    for e0 in edges:
        m = (z <= e0 + 1e-9) & (z > e0 - step + 1e-9)
        if m.sum() == 0:
            rows.append((e0, e0 - step, 0, np.nan, np.nan, np.nan, np.nan)); continue
        rows.append((e0, e0 - step, int(m.sum()), uz[m].mean(), np.median(uz[m]), uz[m].min(), uz[m].max()))
    print("\n  [%s] z-slab (mm)      n   mean_uz    med_uz    min_uz    max_uz   drop_uz  pct_of_1.6" % tag)
    prev = None
    for (a, b, n, mu, md, mn, mx) in rows:
        if n == 0:
            print("   %7.2f..%6.2f %5d      --" % (a, b, n)); prev = None; continue
        dr = "" if prev is None else "%8.4f  %6.2f%%" % (prev - mu, 100 * (prev - mu) / 1.6)
        print("   %7.2f..%6.2f %5d  %+.4f  %+.4f  %+.4f  %+.4f  %s" % (a, b, n, mu, md, mn, mx, dr))
        prev = mu
    return rows

rr = zprofile(tib_nodes, ztop + 1, P[tib_nodes, 2].min() - 1, 1.0, "tibia domains only, 1 mm")
zprofile(np.arange(len(P)), ztop + 1, P[:, 2].min() - 1, 1.0, "ALL nodes, 1 mm")

# 5-mm windows with largest absorption (tibia-only)
rr2 = [r for r in rr if r[2] > 0]
print("\n  --- 5-mm windows ranked by absorbed |dz| (tibia domains) ---")
best = []
for i in range(len(rr2) - 4):
    w = rr2[i:i + 5]
    if any(x[2] == 0 for x in w):
        continue
    drop = w[0][3] - w[4][3]
    best.append((drop, w[0][0], w[-1][1], sum(x[2] for x in w)))
for d, a, b, n in sorted(best, reverse=True)[:6]:
    print("   z %7.2f -> %6.2f  drop=%.4f mm  (%.1f%% of 1.6)  nodes=%d" % (a, b, d, 100 * d / 1.6, n))
print("  sum of drops over all slabs = %.4f mm" % (rr2[0][3] - rr2[-1][3]))
