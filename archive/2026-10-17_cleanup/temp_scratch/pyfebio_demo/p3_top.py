import pickle, numpy as np, h5py
D = pickle.load(open(r"temp/pyfebio_demo/p3_cache.pkl", "rb"))
EL, SIG = D["EL"], D["SIG"]
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)

rows = []
for k in EL:
    conn = EL[k]["conn"]; Z = P[conn]; W = U[conn]
    ne, nn, _ = Z.shape
    dm = np.zeros(ne); span = np.zeros(ne)
    for i in range(nn):
        for j in range(i + 1, nn):
            du = np.linalg.norm(W[:, i] - W[:, j], axis=1)
            dx = np.linalg.norm(Z[:, i] - Z[:, j], axis=1)
            dm = np.maximum(dm, du)
            span = np.maximum(span, dx * (du >= dm - 1e-12))
    Xc = Z.mean(1)
    for e in range(ne):
        rows.append((dm[e], k, Xc[e, 0], Xc[e, 1], Xc[e, 2], span[e]))
rows.sort(key=lambda r: -r[0])
print("=" * 120)
print("TOP 30 ELEMENTS BY INTERNAL ABSOLUTE DISPLACEMENT DIFFERENCE |du|max  (mm)")
print("=" * 120)
print("   |du|max   domain                        centroid(x,y,z)          el.span[mm]")
for r in rows[:30]:
    print("   %7.4f   %-30s (%7.2f,%7.2f,%7.2f)  %6.3f" % (r[0], r[1], r[2], r[3], r[4], r[5]))

print()
print("=" * 120)
print("SUM of |du|max over elements, grouped by domain  (a proxy for 'displacement absorbed here')")
print("=" * 120)
from collections import defaultdict
agg = defaultdict(lambda: [0.0, 0.0, 0])   # sum, max, count>0.05
for dm, k, *_ in rows:
    agg[k][0] += dm; agg[k][1] = max(agg[k][1], dm); agg[k][2] += int(dm > 0.05)
tot = sum(v[0] for v in agg.values())
print("  %-32s %8s %8s %8s %8s" % ("domain", "sum|du|", "share", "max", "n>0.05"))
for k, v in sorted(agg.items(), key=lambda kv: -kv[1][0]):
    print("  %-32s %8.2f %7.1f%% %8.4f %8d" % (k, v[0], 100 * v[0] / tot, v[1], v[2]))
print("  TOTAL sum|du| = %.1f" % tot)

print()
print("=" * 120)
print("TOP 15 ELEMENTS BY NORMALISED GRADIENT (|du|/|dx|)")
print("=" * 120)
rows2 = []
for k in EL:
    conn = EL[k]["conn"]; Z = P[conn]; W = U[conn]
    ne, nn, _ = Z.shape
    g = np.zeros(ne)
    for i in range(nn):
        for j in range(i + 1, nn):
            dx = np.linalg.norm(Z[:, i] - Z[:, j], axis=1)
            du = np.linalg.norm(W[:, i] - W[:, j], axis=1)
            g = np.maximum(g, np.where(dx > 1e-9, du / np.maximum(dx, 1e-12), 0))
    Xc = Z.mean(1)
    for e in range(ne):
        rows2.append((g[e], k, Xc[e, 0], Xc[e, 1], Xc[e, 2]))
rows2.sort(key=lambda r: -r[0])
for r in rows2[:15]:
    print("   %7.4f   %-30s (%7.2f,%7.2f,%7.2f)" % r)

print()
print("=" * 120)
print("MESH CONNECTIVITY BETWEEN BONES (shared nodes?)")
print("=" * 120)
groups = {"tibia": ["tibia_cort__hex8", "tibia_cort__tet4", "tibia_spon__tet4", "tibia_marrow__tet4"],
          "fibula": ["fibula_cort__hex8", "fibula_spon__tet4"],
          "talus": ["talus_cort__hex8", "talus_spon__tet4"],
          "calcaneus": ["calcaneus_cort__hex8", "calcaneus_spon__tet4"],
          "navicular": ["navicular_spon__hex8"],
          "cart_tib": ["tibiotalar_tibia__penta6"], "cart_tal": ["tibiotalar_talus__penta6"]}
N = {g: set(np.concatenate([EL[k]["conn"].ravel() for k in ks]).tolist()) for g, ks in groups.items()}
for g in groups:
    for g2 in groups:
        if g < g2:
            n = len(N[g] & N[g2])
            if n:
                print("  %-10s & %-10s share %4d nodes" % (g, g2, n))
print("  (ligament elements bridge the rest)")
