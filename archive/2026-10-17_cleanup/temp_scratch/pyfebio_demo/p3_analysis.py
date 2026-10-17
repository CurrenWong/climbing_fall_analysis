"""p3_analysis.py - reads p3_main results by re-computing; prints the deliverable tables."""
import h5py, numpy as np, re

H5 = r"temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = r"temp/thums/febio_lig/ankle_lig_WK.feb"
txt = open(FEB, errors="replace").read()
NS = {}
for m in re.finditer(r'<NodeSet name="([^"]+)">([^<]*)</NodeSet>', txt):
    NS[m.group(1)] = np.array([int(x) for x in re.split(r"[,\s]+", m.group(2).strip()) if x])
EL = {}
for m in re.finditer(r'<Elements type="(\w+)" name="([^"]+)">(.*?)</Elements>', txt, re.S):
    typ, name, body = m.group(1), m.group(2), m.group(3)
    ids, conn = [], []
    for e in re.finditer(r'<elem id="(\d+)">([^<]*)</elem>', body):
        ids.append(int(e.group(1)))
        conn.append([int(x) for x in re.split(r"[,\s]+", e.group(2).strip()) if x])
    EL[name] = dict(type=typ, eid=np.array(ids), conn=np.array(conn) - 1)

MAT = {"tibia_cort__hex8": (7300, .3), "tibia_cort__tet4": (7300, .3), "tibia_spon__tet4": (7300, .3),
       "tibia_marrow__tet4": (7300, .3), "fibula_cort__hex8": (7300, .3), "fibula_spon__tet4": (7300, .3),
       "talus_cort__hex8": (7300, .3), "talus_spon__tet4": (7300, .3), "calcaneus_cort__hex8": (7300, .3),
       "calcaneus_spon__tet4": (7300, .3), "navicular_spon__hex8": (7300, .3),
       "tibiotalar_tibia__penta6": (10, .4), "tibiotalar_talus__penta6": (10, .4)}
for k in list(EL):
    if k.startswith("lig_"):
        MAT[k] = (260, .4)
TIBIA = ["tibia_cort__hex8", "tibia_cort__tet4", "tibia_spon__tet4", "tibia_marrow__tet4"]
FOOT = ["talus_cort__hex8", "talus_spon__tet4", "calcaneus_cort__hex8", "calcaneus_spon__tet4",
        "navicular_spon__hex8"]

with h5py.File(H5, "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)
    R = f["states/20/node_data/reaction forces/all"][:].astype(float)
    h5dom = {k: f["meshes/0/domains"][k][:] for k in f["meshes/0/domains"].keys()}
    SIG = {}
    for k in EL:
        g = k if k in h5dom else None
        if g is None:
            want = set(EL[k]["eid"].tolist())
            for u_ in h5dom:
                if u_.startswith("__unnamed_") and set(h5dom[u_][:, 0].tolist()) == want:
                    g = u_; break
        SIG[k] = np.asarray(f["states/20/element_data/stress/%s" % g][:], float)

def _J(dN, Z):
    dNb = np.broadcast_to(dN, (Z.shape[0], dN.shape[0], dN.shape[1]))
    return np.einsum("eni,ena->eia", dNb, Z)

def vol_and_grad(name):
    d = EL[name]; typ = d["type"]; Z = P[d["conn"]]; ne, nn, _ = Z.shape
    if typ == "tet4":
        V = np.abs(np.einsum("ij,ij->i", np.cross(Z[:, 1] - Z[:, 0], Z[:, 2] - Z[:, 0]),
                             Z[:, 3] - Z[:, 0])) / 6.0
    elif typ == "hex8":
        S = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1), (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
        gp, w = np.polynomial.legendre.leggauss(2); V = np.zeros(ne)
        for xi, wx in zip(gp, w):
            for et, wy in zip(gp, w):
                for ze, wz in zip(gp, w):
                    dN = np.array([[.125 * s1 * (1 - s2 * et) * (1 - s3 * ze),
                                    .125 * (1 - s1 * xi) * s2 * (1 - s3 * ze),
                                    .125 * (1 - s1 * xi) * (1 - s2 * et) * s3] for s1, s2, s3 in S])
                    V += np.linalg.det(_J(dN, Z)) * wx * wy * wz
    elif typ == "penta6":
        tri = [(1 / 6, 1 / 6, 2 / 3), (1 / 6, 2 / 3, 1 / 6), (2 / 3, 1 / 6, 1 / 6)]
        gpz, wz = np.polynomial.legendre.leggauss(2); dL = np.array([[1, 0], [0, 1], [-1, -1]], float)
        V = np.zeros(ne)
        for (L1, L2, L3) in tri:
            for ze, w in zip(gpz, wz):
                wb = (1 - ze) / 2; wt = (1 + ze) / 2; dN = np.zeros((nn, 3))
                for i, L in enumerate([L1, L2, L3]):
                    dN[i, :2] = dL[i] * wb; dN[i, 2] = L * (-.5)
                    dN[i + 3, :2] = dL[i] * wt; dN[i + 3, 2] = L * .5
                V += np.linalg.det(_J(dN, Z)) * (1 / 6) * w
    elif typ == "quad4":
        gp, w = np.polynomial.legendre.leggauss(2); A = np.zeros(ne)
        for xi, wx in zip(gp, w):
            for et, wy in zip(gp, w):
                dN = np.array([[-(1 - et), -(1 - xi)], [(1 - et), -(1 + xi)],
                               [(1 + et), (1 + xi)], [-(1 + et), (1 - xi)]]) / 4.0
                A += np.abs(np.linalg.det(_J(dN, Z[:, :, :2]))) * wx * wy
        V = A * 2.0
    Xc = Z.mean(1); Uc = U[d["conn"]].mean(1)
    Xr = Z - Xc[:, None, :]; Ur = U[d["conn"]] - Uc[:, None, :]
    H = np.linalg.solve(np.einsum("nki,nkj->nij", Xr, Xr), np.einsum("nki,nkj->nij", Xr, Ur))
    eps = 0.5 * (H + np.transpose(H, (0, 2, 1)))
    return typ, V, H, eps, Xc

DATA = {}
for k in EL:
    DATA[k] = vol_and_grad(k)

# ---------------- Voigt order check ----------------
print("=== Voigt-order check (reconstruct sigma from eps, E=7300 nu=0.3) ===")
k = "tibia_cort__hex8"
E, nu = MAT[k]; lam = E * nu / ((1 + nu) * (1 - 2 * nu)); mu = E / (2 * (1 + nu))
eps = DATA[k][3]; sig = SIG[k]
tr = np.trace(eps, axis1=1, axis2=2)
pred = lam * tr[:, None, None] * np.eye(3)[None] + 2 * mu * eps
pv = np.stack([pred[:, 0, 0], pred[:, 1, 1], pred[:, 2, 2], pred[:, 0, 1], pred[:, 1, 2], pred[:, 0, 2]], 1)
for i in range(6):
    cc = [np.corrcoef(pv[:, i], sig[:, j])[0, 1] for j in range(6)]
    print("  pred slot %d  corr with reported: %s   best=%d" % (i, np.round(cc, 3), int(np.argmax(np.abs(cc)))))
print("  reported |sigma| mean:", np.round(np.abs(sig).mean(0), 4))
print("  predicted|sigma| mean:", np.round(np.abs(pv).mean(0), 4))

F_TOT = 592.60
u_BC = 1.6
print("\n=== energy calibration ===")
Utot = 0.0; Udom = {}
for k in EL:
    E, nu = MAT[k]; V = DATA[k][1]; eps = DATA[k][3]; sig = SIG[k]
    # Voigt: sxx,syy,szz,sxy,syz,sxz
    w = (sig[:, 0] * eps[:, 0, 0] + sig[:, 1] * eps[:, 1, 1] + sig[:, 2] * eps[:, 2, 2]
         + 2 * sig[:, 3] * eps[:, 0, 1] + 2 * sig[:, 4] * eps[:, 1, 2] + 2 * sig[:, 5] * eps[:, 0, 2])
    Udom[k] = 0.5 * w * V
    Utot += Udom[k].sum()
print("  sum 0.5*sigma:eps*V = %.1f N.mm   (external work 0.5*F*u = %.1f N.mm)" % (
    Utot, 0.5 * F_TOT * u_BC))

np.savez(r"temp/pyfebio_demo/p3_cache.npz", **{"U_" + k.replace("__", "_"): Udom[k] for k in EL})
import pickle
with open(r"temp/pyfebio_demo/p3_cache.pkl", "wb") as fh:
    pickle.dump(dict(DATA={k: (DATA[k][0], DATA[k][1], DATA[k][2], DATA[k][3], DATA[k][4]) for k in EL},
                     SIG=SIG, Udom=Udom, NS=NS, EL={k: EL[k] for k in EL}), fh)
print("  cache written")
