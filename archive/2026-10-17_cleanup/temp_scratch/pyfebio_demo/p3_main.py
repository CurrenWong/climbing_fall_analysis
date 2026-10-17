"""p3_main.py - who absorbs the 1.6 mm?  (read-only measurement)"""
import h5py, numpy as np, re, json

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
    EL[name] = dict(type=typ, eid=np.array(ids), conn=np.array(conn) - 1)  # -> 0-based

with h5py.File(H5, "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(np.float64)
    U = f["states/20/node_data/displacement/all"][:].astype(np.float64)
    R = f["states/20/node_data/reaction forces/all"][:].astype(np.float64)
    est = f["states/20/mesh/element_state"][:]
    h5dom = {k: f["meshes/0/domains"][k][:] for k in f["meshes/0/domains"].keys()}
    SIG = {}
    ligmap = {}
    for k in EL:
        if k in h5dom:
            SIG[k] = np.asarray(f["states/20/element_data/stress/%s" % k][:], float)
        else:
            want = set(EL[k]["eid"].tolist())
            for u in h5dom:
                if u.startswith("__unnamed_") and set(h5dom[u][:, 0].tolist()) == want:
                    ligmap[k] = u
                    SIG[k] = np.asarray(f["states/20/element_data/stress/%s" % u][:], float)
                    break
            if k not in SIG:
                raise KeyError("no stress group for " + k)
print("element_state unique:", np.unique(est), "counts:", np.bincount(est.astype(int)))
print("ligament stress-group map:", {k: v for k, v in ligmap.items()})

# ---------- element volumes (isoparametric, order as in .feb) ----------
def vol_and_grad(name):
    d = EL[name]; typ = d["type"]; Z = P[d["conn"]]           # (ne, nn, 3)
    ne, nn, _ = Z.shape
    def _J(dN, Zloc):
        nc = Zloc.shape[2]
        dNb = np.broadcast_to(dN, (Zloc.shape[0], dN.shape[0], dN.shape[1]))
        return np.einsum("eni,ena->eia", dNb, Zloc)
    if typ == "tet4":
        a = Z[:, 1] - Z[:, 0]; b = Z[:, 2] - Z[:, 0]; c = Z[:, 3] - Z[:, 0]
        V = np.abs(np.einsum("ij,ij->i", np.cross(a, b), c)) / 6.0
    elif typ == "hex8":
        S = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
             (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
        gp, w = np.polynomial.legendre.leggauss(2)
        V = np.zeros(ne)
        for xi, wx in zip(gp, w):
            for et, wy in zip(gp, w):
                for ze, wz in zip(gp, w):
                    dN = np.array([[0.125 * s1 * (1 - s2 * et) * (1 - s3 * ze),
                                    0.125 * (1 - s1 * xi) * s2 * (1 - s3 * ze),
                                    0.125 * (1 - s1 * xi) * (1 - s2 * et) * s3] for s1, s2, s3 in S])
                    V += np.linalg.det(_J(dN, Z)) * wx * wy * wz
    elif typ == "penta6":
        tri = [(1 / 6, 1 / 6, 2 / 3), (1 / 6, 2 / 3, 1 / 6), (2 / 3, 1 / 6, 1 / 6)]
        gpz, wz = np.polynomial.legendre.leggauss(2)
        dL = np.array([[1, 0], [0, 1], [-1, -1]], float)
        V = np.zeros(ne)
        for (L1, L2, L3) in tri:
            for ze, w in zip(gpz, wz):
                wb = (1 - ze) / 2; wt = (1 + ze) / 2
                dN = np.zeros((nn, 3))
                for i, L in enumerate([L1, L2, L3]):
                    dN[i, :2] = dL[i] * wb; dN[i, 2] = L * (-0.5)
                    dN[i + 3, :2] = dL[i] * wt; dN[i + 3, 2] = L * (0.5)
                V += np.linalg.det(_J(dN, Z)) * (1 / 6) * w      # tri area weight = 1/6
    elif typ == "quad4":
        gp, w = np.polynomial.legendre.leggauss(2)
        A = np.zeros(ne)
        for xi, wx in zip(gp, w):
            for et, wy in zip(gp, w):
                dN = np.array([[-(1 - et), -(1 - xi)], [(1 - et), -(1 + xi)],
                               [(1 + et), (1 + xi)], [-(1 + et), (1 - xi)]]) / 4.0
                J2 = _J(dN, Z[:, :, :2])
                A += np.abs(np.linalg.det(J2)) * wx * wy
        V = A * 2.0                                          # shell thickness 2.0 mm
    else:
        raise ValueError(typ)
    # least-squares displacement gradient per element
    Xc = Z.mean(1); Uc = U[d["conn"]].mean(1)
    Xr = Z - Xc[:, None, :]; Ur = U[d["conn"]] - Uc[:, None, :]
    XtX = np.einsum("nki,nkj->nij", Xr, Xr)
    XtU = np.einsum("nki,nkj->nij", Xr, Ur)
    H = np.linalg.solve(XtX, XtU)                            # H[n,i,j] = du_i/dx_j
    eps = 0.5 * (H + np.transpose(H, (0, 2, 1)))
    return typ, V, H, eps

def voigt6(s):
    return np.stack([s[:, 0], s[:, 1], s[:, 2], s[:, 3], s[:, 4], s[:, 5]], 1)

RES = {}
order = ["tibia_cort__hex8", "tibia_cort__tet4", "tibia_spon__tet4", "tibia_marrow__tet4",
         "tibiotalar_tibia__penta6", "tibiotalar_talus__penta6",
         "fibula_cort__hex8", "fibula_spon__tet4",
         "talus_cort__hex8", "talus_spon__tet4",
         "calcaneus_cort__hex8", "calcaneus_spon__tet4", "navicular_spon__hex8"]
ligorder = [k for k in EL if k.startswith("lig_")]

print("\n=== element volume sanity ===")
totV = 0.0
for k in order + ligorder:
    typ, V, H, eps = vol_and_grad(k)
    RES[k] = dict(typ=typ, V=V, H=H, eps=eps, conn=EL[k]["conn"], n=len(V))
    totV += V.sum()
    print("  %-28s %-7s n=%5d  V_tot=%10.2f mm^3  V[min,max]=[%.4f,%.2f]  negV=%d" % (
        k, typ, len(V), V.sum(), V.min(), V.max(), int((V <= 0).sum())))
print("  TOTAL volume = %.1f mm^3" % totV)
