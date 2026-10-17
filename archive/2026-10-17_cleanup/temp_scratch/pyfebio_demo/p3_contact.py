import pickle, numpy as np, h5py
D = pickle.load(open(r"temp/pyfebio_demo/p3_cache.pkl", "rb"))
EL = D["EL"]
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)
    S = {k: f["meshes/0/surfaces"][k][:] for k in f["meshes/0/surfaces"].keys()}
Pd = P + U

def surf_normals(s, want=None):
    f = s[:, 1:][:, :3]                      # use first 3 nodes of each face as a triangle
    p0, p1, p2 = P[f[:, 0]], P[f[:, 1]], P[f[:, 2]]
    n = np.cross(p1 - p0, p2 - p0)
    n = n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    # orient outward = away from the centroid of the mesh region
    c = p0.mean(0)
    flip = ((p0 - c) * n).sum(1) < 0
    n[flip] *= -1
    nodes = np.unique(s[:, 1:][s[:, 1:] >= 0].ravel())
    Nn = {}
    for t in range(len(f)):
        for i in f[t]:
            Nn.setdefault(int(i), []).append(n[t])
    Nn = {i: np.mean(v, 0) / np.linalg.norm(np.mean(v, 0)) for i, v in Nn.items()}
    gm = np.mean(list(Nn.values()), 0); gm /= np.linalg.norm(gm)
    for i in nodes:
        Nn.setdefault(int(i), gm)
    return nodes, Nn

print("=" * 120)
print("CONTACT KINEMATICS: tibial cartilage outer surface vs talar cartilage outer surface")
print("=" * 120)
ia, Na = surf_normals(S["surf_tibiotalar_tibia"])
ib, Nb = surf_normals(S["surf_tibiotalar_talus"])
d2 = ((P[ia][:, None, :] - P[ib][None, :, :]) ** 2).sum(2)
j = d2.argmin(1)
gref = np.sqrt(d2[np.arange(len(ia)), j])
m = gref < 3.0
print("  pairs: %d (contact zone gap_ref<3mm: %d)" % (len(ia), m.sum()))
nrm = np.array([Na[i] for i in ia[m]])
print("  reference normal n (mean over pairs): [%+.3f %+.3f %+.3f]  -> angle to z = %.1f deg" % (
    nrm.mean(0)[0], nrm.mean(0)[1], nrm.mean(0)[2],
    np.degrees(np.arccos(abs(nrm.mean(0)[2])))))
print("  per-node |n_z| median = %.3f  (=> surface obliquity to the load axis)" % np.median(np.abs(nrm[:, 2])))
rows = []
for kk, i in enumerate(ia):
    if not m[kk]:
        continue
    n = Na[i]
    du = U[ib[j[kk]]] - U[i]                 # relative disp of talar surface w.r.t. tibial surface
    dn = du @ n
    dt = np.linalg.norm(du - dn * n)
    rows.append((gref[kk], dn, dt, du[2], du @ (n * np.sign(n[2]))))
rows = np.array(rows)
print("  relative displacement (talar w.r.t. tibial), over contact zone:")
print("    normal component  dn : mean=%+.4f  med=%+.4f  min=%+.4f  max=%+.4f" % (
    rows[:, 1].mean(), np.median(rows[:, 1]), rows[:, 1].min(), rows[:, 1].max()))
print("    tangential  |dt|     : mean=%+.4f  med=%+.4f  min=%+.4f  max=%+.4f" % (
    rows[:, 2].mean(), np.median(rows[:, 2]), rows[:, 2].min(), rows[:, 2].max()))
print("    du_z (total axial)   : mean=%+.4f  med=%+.4f" % (rows[:, 3].mean(), np.median(rows[:, 3])))
print("    |dn| vs |dt| medians: %.4f vs %.4f  -> tangential share = %.0f%%" % (
    abs(np.median(rows[:, 1])), np.median(rows[:, 2]),
    100 * np.median(rows[:, 2]) / (abs(np.median(rows[:, 1])) + np.median(rows[:, 2]))))
print("    gap change (gref - gdef) : med=%+.4f mean=%+.4f" % (
    np.median(rows[:, 0] - np.linalg.norm(Pd[ia[m]] - Pd[ib[j[m]]], axis=1)),
    (rows[:, 0] - np.linalg.norm(Pd[ia[m]] - Pd[ib[j[m]]], axis=1)).mean()))

print()
print("=" * 120)
print("BONE-TO-BONE: tibial articular bone surface vs talar bone surface")
print("=" * 120)
ia2, Na2 = surf_normals(S["bonesurf_tibiotalar_tibia"])
ib2, _ = surf_normals(S["bonesurf_tibiotalar_talus"])
d2 = ((P[ia2][:, None, :] - P[ib2][None, :, :]) ** 2).sum(2)
j2 = d2.argmin(1); gref2 = np.sqrt(d2[np.arange(len(ia2)), j2]); m2 = gref2 < 4.0
dn = np.array([ (U[ib2[j2[k]]] - U[ia2[k]]) @ Na2[ia2[k]] for k in range(len(ia2)) if m2[k]])
dt = np.array([ np.linalg.norm((U[ib2[j2[k]]] - U[ia2[k]]) - ((U[ib2[j2[k]]] - U[ia2[k]]) @ Na2[ia2[k]]) * Na2[ia2[k]])
                for k in range(len(ia2)) if m2[k]])
duz = np.array([ (U[ib2[j2[k]]] - U[ia2[k]])[2] for k in range(len(ia2)) if m2[k]])
print("  n=%d pairs; normal dn med=%+.4f mean=%+.4f | tangential med=%.4f mean=%.4f | du_z med=%+.4f" % (
    m2.sum(), np.median(dn), dn.mean(), np.median(dt), dt.mean(), np.median(duz)))

print()
print("=" * 120)
print("CARTILAGE SURFACE MOTION (inner/bone side vs outer/joint side)")
print("=" * 120)
for k, tag in [("tibiotalar_tibia__penta6", "tibial cartilage"), ("tibiotalar_talus__penta6", "talar cartilage")]:
    conn = EL[k]["conn"]
    inner = np.unique(conn[:, :3]); outer = np.unique(conn[:, 3:])
    print("  %s : inner(bone side) n=%d  uz mean=%+.4f [%+.4f,%+.4f]" % (
        tag, len(inner), U[inner, 2].mean(), U[inner, 2].min(), U[inner, 2].max()))
    print("  %s : outer(joint side) n=%d  uz mean=%+.4f [%+.4f,%+.4f]  |uxy|max=%.4f" % (
        tag, len(outer), U[outer, 2].mean(), U[outer, 2].min(), U[outer, 2].max(),
        np.linalg.norm(U[outer, :2], axis=1).max()))
    print("     inner-outer uz difference: med=%+.4f  -> layer stays intact (small thickness change)" % (
        np.median(U[inner, 2].mean() - U[outer, 2].mean())))
    # theta: distance between inner/outer triangles
    ti = P[conn[:, :3]].mean(1); to = P[conn[:, 3:]].mean(1)
    L0 = np.linalg.norm(to - ti, axis=1)
    L1 = np.linalg.norm((to + U[conn[:, 3:]].mean(1)) - (ti + U[conn[:, :3]].mean(1)), axis=1)
    print("     layer thickness: med %.4f -> %.4f (d=%+.4f, %+.2f%%)" % (
        np.median(L0), np.median(L1), np.median(L1 - L0), 100 * np.median((L1 - L0) / L0)))
