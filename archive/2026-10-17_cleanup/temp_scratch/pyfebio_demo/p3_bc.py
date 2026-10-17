import h5py, numpy as np

H5 = r"temp/thums/febio_lig/ankle_lig_WK_stress.h5"

with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"]
    nid = nd["id"][:]
    P = np.column_stack([nd["x"][:], nd["y"][:], nd["z"][:]]).astype(np.float64)
    nsets = {k: f["meshes/0/nodesets"][k][:] for k in f["meshes/0/nodesets"].keys()}
    doms = {k: f["meshes/0/domains"][k][:] for k in f["meshes/0/domains"].keys()}
    U = {}
    for s in range(21):
        U[s] = f["states/%d/node_data/displacement/all" % s][:].astype(np.float64)
    R = f["states/20/node_data/reaction forces/all"][:].astype(np.float64)

id2idx = {int(i): j for j, i in enumerate(nid)}
def sel(ids):
    return np.array([id2idx[int(i)] for i in ids], dtype=int)

print("=== displacement history (mean over sets) ===")
for s in range(21):
    ix = sel(nsets["tibia_top_ns"])
    jx = sel(nsets["fibula_top_ns"])
    print(" state %2d  tibia_top uz mean=%+.6f  fibula_top uz mean=%+.6f   max|uz| all=%+.6f" % (
        s, U[s][ix, 2].mean(), U[s][jx, 2].mean(), U[s][:, 2].min()))

print()
for name in ["tibia_top_ns", "fibula_top_ns", "talus_low_ns", "calcaneus_bottom_ns"]:
    ix = sel(nsets[name])
    u = U[20][ix]
    print("### %s  n=%d  z_ref span=%.4f (%.3f .. %.3f)" % (
        name, len(ix), P[ix, 2].max() - P[ix, 2].min(), P[ix, 2].min(), P[ix, 2].max()))
    for j, c in enumerate("xyz"):
        print("   du%s : mean=%+.6f std=%.3e min=%+.6f max=%+.6f  ptp=%.4e" % (
            c, u[:, j].mean(), u[:, j].std(), u[:, j].min(), u[:, j].max(), np.ptp(u[:, j])))
    print("   |u| max=%.6f" % np.linalg.norm(u, axis=1).max())
    rr = R[ix]
    print("   sum R = (%.3f, %.3f, %.3f) N" % tuple(rr.sum(0)))
    print()

print("global sum R = (%.3f, %.3f, %.3f) N" % tuple(R.sum(0)))
