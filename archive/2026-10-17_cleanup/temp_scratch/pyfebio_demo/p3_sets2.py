import h5py, numpy as np, re

H5 = r"temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = r"temp/thums/febio_lig/ankle_lig_WK.feb"
txt = open(FEB, errors="replace").read()

NS = {}
for m in re.finditer(r'<NodeSet name="([^"]+)">([^<]*)</NodeSet>', txt):
    name = m.group(1)
    ids = [int(x) for x in re.split(r"[,\s]+", m.group(2).strip()) if x]
    NS[name] = np.array(ids)
print("feb NodeSets:", {k: len(v) for k, v in NS.items()})

with h5py.File(H5, "r") as f:
    nid = f["meshes/0/nodes"]["id"][:].astype(int)
    U = f["states/20/node_data/displacement/all"][:].astype(np.float64)
    R = f["states/20/node_data/reaction forces/all"][:].astype(np.float64)
    ns_h5 = {k: f["meshes/0/nodesets"][k][:] for k in f["meshes/0/nodesets"].keys()}

print("node ids ordered 1..N:", bool(np.array_equal(nid, np.arange(1, len(nid) + 1))))
print("dtypes U,R:", U.dtype, R.dtype)

print("\n=== per-set reaction sums (feb NodeSets) ===")
print("  %-22s %4s %10s %10s %10s | %10s %10s %10s" % ("set", "n", "sumRx", "sumRy", "sumRz", "|Rx|", "|Ry|", "|Rz|"))
for n, ids in NS.items():
    ix = ids - 1
    print("  %-22s %4d %10.4f %10.4f %10.4f | %10.3f %10.3f %10.3f" % (
        n, len(ix), R[ix, 0].sum(), R[ix, 1].sum(), R[ix, 2].sum(),
        np.abs(R[ix, 0]).sum(), np.abs(R[ix, 1]).sum(), np.abs(R[ix, 2]).sum()))
allmask = np.zeros(len(nid), bool)
for ids in NS.values():
    allmask[ids - 1] = True
print("  %-22s %4d %10.4f %10.4f %10.4f" % ("UNION(BC sets)", allmask.sum(),
      R[allmask, 0].sum(), R[allmask, 1].sum(), R[allmask, 2].sum()))
print("  %-22s %4d %10.4f %10.4f %10.4f" % ("ALL nodes", len(nid), R[:, 0].sum(), R[:, 1].sum(), R[:, 2].sum()))
print("  %-22s %4d %10.4f %10.4f %10.4f" % ("FREE nodes (not BC)", (~allmask).sum(),
      R[~allmask, 0].sum(), R[~allmask, 1].sum(), R[~allmask, 2].sum()))
print("  free-node max|R| = %.6f  (residual scale)" % np.abs(R[~allmask]).max())

print("\n=== u on BC sets (max |component|) ===")
for n, ids in NS.items():
    ix = ids - 1
    print("  %-22s |ux|max=%.3e |uy|max=%.3e |uz|max=%.3e  uz[%.5f,%.5f]" % (
        n, np.abs(U[ix, 0]).max(), np.abs(U[ix, 1]).max(), np.abs(U[ix, 2]).max(),
        U[ix, 2].min(), U[ix, 2].max()))

print("\n=== h5 internal nodesets vs feb (sanity) ===")
for k in sorted(ns_h5):
    v = set(int(x) for x in ns_h5[k])
    if k in NS:
        fv = set(int(x) for x in NS[k])
        print("  %-22s h5=%4d feb=%4d  |diff|=%d" % (k, len(v), len(fv), len(v ^ fv)))
    else:
        print("  %-22s h5=%4d (feb: absent)" % (k, len(v)))
