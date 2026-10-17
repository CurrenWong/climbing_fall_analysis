import h5py, numpy as np, re

H5 = r"temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = r"temp/thums/febio_lig/ankle_lig_WK.feb"
txt = open(FEB, errors="replace").read()
def feb_ns(name):
    m = re.search(r'<NodeSet name="%s">(.*?)</NodeSet>' % re.escape(name), txt, re.S)
    return np.array(sorted(int(x) for x in re.split(r"[,\s]+", m.group(1).strip()) if x))

sets = {n: set(feb_ns(n).tolist()) for n in ["tibia_top_ns", "fibula_top_ns", "talus_low_ns",
                                             "calcaneus_bottom_ns", "navicular_all_ns"]}
with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"]
    nid = nd["id"][:].astype(int)
    P = np.column_stack([nd["x"][:], nd["y"][:], nd["z"][:]]).astype(np.float64)
    U = f["states/20/node_data/displacement/all"][:].astype(np.float64)
    R = f["states/20/node_data/reaction forces/all"][:].astype(np.float64)
    dtypes = (f["states/20/node_data/reaction forces/all"].dtype,
              f["states/20/node_data/displacement/all"].dtype)

print("h5 dtypes (R,U):", dtypes)
print("Σ Rz over ALL nodes = %.6f ; Σ|Rz| = %.3f" % (R[:, 2].sum(), np.abs(R[:, 2]).sum()))

print("\n=== Σ Rz per set (feb sets) ===")
tot = 0.0
for n, s in sets.items():
    m = np.array([i - 1 for i in sorted(s) if 1 <= i <= len(nid)])
    if len(m) != len(s):
        print("  (note: %s has %d dangling id(s) > n_nodes)" % (n, len(s) - len(m)))
    tot += R[m, 2].sum()
    print("  %-22s n=%3d ΣRz=%9.3f  Σ|Rz|=%9.3f  max|Rz|=%8.3f" % (
        n, len(m), R[m, 2].sum(), np.abs(R[m, 2]).sum(), np.abs(R[m, 2]).max()))
print("  TOTAL over sets: %.3f" % tot)

big = np.abs(R[:, 2]) > 1.0
print("\n=== nodes with |Rz|>1 N : %d nodes, Σ=%9.3f, Σ|.|=%9.3f" % (
    big.sum(), R[big, 2].sum(), np.abs(R[big, 2]).sum()))
print("   z range %.2f .. %.2f | x %.2f..%.2f | y %.2f..%.2f" % (
    P[big, 2].min(), P[big, 2].max(), P[big, 0].min(), P[big, 0].max(), P[big, 1].min(), P[big, 1].max()))
labels = {}
for n, s in sets.items():
    for i in s:
        labels.setdefault(i, []).append(n)
from collections import Counter
c = Counter()
for i in nid[big]:
    c[",".join(sorted(labels.get(int(i), ["<none>"])))] += 1
for k, v in c.most_common():
    print("   %-60s %d" % (k, v))

neg = R[:, 2] < -1.0
print("\n  nodes with Rz < -1 N: %d  Σ=%.3f" % (neg.sum(), R[neg, 2].sum()))
if neg.sum():
    ii = np.where(neg)[0]
    for j in np.argsort(R[ii, 2])[:10]:
        i = ii[j]
        print("    id=%5d P=(%8.2f,%8.2f,%8.2f) Rz=%9.3f uz=%+.4f sets=%s" % (
            nid[i], P[i, 0], P[i, 1], P[i, 2], R[i, 2], U[i, 2], labels.get(int(nid[i]), [])))

pos = R[:, 2] > 1.0
ii = np.where(pos)[0]
print("\n  top 10 Rz>0:")
for j in np.argsort(-R[ii, 2])[:10]:
    i = ii[j]
    print("    id=%5d P=(%8.2f,%8.2f,%8.2f) Rz=%9.3f uz=%+.4f sets=%s" % (
        nid[i], P[i, 0], P[i, 1], P[i, 2], R[i, 2], U[i, 2], labels.get(int(nid[i]), [])))
