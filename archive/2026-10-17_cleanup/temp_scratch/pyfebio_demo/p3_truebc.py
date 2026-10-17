import h5py, numpy as np, re

H5 = r"temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = r"temp/thums/febio_lig/ankle_lig_WK.feb"
txt = open(FEB, errors="replace").read()

def feb_ns(name):
    m = re.search(r'<NodeSet name="%s">(.*?)</NodeSet>' % re.escape(name), txt, re.S)
    return np.array(sorted(int(x) for x in re.split(r"[,\s]+", m.group(1).strip()) if x)) if m else None

with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"]
    nid = nd["id"][:].astype(int)
    P = np.column_stack([nd["x"][:], nd["y"][:], nd["z"][:]]).astype(np.float64)
    nsets = {k: f["meshes/0/nodesets"][k][:] for k in f["meshes/0/nodesets"].keys()}
    U = f["states/20/node_data/displacement/all"][:].astype(np.float64)
    R = f["states/20/node_data/reaction forces/all"][:].astype(np.float64)

id2idx = {int(i): j for j, i in enumerate(nid)}
def sel(ids): return np.array([id2idx[int(i)] for i in ids], dtype=int)

def report(tag, ids):
    ix = sel(ids)
    u = U[ix]
    print("  %-34s n=%3d  |ux|max=%.2e |uy|max=%.2e |uz|max=%.2e  uz[min,max]=[%.5f,%.5f]  sumRz=%8.3f" % (
        tag, len(ix), np.abs(u[:, 0]).max(), np.abs(u[:, 1]).max(), np.abs(u[:, 2]).max(),
        u[:, 2].min(), u[:, 2].max(), R[ix, 2].sum()))
    return ix

print("=== candidate BC sets: prescribed-z (expect uz = -1.6 for all) ===")
report("h5 nodesets/tibia_top_ns", nsets["tibia_top_ns"])
report("feb NodeSet  tibia_top_ns", feb_ns("tibia_top_ns"))
report("h5 nodesets/fibula_top_ns", nsets["fibula_top_ns"])
report("feb NodeSet  fibula_top_ns", feb_ns("fibula_top_ns"))

print("\n=== candidate fixed sets (expect u = 0 for all) ===")
report("h5 nodesets/calcaneus_bottom_ns", nsets["calcaneus_bottom_ns"])
report("feb NodeSet  calcaneus_bottom_ns", feb_ns("calcaneus_bottom_ns"))
report("h5 nodesets/talus_low_ns", nsets["talus_low_ns"])
report("feb NodeSet  talus_low_ns", feb_ns("talus_low_ns"))

print("\n=== FIND the true prescribed set from the data ===")
for lo, hi, tag in [(-1.61, -1.595, "uz in [-1.61,-1.595]")]:
    m = (U[:, 2] >= lo) & (U[:, 2] <= hi)
    print("  %s : %d nodes ; Rz sum = %.3f N" % (tag, m.sum(), R[m, 2].sum()))
    ids = set(nid[m].tolist())
    for nm in ["tibia_top_ns", "fibula_top_ns"]:
        h5v = set(int(x) for x in nsets[nm]); fv = set(feb_ns(nm).tolist())
        print("     vs h5 %-16s : extra_in_data=%3d missing_from_data=%3d" % (
            nm, len(ids - h5v), len(h5v - ids)))
        print("     vs feb %-15s : extra_in_data=%3d missing_from_data=%3d" % (
            nm, len(ids - fv), len(fv - ids)))

print("\n=== FIND the true zero-displacement sets ===")
zm = (np.abs(U).max(1) < 1e-9)
print("  nodes with u == 0 exactly: %d  ; Rz sum = %.3f N ; |R| sum=%.3f" % (
    zm.sum(), R[zm, 2].sum(), np.abs(R[zm]).sum()))
ids = set(nid[zm].tolist())
for nm in ["calcaneus_bottom_ns", "talus_low_ns", "tibia_top_ns", "fibula_top_ns", "navicular_all_ns"]:
    h5v = set(int(x) for x in nsets[nm]); fv = set(feb_ns(nm).tolist())
    print("  %-20s data\\h5=%3d h5\\data=%3d | data\\feb=%3d feb\\data=%3d" % (
        nm, len(ids - h5v), len(h5v - ids), len(ids - fv), len(fv - ids)))

print("\n  nodes with ux==0 and uy==0 (any uz): %d" % (((np.abs(U[:, 0]) < 1e-12) & (np.abs(U[:, 1]) < 1e-12)).sum()))
print("  nodes with ux==0 & uy==0 & uz==0 : %d" % zm.sum())
print("  nodes with ux==0 & uy==0 & uz!=0 : %d" % (((np.abs(U[:, 0]) < 1e-12) & (np.abs(U[:, 1]) < 1e-12) & ~zm).sum()))
print("  nodes with Rz != 0 : %d, sum=%.3f" % ((np.abs(R[:, 2]) > 1e-9).sum(), R[:, 2].sum()))
