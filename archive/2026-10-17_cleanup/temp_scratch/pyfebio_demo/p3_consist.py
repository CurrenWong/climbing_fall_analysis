import h5py, numpy as np, re

H5 = r"temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = r"temp/thums/febio_lig/ankle_lig_WK.feb"
txt = open(FEB, errors="replace").read()

# --- feb nodes
febnodes = {}
for m in re.finditer(r'<node id="(\d+)">([^<]*)</node>', txt):
    i = int(m.group(1))
    febnodes[i] = tuple(float(v) for v in m.group(2).split(","))
print("feb <node> count:", len(febnodes))

# --- feb elements
febels = {}
for m in re.finditer(r'<Elements type="(\w+)" name="([^"]+)">(.*?)</Elements>', txt, re.S):
    typ, name, body = m.group(1), m.group(2), m.group(3)
    d = {}
    for e in re.finditer(r'<elem id="(\d+)">([^<]*)</elem>', body):
        d[int(e.group(1))] = tuple(int(x) for x in re.split(r"[,\s]+", e.group(2).strip()) if x)
    febels[name] = (typ, d)
print("feb element domains:", len(febels), " total elems:", sum(len(d) for _, d in febels.values()))

with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"]
    nid = nd["id"][:]; P = np.column_stack([nd["x"][:], nd["y"][:], nd["z"][:]]).astype(np.float64)
    doms = {k: f["meshes/0/domains"][k][:] for k in f["meshes/0/domains"].keys()}
    nsets = {k: f["meshes/0/nodesets"][k][:] for k in f["meshes/0/nodesets"].keys()}
    U20 = f["states/20/node_data/displacement/all"][:].astype(np.float64)

# --- coordinate consistency on shared ids
shared = np.array([i for i in nid.tolist() if i in febnodes])
idx = np.array([j for j, i in enumerate(nid.tolist()) if i in febnodes])
fp = np.array([febnodes[i] for i in shared])
d = np.abs(P[idx] - fp)
print("\ncoords: shared ids=%d / %d   max|dP| = %.3e mm  (per axis %s)" % (
    len(shared), len(nid), d.max(), np.round(d.max(0), 8)))
print("  frac |dP|>1e-3:", float((d.max(1) > 1e-3).mean()))

# --- element connectivity consistency
print("\n=== element connectivity: h5 vs feb ===")
h5names = {}
for k in doms:
    if k.startswith("__unnamed_"):
        continue
    h5names[k] = doms[k]
for k in sorted(h5names):
    if k not in febels:
        print("  %-28s NOT IN FEB" % k); continue
    typ, fe = febels[k]
    h5e = h5names[k]
    h5map = {int(r[0]): tuple(int(x) for x in r[1:]) for r in h5e}
    same = ok = 0
    bad = []
    for eid, conn in fe.items():
        if eid in h5map:
            same += 1
            if tuple(h5map[eid]) == tuple(conn):
                ok += 1
            elif len(bad) < 3:
                bad.append((eid, conn, h5map[eid]))
    print("  %-28s feb=%-5s n_feb=%4d n_h5=%4d shared=%4d identical=%4d %s" % (
        k, typ, len(fe), len(h5e), same, ok, ("MISMATCH e.g. " + str(bad[:1])) if ok < same else "OK"))

print("\n=== h5 nodeset vs feb nodeset ===")
def feb_ns(name):
    m = re.search(r'<NodeSet name="%s">(.*?)</NodeSet>' % re.escape(name), txt, re.S)
    return set(int(x) for x in re.split(r"[,\s]+", m.group(1).strip()) if x) if m else None
for nm in sorted(nsets):
    fv = feb_ns(nm)
    h5v = set(int(x) for x in nsets[nm])
    if fv is None:
        print("  %-22s h5=%-4d feb=ABSENT" % (nm, len(h5v))); continue
    print("  %-22s h5=%-4d feb=%-4d  h5\\feb=%d  feb\\h5=%d" % (
        nm, len(h5v), len(fv), len(h5v - fv), len(fv - h5v)))
    if nm == "tibia_top_ns":
        print("     h5\\feb ids:", sorted(h5v - fv))
        print("     feb\\h5 ids:", sorted(fv - h5v))

print("\n=== exact uz of tibia_top_ns ===")
id2idx = {int(i): j for j, i in enumerate(nid)}
ix = np.array([id2idx[int(i)] for i in nsets["tibia_top_ns"]])
uz = U20[ix, 2]
u = np.unique(np.round(uz, 6))
print("  unique uz (rounded 1e-6):", u)
print("  count == -1.6 exactly:", int((uz == -1.6).sum()), "  count==0:", int((uz == 0).sum()))
