import h5py, numpy as np, re, collections

H5 = r"temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = r"temp/thums/febio_lig/ankle_lig_WK.feb"

with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"]
    nid = nd["id"][:]
    P = np.column_stack([nd["x"][:], nd["y"][:], nd["z"][:]]).astype(np.float64)
    nsets = {k: f["meshes/0/nodesets"][k][:] for k in f["meshes/0/nodesets"].keys()}
    U20 = f["states/20/node_data/displacement/all"][:].astype(np.float64)
    R20 = f["states/20/node_data/reaction forces/all"][:].astype(np.float64)

id2idx = {int(i): j for j, i in enumerate(nid)}
def sel(ids):
    return np.array([id2idx[int(i)] for i in ids], dtype=int)

Rz = R20[:, 2]
print("sum Rz = %.3f" % Rz.sum())
nz = np.where(np.abs(Rz) > 0.5)[0]
print("nodes with |Rz|>0.5 N :", len(nz), " sum=%.3f" % Rz[nz].sum(), " min=%.3f max=%.3f" % (Rz[nz].min(), Rz[nz].max()))
print("  their z range: %.3f .. %.3f" % (P[nz, 2].min(), P[nz, 2].max()))
rs = np.sort(np.abs(Rz))[::-1]
print("  top 15 |Rz|:", np.round(rs[:15], 3))
print("  count |Rz|>10N:", int((np.abs(Rz) > 10).sum()), " >50N:", int((np.abs(Rz) > 50).sum()))

print()
print("=== memb of tibia_top_ns ===")
ix = sel(nsets["tibia_top_ns"])
for nm, jj in [("dux", 0), ("duy", 1), ("duz", 2)]:
    v = U20[ix, jj]
    print("  %s  exact0 count=%d   |v+1.6|<0.01 count=%d" % (nm, int((np.abs(v) < 1e-9).sum()), int((np.abs(v + 1.6) < 0.01).sum())))

uz = U20[ix, 2]
hist = collections.Counter(np.round(uz, 1))
print("  duz rounded histogram:", sorted(hist.items()))

print()
print("  nodes in tibia_top_ns with |duz|<0.5: (id, x,y,z, ux,uy,uz, Rz)")
m = np.abs(uz) < 0.5
sub = ix[m]
order = np.argsort(P[sub, 2])
for k in order:
    i = sub[k]
    print("   id=%5d  P=(%9.3f,%9.3f,%9.3f)  u=(%+.4f,%+.4f,%+.4f)  Rz=%+9.3f" % (
        nid[i], P[i, 0], P[i, 1], P[i, 2], U20[i, 0], U20[i, 1], U20[i, 2], R20[i, 2]))

print()
print("=== feb: parse NodeSets + tibia_top Surface node list ===")
txt = open(FEB, errors="replace").read()

def parse_nodeset(name):
    m = re.search(r'<NodeSet name="%s">(.*?)</NodeSet>' % re.escape(name), txt, re.S)
    if not m:
        return None
    return np.array([int(x) for x in re.split(r"[,\s]+", m.group(1).strip()) if x])

for nm in ["tibia_top_ns", "fibula_top_ns", "talus_low_ns", "calcaneus_bottom_ns", "navicular_all_ns"]:
    v = parse_nodeset(nm)
    print("  feb %-20s n=%s" % (nm, len(v) if v is not None else None))
    if nm == "tibia_top_ns":
        h5v = np.sort(nsets["tibia_top_ns"])
        print("     identical to h5 set:", bool(np.array_equal(np.sort(v), h5v)))

m = re.search(r'<Surface name="tibia_top">(.*?)</Surface>', txt, re.S)
tris = re.findall(r"<tri\d?>([^<]*)</tri\d?>", m.group(1))
snodes = set()
for t in tris:
    for x in re.split(r"[,\s]+", t.strip()):
        if x:
            snodes.add(int(x))
print("  Surface tibia_top: %d tris, %d unique nodes" % (len(tris), len(snodes)))
sn = np.array(sorted(snodes))
print("  surface node set == tibia_top_ns set :", set(sn.tolist()) == set(nsets["tibia_top_ns"].tolist()))
print("  surface nodes subset of ns :", set(sn.tolist()) <= set(nsets["tibia_top_ns"].tolist()))
print("  n surface-only (in ns not surface):", len(set(nsets["tibia_top_ns"].tolist()) - set(sn.tolist())))
si = sel(sn)
print("  surface duz: mean=%+.5f std=%.4e min=%+.5f max=%+.5f" % (
    U20[si, 2].mean(), U20[si, 2].std(), U20[si, 2].min(), U20[si, 2].max()))
print("  surface sum Rz = %.3f N" % R20[si, 2].sum())
print("  surface z span %.3f .. %.3f" % (P[si, 2].min(), P[si, 2].max()))
print("  surface y span %.3f .. %.3f   x span %.3f .. %.3f" % (P[si,1].min(), P[si,1].max(), P[si,0].min(), P[si,0].max()))
