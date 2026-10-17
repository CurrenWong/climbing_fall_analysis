import pickle, numpy as np, h5py
D = pickle.load(open(r"temp/pyfebio_demo/p3_cache.pkl", "rb"))
DATA, SIG, EL = D["DATA"], D["SIG"], D["EL"]
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    P = np.column_stack([f["meshes/0/nodes"]["x"][:], f["meshes/0/nodes"]["y"][:],
                         f["meshes/0/nodes"]["z"][:]]).astype(float)
    U = f["states/20/node_data/displacement/all"][:].astype(float)

SOLID = [k for k in EL if not k.startswith("lig_")]
LIG = [k for k in EL if k.startswith("lig_")]

# element z extents
ext = {}
for k in EL:
    Z = P[EL[k]["conn"]]
    ext[k] = (Z[:, :, 2].min(1), Z[:, :, 2].max(1))

DZ = 2.0
z0, z1 = P[:, 2].min(), P[:, 2].max()
levels = np.arange(z1, z0, -DZ)
print("=" * 118)
print("AXIAL FORCE TRANSMITTED ACROSS HORIZONTAL CUT PLANES  (F = sum sigma_zz * V_overlap / DZ)")
print("  negative sigma_zz = compression ; reported F_push = -sum(sigma_zz V)/DZ  [N]")
print("=" * 118)
print("  z-slab                 F_solid[N]   region(s) crossed")
for zz in levels:
    a, b = zz, zz - DZ
    tot = 0.0
    who = {}
    for k in SOLID:
        lo, hi = ext[k]
        span = np.maximum(hi - lo, 1e-9)
        ov = np.clip(np.minimum(hi, a) - np.maximum(lo, b), 0, None)
        frac = ov / span
        m = frac > 0
        if not m.any():
            continue
        v = -((-SIG[k][m, 2]) * DATA[k][1][m] * frac[m]).sum() / DZ
        tot += v
        who[k.split("__")[0]] = who.get(k.split("__")[0], 0.0) + v
    lig = 0.0
    for k in LIG:
        lo, hi = ext[k]
        m = (lo <= a) & (hi >= b)
        if m.any():
            lig += -((-SIG[k][m, 2]) * DATA[k][1][m]).sum() / DZ
    top = sorted(who.items(), key=lambda kv: -abs(kv[1]))[:3]
    print("   %7.2f..%7.2f  %9.2f   lig=%+7.2f | %s" % (
        a, b, tot, lig, "  ".join("%s:%.0f" % kv for kv in top)))

print()
print("=" * 118)
print("WHERE ARE THE REACTION NODES?  (|Rz|>1 N)")
print("=" * 118)
NS = D["NS"]
sets = {n: set((NS[n] - 1).tolist()) for n in NS}
def which(i):
    return ",".join(sorted(n for n, s in sets.items() if i in s)) or "-"
grp = {}
for k in EL:
    g = ("tibia" if k.startswith("tibia_") else "fibula" if k.startswith("fibula_") else
         "cartilage" if k.startswith("tibiotalar") else "foot" if k.split("__")[0] in
         ("talus_cort", "talus_spon", "calcaneus_cort", "calcaneus_spon", "navicular_spon") else "ligament")
    for i in EL[k]["conn"].ravel():
        grp.setdefault(int(i), set()).add(g)
big = np.abs(D["Udom"] and np.zeros(1))  # dummy
with h5py.File(r"temp/thums/febio_lig/ankle_lig_WK_stress.h5", "r") as f:
    R = f["states/20/node_data/reaction forces/all"][:].astype(float)
m = np.abs(R[:, 2]) > 1.0
ii = np.where(m)[0]
from collections import Counter
c = Counter()
for i in ii:
    c[("+".join(sorted(grp.get(int(i), {"?"}))), which(int(i)))] += 1
for k, v in c.most_common(12):
    s = [i for i in ii if ("+".join(sorted(grp.get(int(i), {"?"}))), which(int(i))) == k]
    print("  %-46s n=%4d  sumRz=%+9.2f  z=[%7.2f,%7.2f]  uz=[%+.3f,%+.3f]" % (
        "%s | %s" % k, v, R[s, 2].sum(), P[s, 2].min(), P[s, 2].max(), U[s, 2].min(), U[s, 2].max()))
