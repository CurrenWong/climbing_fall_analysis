import h5py, numpy as np, os, sys

H5 = r"temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = r"temp/thums/febio_lig/ankle_lig_WK.feb"

with h5py.File(H5, "r") as f:
    def show(g, path="", depth=0, maxd=3):
        for k in g.keys():
            o = g[k]
            p = path + "/" + k
            if isinstance(o, h5py.Group):
                print("  " * depth + "G " + p + "  (" + str(len(o.keys())) + " keys)")
                if depth < maxd:
                    show(o, p, depth + 1, maxd)
            else:
                try:
                    sh = o.shape
                except Exception:
                    sh = "?"
                print("  " * depth + "D " + p + "  shape=" + str(sh) + " dtype=" + str(o.dtype))
    print("=== H5 TREE (depth3) ===")
    show(f)
    print()
    print("=== nodesets ===")
    for k in f["meshes/0/nodesets"].keys():
        print("  ", k, f["meshes/0/nodesets"][k].shape)
    print("=== domains ===")
    for k in f["meshes/0/domains"].keys():
        print("  ", k, f["meshes/0/domains"][k].shape)
    print("=== states ===")
    st = f["states"]
    print("  n states:", len(st.keys()))
    s0 = list(st.keys())[0]
    print("  state keys sample:", s0)
    for k in st[s0].keys():
        o = st[s0][k]
        if isinstance(o, h5py.Group):
            print("   G", k, list(o.keys())[:10])
        else:
            print("   D", k, o.shape)
    nd = st[s0]["node_data"]
    print("  node_data fields:", list(nd.keys()))
    for k in nd.keys():
        print("    ", k, list(nd[k].keys()) if isinstance(nd[k], h5py.Group) else nd[k].shape)
    ed = st[s0].get("element_data")
    if ed is not None:
        print("  element_data domains:", len(ed.keys()))
        for k in list(ed.keys()):
            print("    ", k, list(ed[k].keys()), ed[k][list(ed[k].keys())[0]].shape)
