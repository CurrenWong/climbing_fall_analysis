"""读 .xplt：支反力 + 应力 + 与吴恺 600 N / 3.74 MPa 定标对照。

可用 XPLT 环境变量指定输入（默认 temp/thums/febio/ankle_tibiotalar.xplt）。
"""
import os
from pathlib import Path
import numpy as np
import h5py
import pyfebio.xplt as X

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")

_VAR = {"node_data": "PLT_DIC_NODAL", "surface_data": "PLT_DIC_SURFACE",
        "element_data": "PLT_DIC_DOMAIN"}
_SET = {"node_data": "node_sets", "surface_data": "surfaces",
        "element_data": "domains"}


def tolerant_parse_state(buffer, state_cnt, xdictionary, mesh_dict, f):
    i = 0
    while i < len(buffer) - 8:
        tag, offset = X.parse_prefix(buffer[i:i + 8])
        child = buffer[i + 8:i + 8 + offset]
        nm = X.TAG_LUT[tag].name
        if nm == "PLT_STATE_HEADER":
            grp = f"states/{state_cnt}"
            if grp not in f:
                f.create_group(grp)
            for k, v in X._parse_header(child).items():
                f[grp].attrs[k] = v
        elif nm == "PLT_MESH_STATE":
            for k, v in X._parse_mesh_state(child).items():
                f.create_dataset(f"states/{state_cnt}/mesh/{k}", data=v)
        elif nm == "PLT_OBJECTS_STATE":
            X._parse_objects_state(child)
        elif nm == "PLT_STATE_DATA":
            for key, items in X._parse_state_data(child).items():
                tbl = mesh_dict.get(_SET[key])
                for item in items:
                    for set_id, data in item["data"].items():
                        vi = item["id"][0] - 1
                        vn = xdictionary[_VAR[key]][vi].name
                        it = xdictionary[_VAR[key]][vi].itype
                        sid = (set_id + 1) if key == "node_data" else set_id
                        sn = None
                        try:
                            sn = tbl.get(sid) if isinstance(tbl, dict) else tbl[sid]
                            if sn is None and isinstance(tbl, dict):
                                sn = tbl.get(int(sid)) or tbl.get(set_id)
                        except Exception:
                            sn = None
                        f.create_dataset(
                            f"states/{state_cnt}/{key}/{vn}/{sn or f'__set{sid}'}",
                            data=data.reshape(X.VAR_SHAPE_LUT[it]))
        i += offset + 8


X.parse_state = tolerant_parse_state
src = Path(os.environ.get("XPLT", "").strip() or (W / "ankle_tibiotalar.xplt"))
dst = W / "_read.h5"
X.to_hdf5(str(src), str(dst))
print(f"[ok] {src.name} -> {dst.stat().st_size/1024/1024:.1f} MB")

BY_COUNT = {2415: "tibia_cort(hex8)", 6: "tibia_cort(tet4)", 6651: "tibia_end_spon(松质)",
            14482: "tibia_marrow(骨髓)", 394: "talus_cort(hex8)",
            1766: "talus_spon(松质)", 458: "tibiotalar_tibia(软骨)",
            613: "tibiotalar_talus(软骨)"}


def vm(s):
    a, b, c, d, e, g = (s[:, i] for i in range(6))
    return np.sqrt(0.5 * ((a - b) ** 2 + (b - c) ** 2 + (c - a) ** 2)
                   + 3.0 * (d ** 2 + e ** 2 + g ** 2))


AXIS = np.array([0.6766, 0.1363, 0.7237])       # 距骨->胫骨 加载轴
AXIS = AXIS / np.linalg.norm(AXIS)

with h5py.File(dst, "r") as f:
    sts = sorted(f["states"].keys(), key=int)
    st = sts[-1]
    print(f"state {st}（共 {len(sts)} 个）")
    nd = f[f"states/{st}/node_data"]
    print(f"  node_data 变量: {list(nd)}")

    def arr_of(node):
        return np.asarray(node["all"] if isinstance(node, h5py.Group) else node)

    for k in nd:
        a = arr_of(nd[k])
        print(f"    {k}: {a.shape}")
        if "reaction" in k.lower():
            nrm = np.linalg.norm(a, axis=1)
            print(f"      -> |R| 中位 {nrm.mean():.4f}  max {nrm.max():.4f}")
            print(f"      -> 沿加载轴分量之和 = {float((a @ AXIS).sum()):.3f}  "
                  f"（N；负=与轴反向，即受压）")
    print(f"\n  element_data 变量: {list(f[f'states/{st}/element_data'])}")
    grp = f[f"states/{st}/element_data/stress"]
    print(f"\n{'域':26s} {'单元':>7s} {'中位':>9s} {'p99':>9s} {'max':>9s}  (von Mises, MPa)")
    for k in grp:
        arr = np.asarray(grp[k])
        v = vm(arr)
        print(f"{BY_COUNT.get(arr.shape[0], k):26s} {arr.shape[0]:>7,} "
              f"{np.median(v):>9.4f} {np.percentile(v,99):>9.4f} {v.max():>9.4f}")
    d = arr_of(f[f"states/{st}/node_data/displacement"])
    ax = d @ AXIS
    print(f"\n位移: |u| max {np.linalg.norm(d,axis=1).max():.4f} mm；"
          f"沿加载轴 {ax.min():.4f} .. {ax.max():.4f} mm")
