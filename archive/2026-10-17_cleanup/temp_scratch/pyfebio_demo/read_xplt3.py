"""A5: 读 .xplt 取应力。pyfebio 的 parse_state 对自定义 Surface/SolidDomain 集合会
KeyError（mesh_dict 查名字失败），这里打一个容错补丁：查不到就合成名字，不让整份文件读不出来。
"""
from pathlib import Path
import numpy as np
import h5py
import pyfebio.xplt as X

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")

_VAR_LUT = {"node_data": "PLT_DIC_NODAL", "surface_data": "PLT_DIC_SURFACE",
            "element_data": "PLT_DIC_DOMAIN"}
_SET_LUT = {"node_data": "node_sets", "surface_data": "surfaces",
            "element_data": "domains"}


def tolerant_parse_state(buffer, state_cnt, xdictionary, mesh_dict, f):
    print(f"  [patch] mesh_dict keys={list(mesh_dict)} "
          f"domains={mesh_dict.get('domains')}")
    i = 0
    while i < len(buffer) - 8:
        tag, offset = X.parse_prefix(buffer[i:i + 8])
        child = buffer[i + 8:i + 8 + offset]
        name = X.TAG_LUT[tag].name
        if name == "PLT_STATE_HEADER":
            for k, v in X._parse_header(child).items():
                f.create_group(f"states/{state_cnt}") if f"states/{state_cnt}" not in f else None
                f[f"states/{state_cnt}"].attrs[k] = v
        elif name == "PLT_MESH_STATE":
            for k, v in X._parse_mesh_state(child).items():
                f.create_dataset(f"states/{state_cnt}/mesh/{k}", data=v)
        elif name == "PLT_OBJECTS_STATE":
            X._parse_objects_state(child)
        elif name == "PLT_STATE_DATA":
            dd = X._parse_state_data(child)
            for key, items in dd.items():
                tbl = mesh_dict.get(_SET_LUT[key])
                for item in items:
                    for set_id, data in item["data"].items():
                        vi = item["id"][0] - 1
                        var_name = xdictionary[_VAR_LUT[key]][vi].name
                        itype = xdictionary[_VAR_LUT[key]][vi].itype
                        sid = (set_id + 1) if key == "node_data" else set_id
                        set_name = None
                        try:
                            if isinstance(tbl, dict):
                                set_name = tbl.get(sid)
                                if set_name is None:
                                    set_name = tbl.get(int(sid)) or tbl.get(set_id)
                            else:
                                set_name = tbl[sid]
                        except Exception:
                            set_name = None
                        if set_name is None:
                            set_name = f"__set{sid}"
                        dset = f"states/{state_cnt}/{key}/{var_name}/{set_name}"
                        f.create_dataset(dset, data=data.reshape(X.VAR_SHAPE_LUT[itype]))
        i += offset + 8


X.parse_state = tolerant_parse_state

src = W / "ankle_tibiotalar.xplt"
dst = W / "ankle_tibiotalar.h5"
X.to_hdf5(str(src), str(dst))
print(f"[ok] 转换完成 {dst.stat().st_size/1024/1024:.1f} MB")

with h5py.File(dst, "r") as f:
    states = sorted(f["states"].keys(), key=int)
    print(f"states: {len(states)} 个 -> {states[:3]}...{states[-3:]}")
    last = states[-1]
    def walk(g, pre="", d=0):
        for k in g:
            o = g[k]
            if isinstance(o, h5py.Group):
                if d < 3:
                    walk(o, pre + k + "/", d + 1)
            else:
                print(f"    {pre}{k}  shape={o.shape}")
    print(f"\n最后一个 state ({last}) 的数据:")
    walk(f[f"states/{last}"])
