"""读 xplt 每个 state 的 (time, 沿轴支反力) —— 判断模型是否在同一加载状态比较。

★ 为什么必须看 time：FEBio 的 time_stepper 会自适应放大步长，
  "完成 14 步" ≠ "加载到 70%"。之前拿 436.8 / 709 / 969 N 互比而不看 time，
  得出的是「不同加载状态下的支反力」，结论无效。
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X

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

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
AX = np.array([0.6766, 0.1363, 0.7237])
AX = AX / np.linalg.norm(AX)


def look(fp: Path, label: str = ""):
    h5 = W / f"_t_{fp.stem}.h5"
    X.to_hdf5(str(fp), str(h5))
    with h5py.File(h5, "r") as f:
        sts = sorted(f["states"].keys(), key=int)
        print(f"\n=== {label or fp.name}  ({len(sts)} states) ===")
        print(f"{'state':>6s} {'time':>10s} {'R_axis(N)':>12s} {'uz_min':>10s}")
        for s in sts:
            g = f[f"states/{s}"]
            t = float(np.asarray(g.attrs.get("time", np.nan)).ravel()[0])
            if "node_data/reaction forces" not in g:
                continue
            R = np.asarray(g["node_data/reaction forces/all"])
            D = np.asarray(g["node_data/displacement/all"])
            print(f"{s:>6s} {t:>10.5f} {float((R @ AX).sum()):>12.3f} "
                  f"{float((D @ AX).min()):>10.4f}")


if __name__ == "__main__":
    for a in sys.argv[1:]:
        p = Path(a)
        look(p if p.is_absolute() else (W / p))
