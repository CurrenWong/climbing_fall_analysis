"""pyfebio 容错补丁 + 结果读取（复用）。

两处补丁（详见 read_lig.py / MEMORY）：
  ① parse_mesh 对 shell/beam 域没有 "name" 字段 ⇒ KeyError 'name'
  ② parse_state 对自定义 Surface 集合会 KeyError
"""
import textwrap
import inspect
import h5py
import numpy as np
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
AX = np.array([0.8153, 0.135, 0.5631])


def _patch():
    import pyfebio.xplt as X
    _src = textwrap.dedent(inspect.getsource(X.parse_mesh))
    for _o, _n in (('domain["name"]', 'domain.get("name", str(domain["id"][0]))'),
                   ('surface["name"]', 'surface.get("name", str(surface["id"][0]))'),
                   ("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
                   ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
                   ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
                   ("{nodeset['name']}", "{nodeset.get('name', str(nodeset['id'][0]))}"),
                   ("{elementset['name']}",
                    "{elementset.get('name', str(elementset['id'][0]))}")):
        _src = _src.replace(_o, _n)
    _ns = X.__dict__
    exec(compile(_src, "<pm>", "exec"), _ns)
    X.parse_mesh = _ns["parse_mesh"]

    _ps = textwrap.dedent(inspect.getsource(X.parse_state))
    for _a, _b in (
        ('set_name = mesh_dict[set_lut[key]][set_id + 1]',
         'set_name = mesh_dict[set_lut[key]].get(set_id + 1) or '
         'mesh_dict[set_lut[key]].get(int(set_id) + 1) or f"set{set_id}"'),
        ('set_name = mesh_dict[set_lut[key]][set_id]',
         'set_name = mesh_dict[set_lut[key]].get(set_id) or '
         'mesh_dict[set_lut[key]].get(int(set_id)) or f"set{set_id}"'),
    ):
        _ps = _ps.replace(_a, _b)
    _ns2 = X.__dict__
    exec(compile(_ps, "<ps>", "exec"), _ns2)
    X.parse_state = _ns2["parse_state"]
    return X


_X = _patch()


def load_last(h5name="_last.h5"):
    """读最新 .xplt 的末态，返回 dict。"""
    xp = sorted(W.glob("ankle_lig_WK*.xplt"),
                key=lambda p: -p.stat().st_mtime)
    if not xp:
        raise RuntimeError("找不到 .xplt")
    dst = W / h5name
    _X.to_hdf5(str(xp[0]), str(dst))
    out = {"h5": dst}
    with h5py.File(dst, "r") as h:
        st = sorted(int(k) for k in h["states"])
        last = max(st)
        out["nstate"] = len(st)
        out["D"] = np.asarray(
            h[f"states/{last}/node_data/displacement/all"][()],
            float).reshape(-1, 3)
        out["R"] = np.asarray(
            h[f"states/{last}/node_data/reaction forces/all"][()],
            float).reshape(-1, 3)
        out["last"] = last
        out["F"] = float((out["R"] @ AX).sum())
        out["umax"] = float(np.linalg.norm(out["D"], axis=1).max())
        nd = h["meshes/0/nodes"][()]
        out["Xp"] = (np.column_stack([nd["x"], nd["y"], nd["z"]]).astype(float)
                     if nd.dtype.names
                     else np.asarray(nd, float).reshape(-1, 3))
        cart = {}
        for k in h[f"states/{last}/element_data/stress"]:
            if "tibiotalar" in k:
                S = np.asarray(
                    h[f"states/{last}/element_data/stress/{k}"][()], float)
                if S.shape[1] >= 6:
                    sx, sy, sz, sxy, syz, szx = (S[:, 0], S[:, 1], S[:, 2],
                                                 S[:, 3], S[:, 4], S[:, 5])
                    vm = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2
                                        + (sz - sx) ** 2
                                        + 6 * (sxy ** 2 + syz ** 2 + szx ** 2)))
                    cart[k] = (float(np.percentile(vm, 99)), float(vm.max()))
        out["cart"] = cart
    return out
