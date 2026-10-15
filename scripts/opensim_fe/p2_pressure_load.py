# -*- coding: utf-8 -*-
"""修法 ①：把「底部轴向带 59 个节点力」换成**端面分布压强**（加法式注入）。

依据（本仓库 P2_V2 实测）：σ 峰值随载荷带边缘迁移、σ_max 随带宽 −66%
⇒ 峰值是**载荷施加方式**的伪影（集中力奇异 + 带边缘剪切集中），非物理。

修法：在骨 −Z 端**端面**上施加 `<surface_load type="pressure">`，
总力不变（P = load_n / A_face），去掉 59 个集中力与其带边缘。

⚠️ 单一改动：只换载荷施加方式；注册修正保持；其余逐字不动。
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import bone_feb as bf                  # noqa: E402
import g4_fe_p1b as p1b                # noqa: E402
import g4_fe_plate as plate            # noqa: E402
import p2_head_contact as p2           # noqa: E402
import p4_pb_probe as p4               # noqa: E402
import p5_head_contact as p5           # noqa: E402

_ORIG_PG = p5.plate_geom_from_mesh


def regfix(nodes_rot, geom):
    pg = _ORIG_PG(nodes_rot, geom)
    z_top = float(max(nodes_rot[np.asarray(f["poly"], int) - 1, 2].max()
                      for f in geom["cap"]))
    zw = z_top + p5.GAP_MM
    pg["z_wall"] = zw
    pg["z_bot"] = zw - plate.PLATE_THICK_MM
    pg["z_mid"] = zw - 0.5 * plate.PLATE_THICK_MM
    return pg


p5.plate_geom_from_mesh = regfix


def load_face_geom(mesh, spec, with_shell=True, axdiv_max=-0.5):
    """−Z 端**端面**：外表面中法向≈−z 的面片。返回 facets / A / nodes。"""
    nodes_rot = bf._align_pca_rotate(mesh["nodes"])
    axis_idx, a_min, a_max, top_local, bot_local = bf._select_top_bottom(
        nodes_rot, axis=None, top_frac=p2.TOP_FRAC, bot_frac=p2.LOAD_FRAC)
    elems = (p2._combined_elems_all(mesh, spec) if with_shell
             else p5._combined_elems(mesh, spec))
    facets = p1b.oriented_facets(nodes_rot, elems)
    band = p1b.band_facets(facets, bot_local)
    e = np.zeros(3)
    e[axis_idx] = 1.0
    for f in band:
        f["axdiv"] = float(np.dot(f["normal"], e))
    face = [f for f in band if f["axdiv"] < axdiv_max]     # 法向 −z ⇒ 端面
    A = float(sum(f["area"] for f in face))
    nodes = set()
    for f in face:
        nodes.update(int(x) for x in f["poly"])
    return {"facets": face, "area_mm2": A, "n_nodes": len(nodes),
            "axis_idx": axis_idx, "a_min": float(a_min), "a_max": float(a_max),
            "n_band": int(len(band)), "n_face": len(face),
            "bot_band_z": [float(nodes_rot[bot_local - 1, axis_idx].min()),
                           float(nodes_rot[bot_local - 1, axis_idx].max())]}


def surface_xml(face: dict, name="load_face") -> str:
    out = [f'    <Surface name="{name}">']
    for i, f in enumerate(face["facets"]):
        poly = [int(p) for p in f["poly"]]
        tag = "tri3" if len(poly) == 3 else ("quad4" if len(poly) == 4 else None)
        if tag is None:
            raise ValueError("端面含 %d 顶点面片" % len(poly))
        out.append(f'      <{tag} id="{i + 1}">{",".join(map(str, poly))}</{tag}>')
    out.append("    </Surface>")
    out.append('    <!-- load_face: end-face pressure (fix 1), %d facets %.2f mm2 -->'
               % (face["n_face"], face["area_mm2"]))
    return "\n".join(out)


def inject_pressure(text: str, face: dict, load_n: float, lc: str = "2") -> tuple[str, float]:
    """替换 <Loads> 为端面压强。返回 (text, P_MPa)。"""
    P = float(load_n) / face["area_mm2"]        # N/mm² = MPa
    # 1) Surface 插到 <SurfacePair 之前（同一 <Mesh> 内）
    if '<SurfacePair' not in text:
        raise ValueError("未找到 <SurfacePair")
    text = text.replace('<SurfacePair', surface_xml(face) + "\n    <SurfacePair", 1)
    # 2) <Loads> 整段替换
    new = (f'      <Loads>\n'
           f'        <surface_load type="pressure" surface="load_face">\n'
           f'          <pressure lc="{lc}">{P:.10g}</pressure>\n'
           f'        </surface_load>\n'
           f'      </Loads>\n')
    text, n = re.subn(r"\s*<Loads>.*?</Loads>\n", "\n" + new, text,
                      count=1, flags=re.S)
    if n != 1:
        raise ValueError("未找到 <Loads>")
    if '<nodal_load' in text:
        raise ValueError("nodal_load 残留")
    return text, P


def main() -> int:
    ap = __import__("argparse").ArgumentParser()
    ap.add_argument("--heights", nargs="*", type=int, default=[10])
    ap.add_argument("--sign", type=float, default=1.0, help="压强符号（+1/-1）")
    args = ap.parse_args()

    p2.BASE = ROOT / "temp" / "opensim_fe" / "p2_pressure"

    ctx = p4.load_ctx("T6")
    part, mesh, solid_names, shell_specs, dom = ctx[0], ctx[1], ctx[2], ctx[3], ctx[4]
    geom = p2.cap_geom(mesh, part.spec, top_frac=p2.TOP_FRAC, include_shell=True)
    pg = regfix(bf._align_pca_rotate(mesh["nodes"]), geom)
    face = load_face_geom(mesh, part.spec, with_shell=True)
    print("=" * 88)
    print("修法① 端面压强（替换底部轴向带 59 节点力）")
    print("=" * 88)
    print("  端面: %d 面片 / %.2f mm² / %d 节点（带宽内共 %d 面片）"
          % (face["n_face"], face["area_mm2"], face["n_nodes"], face["n_band"]))
    print("  底部带 z ∈ [%.2f, %.2f]（轴 %d），骨 z ∈ [%.2f, %.2f]"
          % (face["bot_band_z"][0], face["bot_band_z"][1], face["axis_idx"],
             face["a_min"], face["a_max"]))

    en = p2.end_nodes(mesh)
    share = p2.LOAD_SHARE["T6"]
    joint_loads = ctx[6]
    rows = []
    for hno in args.heights:
        load = float(joint_loads[hno - 1]) * share
        d = p2.BASE / "T6" / ("h%02d" % hno)
        d.mkdir(parents=True, exist_ok=True)
        feb = d / "T6.feb"
        h5 = feb.with_suffix(".hdf5")
        # 先按常规路径生成（含包板注册修正），再注入压强
        p2.build_one("T6", feb, load, ctx, en, geom, pg, p2.penalty_for(load),
                     with_shell=True)
        text = feb.read_text(encoding="ISO-8859-1")
        text, P = inject_pressure(text, face, load * args.sign, lc="2")
        feb.write_text(text, encoding="ISO-8859-1")
        print("\n  h=%d  载荷 %.1f N  ⇒  P = %.6f MPa  (%d 面片 / %.2f mm²)"
              % (hno, load, P, face["n_face"], face["area_mm2"]), flush=True)
        t = time.time()
        rec = {"h_no": hno, "load_n": load, "P_mpa": P, "sign": args.sign}
        try:
            bf.run_bone_feb(feb, workdir=feb.parent, timeout=1800)
            if feb.with_suffix(".xplt").is_file():
                bf.post_process(feb.with_suffix(".xplt"), h5,
                                solid_names=solid_names,
                                shell_specs=shell_specs, label="T6")
            rec["converged"] = True
        except Exception as e:                       # noqa: BLE001
            rec["converged"] = False
            rec["error"] = "%s: %s" % (type(e).__name__, str(e)[:150])
        rec["seconds"] = round(time.time() - t, 1)
        if h5.is_file():
            with h5py.File(h5, "r") as h:
                sid = sorted(h["states"].keys(), key=int)[-1]
                co = h["meshes/0/nodes"][:]
                xyz = np.stack([co["x"], co["y"], co["z"]], 1).astype(float)
                dk = [k for k in h["meshes/0/domains"].keys()
                      if int(h[f"meshes/0/domains/{k}"].shape[0]) == 712][0]
                conn = h[f"meshes/0/domains/{dk}"][:]
                st = h[f"states/{sid}/element_data/stress/{dk}"][:].astype(float)
                cen = xyz[conn[:, 1:4] - 1].mean(axis=1)
                v = np.sqrt(0.5 * ((st[:, 0] - st[:, 1]) ** 2
                                   + (st[:, 1] - st[:, 2]) ** 2
                                   + (st[:, 2] - st[:, 0]) ** 2)
                            + 3 * (st[:, 3] ** 2 + st[:, 4] ** 2 + st[:, 5] ** 2))
                z0, z1 = xyz[:, 2].min(), xyz[:, 2].max()
                sp = z1 - z0
                pk = np.argsort(-v)[:6]
                rec.update({"sigma_max": float(v.max()),
                            "sigma_median": float(np.median(v)),
                            "peak_zfrac": [float((cen[j, 2] - z0) / sp) for j in pk],
                            "u_mean_z": None})
                u = np.asarray(h[f"states/{sid}/node_data/displacement/T6_nodes"],
                               dtype=float)
                rec["u_mean_z"] = float(u[:, 2].mean())
                rec["u_max"] = float(np.linalg.norm(u, axis=1).max())
            try:
                bs = p2.load_bear_stats(h5, load, "T6")
                rec.update({"bear_frac": bs["bear_frac"],
                            "R_axis_check": bs["R_axis_check"],
                            "contact_force_z": bs["contact_f_bone_n"][2]})
            except Exception as e:                   # noqa: BLE001
                rec["bear_error"] = str(e)[:120]
        rows.append(rec)
        if rec.get("sigma_max") is not None:
            msg = ("σ_max=%.1f σ_med=%.1f 峰值z_frac=%s u_mean_z=%.3f bear=%s Fz=%s"
                   % (rec["sigma_max"], rec["sigma_median"],
                      [round(x, 3) for x in rec["peak_zfrac"]],
                      rec["u_mean_z"], rec.get("bear_frac"),
                      rec.get("contact_force_z")))
        else:
            msg = str(rec.get("error", "?"))
        print("     conv=%s  %s" % (rec["converged"], msg), flush=True)

    out = {"fix": "surface pressure on end face (replace nodal band force)",
           "face": {k: v for k, v in face.items() if k != "facets"},
           "rows": rows}
    dst = ROOT / "results" / "opensim_fe" / ("P2_PRESSURE_LOAD_sign%+d.json" % args.sign)
    dst.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=str),
                   encoding="utf-8")
    print("\n产物:", dst.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
