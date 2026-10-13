# -*- coding: utf-8 -*-
"""FE_BC · W3-a：单骨**逐高度 FE 重解**，替换线性外推。

对照口径（`s4_alignment_bc2.py:173-183`，长骨）：
    现:   sigma_h(i) = sigma_stat * loads[joint][i] / fe_solve_n
    本件: sigma_h(i) = sigma_stat( 在 loads[joint][i] 下**真解** )      # 无外推

判据不变：`m[i] = (ends_h(i) >= sigma_c_ends) | (shaft_h(i) >= sigma_c_shaft)`
（σ_c 取自 `paper_alignment.REGIONS` 的 `_sc`，与 S4 同源）

纪律：冻结输入只读；输出只在 `temp/opensim_fe/fe_bc_perheight/` + 一个 results json。
BC 固定 BC2 条件（`top_frac=0.30, bot_frac=0.30, align_pca=True, time_steps=1`）——
**本件只改「要不要外推」，不改 BC、不改判据。**
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

import bone_feb as bf  # noqa: E402
import fracture_matrix as fm  # noqa: E402
import paper_alignment as pa  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
OUT_JSON = RES / "FE_BC_PERHEIGHT.json"


def statistic(h5: Path, dom: str) -> dict:
    """与 `s4_alignment_bc2.py:159-171` 逐字相同（median 口径）。"""
    import h5py

    with h5py.File(h5, "r") as h:
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
        conn = h["meshes/0/domains"][dom][:]
        sids = sorted(h["states"].keys(), key=int)
        stress = h[f"states/{sids[-1]}/element_data/stress/{dom}"][:].astype(float)
    idx = conn[:, 1:] - 1
    cent = coords[idx].mean(axis=1)
    p = pa.partition_ends_shaft(coords, cent)
    vm = pa.vm_voigt(stress)
    return {"ends_stat_mpa": float(np.median(vm[p["ends_mask"]])),
            "shaft_stat_mpa": float(np.median(vm[p["shaft_mask"]])),
            "vm_max_mpa": float(vm.max())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bone", default="femur_r")
    ap.add_argument("--height-start", type=int, default=0)
    ap.add_argument("--height-end", type=int, default=50)
    args = ap.parse_args()

    bone = args.bone
    out_dir = ROOT / "temp" / "opensim_fe" / "fe_bc_perheight" / bone
    out_dir.mkdir(parents=True, exist_ok=True)

    frozen = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    fmj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    heights = np.asarray(fmj["heights_m"], dtype=float)
    r1 = json.loads((RES / "bc_robust_metric_route1.json").read_text(encoding="utf-8"))
    dom = {r["bone"]: r for r in r1["rows"]}[bone]["solid_dom"]

    part_row = next(p for p in frozen["parts"] if p["bone"] == bone)
    joint = part_row["joint"]
    fe_solve = float(part_row["fe_solve_load_n"])
    loads = np.asarray(fmj["loads_n"][joint], dtype=float)

    # 该部位在 REGIONS 里的强度阈值 + 子区域口径（与 S4 同源）
    reg = next(r for r in pa.REGIONS if r[1] == bone)
    cn = reg[0]
    subs = reg[3]
    is_long = bone in pa.LONG_BONES
    sig_e = pa._sc(subs[0][1]) if is_long else None
    sig_s = pa._sc(subs[1][1]) if is_long else None
    sig_whole = pa._sc(subs[0][1]) if not is_long else None
    # 冻结 S4 的 ref 统计量（线性外推的对照基线）
    s4 = json.loads((RES / "s4_alignment_bc2.json").read_text(encoding="utf-8"))
    REF_YN = s4["detail_yi"][cn]
    print(f"[ph] 冻结 S4 ref: {json.dumps({k: v for k, v in REF_YN.items() if isinstance(v, (int, float, str))}, ensure_ascii=False)}",
          flush=True)

    part = {p.bone: p for p in fm.build_parts()}[bone]
    pids = [p for sd in part.spec.solid_domains for p in sd.pids]
    pids += [p for sh in part.spec.shell_domains for p in sh.pids]
    mats = fm.read_materials(fm.MASTER, pids)
    for sd in part.spec.solid_domains:
        sd.material = mats[sd.pids[0]]
    for sh in part.spec.shell_domains:
        sh.material = mats[sh.pids[0]]
    print(f"[ph] {bone}: joint={joint} fe_solve={fe_solve:g} N "
          f"domain={dom} cn={cn} long={is_long} "
          f"sig_e={sig_e} sig_s={sig_s} sig_whole={sig_whole}", flush=True)
    t0 = time.time()
    mesh = fm.extract_mesh(part.spec, fm.MASTER)
    print(f"[ph] mesh: {mesh['n_nodes']} nodes ({time.time()-t0:.1f}s)", flush=True)

    solid_names = [sd.mesh_name for sd in part.spec.solid_domains]
    shell_specs = [
        (sh.mesh_name, int(mesh["by_shell"][sh.mesh_name]["q4"].shape[0]
                           + mesh["by_shell"][sh.mesh_name]["t3"].shape[0]))
        for sh in part.spec.shell_domains]

    rows: list[dict] = []
    for i in range(args.height_start, min(args.height_end, len(heights))):
        load = float(loads[i])
        d = out_dir / f"h{i+1:02d}"
        d.mkdir(parents=True, exist_ok=True)
        feb = d / f"{bone}.feb"
        h5 = feb.with_suffix(".hdf5")
        rec = {"i": i, "h_m": float(heights[i]), "load_n": load}
        if not h5.is_file():
            t = time.time()
            try:
                bf.build_bone_feb(part.spec, feb, load_n=load, top_frac=0.30,
                                  bot_frac=0.30, axis=None, with_shell=True,
                                  time_steps=1, mesh=mesh, align_pca=True)
                bf.run_bone_feb(feb, workdir=feb.parent, timeout=1800)
                if feb.with_suffix(".xplt").is_file():
                    bf.post_process(feb.with_suffix(".xplt"), h5,
                                    solid_names=solid_names,
                                    shell_specs=shell_specs, label=bone)
            except Exception as exc:  # noqa: BLE001
                rec["error"] = f"{type(exc).__name__}: {exc}"
            rec["seconds"] = round(time.time() - t, 1)
        if not h5.is_file():
            rec.setdefault("error", "hdf5 未产出")
            rows.append(rec)
            print(f"[ph] h={rec['h_m']:5.1f}m load={load:9.1f}N  FAIL {rec['error']}",
                  flush=True)
            continue
        st = statistic(h5, dom)
        rec.update(st)
        # 逐高度：直接用真解统计量，**不再缩放**
        if is_long:
            rec["ends_ratio"] = st["ends_stat_mpa"] / sig_e
            rec["shaft_ratio"] = st["shaft_stat_mpa"] / sig_s
            rec["frac"] = bool(rec["ends_ratio"] >= 1.0 or rec["shaft_ratio"] >= 1.0)
        else:
            rec["whole_ratio"] = st["shaft_stat_mpa"] / sig_whole
            rec["frac"] = bool(rec["whole_ratio"] >= 1.0)
        # 线性外推的旧读数（对照）：用冻结 S4 的 ref 统计量 × load/fe_solve
        sc = load / fe_solve
        if is_long:
            rec["lin_ends_h"] = float(REF_YN["ends_stat_mpa"]) * sc
            rec["lin_shaft_h"] = float(REF_YN["shaft_stat_mpa"]) * sc
            rec["lin_frac"] = bool((rec["lin_ends_h"] >= sig_e)
                                   or (rec["lin_shaft_h"] >= sig_s))
        rows.append(rec)
        print(f"[ph] h={rec['h_m']:5.1f}m load={load:9.1f}N  "
              f"ends={st['ends_stat_mpa']:8.2f} shaft={st['shaft_stat_mpa']:9.2f}  "
              f"frac={rec['frac']}  ({rec.get('seconds','cached')}s)", flush=True)

    ok = [r for r in rows if "frac" in r]
    nfrac = sum(1 for r in ok if r["frac"])
    print(f"\n[ph] {bone}: {len(ok)} 个高度可解，其中判骨折 {nfrac} 个")
    if ok:
        first = next((r["h_m"] for r in ok if r["frac"]), None)
        print(f"[ph] 首次骨折高度 = {first} m "
              f"(全骨折={'是' if nfrac == len(ok) else '否'})")

    OUT_JSON.write_text(json.dumps({
        "bone": bone, "cn": cn, "joint": joint, "domain": dom,
        "is_long": is_long, "fe_solve_n": fe_solve,
        "thresholds": {"ends": sig_e, "shaft": sig_s, "whole": sig_whole},
        "scaling": "per-height re-solve (no extrapolation)",
        "rows": rows, "n_solved": len(ok), "n_fracture": nfrac,
        "first_fracture_h_m": next((r["h_m"] for r in ok if r["frac"]), None),
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[ph] 写出 {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
