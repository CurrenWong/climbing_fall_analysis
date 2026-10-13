# -*- coding: utf-8 -*-
"""FE_BC · W4-a：让股骨在 **8–58 kN 全域可解**（F-1 复工，目标已改写）。

背景（`results/opensim_fe/FE_BC_PERHEIGHT.md` §3.2）：
  现 BC = `bot_frac=0.30` 的**全 3-DOF 刚性夹持带** ⇒ 固定端弯矩 ⇒
  **≥20 kN 起 FEBio 不收敛（37/50 高度失败）** + 高载端解分支非单调。
  ⇒ 逐高度重解（修 F-2 的前置）**不可交付**。

**验收门槛（可判定）**：全域收敛 —— 50/50 高度 `rc=0`。
本脚本先在**最难的高度**试变体（能过最难的就大概率全域可解），再对胜者跑全扫。

变体（全部**加法式**，默认路径不受影响）：
  A0  baseline            : 全 3-DOF 夹持（现状）        + 无侧向支撑
  A1  axial_only          : 仅轴向夹持                   + 无侧向支撑（可能刚体奇异）
  A2  axial_lat_90        : 仅轴向夹持                   + 侧向 [0.05,0.95]
  A3  full_lat_90         : 全 3-DOF 夹持                + 侧向 [0.05,0.95]
  A4  axial_lat_80        : 仅轴向夹持                   + 侧向 [0.10,0.90]
  A5  axial_lat_full      : 仅轴向夹持                   + 侧向 [0.02,0.98]

纪律：冻结输入只读；输出只在 `temp/opensim_fe/fe_bc_stab/` + 一个 results json。
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
OUT_JSON = RES / "FE_BC_STAB.json"

VARIANTS: list[tuple[str, bool, list | None, str]] = [
    ("A0_baseline",      False, None,               "现状：全 3-DOF 夹持 + 无侧向"),
    ("A1_axial_only",    True,  None,               "仅轴向夹持 + 无侧向（可能刚体奇异）"),
    ("A2_axial_lat_90",  True,  [(0.05, 0.95)],     "仅轴向夹持 + 侧向 90%"),
    ("A3_full_lat_90",   False, [(0.05, 0.95)],     "全 3-DOF 夹持 + 侧向 90%"),
    ("A4_axial_lat_80",  True,  [(0.10, 0.90)],     "仅轴向夹持 + 侧向 80%"),
    ("A5_axial_lat_full", True, [(0.02, 0.98)],     "仅轴向夹持 + 侧向 96%"),
]


def statistic(h5: Path, dom: str) -> dict:
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
    ap.add_argument("--heights", default="50,20,8,5,1",
                    help="逗号分隔的 1-based 高度号（先难后易）")
    ap.add_argument("--variants", default="")
    args = ap.parse_args()

    bone = args.bone
    base_dir = ROOT / "temp" / "opensim_fe" / "fe_bc_stab" / bone
    base_dir.mkdir(parents=True, exist_ok=True)

    fmj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    heights = np.asarray(fmj["heights_m"], dtype=float)
    r1 = json.loads((RES / "bc_robust_metric_route1.json").read_text(encoding="utf-8"))
    dom = {r["bone"]: r for r in r1["rows"]}[bone]["solid_dom"]
    part_row = next(p for p in fmj["parts"] if p["bone"] == bone)
    joint = part_row["joint"]
    loads = np.asarray(fmj["loads_n"][joint], dtype=float)

    part = {p.bone: p for p in fm.build_parts()}[bone]
    pids = [p for sd in part.spec.solid_domains for p in sd.pids]
    pids += [p for sh in part.spec.shell_domains for p in sh.pids]
    mats = fm.read_materials(fm.MASTER, pids)
    for sd in part.spec.solid_domains:
        sd.material = mats[sd.pids[0]]
    for sh in part.spec.shell_domains:
        sh.material = mats[sh.pids[0]]
    t0 = time.time()
    mesh = fm.extract_mesh(part.spec, fm.MASTER)
    print(f"[stab] {bone}: mesh {mesh['n_nodes']} nodes ({time.time()-t0:.1f}s), "
          f"domain={dom}, joint={joint}", flush=True)

    solid_names = [sd.mesh_name for sd in part.spec.solid_domains]
    shell_specs = [
        (sh.mesh_name, int(mesh["by_shell"][sh.mesh_name]["q4"].shape[0]
                           + mesh["by_shell"][sh.mesh_name]["t3"].shape[0]))
        for sh in part.spec.shell_domains]

    hlist = [int(x) for x in args.heights.split(",") if x.strip()]
    rows: list[dict] = []
    for hno in hlist:
        i = hno - 1
        load = float(loads[i])
        for name, axial_only, lat, note in VARIANTS:
            if args.variants and args.variants not in name:
                continue
            d = base_dir / f"{name}_h{hno:02d}"
            d.mkdir(parents=True, exist_ok=True)
            feb = d / f"{bone}.feb"
            h5 = feb.with_suffix(".hdf5")
            rec = {"variant": name, "note": note, "h_no": hno,
                   "h_m": float(heights[i]), "load_n": load,
                   "axial_only": axial_only,
                   "lateral_support": [list(b) for b in lat] if lat else None}
            if not h5.is_file():
                t = time.time()
                try:
                    bf.build_bone_feb(part.spec, feb, load_n=load, top_frac=0.30,
                                      bot_frac=0.30, axis=None, with_shell=True,
                                      time_steps=1, mesh=mesh, align_pca=True,
                                      lateral_support=lat,
                                      bottom_fix_axial_only=axial_only)
                    bf.run_bone_feb(feb, workdir=feb.parent, timeout=1800)
                    if feb.with_suffix(".xplt").is_file():
                        bf.post_process(feb.with_suffix(".xplt"), h5,
                                        solid_names=solid_names,
                                        shell_specs=shell_specs, label=bone)
                    rec["converged"] = True
                except Exception as exc:  # noqa: BLE001
                    rec["converged"] = False
                    rec["error"] = f"{type(exc).__name__}: {exc}"
                rec["seconds"] = round(time.time() - t, 1)
            else:
                rec["converged"] = True
            if h5.is_file():
                rec.update(statistic(h5, dom))
            rows.append(rec)
            print(f"[stab] h={hno:2d} ({load:9.1f}N) {name:20s} "
                  f"{'OK ' if rec.get('converged') else 'FAIL'} "
                  + (f"ends={rec['ends_stat_mpa']:8.2f} shaft={rec['shaft_stat_mpa']:9.2f}"
                     if "shaft_stat_mpa" in rec else f"({rec.get('error','')[:60]})")
                  + f"  [{rec.get('seconds','cached')}s]", flush=True)

    OUT_JSON.write_text(json.dumps(
        {"bone": bone, "domain": dom, "joint": joint,
         "gate": "full-domain convergence (rc=0 at all 50 heights)",
         "variants": [{"name": n, "axial_only": a, "lateral_support": l, "note": s}
                      for n, a, l, s in VARIANTS],
         "rows": rows}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n[stab] 写出 {OUT_JSON}")
    # 小结：每个变体在测试过的最高载荷上是否收敛
    print("\n=== 小结（逐变体，测试过的高度） ===")
    for name, *_ in VARIANTS:
        rs = [r for r in rows if r["variant"] == name]
        if not rs:
            continue
        ok = [r["h_no"] for r in rs if r.get("converged")]
        bad = [r["h_no"] for r in rs if not r.get("converged")]
        print(f"  {name:20s} OK h={ok}  FAIL h={bad}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
