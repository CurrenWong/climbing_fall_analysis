# -*- coding: utf-8 -*-
"""FE_BC · F-2 决定性验证：股骨的「解一次 + 线性外推」到底错多少？

动机（`results/opensim_fe/FE_BC_PROBE.md` §4）：
  冻结 `femur_r_fe.json` 自带线性自检 —— 在 `ref_load` 与 `ref_load/2` 各解一次，
  得 `ratio_median=5.79` / `ratio_max=9.80`（**期望 2.0**）。
  而 A1/A2 流水线是**线性外推**（`fracture_matrix.py:606`: `sigma(h)=sigma_ref*load(h)/ref_load`）
  ⇒ 若真非线性，则低载端被**系统性高估** ⇒ 直接解释 A2 股骨 `all_one`。

本脚本：对股骨在同一网格上解一条**载荷阶梯**，把**实测 σ(L)** 与
**线性外推预测 σ(ref)·(L/ref)** 并排，量出误差倍数。

纪律：冻结输入只读；输出只在 `temp/fe_out/fe_bc_linearity/` + 一个 results json。
BC 固定为 BC2 条件（`top_frac=0.30, bot_frac=0.30, align_pca=True, time_steps=1`），
**本脚本不改任何 BC** —— 只测「线性外推」这一个假设。
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

BONE = "femur_r"
OUT_DIR = ROOT / "temp" / "opensim_fe" / "fe_bc_linearity"
OUT_JSON = ROOT / "results" / "opensim_fe" / "FE_BC_LINEARITY.json"


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
            "vm_median_mpa": float(np.median(vm)),
            "vm_p95_mpa": float(np.percentile(vm, 95)),
            "vm_max_mpa": float(vm.max())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fractions", default="1,0.5,0.25,0.125",
                    help="ref_load 的倍数（逗号分隔）")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    frozen = json.loads(
        (ROOT / "results" / "opensim_fe" / "fracture_matrix_bc2.json")
        .read_text(encoding="utf-8"))
    ref = next(p for p in frozen["parts"] if p["bone"] == BONE)
    rl = float(ref["fe_solve_load_n"])
    r1 = json.loads((ROOT / "results" / "opensim_fe" / "bc_robust_metric_route1.json")
                    .read_text(encoding="utf-8"))
    dom = {r["bone"]: r for r in r1["rows"]}[BONE]["solid_dom"]

    part = {p.bone: p for p in fm.build_parts()}[BONE]
    pids = [p for sd in part.spec.solid_domains for p in sd.pids]
    pids += [p for sh in part.spec.shell_domains for p in sh.pids]
    print(f"[lin] materials for {len(pids)} pids ...", flush=True)
    mats = fm.read_materials(fm.MASTER, pids)
    for sd in part.spec.solid_domains:
        sd.material = mats[sd.pids[0]]
    for sh in part.spec.shell_domains:
        sh.material = mats[sh.pids[0]]
    print("[lin] extracting mesh (once) ...", flush=True)
    t0 = time.time()
    mesh = fm.extract_mesh(part.spec, fm.MASTER)
    print(f"[lin] mesh: {mesh['n_nodes']} nodes ({time.time()-t0:.1f}s)", flush=True)

    fracs = [float(x) for x in args.fractions.split(",") if x.strip()]
    rows: list[dict] = []
    for fr in fracs:
        load = rl * fr
        d = OUT_DIR / f"f{fr:g}"
        d.mkdir(parents=True, exist_ok=True)
        feb = d / f"{BONE}.feb"
        h5 = feb.with_suffix(".hdf5")
        rec = {"fraction": fr, "load_n": load}
        if not h5.is_file():
            t = time.time()
            bf.build_bone_feb(part.spec, feb, load_n=load, top_frac=0.30,
                              bot_frac=0.30, axis=None, with_shell=True,
                              time_steps=1, mesh=mesh, align_pca=True)
            bf.run_bone_feb(feb, workdir=feb.parent, timeout=1800)
            # .hdf5 是 post_process 写出的（xplt -> hdf5），不是 run 直接产
            if feb.with_suffix(".xplt").is_file():
                bf.post_process(
                    feb.with_suffix(".xplt"), h5,
                    solid_names=[sd.mesh_name for sd in part.spec.solid_domains],
                    shell_specs=[
                        (sh.mesh_name,
                         int(mesh["by_shell"][sh.mesh_name]["q4"].shape[0]
                             + mesh["by_shell"][sh.mesh_name]["t3"].shape[0]))
                        for sh in part.spec.shell_domains],
                    label=BONE)
            rec["seconds"] = round(time.time() - t, 1)
            if not h5.is_file():
                rec["error"] = "hdf5 未产出"
                print(f"[lin] f={fr:g} load={load:.0f}N FAIL", flush=True)
                rows.append(rec)
                continue
        rec.update(statistic(h5, dom))
        rows.append(rec)
        print(f"[lin] f={fr:<6g} load={load:9.1f}N  ends={rec['ends_stat_mpa']:8.2f}"
              f"  shaft={rec['shaft_stat_mpa']:8.2f}  vm_max={rec['vm_max_mpa']:9.2f}"
              f"  ({rec.get('seconds','cached')}s)", flush=True)

    base = next((r for r in rows if r["fraction"] == 1.0), None)
    if base:
        print("\n=== 线性外推误差（实测 vs σ(ref)·f） ===")
        print(f"{'f':>7} {'load_N':>10} {'shaft_meas':>11} {'shaft_lin':>11} "
              f"{'lin/meas':>9} {'ends_meas':>10} {'ends_lin':>10}")
        for r in rows:
            if "shaft_stat_mpa" not in r:
                continue
            sl = base["shaft_stat_mpa"] * r["fraction"]
            el = base["ends_stat_mpa"] * r["fraction"]
            r["shaft_linear_pred_mpa"] = sl
            r["ends_linear_pred_mpa"] = el
            r["linear_overprediction_shaft"] = sl / r["shaft_stat_mpa"]
            r["linear_overprediction_ends"] = el / r["ends_stat_mpa"]
            print(f"{r['fraction']:>7g} {r['load_n']:>10.1f} "
                  f"{r['shaft_stat_mpa']:>11.2f} {sl:>11.2f} "
                  f"{sl/r['shaft_stat_mpa']:>9.2f} "
                  f"{r['ends_stat_mpa']:>10.2f} {el:>10.2f}")

    OUT_JSON.write_text(json.dumps({
        "bone": BONE, "domain": dom, "ref_load_n": rl,
        "bc": {"top_frac": 0.30, "bot_frac": 0.30, "align_pca": True,
               "time_steps": 1, "rules": "linear scaling: sigma(h)=sigma_ref*load(h)/ref_load"},
        "rows": rows, "reads_frozen_only": True,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n[lin] 写出 {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
