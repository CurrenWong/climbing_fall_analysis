"""g9_sigma_level_sens.py -- G9 腓骨 σ 水平口径敏感性（同解、多统计量、同一 G3 链）。

**动机**（`docs/methods/G9_腓骨FE侧向支撑设计.md` §11.6/§11.7）：
`fib_cort` 的应力分布极偏斜（max≈1 100 MPa 而带符号 median≈0）⇒ `median` 落在噪声地板、
**刀锋敏感**（同载荷、deck 逐位相同的两次解 `median|σ_zz|` 差 23%，2× 载荷只升 1.15×）。
载荷流正对照证明**两个解都满足平衡**（自由区载流 92%/91% ≈ 载荷、`A_eff` 76.01 mm² 两解一致）
⇒ 差异来自统计量，不是解。

**本脚本回答**：换哪个统计量/口径，才能既可复现、又不让腓骨退化到 `all_zero`/`all_one`？
候选（同一 FE 解）：

| 键 | 口径 | 性质 |
|---|---|---|
| `vm_median` | `median` σ_vm | **G3 现冻结口径**（基准，不覆盖） |
| `abs_szz_median` | `median` \\|σ_zz\\| | §11.6 证伪（不可复现） |
| `abs_szz_p95` | `p95` \\|σ_zz\\| | 偏斜稳健，但绝对水平高 |
| `szz_mean_signed` | `mean` σ_zz（体积平均） | 最稳健（两解差 0.5%） |
| `F_over_Aeff` | `F/A_eff`（**平衡锚定**） | `A_eff` 实测且可由 `Σσ·A = F` 独立验证 |

每档都接**同一 G3 机制**（`g8.build_mixed_inputs` + `g3.run_config`，种子/COV/判定零改动），
输出 `util_max`、`n_fractured_h`、`first_fracture_h`、`status`、`OR`、`ci95`。

复现（仓库根, PYTHONPATH=src）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g9_sigma_level_sens.py
    # 可选: --bands long95 --shares 0.050,0.100 --out-suffix _x

产物（纯增量）::

    results/opensim_fe/g9_sigma_level_sens.json
    results/opensim_fe/G9_SIGMA_LEVEL_SENS_REPORT.md
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import bone_feb as bf  # noqa: E402
import fracture_matrix as fm  # noqa: E402
from scripts.opensim_fe import fibula_fe_audit as audit  # noqa: E402
from scripts.opensim_fe import g3_subject_variability as g3  # noqa: E402
from scripts.opensim_fe import g8_fibula_or_share as g8  # noqa: E402
from scripts.opensim_fe import g9_fibula_fe_support as g9  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
TEMP = ROOT / "temp" / "opensim_fe" / "g9_sigma_level_sens"
FIBULA, FIBULA_SITE = g9.FIBULA, g9.FIBULA_SITE
SIGMA_C = g9.FIBULA_SIGMA_C_WEAK          # 70 MPa（g9 冻结；沿用）
DOM = "fib_cort"
# 自由区/固定区/载荷带（deck 用 top_frac=bot_frac=0.30，见 g9._build_run）
FREE = (0.30, 0.70)
LEVELS = ["vm_median", "abs_szz_median", "abs_szz_p95", "szz_elem_mean_all",
          "szz_free_mean_vol", "F_over_Aeff"]


# ---------------------------------------------------------------------------
def _deck_geometry(feb: Path) -> dict:
    """从 .feb 取节点与 `fib_cort` 的 hex8 体积/形心（用仓库 `_hex8_signed_volumes`）。"""
    t = feb.read_text(encoding="ISO-8859-1")
    nodes = {}
    for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', t):
        nodes[int(m.group(1))] = [float(x) for x in m.group(2).split(",")]
    arr = np.array([nodes[i] for i in sorted(nodes)])
    conn = None
    for m in re.finditer(r'<Elements\s+([^>]*)>(.*?)</Elements>', t, re.S):
        at = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        if at.get("name") != DOM:
            continue
        assert at.get("type") == "hex8", f"{DOM} 不是 hex8？type={at.get('type')}"
        conn = np.array([[int(x) for x in e.group(1).split(",")]
                         for e in re.finditer(r'<elem id="\d+">([^<]+)</elem>',
                                               m.group(2))])
    assert conn is not None, f"deck 里没找到域 {DOM}"
    vol = np.abs(bf._hex8_signed_volumes(arr, conn))      # 1-based 局部索引
    cen = arr[conn - 1].mean(axis=1)
    return {"vol": vol, "z": cen[:, 2], "n_elem": int(len(vol))}


def _section(feb: Path, hdf5: Path) -> dict:
    """自由区载流与 A_eff（正对照：Σσ_zz·V/Δz 应 ≈ 施加载荷）。"""
    g = _deck_geometry(feb)
    with h5py.File(hdf5, "r") as h5:
        sid = sorted(h5["states"].keys(), key=int)[-1]
        s = np.asarray(h5[f"states/{sid}/element_data/stress/{DOM}"][:], float)
    assert s.shape[0] == g["n_elem"], f"应力/几何元素数不符 {s.shape[0]} vs {g['n_elem']}"
    z, vol = g["z"], g["vol"]
    z0, z1 = z.min(), z.max()
    lo, hi = z0 + FREE[0] * (z1 - z0), z0 + FREE[1] * (z1 - z0)
    m = (z >= lo) & (z < hi)
    dz = hi - lo
    szz = s[:, 2]
    return {"free_force_n": float((szz[m] * vol[m]).sum() / dz),
            "A_eff_mm2": float(vol[m].sum() / dz),
            "z_span_mm": [float(z0), float(z1)], "dz_free_mm": float(dz)}


def _levels(feb: Path, hdf5: Path, rec: dict, alloc: float, sec: dict) -> dict:
    """同一 FE 解上取各档 σ 水平（MPa）。

    ⚠️ `szz_elem_mean_all` = **未加权**全体单元均值（含固定区的 0 值）⇒ 会被稀释，
    **不作口径用**，保留只为暴露"用错加权会得到 ~10 MPa 的假水平"。
    `szz_free_mean_vol` = **自由区体积加权**均值（符号量，压缩为负），
    按定义应 = `−F/A_eff` ⇒ 是平衡锚定的**交叉验证**，二者应几乎相等。
    """
    g = _deck_geometry(feb)
    with h5py.File(hdf5, "r") as h5:
        sid = sorted(h5["states"].keys(), key=int)[-1]
        s = np.asarray(h5[f"states/{sid}/element_data/stress/{DOM}"][:], float)
    szz, az = s[:, 2], np.abs(s[:, 2])
    z, vol = g["z"], g["vol"]
    z0, z1 = z.min(), z.max()
    lo, hi = z0 + FREE[0] * (z1 - z0), z0 + FREE[1] * (z1 - z0)
    m = (z >= lo) & (z < hi)
    return {
        "vm_median": float(rec["vm_median_mpa"]),
        "abs_szz_median": float(np.median(az)),
        "abs_szz_p95": float(np.percentile(az, 95)),
        "szz_elem_mean_all": float(np.mean(szz)),
        "szz_free_mean_vol": float((szz[m] * vol[m]).sum() / vol[m].sum()),
        # 平衡锚定：自由区实测截面上的名义压缩应力（正值）
        "F_over_Aeff": float(alloc / sec["A_eff_mm2"]),
    }


def _row_stats(sigma_h: np.ndarray) -> dict:
    return g9._fibula_row_stats(sigma_h, SIGMA_C)


def _or_for(row: np.ndarray, inputs: dict) -> dict:
    mixed, _ = g8.build_mixed_inputs(inputs, np.asarray(row, float))
    configs = {sd: g3.run_config(sd, g3.COV_B_PRIMARY, g3.COV_E_PRIMARY,
                                 g3.COV_A_PRIMARY, mixed, g3.N_SUBJ, True)
               for sd in g3.SEEDS}
    prim = configs[g3.OFFICIAL_SEED]
    fr = g3._site_or_record(prim["or_map"][FIBULA_SITE])
    return {"status": fr.get("status"), "OR": fr.get("OR"), "ci95": fr.get("ci95"),
            "n_fracture": fr.get("n_fracture"), "n": fr.get("n"),
            "a3_n_within": (prim.get("a3") or {}).get("n_within")}


# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bands", default="long95", help="支撑带策略（D1a′ 默认 long95）")
    ap.add_argument("--shares", default="0.050,0.100")
    ap.add_argument("--out-suffix", default="")
    args = ap.parse_args(argv)

    band = g9.BAND_ALL[args.bands]
    shares = [float(x) for x in args.shares.split(",") if x.strip()]
    t0 = time.time()

    frozen = audit.load_frozen()
    heights = np.asarray(frozen["heights"], float)
    ref_idx = int(np.argmin(np.abs(heights - fm.REF_H)))
    ankle_ref = float(frozen["loads"]["ankle_r__vert"][ref_idx])
    ratio = np.asarray(frozen["loads"]["ankle_r__vert"], float) / ankle_ref

    parts = fm.build_parts()
    fib = {p.bone: p for p in parts}[FIBULA]
    pids = [p for sd in fib.spec.solid_domains for p in sd.pids] + \
           [p for sh in fib.spec.shell_domains for p in sh.pids]
    mats = fm.read_materials(fm.MASTER, pids)

    TEMP.mkdir(parents=True, exist_ok=True)
    print(f"[sens] ankle_r@ref={ankle_ref:.3f} N  bands={band}", flush=True)
    mesh = fm.extract_mesh(fib.spec, fm.MASTER)
    print(f"[sens] mesh {mesh['n_nodes']} nodes ({time.time()-t0:.0f}s)", flush=True)

    inputs = g3.load_inputs()
    out: dict = {}
    for sh in shares:
        key = f"share_{sh:.3f}__{args.bands}"
        d = TEMP / key
        d.mkdir(parents=True, exist_ok=True)
        alloc = ankle_ref * sh
        rec = g9._solve(fib, mesh, mats, d / "fib_sens.feb", alloc, band, szz=True)
        if not rec.get("ok"):
            out[key] = {"ok": False, "share": sh, "alloc_n": alloc,
                        "error": rec.get("error"), "stage": rec.get("stage")}
            print(f"[sens] {key}: FE FAIL ({rec.get('stage')}) {rec.get('error')}",
                  flush=True)
            continue
        sec = _section(d / "fib_sens.feb", d / "fib_sens.hdf5")
        lev = _levels(d / "fib_sens.feb", d / "fib_sens.hdf5", rec, alloc, sec)
        per: dict = {}
        for name in LEVELS:
            if not np.isfinite(lev[name]):
                continue
            row = lev[name] * ratio
            st = _row_stats(row)
            per[name] = {"level_mpa": lev[name],
                         "sigma_h_mpa": [float(x) for x in row],
                         "row": st, "g3": _or_for(row, inputs)}
        out[key] = {"ok": True, "share": sh, "alloc_n": float(alloc),
                    "bands": [list(b) for b in band],
                    "equilibrium": {**sec,
                                    "load_n": float(alloc),
                                    "force_over_load": sec["free_force_n"] / alloc},
                    "vm_max_mpa": rec["vm_max_mpa"], "levels": per}
        print(f"[sens] {key}: A_eff={sec['A_eff_mm2']:.2f} mm²  "
              f"载流={sec['free_force_n']:.1f} N ({sec['free_force_n']/alloc*100:.1f}% of "
              f"{alloc:.1f} N)", flush=True)
        for name in LEVELS:
            if name not in per:
                continue
            p, g = per[name], per[name]["g3"]
            print(f"        {name:17s} level={p['level_mpa']:9.3f} MPa  "
                  f"util_max={p['row']['ratio_max']:7.3f}  "
                  f"骨折 h={p['row']['n_fractured_h']:3d} "
                  f"(首折 {p['row']['first_fracture_h']})  "
                  f"status={g['status']:8s} OR={g['OR']}", flush=True)

    res = {"meta": {"title": "G9 fibula sigma-level sensitivity (same FE solution)",
                    "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "bands": [list(b) for b in band],
                    "sigma_c_mpa": SIGMA_C, "levels": LEVELS,
                    "free_zone_zfrac": list(FREE),
                    "seeds": g3.SEEDS, "official_seed": g3.OFFICIAL_SEED,
                    "n_subj": g3.N_SUBJ,
                    "design_doc": "docs/methods/G9_腓骨FE侧向支撑设计.md §11.6/§11.7",
                    "frozen_measure_untouched": "vm_median（G3 现冻结口径）仅作基准，本脚本不覆盖它"},
           "combos": out, "elapsed_s": time.time() - t0}
    (RES / f"g9_sigma_level_sens{args.out_suffix}.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_report(res, RES / f"G9_SIGMA_LEVEL_SENS_REPORT{args.out_suffix}.md")
    print("=" * 100)
    print(f"wrote: {RES / f'g9_sigma_level_sens{args.out_suffix}.json'}")
    print(f"wrote: {RES / f'G9_SIGMA_LEVEL_SENS_REPORT{args.out_suffix}.md'}")
    return 0


def _write_report(res: dict, path: Path) -> None:
    L = ["# G9 · 腓骨 σ 水平口径敏感性（同一 FE 解 × 5 档统计量 × 同一 G3 链）", "",
         f"生成 {res['meta']['generated_at']}｜支撑带 {res['meta']['bands']}｜"
         f"σ_c={res['meta']['sigma_c_mpa']} MPa｜种子 {res['meta']['seeds']}"
         f"（official={res['meta']['official_seed']}，n={res['meta']['n_subj']}）", "",
         "## 1. 正对照（载荷流 / 截面）", "",
         "| 组合 | 分配载荷 | 自由区 A_eff (mm²) | 自由区载流 (N) | 载流/载荷 |",
         "|---|---|---|---|---|"]
    for k, v in res["combos"].items():
        if not v.get("ok"):
            L.append(f"| {k} | — | FE FAIL | — | — |")
            continue
        e = v["equilibrium"]
        L.append(f"| {k} | {e['load_n']:.1f} | {e['A_eff_mm2']:.2f} | "
                 f"{e['free_force_n']:.1f} | **{e['force_over_load']*100:.1f}%** |")
    L += ["", "> 判据：自由区任意截面载流 ≈ 分配载荷，且 `A_eff` 在各解间重现。", "",
          "## 2. 各口径的判定", "",
          "| 组合 | 口径 | level@alloc (MPa) | util_max | 骨折 h | 首折 h | status | OR |",
          "|---|---|---|---|---|---|---|---|"]
    for k, v in res["combos"].items():
        if not v.get("ok"):
            continue
        for name, p in v["levels"].items():
            r, g = p["row"], p["g3"]
            L.append(f"| {k} | `{name}` | {p['level_mpa']:.3f} | {r['ratio_max']:.3f} | "
                     f"{r['n_fractured_h']} | {r['first_fracture_h']} | {g['status']} | "
                     f"{g['OR']} |")
    L += ["", "## 3. 读法", "",
          "- **`vm_median` 是 G3 现冻结口径**，此处仅作基准，**不覆盖**。",
          "- `abs_szz_median` 在 §11.6 已证伪（同载荷两解差 23%、2× 载荷只升 1.15×）。",
          "- `abs_szz_p95` / `szz_free_mean_vol` 可复现，但 p95 **不满足 DoD-1 的线性要求**"
          "（2× 载荷时 p95 升 11.2× —— 分位数随【受载区扩大】而跳），故不作 σ(h) 的水平。",
          "- `szz_elem_mean_all` 是**未加权**全体均值（含固定区 0 值）⇒ 被稀释约 1.7×，"
          "**不作口径**，只用于暴露加权错误。",
          "- `szz_free_mean_vol` = 自由区体积加权轴向均值（压缩为负）⇒ 按定义 = `−F/A_eff` ⇒ "
          "与 `F_over_Aeff` 互为交叉验证。",
          "- `F_over_Aeff` 由**实测截面**与**平衡**锚定，A_eff 各解重现；水平可被载流独立验证。",
          "- OR 由 `load(h)` 形状 + 交叉高度决定 ⇒ 换统计量改的是 **util 的绝对水平**；"
          "口径选择须避免把腓骨推成 `all_zero` 或 `all_one`。", ""]
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())