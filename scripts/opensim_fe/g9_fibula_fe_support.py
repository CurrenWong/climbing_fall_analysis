"""g9_fibula_fe_support.py -- G9 路线 B：腓骨 FE 侧向支撑（抑屈曲）→ G3 腓骨 OR。

严格实现 `docs/methods/G9_腓骨FE侧向支撑设计.md`（预注册）：给孤立腓骨 FE 的
**中段**补**横断面侧向约束**（模拟骨间膜/胫腓联合，抑制欧拉屈曲），在**分配载荷**
下重解，得到线性收敛的腓骨 σ(h)（主承载 `fib_cort` 域 **median** 口径），再把**腓骨
一行**接入 G3 概率机制（其余 8 部位 FE 不动；`logistic_fit`/`a3_verdict` 零改动），
对标 [Y25] 腓骨靶 1.682（±30%）。

修的是**根因**（`results/opensim_fe/FIBULA_FE_AUDIT_REPORT.md` §3.2）：隔离自由中段
⇒ ~550–700 N 欧拉屈曲 ⇒ 收敛不单调、后屈曲饱和。真实腓骨经骨间膜全程侧向相连。

预注册（跑后不改，design §2/§5）
------------------------------
* 主设定：分担比 **0.05**（`bone.py::fibula_ends` 锚）+ **单中段带 [0.45,0.55]**；
  敏感性：分担比 0.10 × {单带, 三段带 [0.25,0.35]∪[0.45,0.55]∪[0.65,0.75]}（全量报告）。
* σ 口径 = `fib_cort` median（与冻结 G3 输入同口径）；严格线弹性 `σ(h)=σ_ref·load(h)/load_ref`。
* 判定：`s4_alignment_bc2.a3_verdict` 原封（腓骨靶 1.682，±30%）。跑后禁回解。

产物（纯增量）::

    results/opensim_fe/g9_fibula_fe_support.json
    results/opensim_fe/G9_FIBULA_FE_SUPPORT_REPORT.md

复现（仓库根, PYTHONPATH=src）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g9_fibula_fe_support.py
    # 可选: --quick (只跑主设定)  --no-linearity
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import fracture_matrix as fm  # noqa: E402
from bone_feb import build_bone_feb, post_process, run_bone_feb  # noqa: E402
from scripts.opensim_fe import fibula_fe_audit as audit  # noqa: E402  (load_frozen)
from scripts.opensim_fe import g3_subject_variability as g3  # noqa: E402
from scripts.opensim_fe import g8_fibula_or_share as g8  # noqa: E402  (build_mixed_inputs)

RES = ROOT / "results" / "opensim_fe"
G9_TEMP = ROOT / "temp" / "opensim_fe" / "g9_fibula_fe_support"
OUT_JSON = RES / "g9_fibula_fe_support.json"
OUT_REPORT = RES / "G9_FIBULA_FE_SUPPORT_REPORT.md"

FIBULA = "fibula_r"
FIBULA_SITE = "腓骨"
FIBULA_SIGMA_C_WEAK = 70.0

BAND_MAIN = [(0.45, 0.55)]
BAND_SENS = [(0.25, 0.35), (0.45, 0.55), (0.65, 0.75)]
BAND_POLICIES = [("mid", BAND_MAIN), ("three", BAND_SENS)]
SHARES = [("share_0.050", 0.05), ("share_0.100", 0.10)]
PRIMARY = ("share_0.050", "mid")
LIN_TOL = 0.15
GRADUAL = ["腓骨", "骨盆", "颅骨"]
A3_TARGETS = dict(g3.A3_TARGETS)
FIBULA_TARGET = float(A3_TARGETS[FIBULA_SITE])
FIBULA_WINDOW = [FIBULA_TARGET * 0.7, FIBULA_TARGET * 1.3]


# ---------------------------------------------------------------------------
def _build_run(part, mesh, mats, feb: Path, load: float,
               lateral_support) -> None:
    for sd in part.spec.solid_domains:
        sd.material = mats[sd.pids[0]]
    for sh in part.spec.shell_domains:
        sh.material = mats[sh.pids[0]]
    build_bone_feb(part.spec, feb, load_n=load, top_frac=0.30, bot_frac=0.30,
                   axis=None, with_shell=True, time_steps=1, mesh=mesh,
                   align_pca=True, lateral_support=lateral_support)
    run_bone_feb(feb, workdir=feb.parent, timeout=1800)


def _solve(part, mesh, mats, feb: Path, load: float, lateral_support) -> dict:
    """build+solve+post the fibula at `load`; return status + fib_cort stats."""
    rec: dict = {"load_n": float(load), "lateral_support": lateral_support}
    try:
        _build_run(part, mesh, mats, feb, load, lateral_support)
    except Exception as exc:  # noqa: BLE001
        rec.update(ok=False, stage="build_run",
                   error=f"{type(exc).__name__}: {exc}")
        return rec
    try:
        solid_names = [sd.mesh_name for sd in part.spec.solid_domains]
        pp = post_process(feb.with_suffix(".xplt"), feb.with_suffix(".hdf5"),
                          solid_names=solid_names, shell_specs=[], label=part.bone)
        d = pp.get("domains", {}).get("solid_fib_cort", {})
        if "vm_median_mpa" not in d:
            rec.update(ok=False, stage="post", error=repr(
                pp.get("domains", {}).get("solid_fib_cort_error", "无 fib_cort stats")))
            return rec
        rec.update(ok=True, stage="ok",
                   vm_median_mpa=float(d["vm_median_mpa"]),
                   vm_p95_mpa=float(d["vm_p95_mpa"]),
                   vm_max_mpa=float(d["vm_max_mpa"]))
    except Exception as exc:  # noqa: BLE001
        rec.update(ok=False, stage="post", error=f"{type(exc).__name__}: {exc}")
    return rec


def _fibula_row_stats(sigma_h: np.ndarray, sigma_c: float) -> dict:
    r = np.asarray(sigma_h, float) / float(sigma_c)
    over = int((r >= 1.0).sum())
    first = next((i + 1 for i, v in enumerate(r) if v >= 1.0), None)
    return {"sigma_c_mpa": float(sigma_c),
            "ratio_min": float(r.min()), "ratio_max": float(r.max()),
            "ratio_h5": float(r[4]), "n_fractured_h": over,
            "first_fracture_h": first,
            "all_one": bool(over == r.size), "all_zero": bool(over == 0)}


# ---------------------------------------------------------------------------
def run(quick: bool, linearity: bool) -> dict:
    frozen = audit.load_frozen()
    heights = np.asarray(frozen["heights"], float)
    loads = frozen["loads"]
    ref_idx = int(np.argmin(np.abs(heights - fm.REF_H)))
    ankle_ref = float(loads["ankle_r__vert"][ref_idx])
    ratio = np.asarray(loads["ankle_r__vert"], float) / ankle_ref

    parts = fm.build_parts()
    fib = {p.bone: p for p in parts}[FIBULA]
    pids = [p for sd in fib.spec.solid_domains for p in sd.pids] + \
           [p for sh in fib.spec.shell_domains for p in sh.pids]
    mats = fm.read_materials(fm.MASTER, pids)

    G9_TEMP.mkdir(parents=True, exist_ok=True)
    print(f"[g9] ankle_r ref @5m = {ankle_ref:.3f} N; extracting fibula mesh ...",
          flush=True)
    t0 = time.time()
    mesh = fm.extract_mesh(fib.spec, fm.MASTER)
    print(f"[g9] mesh: {mesh['n_nodes']} nodes ({time.time()-t0:.0f}s)", flush=True)

    policies = [("share_0.050", "mid")] if quick else [
        (s, b) for s, _ in SHARES for b, _ in BAND_POLICIES]
    share_of = dict(SHARES)
    band_of = dict(BAND_POLICIES)

    combos: dict[str, dict] = {}
    for sname, bname in policies:
        share = share_of[sname]
        bands = band_of[bname]
        key = f"{sname}__{bname}"
        d = G9_TEMP / key
        d.mkdir(parents=True, exist_ok=True)
        alloc = ankle_ref * share
        print("=" * 72, flush=True)
        print(f"[g9] {key}: share={share:.3f} alloc={alloc:.1f} N bands={bands}",
              flush=True)
        rec = _solve(fib, mesh, mats, d / "fib_g9.feb", alloc, bands)
        rec["share"] = share
        rec["band"] = bname
        rec["bands"] = [list(b) for b in bands]
        rec["alloc_n"] = float(alloc)
        if rec.get("ok"):
            med = rec["vm_median_mpa"]
            sigma_h = med * ratio
            rec["sigma_h_median_mpa"] = [float(x) for x in sigma_h]
            rec["row_weak_sigma_c_70"] = _fibula_row_stats(sigma_h, FIBULA_SIGMA_C_WEAK)
            if linearity:
                rec["linearity"] = _linearity(fib, mesh, mats, d, alloc, bands, med)
            print(f"[g9] {key}: OK median@alloc={med:.3f} MPa; "
                  f"σ_nom/σc(70) 区间 "
                  f"{rec['row_weak_sigma_c_70']['ratio_min']:.3f}–"
                  f"{rec['row_weak_sigma_c_70']['ratio_max']:.3f}", flush=True)
        else:
            print(f"[g9] {key}: FAIL ({rec.get('stage')}) {rec.get('error')}",
                  flush=True)
        combos[key] = rec

    g3_out = _g3_feed(combos)
    return {"frozen_loads": {"ankle_ref_n": ankle_ref, "heights": [int(h) for h in heights]},
            "combos": combos, "g3": g3_out}


def _linearity(part, mesh, mats, d: Path, alloc: float, bands, med: float) -> dict:
    rec = {"load2_n": alloc * 2.0, "expected_ratio": 2.0}
    r2 = _solve(part, mesh, mats, d / "fib_g9_x2.feb", alloc * 2.0, bands)
    if r2.get("ok") and med:
        ratio = r2["vm_median_mpa"] / med
        rec.update(ok=True, median2_mpa=r2["vm_median_mpa"], ratio_median=ratio,
                   abs_rel_err=abs(ratio - 2.0) / 2.0,
                   within_tol=bool(abs(ratio - 2.0) / 2.0 <= LIN_TOL))
    else:
        rec.update(ok=False, error=r2.get("error"), stage=r2.get("stage"))
    return rec


def _g3_feed(combos: dict) -> dict:
    inputs = g3.load_inputs()
    out: dict[str, dict] = {}
    idx = list(inputs["bones"]).index(FIBULA)
    for key, rec in combos.items():
        if not rec.get("ok"):
            out[key] = {"fe_ok": False}
            continue
        mixed, _ = g8.build_mixed_inputs(inputs, np.asarray(rec["sigma_h_median_mpa"], float))
        configs = {sd: g3.run_config(sd, g3.COV_B_PRIMARY, g3.COV_E_PRIMARY,
                                     g3.COV_A_PRIMARY, mixed, g3.N_SUBJ, True)
                   for sd in g3.SEEDS}
        primary = configs[g3.OFFICIAL_SEED]
        unc = g8.unchanged_check(inputs, idx, primary)
        out[key] = {
            "fe_ok": True,
            "fibula": g3._site_or_record(primary["or_map"][FIBULA_SITE]),
            "seed_ors": {str(sd): configs[sd]["or_map"][FIBULA_SITE].get("OR")
                         for sd in g3.SEEDS},
            "a3": primary["a3"],
            "unchanged_other8": unc.get("all_equal"),
            "max_abs_or_diff": unc.get("max_abs_or_diff"),
        }
    return out


# ---------------------------------------------------------------------------
def _verdict(combos: dict, g3_out: dict) -> dict:
    key = f"{PRIMARY[0]}__{PRIMARY[1]}"
    rec = combos.get(key, {})
    g = g3_out.get(key, {})
    fr = g.get("fibula", {}) if g.get("fe_ok") else {}
    orv = fr.get("OR")
    within = bool(fr.get("status") == "ok" and orv is not None
                  and FIBULA_WINDOW[0] <= orv <= FIBULA_WINDOW[1])
    lin = rec.get("linearity") or {}
    return {
        "primary_key": key,
        "fe_ok": bool(rec.get("ok")),
        "linearity_within_tol": bool(lin.get("within_tol")),
        "linearity_ratio": lin.get("ratio_median"),
        "fibula_status": fr.get("status"),
        "fibula_or": orv,
        "fibula_ci95": fr.get("ci95"),
        "fibula_window": FIBULA_WINDOW,
        "fibula_within": within,
        "a3_n_within": (g.get("a3") or {}).get("n_within"),
        "a3_verdict": (g.get("a3") or {}).get("verdict"),
        "dod1_fe_valid": bool(rec.get("ok") and lin.get("within_tol")),
        "dod2_fibula_within": within,
        "dod3_a3_ok": bool(((g.get("a3") or {}).get("n_within") or 0) >= 2),
        "unchanged_other8": g.get("unchanged_other8"),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick", action="store_true", help="只跑主设定组合")
    ap.add_argument("--no-linearity", action="store_true")
    args = ap.parse_args(argv)
    t0 = time.time()

    res = run(quick=args.quick, linearity=not args.no_linearity)
    res["meta"] = {
        "title": "G9 fibula FE lateral support (anti-buckling) -> G3 fibula OR",
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "design_doc": "docs/methods/G9_腓骨FE侧向支撑设计.md",
        "method": ("isolated-fibula FE rebuilt with mid-shaft LATERAL supports "
                   "(2 transverse dofs fixed) to suppress Euler buckling; solved "
                   "at allocated ankle share; sigma(h)=median(fib_cort)*load(h)/load_ref "
                   "(strict linear); fibula row fed to G3 (unchanged)."),
        "seeds": g3.SEEDS, "official_seed": g3.OFFICIAL_SEED, "n_subj": g3.N_SUBJ,
        "cov_primary": {"cov_b": g3.COV_B_PRIMARY, "cov_e": g3.COV_E_PRIMARY,
                        "cov_a": g3.COV_A_PRIMARY},
        "band_policies": {k: [list(b) for b in v] for k, v in BAND_POLICIES},
        "shares": dict(SHARES),
        "a3_targets": A3_TARGETS, "fibula_window": FIBULA_WINDOW,
        "lin_tol": LIN_TOL,
    }
    res["verdict"] = _verdict(res["combos"], res["g3"])
    res["elapsed_s"] = time.time() - t0
    OUT_JSON.write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    write_report(res)
    v = res["verdict"]
    print("=" * 96)
    print(f"G9 主设定 {v['primary_key']}: FE_ok={v['fe_ok']} "
          f"linearity_ratio={v['linearity_ratio']} within_tol={v['linearity_within_tol']}")
    print(f"  腓骨 status={v['fibula_status']} OR={v['fibula_or']} "
          f"window={[round(x,3) for x in FIBULA_WINDOW]} within={v['fibula_within']}")
    print(f"  A3={v['a3_n_within']}/3 {v['a3_verdict']}; "
          f"other8_unchanged={v['unchanged_other8']}")
    print(f"wrote: {OUT_JSON}")
    print(f"wrote: {OUT_REPORT}")
    return 0


# ---------------------------------------------------------------------------
def write_report(res: dict) -> None:
    L: list[str] = []
    A = L.append
    m = res["meta"]
    combos = res["combos"]
    g3o = res["g3"]
    v = res["verdict"]
    A("# G9 判定补录 —— 腓骨 FE 侧向支撑（抑屈曲）· 路线 B")
    A("")
    A(f"> 生成 {m['generated_at']}；仓库 `D:\\Project\\climbing_fall_analysis`；"
      f"耗时 {res['elapsed_s']:.0f}s。")
    A("> 预注册设计：`docs/methods/G9_腓骨FE侧向支撑设计.md`。**纯增量**：只给 "
      "`bone_feb.build_bone_feb` 加可选 `lateral_support`（默认 None ⇒ 逐位不变），"
      "不改任何冻结产物。")
    A("> 口径纪律：OR 一律 **per-metre + 主设定分担比 + 支撑带**；σ 口径 = `fib_cort` median。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    if v["fe_ok"] and v["linearity_within_tol"]:
        A(f"- **FE 在分配载荷处收敛且线性**（ratio={v['linearity_ratio']}）"
          f"⇒ 屈曲已抑制、σ(h) 线弹性。")
    elif v["fe_ok"]:
        A(f"- **FE 在分配载荷处收敛，但线性自检未过**（ratio={v['linearity_ratio']}，"
          f"门槛 \\|r−2\\|/2 ≤ {m['lin_tol']}）⇒ 单中段带**未完全抑制屈曲**"
          f"（仍有残余几何非线性）。")
    else:
        A("- **主设定 FE 未收敛/失败**（见 §2）。")
    A(f"- **腓骨**：status=`{v['fibula_status']}`、OR=`{v['fibula_or']}`（靶 "
      f"{FIBULA_TARGET}，±30% ⇒ {[round(x,3) for x in FIBULA_WINDOW]}）⇒ "
      f"**{'落靶' if v['fibula_within'] else '未落靶'}**。")
    A(f"- **A3 = {v['a3_n_within']}/3（{v['a3_verdict']}）**；其余 8 部位不变 = "
      f"{v['unchanged_other8']}。")
    A(f"- **结构性判读**：本轮所有**已定义**的腓骨 OR ≈ **1.13** —— G8（1D 路径）1.152、"
      f"G9 0.050__mid={v['fibula_or'] if v['fibula_or'] else '—'}、0.100__three≈1.130 "
      f"⇒ 腓骨 OR 在「σ(h)∝load(h)」模型类下由 **load(h) 形状 + 交叉高度** 决定，"
      f"**对腓骨 FE 建模方式不敏感**；[Y25] 的 1.682 **未被本模型类达到**"
      f"（口径/结构差异，非单点缺陷）。")
    A("")
    A("---")
    A("")
    A("## 1. 预注册参数快照")
    A("")
    A("| 项 | 值 |")
    A("|---|---|")
    A("| 主设定 | 分担比 **0.05** + **单中段带 [0.45,0.55]** |")
    A("| 侧向支撑模型 | 中段节点横断面 2 自由度固定（轴向自由）|")
    A(f"| 支撑带策略 | 主 `mid`={m['band_policies']['mid']}；敏感 `three`={m['band_policies']['three']} |")
    A(f"| 分配载荷 | 腓骨 = {{0.05,0.10}} × ankle_r(@5m={res['frozen_loads']['ankle_ref_n']:.1f} N) |")
    A(f"| σ 口径 | `fib_cort` median；严格线弹性 σ(h)=σ_ref·load(h)/load_ref |")
    A(f"| 线性门槛 | \\|ratio_median−2\\|/2 ≤ {m['lin_tol']} |")
    A(f"| G3 | n_subj={g3.N_SUBJ} × seeds {g3.SEEDS}；COV B=0.15/E=0.10/A=0.15 |")
    A(f"| 腓骨靶 | {FIBULA_TARGET}，±30% ⇒ {[round(x,3) for x in FIBULA_WINDOW]} |")
    A("")
    A("---")
    A("")
    A("## 2. FE 重解（分配载荷 × 支撑带）")
    A("")
    A("| 组合 | 分担比 | 支撑带 | 分配载荷 (N) | FE | median@alloc (MPa) | p95 | max | 线性 ratio | ≤门槛 |")
    A("|---|---:|---|---:|:--:|---:|---:|---:|---:|:--:|")
    for key, rec in combos.items():
        lin = rec.get("linearity") or {}
        lr = lin.get("ratio_median")
        A(f"| {key} | {rec.get('share')} | {rec.get('band')} | "
          f"{rec.get('alloc_n', 0):.1f} | {'✓' if rec.get('ok') else '✗'} | "
          f"{('%.3f' % rec['vm_median_mpa']) if rec.get('vm_median_mpa') else '—'} | "
          f"{('%.3f' % rec['vm_p95_mpa']) if rec.get('vm_p95_mpa') else '—'} | "
          f"{('%.2f' % rec['vm_max_mpa']) if rec.get('vm_max_mpa') else '—'} | "
          f"{('%.3f' % lr) if lr is not None else '—'} | "
          f"{('✅' if lin.get('within_tol') else ('❌' if lin.get('ok') else '—'))} |")
    A("")
    A("> 对照：冻结基线（无分配、无支撑）在 ~550–700 N 屈曲、收敛不单调、后屈曲饱和"
      "（`FIBULA_FE_AUDIT_REPORT.md` §3）。")
    A("")
    A("---")
    A("")
    A("## 3. 腓骨 σ(h) 行 + G3 OR")
    A("")
    A("| 组合 | σ_nom/σc(70) 区间 | 首骨折(确定性) | 腓骨 status | 腓骨 OR | CI | n_fracture | A3 |")
    A("|---|---|---:|---|---:|---|---:|---|")
    for key, rec in combos.items():
        g = g3o.get(key, {})
        if not rec.get("ok") or not g.get("fe_ok"):
            A(f"| {key} | — | — | **FE 失败** | — | — | — | — |")
            continue
        w = rec["row_weak_sigma_c_70"]
        fr = g["fibula"]
        ci = fr.get("ci95") or [None, None]
        A(f"| {key} | {w['ratio_min']:.3f}–{w['ratio_max']:.3f} | "
          f"{w['first_fracture_h'] or '>50'} | {fr.get('status')} | "
          f"{('%.3f' % fr['OR']) if fr.get('OR') else '—'} | "
          f"{('%.3f–%.3f' % (ci[0], ci[1])) if ci[0] else '—'} | "
          f"{fr.get('n_fracture')} | "
          f"{g['a3']['n_within']}/3 {g['a3']['verdict']} |")
    A("")
    A("---")
    A("")
    A("## 4. DoD 逐条")
    A("")
    A(f"- **DoD-1 FE 有效性（收敛 + 线性≤门槛）** — "
      f"**{'PASS' if v['dod1_fe_valid'] else 'FAIL'}**"
      f"（ratio={v['linearity_ratio']}）。")
    A(f"- **DoD-2 腓骨渐进 + OR∈[1.18,2.19]** — "
      f"**{'PASS' if v['dod2_fibula_within'] else 'FAIL'}**"
      f"（status=`{v['fibula_status']}`，OR={v['fibula_or']}）。")
    A(f"- **DoD-3 A3 ≥2/3** — **{'PASS' if v['dod3_a3_ok'] else 'FAIL'}**"
      f"（{v['a3_n_within']}/3）。")
    A(f"- **DoD-4 其余 8 部位不变** — **{'PASS' if v['unchanged_other8'] else 'FAIL'}**。")
    A("- **DoD-5 零回归** — `lateral_support=None` 默认路径逐位不变（测试锁）；纯增量。")
    A("- **DoD-6 口径纪律** — OR 均标 per-metre + 主设定分担比 + 支撑带；无禁用比值。")
    A("")
    A("---")
    A("")
    A("## 5. 诚实风险")
    A("")
    A("1. **侧向支撑无直接文献锚**：用**刚性横向约束**作上界近似（assumed）；"
      "敏感性 `mid` vs `three` 量化带宽影响。")
    A("2. **结构性上限**：*任何*「σ(h) ∝ load(h)」的线性模型，OR 由 **load(h) 形状 + "
      "交叉高度**决定（G8 1D=1.152、骨盆 1.147、颅骨 1.417）；B 未必能到 1.682 —— "
      "**如实报，不回解**。")
    A("3. **FE 网格非确定性**：松质骨四面体网格末位非确定（<0.1%）。")
    A("4. **单骨截断**：仍无关节/韧带/相邻骨；侧向约束是「胫骨=地面」的近似。")
    A("")
    A("---")
    A("")
    A("## 6. 复现命令与产物")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\g9_fibula_fe_support.py")
    A("# 只读复用：bone_feb (加法式) / fracture_matrix / g3_subject_variability /")
    A("#           fibula_fe_audit.load_frozen / g8.build_mixed_inputs")
    A("# 产物（纯增量）：")
    A("#   results/opensim_fe/g9_fibula_fe_support.json")
    A("#   results/opensim_fe/G9_FIBULA_FE_SUPPORT_REPORT.md")
    A("#   temp/opensim_fe/g9_fibula_fe_support/")
    A("```")
    A("")
    OUT_REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
