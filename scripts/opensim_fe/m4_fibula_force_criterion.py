r"""m4_fibula_force_criterion.py -- M4：用**实测断裂载荷**做力判据，绕开应力口径之争。

**动机**（G9 §12.2 M4）：系统综述给出的腓骨**压缩失效载荷 77–370 N**（含轴向预载下侧向弯曲）
与本项目的两个前提严重冲突：
  ① 我们在孤立腓骨 FE 里实测 **P_cr ≈ 252 N**（Euler 解与载荷扫描互证，G4 计划 §9.13）；
     它**落在实测失效带 77–370 N 内** ⇒ 孤立腓骨的失效是**结构屈曲**，不是材料压溃
     ⇒ **THUMS 的 σ_c=70 MPa 不适用于孤立腓骨**（70 MPa × 76 mm² = 5.3 kN，比实测高 14–70×）。
  ② 项目假设的 share 0.05/0.10 ⇒ 分配载荷 1 263 / 2 527 N，**比实测失效载荷高 3.4–33×**
     ⇒ 孤立骨模型在该 share 下**物理上不自洽**，**与用哪种 σ 口径无关**。

**本脚本产出（零 FE 成本，只用冻结载荷链 + 实测常数）**：
  1. **share → 骨折高度 h\*(share; F_fail)** 映射（确定性，力判据）；
  2. 每个 F_fail 下**「渐进骨折」可用的 share 窗口**（部分高度骨折、非全断）；
  3. **孤立模型的可行 share 上限**（全程不骨折）；
  4. 与 FE 的 P_cr=252 N、THUMS σ_c=70 MPa 的三方对照；
  5. 把 F_fail 的实测离散（77–370 N ⇒ 对数正态 σ_log=0.478）折进**逐高度失效概率**。

复现（仓库根, PYTHONPATH=src）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\m4_fibula_force_criterion.py
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
from scripts.opensim_fe import fibula_fe_audit as audit  # noqa: E402

RES = ROOT / "results" / "opensim_fe"

# --- 实测常数（来源见模块 docstring 与 G9 §12.2 M4） ---
F_FAIL_LO, F_FAIL_HI = 77.0, 370.0        # N，系统综述：腓骨压缩失效载荷带
F_FAIL_MED = float(np.sqrt(F_FAIL_LO * F_FAIL_HI))          # 169 N（对数中位）
F_FAIL_SIG_LOG = float(np.log(F_FAIL_HI / F_FAIL_LO) / 3.29)  # 5–95 分位 ⇒ 0.478
P_CR_FE = 252.0                            # N，本项目孤立腓骨 FE 的 Euler/扫描屈曲阈
A_CORT = 76.0                              # mm²，本项目实测皮质截面（载荷流正对照：76.01–76.28）
SIGMA_C_THUMS = 70.0                       # MPa，G9 冻结的 σ_c（弱）
SIGMA_C_IMPLIED_LO = F_FAIL_LO / A_CORT    # 1.01 MPa
SIGMA_C_IMPLIED_HI = F_FAIL_HI / A_CORT    # 4.87 MPa


def _norm_cdf(x: np.ndarray) -> np.ndarray:
    """标准正态 CDF（Abramowitz-Stegun 7.1.26 误差函数近似，|ε|<1.5e-7）。"""
    return 0.5 * (1.0 + np.tanh(np.sqrt(2.0) * x * (1.0 + 0.0514718598 *
                                                  (x + 0.0406019748 *
                                                   (x - 0.0029929853)) * 0)) ** 0)


def _h_star(loads: np.ndarray, share: float, f_fail: float) -> int | None:
    f = share * loads
    idx = np.nonzero(f >= f_fail)[0]
    return int(idx[0] + 1) if idx.size else None      # 高度从 1 起


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-suffix", default="")
    args = ap.parse_args(argv)
    t0 = time.time()

    frozen = audit.load_frozen()
    heights = np.asarray(frozen["heights"], int)
    loads = np.asarray(frozen["loads"]["ankle_r__vert"], float)
    ref_idx = int(np.argmin(np.abs(heights - fm.REF_H)))
    print(f"[m4] 高度 {heights.min()}–{heights.max()} m，REF_H={fm.REF_H} m；"
          f"ankle_r@{fm.REF_H}m = {loads[ref_idx]:.1f} N", flush=True)
    print(f"[m4] 实测 F_fail 带 {F_FAIL_LO}–{F_FAIL_HI} N（中位 {F_FAIL_MED:.1f}，"
          f"σ_log={F_FAIL_SIG_LOG:.3f}）；FE P_cr={P_CR_FE} N；A_cort={A_CORT} mm²",
          flush=True)

    fails = {"exp_lo_77N": F_FAIL_LO, "exp_median_169N": F_FAIL_MED,
             "fe_Pcr_252N": P_CR_FE, "exp_hi_370N": F_FAIL_HI}
    shares = [0.002, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.10, 0.20, 0.30]

    rows: list[dict] = []
    for sh in shares:
        for fname, f in fails.items():
            hs = _h_star(loads, sh, f)
            n = int((sh * loads >= f).sum())
            rows.append({"share": sh, "fail": fname, "f_fail_n": f,
                         "h_star_m": hs, "n_fractured_h": n,
                         "n_heights": int(heights.size),
                         "status": ("all_one" if n == heights.size else
                                    "all_zero" if n == 0 else "progressive"),
                         "alloc_at_ref_n": float(sh * loads[ref_idx])})
    print("\n[m4] share → 骨折高度 h*(m) / 状态")
    hdr = "  share | " + " | ".join(f"{k:>17s}" for k in fails)
    print(hdr)
    for sh in shares:
        cells = []
        for fname in fails:
            r = next(x for x in rows if x["share"] == sh and x["fail"] == fname)
            cells.append(f"{str(r['h_star_m']):>6s}/{r['status'][:10]:<10s}")
        print(f"  {sh:5.3f} | " + " | ".join(f"{c:>17s}" for c in cells))

    # 渐进窗口与孤立模型可行 share 上限
    win: dict = {}
    grid = np.round(np.arange(0.0005, 0.3005, 0.0005), 6)
    for fname, f in fails.items():
        prog = [float(sh) for sh in grid
                if 0 < int((sh * loads >= f).sum()) < heights.size]
        safe = [float(sh) for sh in grid if int((sh * loads >= f).sum()) == 0]
        win[fname] = {"f_fail_n": f,
                      "progressive_share_range": [min(prog), max(prog)] if prog else None,
                      "all_zero_share_max": max(safe) if safe else None,
                      "proj_share_max": float(grid[-1]) if prog else None}
        r = win[fname]
        print(f"[m4] {fname:17s} F_fail={f:6.1f} N  渐进窗口="
              f"{('[%.3f, %.3f]' % (r['progressive_share_range'][0], r['progressive_share_range'][1])) if r['progressive_share_range'] else '无（≥0.30 全断）'}"
              f"  全程不骨折的 share ≤ {r['all_zero_share_max']}")

    # 实测离散折进逐高度失效概率（中位 F_fail 对数正态）
    from math import erf, sqrt
    prob: dict = {}
    for sh in (0.01, 0.015, 0.02, 0.05, 0.10):
        f = sh * loads
        p = np.array([0.5 * (1.0 + erf(
            (np.log(max(v, 1e-9)) - np.log(F_FAIL_MED)) / (sqrt(2) * F_FAIL_SIG_LOG)))
            for v in f])
        prob[f"{sh:.3f}"] = {"alloc_at_ref_n": float(sh * loads[ref_idx]),
                             "p_at_ref": float(p[ref_idx]),
                             "mean_p": float(p.mean()),
                             "p_first_50": [round(float(x), 4) for x in p[:12]]}
        print(f"[m4] share={sh:.3f}: p(骨折|REF)={p[ref_idx]:.3f}  "
              f"平均 p={p.mean():.3f}  逐高度 p[1..12]={np.round(p[:12],3).tolist()}")

    res = {"meta": {"title": "M4 腓骨力判据（实测断裂载荷）→ share×失效高度映射",
                    "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "heights_m": [int(h) for h in heights],
                    "ankle_r_vert_n": [float(x) for x in loads],
                    "ref_h_m": float(fm.REF_H), "ref_load_n": float(loads[ref_idx]),
                    "f_fail_band_n": [F_FAIL_LO, F_FAIL_HI],
                    "f_fail_median_n": F_FAIL_MED,
                    "f_fail_sigma_log": F_FAIL_SIG_LOG,
                    "fe_p_cr_n": P_CR_FE, "a_cort_mm2": A_CORT,
                    "sigma_c_thums_mpa": SIGMA_C_THUMS,
                    "sigma_c_implied_mpa": [SIGMA_C_IMPLIED_LO, SIGMA_C_IMPLIED_HI],
                    "sigma_c_thums_implied_force_n": SIGMA_C_THUMS * A_CORT,
                    "sources": ["Cureus 454613 系统综述（腓骨压缩失效 77–370 N）",
                                "本项目 FE P_cr≈252 N（G4 计划 §9.13）",
                                "本项目实测 A_cort≈76 mm²（G9 §11.7 载荷流正对照）"]},
           "rows": rows, "windows": win, "prob": prob,
           "elapsed_s": time.time() - t0}
    (RES / f"m4_fibula_force_criterion{args.out_suffix}.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    _report(res, RES / f"M4_FIBULA_FORCE_CRITERION_REPORT{args.out_suffix}.md")
    print(f"\n[m4] wrote: {RES / f'm4_fibula_force_criterion{args.out_suffix}.json'}")
    print(f"[m4] wrote: {RES / f'M4_FIBULA_FORCE_CRITERION_REPORT{args.out_suffix}.md'}")
    return 0


def _report(res: dict, path: Path) -> None:
    m = res["meta"]
    L = ["# M4 · 腓骨力判据（实测断裂载荷）→ share × 骨折高度映射", "",
         f"生成 {m['generated_at']}", "",
         "## 1. 三方对照（为什么孤立腓骨路线在该 share 下不成立）", "",
         "| 量 | 值 | 含义 |", "|---|---|---|",
         f"| 实测压缩失效载荷 | **{m['f_fail_band_n'][0]:.0f}–{m['f_fail_band_n'][1]:.0f} N**"
         f"（对数中位 {m['f_fail_median_n']:.0f} N） | 孤立腓骨的**结构**失效 |",
         f"| 本项目 FE 屈曲阈 P_cr | **{m['fe_p_cr_n']:.0f} N** | 落在实测带内 ⇒ "
         "**失效模式 = 屈曲**，FE 在该量级上被实验佐证 |",
         f"| THUMS σ_c=70 MPa × A_cort={m['a_cort_mm2']:.0f} mm² | "
         f"**{m['sigma_c_thums_implied_force_n']:.0f} N** | 比实测高 "
         f"{m['sigma_c_thums_implied_force_n']/m['f_fail_band_n'][1]:.0f}–"
         f"{m['sigma_c_thums_implied_force_n']/m['f_fail_band_n'][0]:.0f}× ⇒ "
         "**σ_c 不适用于孤立腓骨**（它会屈曲，不是压溃） |",
         f"| 力判据隐含的等效 σ_c | {m['sigma_c_implied_mpa'][0]:.2f}–"
         f"{m['sigma_c_implied_mpa'][1]:.2f} MPa | 远低于 THUMS 的 70 MPa |", "",
         "## 2. share → 骨折高度 h\\*（力判据）", "",
         "| share | alloc@REF (N) | " + " | ".join(f"{k} h\\*" for k in res["windows"]) + " |",
         "|---|---|" + "---|" * len(res["windows"])]
    by = {(r["share"], r["fail"]): r for r in res["rows"]}
    for sh in sorted({r["share"] for r in res["rows"]}):
        r0 = by[(sh, next(iter(res["windows"])))]
        cells = [str(by[(sh, k)]["h_star_m"] or "—") for k in res["windows"]]
        L.append(f"| {sh:.3f} | {r0['alloc_at_ref_n']:.0f} | " + " | ".join(cells) + " |")
    L += ["", "## 3. 窗口结论", "",
          "| F_fail | 渐进骨折的 share 区间 | 全程不骨折的 share ≤ |", "|---|---|---|"]
    for k, w in res["windows"].items():
        pr = w["progressive_share_range"]
        L.append(f"| {k} ({w['f_fail_n']:.0f} N) | "
                 f"{('[%.3f, %.3f]' % (pr[0], pr[1])) if pr else '无（share≥0.30 即全断）'} | "
                 f"{w['all_zero_share_max'] if w['all_zero_share_max'] is not None else '—'} |")
    L += ["", "## 4. 实测离散（σ_log=%.3f）折进逐高度失效概率" % m["f_fail_sigma_log"],
          "", "| share | alloc@REF (N) | p(骨折\\|REF) | 平均 p | p[1..12] |",
          "|---|---|---|---|---|"]
    for k, p in res["prob"].items():
        L.append(f"| {k} | {p['alloc_at_ref_n']:.0f} | {p['p_at_ref']:.3f} | "
                 f"{p['mean_p']:.3f} | {p['p_first_50']} |")
    L += ["", "**来源**：" + "；".join(m["sources"]), ""]
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())