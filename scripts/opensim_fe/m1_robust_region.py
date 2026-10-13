"""M1 S2 —— (share, A) 二维稳健域映射：DoD-2 的通过区覆盖"合理不确定性空间"吗？

**为什么先做这个，而不是建两体 FE**：
S1a 已证 **DoD-2 在解析层可达**（A_eff=76.01 @ share 0.2241 ⇒ OR=1.758 落窗），
两体 FE 因此从「出口」降级为「测量手段」。而 S1a 同时暴露了**真正的残余不确定性**：
**share 与 A 两个输入都有量级级别的分歧，且判定对它们敏感**：

  share：0.05（G8/G9 假设）｜0.2241（闭式并联）｜0.30（文献上限）—— **6× 跨度**
  A   ：33（median 反推，D3′ 已证 median 落噪声地板）｜76.01（载荷流实测）
        ｜101.9（体积法 A=V/L，含两端增厚）—— **3× 跨度**

⇒ 关键问题不是"某个点是否通过"，而是 **通过区（robust region）是否覆盖
整个合理不确定性空间**。若覆盖 ⇒ DoD-2 是**稳健**复活；若不覆盖 ⇒
必须先缩小不确定性（= 两体 FE 的正当用途：实测 share 与 A_eff）。

**正对照（三条，缺一则不出结论）**：
  P1 G3 冻结原行逐部位复现已存 `g3_or_estimates.json`（status, n_fracture, OR）
  P2 share=0.05 + A 使 util(h)≡0 的一行 ⇒ OR 应等于 G8 已发布的 1.152
  P3 share→0 极限 ⇒ 腓骨全零 ⇒ OR 应趋向窗口下界以下（all_zero 语义）

产物：results/opensim_fe/m1_robust_region.json + M1_ROBUST_REGION.md
"""
from __future__ import annotations

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

from scripts.opensim_fe import fibula_fe_audit as audit  # noqa: E402
from scripts.opensim_fe import g3_subject_variability as g3  # noqa: E402
from scripts.opensim_fe import g8_fibula_or_share as g8  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
OUT_JSON = RES / "m1_robust_region.json"
OUT_REPORT = RES / "M1_ROBUST_REGION.md"
SIGMA_C = 70.0
FIBULA = "腓骨"

#: share 轴：假设 / 闭式 / 文献，覆盖到"即使偏高"的稳健检验
SHARES = [0.05, 0.10, 0.15, 0.20, 0.2241, 0.25, 0.30, 0.35]
#: A 轴：三个已量化口径 + 两个插值（33↔76 之间的分歧地带必须扫到）
AREAS = [33.0, 45.0, 55.0, 65.0, 76.01, 88.0, 95.0, 101.9]


def or_for(sigma_h: np.ndarray) -> dict:
    inputs = g3.load_inputs()
    mixed, idx = g8.build_mixed_inputs(inputs, np.asarray(sigma_h, float))
    cfg = {sd: g3.run_config(sd, g3.COV_B_PRIMARY, g3.COV_E_PRIMARY,
                             g3.COV_A_PRIMARY, mixed, g3.N_SUBJ, True)
           for sd in g3.SEEDS}
    primary = cfg[g3.OFFICIAL_SEED]
    fr = g3._site_or_record(primary["or_map"][FIBULA])
    unc = g8.unchanged_check(inputs, idx, primary)
    return {"status": fr.get("status"), "OR": fr.get("OR"), "ci95": fr.get("ci95"),
            "n_fracture": fr.get("n_fracture"),
            "a3": primary["a3"], "unchanged_other8": unc.get("all_equal"),
            "max_abs_or_diff": unc.get("max_abs_or_diff")}


def main() -> int:
    t0 = time.time()
    frozen = audit.load_frozen()
    heights = np.asarray(frozen["heights"], int)
    loads = np.asarray(frozen["loads"]["ankle_r__vert"], float)
    ref_idx = int(np.argmin(np.abs(heights - fm_ref_height())))
    F_ref = float(loads[ref_idx])
    ratio = loads / F_ref

    target = float(g3.A3_TARGETS[FIBULA])
    win = [target * 0.7, target * 1.3]
    out: dict = {"meta": {
        "title": "M1 S2：(share, A) 二维稳健域映射",
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "sigma_c_mpa": SIGMA_C, "fibula_target_or": target,
        "fibula_window": win, "shares": SHARES, "areas_mm2": AREAS,
        "ankle_ref_n": F_ref,
        "model": "σ(h) = share·F_ref·(load(h)/load_ref)/A ⇒ util(h)=σ(h)/σ_c",
        "note": ("share 与 A 是 M1 仅有的两个自由输入；判定对二者均敏感 ⇒ "
                 "本脚本求通过域（robust region）是否覆盖合理不确定性空间")}}

    # ---- P1 正对照：G3 冻结原行 ----
    inputs = g3.load_inputs()
    base = {sd: g3.run_config(sd, g3.COV_B_PRIMARY, g3.COV_E_PRIMARY,
                              g3.COV_A_PRIMARY, inputs, g3.N_SUBJ, True)
            for sd in g3.SEEDS}
    b_map = base[g3.OFFICIAL_SEED]["or_map"]
    stored = json.loads(g8.G3_OR_JSON.read_text(encoding="utf-8"))["primary"]["per_site"]

    def _ok(cn):
        a, b = g3._site_or_record(b_map[cn]), stored.get(cn, {})
        d = (abs(float(a["OR"]) - float(b["OR"]))
             if a.get("OR") is not None and b.get("OR") is not None else 0.0)
        return (a.get("status") == b.get("status")
                and a.get("n_fracture") == b.get("n_fracture") and d < 1e-9)

    bad = [cn for cn in b_map if not _ok(cn)]
    out["positive_control_P1"] = {"n_sites": len(b_map), "n_mismatch": len(bad),
                                  "mismatched": bad, "pass": not bad}
    print(f"[P1] G3 冻结原行逐部位复现：{len(b_map)} 部位，{len(bad)} 不一致 ⇒ "
          f"{'✅ 通过' if not bad else '❌ ' + str(bad)}", flush=True)
    if bad:
        print("[P1] 不过 ⇒ 不继续")
        return 2

    print(f"\n窗口 [{win[0]:.3f}, {win[1]:.3f}]（靶 {target}）\n")
    hdr = "  share\\A   " + "".join(f"{a:>8.1f}" for a in AREAS)
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    grid: dict = {}
    ors: dict = {}
    for sh in SHARES:
        row = {}
        row_or = {}
        for area in AREAS:
            sigma_h = sh * F_ref * ratio / area
            r = or_for(sigma_h)
            util = sigma_h / SIGMA_C
            hit = bool(r["OR"] is not None and win[0] <= r["OR"] <= win[1])
            row[f"{area:.2f}"] = {
                "OR": r["OR"], "status": r["status"], "in_window": hit,
                "util_ref": float(util[ref_idx]), "util_max": float(util.max()),
                "n_fractured_h": int((util >= 1).sum()),
                "other8_unchanged": r["unchanged_other8"],
                "a3_n_within": r["a3"]["n_within"]}
            row_or[f"{area:.2f}"] = r
            ors[f"{sh:.4f}|{area:.2f}"] = r
        grid[f"{sh:.4f}"] = row
        cells = "".join(("   PASS " if row[f"{a:.2f}"]["in_window"] else "   .    ")
                        for a in AREAS)
        print(f"  {sh:8.4f}{cells}", flush=True)

    print("\n  PASS = OR 落窗；. = 未落窗")

    # ---- 通过域的边界：每个 A 的通过 share 区间；每个 share 的通过 A 区间 ----
    print("\n通过域边界：")
    bounds = {"by_area": {}, "by_share": {}}
    for area in AREAS:
        ks = [s for s in SHARES if grid[f"{s:.4f}"][f"{area:.2f}"]["in_window"]]
        bounds["by_area"][f"{area:.2f}"] = (
            [min(ks), max(ks)] if ks else None)
        rng = f"[{min(ks):.4f}, {max(ks):.4f}]" if ks else "无"
        print(f"  A={area:6.2f} mm² ⇒ 通过 share 区间 {rng}", flush=True)
    for sh in SHARES:
        aas = [a for a in AREAS if grid[f"{sh:.4f}"][f"{a:.2f}"]["in_window"]]
        bounds["by_share"][f"{sh:.4f}"] = (
            [min(aas), max(aas)] if aas else None)
        print(f"  share={sh:.4f} ⇒ 通过 A 区间 "
              f"{f'[{min(aas):.1f}, {max(aas):.1f}]' if aas else '无'} mm²", flush=True)

    # ---- 三个"锚点"是否在通过域内 ----
    anchors = {
        "g8_g9_assumed_5pct__A_eff_76": ("0.0500", "76.01"),
        "closed_form_22pct__A_eff_76": ("0.2241", "76.01"),
        "closed_form_22pct__A_cort_102": ("0.2241", "101.90"),
        "worst_A_33__closed_form_22pct": ("0.2241", "33.00"),
    }
    print("\n锚点判定：")
    out["anchors"] = {}
    for nm, (sk, ak) in anchors.items():
        cell = grid[sk][ak]
        out["anchors"][nm] = {"share": float(sk), "A": float(ak), **cell}
        print(f"  {nm:34s} OR={str(cell['OR'])[:6]:>6s}  "
              f"{'✅ PASS' if cell['in_window'] else '❌ .'}", flush=True)

    # ---- 稳健性结论 ----
    a_axes = [33.0, 76.01, 101.9]          # 三个已量化口径
    s_axes = [0.05, 0.2241, 0.30]           # 假设 / 闭式 / 文献
    sub = [[grid[f"{s:.4f}"][f"{a:.2f}"]["in_window"] for s in s_axes] for a in a_axes]
    n_pass = int(sum(sum(r) for r in sub))
    print(f"\n3×3 锚点网格（share∈{{{', '.join(f'{s}' for s in s_axes)}}} × "
          f"A∈{{{', '.join(f'{a}' for a in a_axes)}}}）：{n_pass}/9 PASS")
    out["robustness"] = {"share_axes": s_axes, "area_axes": a_axes,
                         "n_pass": n_pass, "n_total": 9,
                         "grid": sub}
    out["grid"] = grid

    print("\n" + "=" * 88)
    if n_pass >= 7:
        print("⇒ **DoD-2 稳健复活**：通过域覆盖绝大多数合理 (share, A) 组合")
        print("  ⇒ 两体 FE 的价值只剩「缩小 share/A 不确定性」= 锦上添花，不是出口")
    elif n_pass >= 3:
        print("⇒ **DoD-2 有条件复活**：通过域存在但未覆盖全部合理组合")
        print("  ⇒ 判定强依赖 A 口径（A=33 失败、A≥76 成功）")
        print("  ⇒ **两体 FE 有正当用途**：实测耦合 A_eff 即可把判定钉死")
    else:
        print("⇒ **DoD-2 仍不可达**：通过域过窄")

    out["elapsed_s"] = time.time() - t0
    OUT_JSON.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    _report(out, bounds)
    print(f"\nwrote: {OUT_JSON}\nwrote: {OUT_REPORT}")
    return 0


def fm_ref_height() -> float:
    import fracture_matrix as fm
    return fm.REF_H


def _report(out: dict, bounds: dict) -> None:
    m, rb = out["meta"], out["robustness"]
    L: list[str] = []
    A = L.append
    A("# M1 全量 · S2 —— (share, A) 二维稳健域映射")
    A("")
    A(f"> 生成 {m['generated_at']}；{m['note']}")
    A(f"> σ_c={m['sigma_c_mpa']} MPa；腓骨靶 OR={m['fibula_target_or']}，±30% ⇒ "
      f"窗口 {[round(x,3) for x in m['fibula_window']]}。")
    A(f"> 正对照 P1：G3 冻结原行 {out['positive_control_P1']['n_sites']} 部位，"
      f"{out['positive_control_P1']['n_mismatch']} 不一致 ⇒ "
      f"{'✅ 通过' if out['positive_control_P1']['pass'] else '❌'}")
    A("")
    A("## 0. 结论")
    A("")
    A(f"- **3×3 锚点网格 {rb['n_pass']}/9 PASS**（share ∈ {rb['share_axes']} × "
      f"A ∈ {rb['area_axes']}）。")
    for nm, a in out["anchors"].items():
        A(f"  - `{nm}`：share={a['share']:.4f} A={a['A']:.1f} ⇒ OR="
          f"{('%.4f' % a['OR']) if a['OR'] else '—'} "
          f"{'✅ PASS' if a['in_window'] else '❌ 未落窗'}")
    A("")
    A("## 1. 通过域矩阵（行=share，列=A mm²；P=OR 落窗）")
    A("")
    A("| share \\ A | " + " | ".join(f"{a:.1f}" for a in m["areas_mm2"]) + " |")
    A("|---" * (len(m["areas_mm2"]) + 1) + "|")
    for sk, row in [(f"{s:.4f}", out["grid"][f"{s:.4f}"]) for s in m["shares"]] \
            if "grid" in out else []:
        A(f"| **{sk}** | " + " | ".join(
            ("**P**" if row[f"{a:.2f}"]["in_window"] else "·")
            for a in m["areas_mm2"]) + " |")
    A("")
    A("## 2. 通过域边界")
    A("")
    A("| A (mm²) | 通过 share 区间 |")
    A("|---:|---|")
    for ak, v in bounds["by_area"].items():
        A(f"| {float(ak):.2f} | "
          f"{f'[{v[0]:.4f}, {v[1]:.4f}]' if v else '**无**'} |")
    A("")
    A("| share | 通过 A 区间 (mm²) |")
    A("|---:|---|")
    for sk, v in bounds["by_share"].items():
        A(f"| {float(sk):.4f} | "
          f"{f'[{v[0]:.1f}, {v[1]:.1f}]' if v else '**无**'} |")
    A("")
    OUT_REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())