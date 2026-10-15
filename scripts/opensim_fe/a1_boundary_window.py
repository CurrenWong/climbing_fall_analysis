# -*- coding: utf-8 -*-
"""A1 聚类边界 —— 混 BC 复核 + 边界窗口扫描（零 FE，只读）

背景
----
`S4_ALIGNMENT_BC2_REPORT.md` §1 报 A1 = FAIL，边界 [5.0, 12.0]、Δ = [2.0, 4.0]。
`A1_BOUNDARY_CONVENTION_REPORT.md`（10-13）已证：#1 的 Δ=2.0 是**报边界口径产物**
（换约定即 PASS）、#2 的 Δ∈[3.0,4.0] 是**真差异**。
但它**没有回答**两个自然跟进问题，本脚本补上：

  Q1  A2 那次「混 BC」复核（股骨换 P5 接触 BC）跑的是同一个 `route_yi`，
      会一并算出 A1 —— 但报告只写了 A2。**A1 在混 BC 下到底动没动？**
  Q2  #2 的「真差异」是**结构不可达**，还是只是某个部位的值不对？**有没有窗口？**

结论（本脚本产出，详见 `results/opensim_fe/A1_BOUNDARY_WINDOW_REPORT.md`）
--------------------------------------------------------------
  A1  混 BC 对 A1 **逐位无影响**（[5.0,12.0]/[2.0,4.0] 不变）⇒ 修股骨不移动 A1。
  A2  #2 **跟得动**胸椎首折高度 h*，**窗口 = h* ∈ [16, 18] m**（Δ₂=1.0）；
      h*=20⇒3.0、h*=25⇒8.0、**h*=33⇒16.0（大幅恶化）**。
  A3  **A1 与 A2 修法互相拉扯**：A2 只要胸椎变「阶跃」即命中（任何高度），
      A1 要求**恰落 [16,18] m** ⇒ **禁照搬股骨修法**（见报告 §3）。
  A4  #1 在**所有变体下恒为 5.0**（Δ=2.0，判据严格 < 2.0 ⇒ 永不过）—— 纯口径问题。

性质
----
纯只读、零 FEBio、**不改判据、不改冻结输入、不回溯改判**。
唯一写盘：`results/opensim_fe/a1_boundary_window.json`（新产物）。
⚠️ Q1 依赖 A2 修订时落盘的 patch 产物 `temp/b1_official/…`（temp 为易失目录）；
   若缺失，请先跑 `temp/b1_official_rerun.py` 生成，再重跑本脚本（A 段会跳过并提示）。

复现
----
    cd D:\\Project\\climbing_fall_analysis
    .venv\\Scripts\\python.exe scripts\\opensim_fe\\a1_boundary_window.py
"""
from __future__ import annotations

import datetime as _dt
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts" / "opensim_fe"
RES = ROOT / "results" / "opensim_fe"
sys.path.insert(0, str(SCRIPTS))

import s4_alignment_bc2 as s4   # noqa: E402
import paper_alignment as pa    # noqa: E402

# 判据与目标：**直接引用官方模块**，避免本地副本与官方漂移
TARGETS = list(s4.A1_TARGETS)
TOL = float(s4.A1_TOL)

# A2 修订（混 BC）时落盘的 patch 产物
PATCH_JSON = ROOT / "temp" / "b1_official" / "fracture_matrix_bc2_p5femur.json"
PATCH_DIR = ROOT / "temp" / "b1_official" / "bc2_p5femur"

# 官方默认（用于跑完 A 段后还原）
_ORIG_JSON = s4.BC2_JSON
_ORIG_DIR = s4.BC2_DIR


# ------------------------------------------------------------------ helpers
def _row_map(yi: dict) -> dict:
    """逐部位 0/1 行（1 = 该高度判骨折）。"""
    return {cn: np.asarray(yi["matrix_by_part"][cn], dtype=int)
            for cn in pa.REGION_NAMES}


def _step_at(h: int) -> np.ndarray:
    """合成「在 h 处阶跃」的行（h 之前 0、h 起 1），仅用于 what-if。"""
    return np.array([1 if (i + 1) >= h else 0 for i in range(50)], dtype=int)


def _a1_of(rows: dict, override: dict) -> tuple[dict, list, dict]:
    """把 override 里的行替换掉，重跑官方 聚类 → 相段 → A1 判定。"""
    r2 = {cn: (override[cn] if cn in override else rows[cn])
          for cn in pa.REGION_NAMES}
    M = np.array([[int(r2[cn][i]) for cn in pa.REGION_NAMES] for i in range(50)],
                 dtype=int)
    clus = pa.jaccard_clusters(M, 5)
    phases = s4.three_phases_seq(clus["labels"], M)
    return s4.a1_verdict(phases), phases, clus


def _first(row: np.ndarray):
    for i, v in enumerate(row):
        if v:
            return i + 1
    return None


def _kind(row: np.ndarray) -> str:
    s = int(row.sum())
    return "常值1" if s == 50 else ("常值0" if s == 0 else "阶跃")


def _table(rows: dict, title: str) -> None:
    print("  %s" % title)
    for cn in pa.REGION_NAMES:
        print("    %-8s first=%-6s n_frac=%-3d %s"
              % (cn, _first(rows[cn]) or ">50", int(rows[cn].sum()),
                 _kind(rows[cn])))


# ------------------------------------------------------------------ run A
def run_a_mixed_bc(out: dict) -> dict | None:
    """Q1：混 BC（股骨 → P5 接触）下 A1 是否改变。

    返回混 BC 基线的逐部位行（供 B/C/D 用），失败返回 None。
    """
    print("=" * 74)
    print("A. 混 BC 复核 —— Q1：修股骨会不会移动 A1 边界？")
    print("=" * 74)

    a1_0 = s4.a1_verdict(s4.route_yi("median")["phases"])
    print("  [正对照·未 patch]   A1 边界=%-14s Δ=%-12s -> %s"
          % (a1_0["ours_boundaries_m"], a1_0["abs_delta"], a1_0["verdict"]))

    if not PATCH_JSON.exists():
        print("  ⚠️ 缺 patch 产物：%s" % PATCH_JSON.relative_to(ROOT))
        print("     请先跑 temp/b1_official_rerun.py 生成，再重跑本脚本。")
        print("     （A 段跳过；B/C 段回落到 BC2 基线并明确标注。）")
        out["mixed_bc"] = None
        print()
        return None

    try:
        s4.BC2_JSON = PATCH_JSON
        s4.BC2_DIR = PATCH_DIR
        yi1 = s4.route_yi("median")
    finally:
        s4.BC2_JSON = _ORIG_JSON      # 还原官方默认
        s4.BC2_DIR = _ORIG_DIR

    a1_1 = s4.a1_verdict(yi1["phases"])
    print("  [混 BC·股骨→P5]     A1 边界=%-14s Δ=%-12s -> %s"
          % (a1_1["ours_boundaries_m"], a1_1["abs_delta"], a1_1["verdict"]))

    same = a1_0["ours_boundaries_m"] == a1_1["ours_boundaries_m"]
    print("  ⇒ A1 %s —— %s"
          % ("**逐位不变**" if same else "**有变化**",
             "修股骨不移动 A1 边界" if same else "需重新评估"))

    # 交叉校验：同 run 下 A2 确实翻了（证明确实跑到 patch，不是空跑）
    v2 = s4.a2_verdict(yi1["logistic"])
    print("  [同 run 交叉校验]   A2 n_hit=%d/%d；股骨=%s；股骨首折=%s m"
          % (v2["n_hit"], v2["n_target"], yi1["logistic"]["股骨"]["status"],
             yi1["first"]["股骨"]))
    print()

    out["mixed_bc"] = {
        "baseline": {"boundaries": a1_0["ours_boundaries_m"],
                     "delta": a1_0["abs_delta"], "verdict": a1_0["verdict"]},
        "p5femur": {"boundaries": a1_1["ours_boundaries_m"],
                    "delta": a1_1["abs_delta"], "verdict": a1_1["verdict"]},
        "a1_unchanged": bool(same),
        "a2_cross_check": {"n_hit": v2["n_hit"], "n_target": v2["n_target"],
                           "femur_status": yi1["logistic"]["股骨"]["status"],
                           "femur_first_m": yi1["first"]["股骨"]},
    }
    return _row_map(yi1)


# ------------------------------------------------------------------ run B
def run_b_whatif(rows: dict, out: dict, tag: str) -> None:
    """Q2 前置：修哪些部位能移动 #2？各自修到「股骨的值」(33 m) 与「窗口内」(16 m)。"""
    print("=" * 74)
    print("B. 单变量/组合 what-if（基线 = %s）—— 哪些部位的修复能移动 #2？" % tag)
    print("=" * 74)
    _table(rows, "基线逐部位：")
    print()

    variants = [
        ("基线（无追加修复）", {}),
        ("只修 胸椎 @33（照搬股骨值）", {"胸椎": _step_at(33)}),
        ("只修 腓骨 @33", {"腓骨": _step_at(33)}),
        ("修 胸椎+腓骨 @33", {"胸椎": _step_at(33), "腓骨": _step_at(33)}),
        ("修 胸椎+腓骨 @16（窗口内）",
         {"胸椎": _step_at(16), "腓骨": _step_at(16)}),
    ]
    rec = {}
    for label, over in variants:
        v, phases, _ = _a1_of(rows, over)
        print("  %-30s 边界=%-14s Δ=%-12s -> %s"
              % (label, v["ours_boundaries_m"], v["abs_delta"], v["verdict"]))
        print("      相段: %s" % [p["range_m"] for p in phases])
        rec[label] = {"boundaries": v["ours_boundaries_m"],
                      "delta": v["abs_delta"], "verdict": v["verdict"],
                      "phases": [p["range_m"] for p in phases]}
    out.setdefault("whatif", {})[tag] = rec
    print()


# ------------------------------------------------------------------ run C
def run_c_window(rows: dict, out: dict, tag: str) -> None:
    """Q2：胸椎首折 h* 与 #2 边界的关系，找通过窗口（基线 = tag）。"""
    print("=" * 74)
    print("C. 边界窗口扫描（基线 = %s）—— 胸椎首折 h* 落在哪，#2 才进 (14,18)？" % tag)
    print("=" * 74)
    print("  目标 #1=%.0f m / #2=%.0f m，容差 |Δ|<%.0f ⇒ #1∈(5,9)、#2∈(14,18)"
          % (TARGETS[0], TARGETS[1], TOL))
    print("  （#1 在所有变体下恒 5.0 ⇒ Δ₁=2.0，末格约定下永不过；见 D 段）\n")

    for only, label in ((True, "只改胸椎（腓骨维持原状）"),
                        (False, "胸椎+腓骨同时改")):
        print("  【%s】" % label)
        print("  %-14s %-12s %-10s %s" % ("胸椎首折 h*", "#2 边界", "Δ₂", "判定(#2)"))
        print("  " + "-" * 58)
        rec = {}
        for h in list(range(10, 23)) + [25, 30, 33]:
            over = {"胸椎": _step_at(h)}
            if not only:
                over["腓骨"] = _step_at(h)
            v, _, _ = _a1_of(rows, over)
            b2 = v["ours_boundaries_m"][1]
            d2 = v["abs_delta"][1]
            ok = d2 < TOL
            print("  %-14s %-12s %-10s %s%s"
                  % ("%d m" % h, b2, d2, "PASS" if ok else "FAIL",
                     "   ← ✅ #2 可过" if ok else ""))
            rec[str(h)] = {"b2": b2, "d2": d2, "pass2": bool(ok)}

        window = [int(k) for k, v in sorted(rec.items(), key=lambda kv: int(kv[0]))
                  if v["pass2"]]
        print("  ⇒ #2 通过窗口 = h* ∈ %s" % (window if window else "**空集**"))
        print("     注：h*=33（照搬股骨修法）⇒ Δ₂ = %.1f" % rec["33"]["d2"])
        print()
        out.setdefault("window", {})[tag + ("|只胸椎" if only else "|胸+腓")] = {
            "per_h": rec, "pass_window_h": window,
            "note": "h*=33 ⇒ Δ₂=%.1f" % rec["33"]["d2"]}


# ------------------------------------------------------------------ run D
def run_d_boundary1(rows: dict, out: dict) -> None:
    """#1 恒 5.0 ⇒ 纯口径问题（与 A1_BOUNDARY_CONVENTION_REPORT 呼应）。"""
    print("=" * 74)
    print("D. #1 恒定性 —— 是不是纯口径问题？")
    print("=" * 74)
    seen = set()
    for over in ({}, {"胸椎": _step_at(16)}, {"胸椎": _step_at(33)},
                 {"胸椎": _step_at(16), "腓骨": _step_at(16)},
                 {"腓骨": _step_at(16)}):
        v, _, _ = _a1_of(rows, over)
        seen.add((v["ours_boundaries_m"][0], v["abs_delta"][0]))
    vals = sorted(seen)
    print("  5 个变体扫描下 #1 取值集合 = %s  (报出值, Δ₁)" % vals)
    print("  ⇒ #1 恒为 5.0（Δ₁ = 2.0；判据 %.1f 为**严格小于**）" % TOL)
    print("     真实相变区间 [5, 6] m（由胫骨首折 5 m 锚定，1 m 网格 ±0.5 m）")
    print("     末格 5.0 ⇒ Δ=2.0 FAIL ｜ 中点 5.5 ⇒ 1.5 PASS ｜ 次相首格 6.0 ⇒ 1.0 PASS")
    print("     ⇒ **#1 的白拿条件 = 在预注册里声明报边界口径**（禁赛后改规则）")
    out["boundary1"] = {"values_seen": [list(v) for v in vals],
                        "always": 5.0,
                        "conventions": {"末格": 2.0, "中点": 1.5, "次相首格": 1.0}}
    print()


# ------------------------------------------------------------------ main
def main() -> int:
    print(__doc__)
    print()

    out: dict = {"meta": {
        "generated_at": _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "script": "scripts/opensim_fe/a1_boundary_window.py",
        "targets_m": TARGETS, "tol_m": TOL,
        "caliber": "BC2 正对照（median, align_pca=True）+ 官方 "
                   "three_phases_seq / a1_verdict",
        "note": "零 FE 只读复核；B/C/D 的 override 行是**范围界定探针**，"
                "非结果、不入判",
    }}

    rows_bc2 = _row_map(s4.route_yi("median"))
    rows_mixed = run_a_mixed_bc(out)

    # 基线选择：优先混 BC（当前最优状态）；缺 patch 产物则回落 BC2 并标注
    if rows_mixed is not None:
        rows, tag = rows_mixed, "混 BC（股骨已阶跃@33）"
    else:
        rows, tag = rows_bc2, "BC2 正对照（股骨仍 常值1）"
        print("=" * 74)
        print("A. 混 BC 复核 —— 未测（缺 patch 产物），基线回落 BC2")
        print("=" * 74)
        print()

    print("=" * 74)
    print("基线（%s）逐部位：" % tag)
    print("=" * 74)
    _table(rows, "")
    print()

    run_b_whatif(rows, out, tag)
    run_c_window(rows, out, tag)

    # 对照：BC2 基线上同样扫一遍，凸显「先修谁」的顺序效应
    if rows_mixed is not None:
        r2 = _row_map(s4.route_yi("median"))
        run_b_whatif(r2, out, "BC2 正对照（股骨仍 常值1）")

    run_d_boundary1(rows, out)

    print("=" * 74)
    print("汇总")
    print("=" * 74)
    mb = out.get("mixed_bc")
    if mb:
        print("  ① 混 BC 是否移动 A1：%s（股骨首折 %.0f m；A2 同 run n_hit=%d/%d）"
              % ("**否，逐位不变**" if mb["a1_unchanged"] else "是",
                 mb["a2_cross_check"]["femur_first_m"],
                 mb["a2_cross_check"]["n_hit"], mb["a2_cross_check"]["n_target"]))
    else:
        print("  ① 混 BC：未测（缺 temp/b1_official patch 产物）")
    for k, v in out.get("window", {}).items():
        print("  ② #2 通过窗口 [%s]：h* ∈ %s"
              % (k, v["pass_window_h"] if v["pass_window_h"] else "空集"))
    print("  ③ #1 恒 5.0 ⇒ 需在预注册声明报边界口径（末格→中点）")
    print("  基线 = %s" % tag)

    dst = RES / "a1_boundary_window.json"
    dst.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print("\n  产物：%s" % dst.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
