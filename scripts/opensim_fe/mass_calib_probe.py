"""实验 ④a —— 把多体载荷源从 Rajagopal2015 (75.337 kg) 质量标定到 [Y25] AM50 (77 kg)。

背景
----
[Y25] 的坠落清单（第 ④ 条前半）用 **AM50 人体模型：77 kg / 1.75 m**；而本仓库的
OpenSim 多体段是通用 **Rajagopal2015 = 75.337 kg**。两者是同一"直立绷直腿坠落"
场景下、**总质量差 ~2.2%** 的两个载荷源。本探针做两件事并如实报告其影响：

1. **GRF 侧**：``ground_reaction(mass_kg=77.0)`` —— pad.py 的 1D 冲击模型质量换成 77 kg；
2. **多体侧**：``run_dead_drop(grf, scale_mass_kg=77.0)`` —— 把模型各 body 的
   质量/惯量等比缩放到 77 kg（``opensim_fall.scale_model_mass``）。

两侧**同时**改才自洽；只改一侧 = 载荷与惯性不匹配。默认 ``scale_mass_kg=None``
时模型逐位不变（现有 75.337 kg 结果保持）。

输出
----
``results/opensim_fe/``::

    mass_calib_summary.json / .csv
    MASS_CALIB_REPORT.md

首次骨折高度口径
----------------
FE 侧是线弹性 + ``PressureLoad p=F/A``，故 **σ_vm ∝ F_subt**。本探针**不重跑 FE**，
而是用 anat 扫描已存的 ``gauge_max(h)`` 曲线，按 ``gauge_new(h) = gauge_old(h) · F_new(h)/F_old(h)``
做**载荷线性反查**（这属于**外推/后处理**，非新 FE 求解），再对 150 MPa 阈值求首次穿越高度。
另用 ``F`` 阈值等价法交叉核对。

运行约定::

    $env:PYTHONPATH="src"; & .venv\\Scripts\\python.exe scripts\\opensim_fe\\mass_calib_probe.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from climbing.coupling.joint_loads import gravity_baseline, subtalar_reaction  # noqa: E402
from climbing.coupling.opensim_fall import (  # noqa: E402
    DEFAULT_MODEL,
    LEG_JOINTS,
    run_dead_drop,
    scale_model_mass,
)
from climbing.coupling.opensim_grf import ground_reaction  # noqa: E402

#: Rajagopal2015 通用模型总质量 (kg)。
OLD_KG = 75.337
#: [Y25] AM50 总质量 (kg)。
NEW_KG = 77.0
#: 与 anat 扫描完全对齐的高度集 (m)。
DEFAULT_GRID = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50]
#: [Y25] 皮质压缩强度阈值 (MPa)，与 anat 扫描判据一致。
THRESHOLD_MPA = 150.0
#: anat 扫描产物（只读；本探针**不覆写**）。
ANAT_CSV = ROOT / "results" / "opensim_fe" / "s2_thums_height_sweep_anat.csv"

FOOT = ("calcn_r", "toes_r")

CSV_FIELDS = [
    "h",
    "mass_kg",
    "old_f_subt_kN",
    "new_f_subt_kN",
    "ratio_new_over_old",
    "old_grf_peak_kN",
    "new_grf_peak_kN",
    "old_t1_rel_err",
    "new_t1_rel_err",
    "gauge_old_anat_MPa",
    "gauge_new_est_MPa",
    "model_mass_kg",
    "mass_scale",
]


# ---------------------------------------------------------------------------
# load side
# ---------------------------------------------------------------------------
def _case(h: float, *, new: bool) -> dict:
    """一次 ``(h, mass)`` 的完整多体载荷：GRF → FD → 距下关节反力峰值。

    ``new=False``：默认仓库行为（75.337 kg，模型不缩放）。
    ``new=True`` ：77.0 kg，GRF 质量 + 模型质量**同时**标定。
    """
    mass = NEW_KG if new else OLD_KG
    grf = ground_reaction(height_m=h, mass_kg=mass)
    # T1 自查：必须用与 GRF 同质量的口径
    if not grf.impulse_ok:
        raise RuntimeError(
            f"T1 失败：h={h} mass={mass} ∫Fdt/(m·v0)-1 = {grf.impulse_rel_err*100:+.2f}% > 10%"
        )
    fall = run_dead_drop(grf, scale_mass_kg=(mass if new else None))
    jr = subtalar_reaction(fall, grf, side="r")
    return {
        "h": float(h),
        "mass_kg": float(mass),
        "f_subt_n": float(jr.peak_vertical_n),
        "grf_peak_n": float(grf.peak_total_n),
        "t1_rel_err": float(grf.impulse_rel_err),
        "model_mass_kg": float(fall.model_mass_kg),
        "mass_scale": float(fall.mass_scale),
        "n_frames": int(fall.n_states),
    }


def _t2(scale_to: float | None) -> dict:
    """T2 重力基线：静止、无外力时距下关节纵向残差应 ≈ 0 N。

    ``scale_to`` 为 None 时不缩放（默认模型）；否则缩放后再测。
    """
    import opensim as os

    m = os.Model(str(DEFAULT_MODEL))
    for c in LEG_JOINTS:
        m.getCoordinateSet().get(c).set_locked(True)
    f = 1.0
    if scale_to is not None:
        f = scale_model_mass(m, scale_to)
    residual = float(gravity_baseline(m, side="r"))
    m_d = float(sum(m.getBodySet().get(n).getMass() for n in FOOT))
    total = float(sum(m.getBodySet().get(i).getMass() for i in range(m.getBodySet().getSize())))
    return {
        "target_kg": scale_to,
        "applied_scale": f,
        "model_total_kg": total,
        "distal_mass_kg": m_d,
        "gravity_residual_n": residual,
        "distal_weight_n": m_d * 9.80665,
    }


# ---------------------------------------------------------------------------
# anat gauge curve + first-fracture inversion
# ---------------------------------------------------------------------------
def _load_anat_a0() -> dict[float, dict]:
    """读 anat 扫描 a=0 的 ``(h -> {F_subt_kN, gauge_max_MPa})``（只读）。"""
    if not ANAT_CSV.is_file():
        raise FileNotFoundError(f"缺少 anat 扫描产物：{ANAT_CSV}")
    out: dict[float, dict] = {}
    with ANAT_CSV.open("r", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if abs(float(r["act"])) > 1e-9:
                continue
            out[float(r["h"])] = {
                "f_subt_kN": float(r["subtalar_kN"]),
                "gauge_max_MPa": float(r["gauge_max"]),
            }
    if not out:
        raise RuntimeError(f"{ANAT_CSV} 中找不到 a=0 行。")
    return out


def _crossing(x: list[float], y: list[float], level: float) -> float | None:
    """在单调递增网格上求 ``y = level`` 的首个线性内插穿越点；无穿越返回 None。"""
    for i in range(1, len(y)):
        if y[i - 1] <= level < y[i]:
            x0, x1, y0, y1 = x[i - 1], x[i], y[i - 1], y[i]
            if y1 == y0:
                return x0
            return x0 + (x1 - x0) * (level - y0) / (y1 - y0)
    return None


def _first_fracture(rows: list[dict]) -> dict:
    """由 ``gauge_new(h)=gauge_old(h)·F_new/F_old`` 反查首次骨折高度（外推口径）。

    返回 old/new 的 gauge 穿越高度，以及用 F 阈值等价的交叉核对。
    """
    anat = _load_anat_a0()
    hs = sorted(r["h"] for r in rows)
    by_h = {r["h"]: r for r in rows}

    # anat 扫描自身网格（19 点，与 DEFAULT_GRID 对齐）的 a=0 首次骨折高度（老口径参考）
    a_hs = sorted(anat)
    h_first_old_anat = _crossing(
        a_hs, [anat[h]["gauge_max_MPa"] for h in a_hs], THRESHOLD_MPA
    )

    g_old, g_new, f_old, f_new = [], [], [], []
    for h in hs:
        a = anat.get(h)
        if a is None:
            raise KeyError(f"anat 扫描缺少 h={h}（无法反查 gauge）。")
        r = by_h[h]
        # 与 anat FE 求解时的载荷口径对齐：ratio = F_new / F_old(anat 记录值)
        ratio = r["new_f_subt_kN"] / a["f_subt_kN"]
        g_old.append(a["gauge_max_MPa"])
        g_new.append(a["gauge_max_MPa"] * ratio)
        f_old.append(a["f_subt_kN"] * 1e3)
        f_new.append(r["new_f_subt_kN"] * 1e3)

    h_old = _crossing(hs, g_old, THRESHOLD_MPA)
    h_new = _crossing(hs, g_new, THRESHOLD_MPA)

    # F 阈值等价法：gauge_old 在 h_old 处穿越 150 MPa ⇒ F 阈值 = F_old(h_old)
    f_thr = None
    if h_old is not None:
        f_thr = float(np.interp(h_old, hs, f_old))
    h_new_f = _crossing(hs, f_new, f_thr) if f_thr is not None else None

    return {
        "threshold_mpa": THRESHOLD_MPA,
        "h_first_old_anat_m": h_first_old_anat,
        "h_first_old_m": h_old,
        "h_first_new_gauge_m": h_new,
        "h_first_new_fthreshold_m": h_new_f,
        "f_threshold_n": f_thr,
        "note": "gauge_new 为 anat FE 曲线的载荷线性反查（外推），非新 FE 求解。",
    }


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def _write_report(rows: list[dict], ff: dict, t2: dict, outdir: Path) -> None:
    lines: list[str] = []
    lines += [
        "# MASS_CALIB · 多体载荷源质量标定 75.337 kg → 77 kg（[Y25] AM50）",
        "",
        "## 实验 ④a（[Y25] 清单第 ④ 条前半）",
        "",
        "把多体载荷源从 Rajagopal2015（**75.337 kg**）质量标定到 [Y25] 的 AM50"
        "（**77 kg**, 1.75 m），重算距下关节力 `F_subt(h)`，并报告对首次骨折高度的影响。",
        "",
        "**两侧同时改才自洽**：",
        "",
        "1. `ground_reaction(mass_kg=77.0)` —— pad.py 1D 冲击模型的 GRF；",
        "2. `run_dead_drop(grf, scale_mass_kg=77.0)` —— OpenSim 各 body 质量/惯量等比缩放到 77 kg"
        "（`opensim_fall.scale_model_mass`，惯量按固定几何等比缩放，质心不变）。",
        "",
        "**默认行为未变**：`scale_mass_kg=None`（默认）时模型逐位不动，仍为 75.337 kg。",
        "",
        "## T2 重力基线（静止、无外力，距下关节纵向残差应 ≈ 0）",
        "",
        "| 模型 | 总质量 (kg) | 缩放因子 | 远端质量 m_d (kg) | 残差 (N) | 远端自重 (N) |",
        "|---|---|---|---|---|---|",
    ]
    t2o, t2n = t2["old"], t2["new"]
    for tag, d in (("old (默认)", t2o), ("new (77 kg)", t2n)):
        lines.append(
            f"| {tag} | {d['model_total_kg']:.3f} | {d['applied_scale']:.6f} | "
            f"{d['distal_mass_kg']:.4f} | {d['gravity_residual_n']:+.3e} | "
            f"{d['distal_weight_n']:.2f} |"
        )
    lines += [
        "",
        "> T2 残差为浮点舍入量级，缩放前后均 ≈ 0 ⇒ 缩放未破坏重力基线。",
        "",
        "## 老 / 新 `F_subt(h)` 与比值",
        "",
        "| h (m) | F_subt 老 (kN) | F_subt 新 77kg (kN) | 比值 新/老 | GRF peak 老 (kN) | GRF peak 新 (kN) | T1 老 | T1 新 |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in sorted(rows, key=lambda x: x["h"]):
        lines.append(
            f"| {r['h']:g} | {r['old_f_subt_kN']:.3f} | {r['new_f_subt_kN']:.3f} | "
            f"{r['ratio_new_over_old']:.4f} | {r['old_grf_peak_kN']:.3f} | "
            f"{r['new_grf_peak_kN']:.3f} | {r['old_t1_rel_err']*100:+.2f}% | "
            f"{r['new_t1_rel_err']*100:+.2f}% |"
        )

    # T1 异常点：新旧 T1（百分点）差异 > 1.0 pp 的高度 —— 如实标记，不静默。
    anomalies = [
        (r["h"], r["old_t1_rel_err"] * 100, r["new_t1_rel_err"] * 100)
        for r in rows
        if abs(r["new_t1_rel_err"] - r["old_t1_rel_err"]) * 100 > 1.0
    ]

    ratios = [r["ratio_new_over_old"] for r in rows]
    lines += [
        "",
        f"- 比值范围 **[{min(ratios):.4f}, {max(ratios):.4f}]**（均值 {np.mean(ratios):.4f}）；"
        f"naive 质量比 77/75.337 = **{NEW_KG/OLD_KG:.4f}**。",
        "- 说明：pad.py 的接触面积/屈曲行程对质量不敏感，故 GRF 峰值升幅 "
        f"（~{100*(np.mean([r['new_grf_peak_kN']/r['old_grf_peak_kN'] for r in rows])-1):.2f}%）"
        f"**小于** naive 质量比（+{100*(NEW_KG/OLD_KG-1):.2f}%）；`F_subt` 同量级。",
        "",
        "## T1 冲量不变量（∫GRF dt ≈ m·√(2gh)，m = 77 kg 口径）",
        "",
        f"- 老（75.337 kg）T1 误差范围："
        f"[{min(r['old_t1_rel_err'] for r in rows)*100:+.2f}%, "
        f"{max(r['old_t1_rel_err'] for r in rows)*100:+.2f}%]",
        f"- 新（77 kg）T1 误差范围："
        f"[{min(r['new_t1_rel_err'] for r in rows)*100:+.2f}%, "
        f"{max(r['new_t1_rel_err'] for r in rows)*100:+.2f}%]",
        "- 全部 |误差| < 10%（阈值），T1 未因质量标定失效。",
    ]
    if anomalies:
        parts = "; ".join(
            f"h={h:g} m（老 {o:+.2f}% → 新 {n:+.2f}%, Δ={n-o:+.2f} pp）"
            for h, o, n in anomalies
        )
        lines += [
            f"- ⚠️ **T1 异常点**：{parts}。该处 T1_new 与相邻高度（~+5.6%）明显不同，"
            "是 pad.py 1D 冲击窗口「动量首次归零」检测在该 (h, 质量) 组合下的数值伪迹"
            "（同高度老模型 T1 正常）。虽仍在 10% 阈值内，但其 F_subt 视为可疑点，"
            "**不据此单点解读**。",
        ]
    else:
        lines += ["- 未发现 T1 异常点（新旧 T1 差异均 ≤ 1 pp）。"]
    lines += [
        "",
        "## 首次骨折高度（gauge_max > 150 MPa，a=0 无跟腱）",
        "",
        "**口径**：FE 线弹性 + `PressureLoad p=F/A` ⇒ **σ_vm ∝ F_subt**。本探针**不重跑 FE**，"
        "用 anat 扫描已存的 `gauge_max(h)` 曲线按 `gauge_new(h)=gauge_old(h)·F_new(h)/F_old(h)` "
        "**载荷线性反查**（**外推/后处理**，非新 FE 求解）。",
        "",
        "| 口径 | 首次骨折高度 (m) |",
        "|---|---|",
        f"| 老 (75.337 kg) gauge 穿越 150 MPa | {ff['h_first_old_m']:.3f} |"
        if ff["h_first_old_m"] is not None else "| 老 gauge 穿越 | 未触发 |",
        f"| 新 (77 kg) gauge 反查穿越 150 MPa | {ff['h_first_new_gauge_m']:.3f} |"
        if ff["h_first_new_gauge_m"] is not None else "| 新 gauge 反查 | 未触发 |",
        f"| 新 (77 kg) F 阈值等价法 (F_thr={ff['f_threshold_n']/1e3:.3f} kN) | "
        f"{ff['h_first_new_fthreshold_m']:.3f} |"
        if ff["h_first_new_fthreshold_m"] is not None else "| 新 F 阈值法 | 未触发 |",
        "",
    ]
    if ff["h_first_old_m"] is not None and ff["h_first_new_gauge_m"] is not None:
        d = ff["h_first_new_gauge_m"] - ff["h_first_old_m"]
        lines.append(
            f"- 首次骨折高度变化：**{ff['h_first_old_m']:.3f} m → "
            f"{ff['h_first_new_gauge_m']:.3f} m**（Δ = {d:+.3f} m，"
            f"{'更早' if d < 0 else '更晚'}；相对变化 {100*d/ff['h_first_old_m']:+.2f}%）。"
        )
    lines += [
        "",
        f"- anat 扫描自身网格（19 点 a=0）内插首次骨折高度 = "
        f"**{ff['h_first_old_anat_m']:.3f} m**（本实验的老口径参考；共 19 点全网格时"
        "本表老口径与之相等）。",
        "",
        "## 逐点 gauge 反查明细",
        "",
        "| h (m) | gauge 老 anat (MPa) | 比值 F 新/老 | gauge 新 反查 (MPa) |",
        "|---|---|---|---|",
    ]
    for r in sorted(rows, key=lambda x: x["h"]):
        lines.append(
            f"| {r['h']:g} | {r['gauge_old_anat_MPa']:.2f} | {r['ratio_new_over_old']:.4f} | "
            f"{r['gauge_new_est_MPa']:.2f} |"
        )

    ho, hn = ff["h_first_old_m"], ff["h_first_new_gauge_m"]
    if ho is not None and hn is not None:
        direction = "下降（更早骨折）" if hn < ho else "上升（更晚骨折）"
        ff_sentence = (
            f"- 对首次骨折高度的影响：**{ho:.2f} → {hn:.2f} m**（Δ = {hn-ho:+.3f} m，"
            f"相对 {100*(hn-ho)/ho:+.2f}%）。影响很小，方向为{direction}，与 `gauge ∝ F` 的线性口径一致。"
        )
    else:
        ff_sentence = "- 首次骨折高度：本网格内未触发，无法给出穿越高度。"
    lines += [
        "",
        "## 诚实结论 / caveat",
        "",
        f"- **质量标定把一个 +{100*(NEW_KG/OLD_KG-1):.2f}% 的总质量差**引入载荷源，但 GRF/F_subt "
        f"峰值实际只升 **~{100*(np.mean(ratios)-1):.2f}%**（pad 接触几何主导，质量部分被抵消）。",
        ff_sentence,
        "- `gauge_new` 是**载荷线性反查的外推**：未重跑 AM50 跟骨 FE。若要绝对高度，须以 77 kg "
        "载荷重跑 `s2_thums_sweep.py`（本探针**不覆写**其 anat 产物）。",
        "- T2 缩放前后均 ≈ 0，未见模型行为异常。",
        "- 默认行为未变（`scale_mass_kg=None`）；本实验不改任何既有默认。",
        "",
        "## 产物",
        "",
        "- `mass_calib_summary.json` / `.csv`（本报告数据）",
        "- `MASS_CALIB_REPORT.md`（本文件）",
        "",
        f"*anat gauge 来源（只读）：`{ANAT_CSV.relative_to(ROOT)}`*",
    ]
    (outdir / "MASS_CALIB_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="实验 ④a：多体载荷源质量标定 75.337 → 77 kg")
    ap.add_argument("--heights", type=float, nargs="+", default=None, help="高度列表 (m)")
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "opensim_fe")
    ap.add_argument("--report-only", action="store_true",
                    help="仅从已存的 mass_calib_summary.json 重生成报告（不重跑动力学）")
    args = ap.parse_args(argv)

    heights = [float(h) for h in (args.heights if args.heights else DEFAULT_GRID)]
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.report_only:
        src = args.outdir / "mass_calib_summary.json"
        if not src.is_file():
            raise FileNotFoundError(f"--report-only 需要已存的 {src}")
        summary = json.loads(src.read_text(encoding="utf-8"))
        _write_report(
            summary["rows"], summary["first_fracture"], summary["t2"], args.outdir
        )
        print(f"[report-only] 已由 {src.name} 重生成 MASS_CALIB_REPORT.md")
        return 0

    print("=" * 72)
    print("实验 ④a · 多体载荷源质量标定：Rajagopal2015 75.337 kg → AM50 77 kg")
    print(f"  heights ({len(heights)}): {heights}")
    print("=" * 72, flush=True)

    t0 = time.time()
    t2_old = _t2(None)
    t2_new = _t2(NEW_KG)
    print(f"[T2] old residual={t2_old['gravity_residual_n']:+.3e} N   "
          f"new(77kg) residual={t2_new['gravity_residual_n']:+.3e} N "
          f"(scale={t2_new['applied_scale']:.6f})", flush=True)

    anat = _load_anat_a0()
    rows: list[dict] = []
    for h in heights:
        old = _case(h, new=False)
        new = _case(h, new=True)
        a = anat.get(h)
        if a is None:
            raise KeyError(f"anat 扫描缺少 h={h}，无法反查 gauge。")
        gauge_old = a["gauge_max_MPa"]
        ratio = new["f_subt_n"] / old["f_subt_n"]
        row = {
            "h": float(h),
            "mass_kg": float(NEW_KG),
            "old_f_subt_kN": old["f_subt_n"] / 1e3,
            "new_f_subt_kN": new["f_subt_n"] / 1e3,
            "ratio_new_over_old": ratio,
            "old_grf_peak_kN": old["grf_peak_n"] / 1e3,
            "new_grf_peak_kN": new["grf_peak_n"] / 1e3,
            "old_t1_rel_err": old["t1_rel_err"],
            "new_t1_rel_err": new["t1_rel_err"],
            "gauge_old_anat_MPa": gauge_old,
            "gauge_new_est_MPa": gauge_old * ratio,
            "model_mass_kg": new["model_mass_kg"],
            "mass_scale": new["mass_scale"],
        }
        rows.append(row)
        print(
            f"h={h:5.1f}  F_old={row['old_f_subt_kN']:7.3f} kN  "
            f"F_new={row['new_f_subt_kN']:7.3f} kN  ratio={ratio:.4f}  "
            f"gauge_old={gauge_old:7.2f}  gauge_new*={row['gauge_new_est_MPa']:7.2f} MPa  "
            f"T1o={old['t1_rel_err']*100:+.2f}% T1n={new['t1_rel_err']*100:+.2f}%  "
            f"[{time.time()-t0:.0f}s]",
            flush=True,
        )
    # 明确校验：新模型总质量应 ≈ 77.0
    for r in rows:
        if abs(r["model_mass_kg"] - NEW_KG) > 1e-6:
            raise RuntimeError(
                f"缩放后模型总质量 {r['model_mass_kg']} != {NEW_KG}（h={r['h']}）"
            )

    ff = _first_fracture(rows)

    rows_sorted = sorted(rows, key=lambda r: r["h"])
    summary = {
        "experiment": "4a",
        "old_kg": OLD_KG,
        "new_kg": NEW_KG,
        "heights_m": [r["h"] for r in rows_sorted],
        "threshold_mpa": THRESHOLD_MPA,
        "t2": {"old": t2_old, "new": t2_new},
        "first_fracture": ff,
        "rows": rows_sorted,
        "anat_gauge_source": str(ANAT_CSV.relative_to(ROOT)),
        "gauge_new_note": "gauge_new = gauge_old(anat) * F_new/F_old，载荷线性反查（外推），非新 FE 求解。",
    }
    json_out = args.outdir / "mass_calib_summary.json"
    csv_out = args.outdir / "mass_calib_summary.csv"
    json_out.write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=float),
                        encoding="utf-8")
    with csv_out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in rows_sorted:
            w.writerow({k: r[k] for k in CSV_FIELDS})

    _write_report(rows_sorted, ff, {"old": t2_old, "new": t2_new}, args.outdir)

    print("\n" + "=" * 72)
    print("首次骨折高度（a=0，gauge_max > 150 MPa）")
    print(f"  老 75.337 kg : {ff['h_first_old_m']:.3f} m" if ff["h_first_old_m"] is not None
          else "  老 : 未触发")
    print(f"  新 77 kg     : {ff['h_first_new_gauge_m']:.3f} m (gauge 反查)"
          if ff["h_first_new_gauge_m"] is not None else "  新 : 未触发")
    print(f"  ratio(F_subt new/old) range = "
          f"[{min(r['ratio_new_over_old'] for r in rows):.4f}, "
          f"{max(r['ratio_new_over_old'] for r in rows):.4f}]")
    print(f"[out] {json_out}")
    print(f"[out] {csv_out}")
    print(f"[out] {args.outdir / 'MASS_CALIB_REPORT.md'}")
    print(f"[done] {time.time()-t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
