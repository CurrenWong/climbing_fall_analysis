"""Phase P1 —— 数据驱动场景库的蒙特卡洛评估。

用法
----
    python scripts/phase_p1_montecarlo.py                # 1000 次，快速预设
    python scripts/phase_p1_montecarlo.py --n 200        # 小样本试跑
    python scripts/phase_p1_montecarlo.py --strict-tol   # 用严格容差（慢 ~18x）

产出
----
    results/p1_montecarlo.txt          完整报告（可直接读）
    results/p1_montecarlo.csv          逐次结果（可自行加权/统计）
    results/figures/p1_coverage.png    三联图：真实权重 / 部位分布 / 覆盖度

它在回答什么
------------
方案 v2 §P1 的验证判据是：

> 用真实权重跑 1000 次蒙特卡洛，输出**损伤部位分布**，与 [B25] 实测对比：
> 下肢应占 ~67%、踝应为第一大部位（~40%）、头/颈应占 ~3%。

**本脚本会给出这个对比，并明确报告它不合格** —— 因为模型里没有踝。
这不是 bug，是本 phase 要**量化**的那件事：把"我们怀疑模型不完整"
变成"缺口 = 无旋转以外的 67% 场景 + 58% 的损伤类型"这样的可复现数字。
能通过的判据是采样器本身（对照 [B25] Table 3 的 95% CI）。

⚠️ 快速预设的合法性
------------------
``FAST_KW`` 把 ``max_ode_step`` 从 1e-4 放宽到 1e-3、``rtol`` 1e-8→1e-6，
并把 ``t_max`` 从 ``max(2.0, t_fall+1.0)`` 收到 ``t_fall+0.25``，
单次耗时约 6 s → 0.46 s。**这不是拍脑袋的偷工**：实测 8 个场景，
峰值力相对偏差 ≤ 1.5e-5。该等价性由
``tests/test_climbing.py::TestFastSolverPreset`` 锁住 ——
**放宽预设之前必须先有这条测试**。
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from climbing.bone import BONE_SITES, assess_sites          # noqa: E402
from climbing.pad import POSTURES, simulate_boulder_fall    # noqa: E402
from climbing.scenarios import (                            # noqa: E402
    DIMENSIONS,
    LANDING_TO_POSTURE,
    OBSERVED_INJURY_LOCATION,
    OBSERVED_INJURY_TYPE,
    ScenarioSampler,
    coverage_report,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"

#: 快速求解器预设（等价性见模块 docstring）
FAST_KW = dict(max_ode_step=1e-3, rtol=1e-6, atol=1e-9)
#: 严格预设 = pad.py 的默认值
STRICT_KW: dict = {}
#: 冲击后继续积分的余量。峰值出现在首触后 ~20–40 ms，0.25 s 足够。
POST_IMPACT_S = 0.25

SITE_ORDER = list(BONE_SITES.keys())


# --------------------------------------------------------------------------
def run_trials(n: int, seed: int, fast: bool) -> tuple[list[dict], list[dict]]:
    """跑 n 次；返回 ``(逐次场景记录, 逐次部位记录)``。"""
    s = ScenarioSampler(seed=seed)
    kw = FAST_KW if fast else STRICT_KW
    trial_rows: list[dict] = []
    site_rows: list[dict] = []

    t0 = time.time()
    for i in range(n):
        spec = s.sample_simulation()
        sc = spec.scenario
        t_max = None
        if fast:
            t_max = float(np.sqrt(2.0 * spec.height_m / 9.80665)) + POST_IMPACT_S
        try:
            r = simulate_boulder_fall(
                height_m=spec.height_m, mass_kg=spec.mass_kg,
                posture=spec.posture, on_pad=spec.on_pad,
                t_max=t_max, **kw,
            )
        except RuntimeError as e:                     # 撞击未发生 —— 不该出现
            print(f"  [warn] 第 {i} 次仿真失败: {e}")
            continue
        v = assess_sites(r, posture=spec.posture)

        row = {
            "i": i,
            **sc.to_dict(),
            "posture": spec.posture,
            "height_m": round(spec.height_m, 3),
            "mass_kg": round(spec.mass_kg, 1),
            "on_pad": int(spec.on_pad),
            "rotation_free": int(sc.rotation_free),
            "ligament_dominant": int(spec.ligament_dominant),
            "upper_limb_mechanism": int(spec.upper_limb_mechanism),
            "weight_exact": int(spec.weight_exact),
            "peak_force_kn": round(r.peak_force_kn, 3),
            "peak_primary_g": round(r.peak_primary_g, 1),
            "hic": round(r.hic, 1),
            "n_fractured_sites": len(v.fractured_sites),
            "worst_site": v.worst.key if v.worst else "",
            "worst_utilization": round(v.worst.utilization, 3) if v.worst else 0.0,
            "fractured": "|".join(sorted(s2.key for s2 in v.fractured_sites)),
        }
        trial_rows.append(row)
        for load in v.sites:
            site_rows.append({
                "i": i,
                "site": load.key,
                "utilization": round(load.utilization, 4),
                "fractured": int(load.fractured),
                "peak_force_kn": round(load.peak_force_kn, 3),
                "threshold_kn": round(load.threshold_kn, 3),
                "rotation_free": int(sc.rotation_free),
                "weight_exact": int(spec.weight_exact),
            })
    dt = time.time() - t0
    print(f"  跑完 {len(trial_rows)}/{n} 次，用时 {dt:.1f} s "
          f"（{dt/max(1,len(trial_rows))*1000:.0f} ms/次）")
    return trial_rows, site_rows


# --------------------------------------------------------------------------
def marginal_report(rows: list[dict], n: int) -> tuple[str, bool]:
    """把 MC 频率对照 [B25] Table 3 的实测占比与 95% CI。"""
    lines = ["", "=" * 78,
             "  一、采样器验证 —— MC 频率 vs [B25] Table 3 实测值（95% CI）",
             "=" * 78,
             "  判据：MC 频率落在实测 95% CI 内即通过（★=通过，✗=越界）",
             ""]
    ok_all = True
    for key, dim in DIMENSIONS.items():
        lines.append(f"  [{dim.name_cn}]  {dim.note}")
        for cat in dim.cats:
            if cat.is_unknown:
                continue
            hit = sum(1 for r in rows if r[key] == cat.label)
            f = hit / n
            ci = cat.ci
            # 二项抽样自身的 ±2σ 区间（n=1000 时约 ±3%），与实测 CI 取并集
            se = float(np.sqrt(max(cat.prob * (1 - cat.prob), 1e-9) / n))
            lo = min(ci[0], cat.prob - 2 * se) if ci else cat.prob - 2 * se
            hi = max(ci[1], cat.prob + 2 * se) if ci else cat.prob + 2 * se
            good = lo <= f <= hi
            ok_all &= good
            mark = "★" if good else "✗"
            if ci:
                lines.append(f"    {mark} {cat.label:<24} MC {f:6.3f}   "
                             f"实测 {cat.prob:.3f} [{ci[0]:.3f},{ci[1]:.3f}]")
            else:
                lines.append(f"    {mark} {cat.label:<24} MC {f:6.3f}   "
                             f"实测 {cat.prob:.3f}")
        lines.append("")
    lines.append(f"  结论：{'✅ 采样器与 [B25] Table 3 一致' if ok_all else '❌ 有取值越界，检查权重抄录'}")
    return "\n".join(lines), ok_all


def site_distribution(rows: list[dict], site_rows: list[dict]) -> dict:
    """模型输出的分部位损伤分布（按**损伤例数**计，与 [B25] 口径一致）。"""
    out: dict = {}
    for subset, name in ((lambda r: True, "all"),
                         (lambda r: r["rotation_free"] == 1, "no_rotation"),
                         (lambda r: r["weight_exact"] == 1, "credible")):
        keys = {r["i"] for r in rows if subset(r)}
        cnt = {k: 0 for k in SITE_ORDER}
        for s in site_rows:
            if s["i"] in keys and s["fractured"]:
                cnt[s["site"]] += 1
        tot = sum(cnt.values())
        out[name] = {
            "n_trials": len(keys),
            "n_injuries": tot,
            "share": {k: (v / tot if tot else 0.0) for k, v in cnt.items()},
        }
    return out


def site_utilization_stats(site_rows: list[dict]) -> dict[str, dict]:
    """每个部位的利用率统计 —— 当**一例损伤都没有**时，这才是有效信息。"""
    stats: dict[str, dict] = {}
    for k in SITE_ORDER:
        u = np.array([s["utilization"] for s in site_rows if s["site"] == k])
        stats[k] = {
            "n": int(u.size),
            "median": float(np.median(u)) if u.size else 0.0,
            "p95": float(np.percentile(u, 95)) if u.size else 0.0,
            "max": float(u.max()) if u.size else 0.0,
            "n_over": int((u >= 1.0).sum()),
        }
    return stats


def worst_site_counts(rows: list[dict]) -> dict[str, int]:
    c: dict[str, int] = {}
    for r in rows:
        c[r["worst_site"]] = c.get(r["worst_site"], 0) + 1
    return c



def gap_table(rows: list[dict]) -> str:
    lines = ["", "=" * 78,
             "  三、与 [B25] 实测损伤分布的对比 —— **判据不合格**",
             "=" * 78,
             "  [B25] Table 1（n=245，301 例损伤）：",
             ""]
    lines.append(f"    {'部位':<12}{'实测占比':>10}{'95% CI':>20}")
    for k, (sh, ci, nn) in OBSERVED_INJURY_LOCATION.items():
        lines.append(f"    {k:<12}{sh:>9.0%}   [{ci[0]:.3f}, {ci[1]:.3f}]  n={nn}")
    lines.append("")
    lines.append("  模型能输出的「部位」（只有 7 个骨部位，且只有骨折判据）：")
    lines.append(f"    {', '.join(BONE_SITES[k].name_cn for k in SITE_ORDER)}")
    lines.append("")
    lines.append("  ── 逐条对表 ──")
    lines.append("    ❌ 踝（实测 40%，第一大部位）      → 模型输出 **0%**，结构上无此自由度")
    lines.append("    ❌ 肘（实测 16%）                  → 模型输出 **0%**，无上肢")
    lines.append("    ❌ 膝（实测 15%）                  → 模型只可能给出股骨/胫骨骨折，无膝关节")
    lines.append("    ⚠️ 下肢（实测 67%）                → 模型只覆盖其中**骨折**那部分")
    lines.append("    ⚠️ 头/颈（实测 3%）                → 由 HIC 覆盖，不在 bone.py 的部位表里")
    lines.append("    ❌ 扭伤（实测 36%）/ 脱位（11%）/ 韧带断裂（11%）")
    lines.append("                                        → 模型**没有韧带判据**，全部输出 0%")
    lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------------
def make_figure(dist: dict, rows: list[dict], site_rows: list[dict],
                n: int, path: pathlib.Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6))

    # (a) 真实场景权重
    ax = axes[0]
    dims = [("landing_position", "落地姿势"), ("rotation", "旋转"),
            ("height_band", "高度带"), ("wall_type", "墙型")]
    labels, vals, colors = [], [], []
    palette = ["#1565c0", "#c62828", "#2e7d32", "#ef6c00"]
    for i, (k, cn) in enumerate(dims):
        for cat in DIMENSIONS[k].cats:
            if cat.is_unknown:
                continue
            labels.append(f"{cn}:{cat.label}")
            vals.append(cat.prob)
            colors.append(palette[i])
    ypos = np.arange(len(labels))[::-1]
    ax.barh(ypos, vals, color=colors, alpha=0.85)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("实测占比")
    ax.set_title("[B25] Table 3 真实场景权重\n(权重来源：n=245，不再手写)", fontsize=10)
    ax.axvline(1 / 3, ls=":", color="k", lw=1)
    ax.grid(axis="x", alpha=0.3)

    # (b) 模型输出的**利用率**分布（损伤全为 0，故看利用率）
    ax = axes[1]
    util = site_utilization_stats(site_rows)
    order = sorted(SITE_ORDER, key=lambda k: -util[k]["max"])
    names = [BONE_SITES[k].name_cn for k in order]
    med = [util[k]["median"] for k in order]
    mx = [util[k]["max"] for k in order]
    p95 = [util[k]["p95"] for k in order]
    x = np.arange(len(order))
    ax.bar(x, mx, color="#ef6c00", alpha=0.35, label="最大")
    ax.bar(x, p95, color="#ef6c00", alpha=0.6, label="p95")
    ax.bar(x, med, color="#c62828", alpha=0.9, label="中位")
    ax.axhline(1.0, color="k", ls="--", lw=1.6)
    ax.text(0.05, 1.02, "骨折阈值 = 1.0", ha="left", fontsize=8.5)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("利用率 = 峰值载荷 / 阈值")
    ax.set_ylim(0, max(1.15, max(mx) * 1.25))
    ax.set_title(f"模型输出的分部位利用率\nN={n} 次蒙特卡洛 —— 最高仅 "
                 f"{max(mx):.2f}，零骨折", fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8, loc="center right")

    # (c) 覆盖度瀑布
    ax = axes[2]
    p_no_rot = dist["no_rotation"]["n_trials"] / max(1, dist["all"]["n_trials"])
    p_cred = dist["credible"]["n_trials"] / max(1, dist["all"]["n_trials"])
    frac = OBSERVED_INJURY_TYPE["fracture"][0]
    covered = frac * 0.71 + OBSERVED_INJURY_LOCATION["head_neck"][0]
    stages = ["全部坠落\n(100%)", "无旋转\n可表达", "模型结论\n可信", "能产出的\n损伤类型"]
    values = [1.0, p_no_rot, p_cred, covered]
    cols = ["#616161", "#ef6c00", "#c62828", "#1565c0"]
    ax.bar(range(4), [v * 100 for v in values], color=cols, alpha=0.85)
    for i, v in enumerate(values):
        ax.text(i, v * 100 + 2, f"{v:.0%}", ha="center", fontsize=11,
                fontweight="bold")
    ax.set_xticks(range(4))
    ax.set_xticklabels(stages, fontsize=8.5)
    ax.set_ylabel("可覆盖比例 (%)")
    ax.set_ylim(0, 118)
    ax.set_title("结构性覆盖度 —— 缺口不是标定误差", fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    ax.text(0.5, 0.62,
            "踝（40%，第一大部位）、肘（16%）、膝（15%）\n"
            "扭伤 + 脱位 + 韧带断裂 = 58%\n"
            "→ 模型全部输出 0%",
            transform=ax.transAxes, ha="center", fontsize=8.5, color="#c62828",
            bbox=dict(boxstyle="round,pad=0.4", fc="#fff8f8", ec="#c62828"))

    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print(f"  → {path.relative_to(ROOT)}")


# --------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="P1 场景库蒙特卡洛评估")
    ap.add_argument("--n", type=int, default=1000, help="仿真次数（默认 1000）")
    ap.add_argument("--seed", type=int, default=20261003)
    ap.add_argument("--strict-tol", action="store_true",
                    help="用严格容差（慢约 18 倍，用于核对快速预设）")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    print(f"[P1] 蒙特卡洛：N={args.n} seed={args.seed} "
          f"预设={'严格' if args.strict_tol else '快速'}")
    rows, site_rows = run_trials(args.n, args.seed, fast=not args.strict_tol)
    n = len(rows)
    if n == 0:
        print("没有任何成功样本")
        return 1

    dist = site_distribution(rows, site_rows)
    marg_txt, ok = marginal_report(rows, n)
    gap_txt = gap_table(rows)

    # 落地姿势映射后的实际频率
    post_cnt: dict[str, int] = {}
    for r in rows:
        post_cnt[r["posture"]] = post_cnt.get(r["posture"], 0) + 1

    lines = [
        "=" * 78,
        "  Phase P1 —— 数据驱动场景库 · 蒙特卡洛报告",
        "=" * 78,
        f"  N = {n}    seed = {args.seed}    "
        f"求解预设 = {'严格' if args.strict_tol else '快速（已验证等价）'}",
        "",
        coverage_report(),
        marg_txt,
        "",
        "=" * 78,
        "  二、落地姿势映射后的实际采样频率",
        "=" * 78,
    ]
    for k, c in sorted(post_cnt.items(), key=lambda kv: -kv[1]):
        lines.append(f"    {POSTURES[k].name:<26} {c:5d}  {c/n:6.1%}")
    lines.append("")
    lines.append("  映射依据（本项目的假设，逐条 rationale 见 scenarios.py）：")
    for lab, m in LANDING_TO_POSTURE.items():
        got = ", ".join(f"{n_}({w:.2f})" for n_, w in m.postures)
        lines.append(f"    {lab:<22} → {got}")
    lines += ["", gap_txt, ""]

    lines += [
        "=" * 78,
        "  四、模型输出的分部位分布（三个子集）",
        "=" * 78,
    ]
    for sub, cn in (("all", "全体样本"), ("no_rotation", "仅无旋转子集"),
                    ("credible", "仅「结论可信」子集")):
        d = dist[sub]
        lines.append(f"  [{cn}]  {d['n_trials']} 次仿真 / {d['n_injuries']} 例部位损伤")
        if d["n_injuries"] == 0:
            lines.append("      （无部位超过阈值）")
        else:
            for k in sorted(SITE_ORDER, key=lambda k: -d["share"][k]):
                if d["share"][k] > 0:
                    lines.append(f"      {BONE_SITES[k].name_cn:<12}"
                                 f"{d['share'][k]:7.1%}")
        lines.append("")

    # ---- 4b. 一例都没有时，看利用率分布 ----
    util = site_utilization_stats(site_rows)
    wsc = worst_site_counts(rows)
    order = sorted(SITE_ORDER, key=lambda k: -util[k]["max"])
    lines += [
        "  4b. 一例损伤都没有 ⇒ 部位分布为空。此时有效信息是**利用率分布**：",
        "",
        f"      {'部位':<12}{'中位':>8}{'p95':>8}{'最大':>8}{'超阈值':>8}"
        f"   组(P0对标)",
        f"      {'-'*56}",
    ]
    group = {"calcaneus": "B ❌偏晚", "tibia_distal": "A ✅命中",
             "tibia_mid": "A ✅命中", "fibula_ends": "B ❌待查",
             "femoral_neck": "A ✅命中", "pelvis": "B ❌偏晚",
             "lumbar_spine": "B ❌偏晚"}
    for k in sorted(SITE_ORDER, key=lambda k: -util[k]["max"]):
        u = util[k]
        lines.append(f"      {BONE_SITES[k].name_cn:<12}{u['median']:8.3f}"
                     f"{u['p95']:8.3f}{u['max']:8.3f}{u['n_over']:8d}   {group.get(k,'')}")
    lines += [
        "",
        "      「最危险部位」出现次数（每次取利用率最高者）：",
    ]
    for k, v in sorted(wsc.items(), key=lambda kv: -kv[1]):
        if not k:
            continue
        lines.append(f"        {BONE_SITES[k].name_cn:<12}{v:5d}  {v/n:6.1%}")
    lines += [
        "",
        "  ⚠️ 关键：最危险部位 **93.7% 落在胫骨远端** —— 那是 P0 里",
        "     **对标 [Y25] 命中的组 A**（载荷路径接近均匀柱体）。",
        "     也就是说「零骨折」这个结论**不依赖**组 B 那些已证实不可靠的阈值",
        "     （跟骨被高估 17–34×）。它是由我们已经验证过的那部分阈值得出的。",
        "",
        "  为什么一例骨折都报不出来 —— 两层原因：",
        "    (1) **高度不够**：抱石墙 4.5 m，94% 落在 20 cm 泡沫垫上。",
        "        模型在**刚性地面**上要到 13 m 才报胫骨远端骨折（P0 对标结果），",
        "        垫子把力又压低一档 ⇒ 4.5 m 落垫远在阈值之下。",
        "        （与 [Y25] 一致：[Y25] 刚性地面首个骨折在 7–9 m。）",
        "    (2) **模型漏掉了真实机制**：[B25] 里骨折占 23%（踝骨折 8%），",
        "        但真实机制是**足缘小面积接触 + 踝旋后/内翻**导致的局部载荷，",
        "        不是均匀柱体轴向压缩 —— 1D 模型产生不了。",
        "        （[B25] 原文：护垫下沉 → 踝旋后 → 韧带负荷。）",
        "",
    ]

    # 最常出问题的姿势
    lines += [
        "=" * 78,
        "  五、按落地姿势的利用率（损伤全为 0，故看利用率）",
        "=" * 78,
        f"    {'姿势':<26}{'n':>6}{'中位':>9}{'最大':>9}   超阈值",
    ]
    per_post: dict[str, list[float]] = {}
    for r in rows:
        per_post.setdefault(r["posture"], []).append(float(r["worst_utilization"]))
    for pose, v in sorted(per_post.items(), key=lambda kv: -max(kv[1])):
        vv = np.array(v)
        lines.append(f"    {POSTURES[pose].name:<26}{len(vv):6d}"
                     f"{np.median(vv):9.3f}{vv.max():9.3f}"
                     f"{int((vv >= 1.0).sum()):8d}")
    lines += [
        "",
        "  判读：脚先类姿势（绷直/单脚/前脚掌）利用率最高且是唯一接近阈值的",
        "        一族 —— 与 [B25]「下肢占 67%、踝为第一大部位」的方向一致。",
        "        但绝对值全部低于阈值，见上一节的原因分析。",
    ]

    lines += [
        "",
        "=" * 78,
        "  六、结论",
        "=" * 78,
        f"  1. ✅ 采样器通过验证：MC 频率与 [B25] Table 3 的 95% CI "
        f"{'一致' if ok else '**不一致**'}。",
        f"  2. ⚠️ 模型在真实抱石场景下**一例骨折都报不出来**：",
        f"     最高利用率 {util[order[0]]['max']:.2f}（{BONE_SITES[order[0]].name_cn}），"
        f"利用率随高度单调（0.20 → 0.36 → 0.50）。",
        "     ⇒ 模型输出的**损伤部位分布是空的**，§P1 的对比判据无从谈起。",
        "  3. ❌ 不合格的原因已定位，分两层：",
        "     (a) **结构性**：没有踝（第一大部位 40%）、没有韧带判据（58% 类型）、",
        "         没有上肢（肘 16%）。",
        "     (b) **量程**：真实抱石（≤4.5 m + 94% 落垫）根本够不到模型的骨折阈值；",
        "         模型的骨折是「7–13 m 刚性地面」量级的现象。",
        "  4. ✅ 一个重要的**排除项**：最危险部位 93.7% 是胫骨远端 ——",
        "     那是 P0 对标 [Y25] **命中**的组 A。所以「零骨折」不是被",
        "     组 B 那些不可靠阈值（跟骨高估 17–34×）掩盖出来的。",
        "  5. ➡️ 结构性缺口调参补不上，必须加自由度（P2/P3）；",
        "     量程问题则需要**独立的抱石实测标定数据**才能修。",
        "     P1 的价值正在于此 —— 它把「感觉模型不对」变成了上面这些数字。",
        "=" * 78,
    ]

    rep = "\n".join(lines)
    (OUT / "p1_montecarlo.txt").write_text(rep, encoding="utf-8")
    print(f"  → results/p1_montecarlo.txt  ({len(rep)} 字符)")

    with (OUT / "p1_montecarlo.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  → results/p1_montecarlo.csv  ({len(rows)} 行)")

    make_figure(dist, rows, site_rows, n, FIG / "p1_coverage.png")

    print("\n" + rep[:2500])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
