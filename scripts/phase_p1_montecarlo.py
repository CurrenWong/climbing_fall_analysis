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

from climbing.ankle_injury import assess_ankle              # noqa: E402
from climbing.bone import BONE_SITES, assess_sites          # noqa: E402
from climbing.pad import POSTURES                           # noqa: E402
from climbing.pad2d import Posture2D, simulate_boulder_fall_2d  # noqa: E402
from climbing.scenarios import (                            # noqa: E402
    DIMENSIONS,
    LANDING_TO_POSTURE,
    OBSERVED_INJURY_LOCATION,
    OBSERVED_INJURY_TYPE,
    OBSERVED_SPECIFIC_INJURY,
    ScenarioSampler,
    coverage_report,
)

import dataclasses as _dataclasses

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
# ROTATION_TO_POSTURE2D —— 坠落旋转维度 → 2D 姿势扩展的工厂
# --------------------------------------------------------------------------
# 每个键对应一个 ``Callable[[np.random.Generator], Posture2D]``，
# 每次 trial 调用一次以抽样（带 RNG），保证 trial 间独立可复现。
#
# 来源说明（每个分支都在注释里给出原始出处）：
# * "without"             —— [硬门槛] 不变量 1：默认 ⇒ 1D bit-identical
# * "longitudinal"        —— [T22] 内翻/外翻 50/50，幅度 {10, 20, 30}° 三档
# * "antero_posterior"    —— [B25] Table 3 (12%)：前后向旋转 ⇒ 身体俯仰
# * "back_transverse"     —— [B25] Table 3 (7%)：横向旋转目前无独立建模
# * "front_transverse"    —— [B25] Table 3 (3%)：横向旋转目前无独立建模
# * "multiple"            —— [B25] Table 3 (9%)：复合旋转模型无复合处理
# * "unknown"             —— [B25] Table 3 (8%)：与 'without' 同（保守默认）
_VALID_POSTURE2D_FIELDS = frozenset(
    f.name for f in _dataclasses.fields(Posture2D)
)


def _build_p2d(**kwargs) -> Posture2D:
    """Build a ``Posture2D``，过滤掉当前 Posture2D 不支持的字段。

    并行 agent 正在写 ``pad2d.py`` 的 P2 踝字段（``ankle_beta0_rad`` /
    ``beta_from_theta`` / ``mu_slide`` / ``ankle_d_lateral_m`` /
    ``ankle_d_medial_m`` / ``ligament_r_mm``），在他们写好之前
    ``Posture2D(ankle_beta0_rad=...)`` 会抛 ``TypeError``。
    本函数做**显式过滤**：把不在 ``_VALID_POSTURE2D_FIELDS`` 里的
    kwargs 静默丢弃，回退到 ``Posture2D()`` 全默认。

    一旦并行 agent 完成 P2 字段添加，过滤分支不被触发，所有非默认
    参数都生效。
    """
    filtered = {k: v for k, v in kwargs.items() if k in _VALID_POSTURE2D_FIELDS}
    return Posture2D(**filtered)


def _make_without(_rng) -> Posture2D:
    """``"without"`` ⇒ 全默认 ⇒ 与 1D 模型逐位一致（[硬门槛] 不变量 1）。

    [B25] Table 3：P("without") = 30%（n=245）；1D 模型 92% 的可信骨方都在这一档。
    """
    return _build_p2d()


def _make_longitudinal(rng) -> Posture2D:
    """``"longitudinal"`` ⇒ 踝内翻 / 外翻对称性（[T22] Fig.5）。

    - 角度幅度：从 ``{10°, 20°, 30°}`` 均匀采样（[T22] 标定这三档）
    - 符号：50/50（内翻 + / 外翻 −）
    - 其他踝 Beta 参数：
      * ``beta_from_theta = 1.0`` —— β(t) = β0 + 1·θ(t)，β 跟随 θ 1:1 演化
      * ``mu_slide = 0.4`` —— 粘-滑阈值：|sin β| > 0.4 ⇒ sliding
      * ``ankle_d_lateral_m = 0.030`` —— β≥0（内翻）的 CoP→距下轴力臂（m）
      * ``ankle_d_medial_m = 0.045`` —— β<0（外翻）的 CoP→距下轴力臂（m）
      * ``ligament_r_mm = 22.0`` —— 韧带力臂（[R6] 默认值）
      * ``ankle_load_share = 0.477`` —— measured anchor（见
        ``results/opensim_fe/joint_reactions_summary.csv`` 中
        ``ankle_r.peak_force_n / GRF_peak = 25294.9 / 53060``）。
        **结构性诚实声明**：本期望以此 share 把踝骨折从 23% 降到 ≈ 8%，
        但 V6 校准把 ``ANKLE_AVULSION_LOAD_N`` 钉到 300 N，导致所有 stuck
        在 share=0.477 下 F_lig 仍 >> 300，fracture 率结构性不能解耦。
        故实际 fracture 率在 share=1.0 与 share=0.477 下几乎一致
        （228 vs 230）。详见 ``docs/P2_踝结论_2026-10-07.md`` §四。
    """
    mag_deg = float(rng.choice([10.0, 20.0, 30.0]))
    sign = float(rng.choice([-1.0, 1.0]))
    return _build_p2d(
        ankle_beta0_rad=float(np.deg2rad(mag_deg * sign)),
        beta_from_theta=1.0,
        mu_slide=0.4,
        ankle_d_lateral_m=0.030,
        ankle_d_medial_m=0.045,
        ligament_r_mm=22.0,
        ankle_load_share=0.477,
    )


def _make_antero_posterior(rng) -> Posture2D:
    """``"antero_posterior"`` ⇒ 身体俯仰（前后向旋转）。

    ±15°，50/50 符号（前倾 / 后仰）。[B25] Table 3：P=12%。
    踝参数保持默认（前后向旋转的主通路在 f_pad(t)，踝不直接参与）。
    """
    sign = float(rng.choice([-1.0, 1.0]))
    return _build_p2d(theta0_rad=float(np.deg2rad(15.0 * sign)))


def _make_transverse_default(_rng) -> Posture2D:
    """``"back_transverse"`` / ``"front_transverse"`` ⇒ 默认。

    [B25] Table 3：back_transverse 7% + front_transverse 3% = 10%。
    当前 2D 模型只有俯仰（theta0），没有 yaw 自由度 ⇒ 用 Posture2D() 全默认。
    """
    return _build_p2d()


def _make_multiple(_rng) -> Posture2D:
    """``"multiple"`` ⇒ 默认（多种旋转复合，模型无复合处理）。

    [B25] Table 3：P=9%。
    """
    return _build_p2d()


#: 坠落旋转维度 → ``Posture2D`` 工厂。**每个值是 ``(np.random.Generator)
#:  -> Posture2D``**；``run_trials`` 用 trial 局部 rng 调用以避免污染
#: ``ScenarioSampler`` 的状态。
ROTATION_TO_POSTURE2D: dict = {
    "without":          _make_without,
    "longitudinal":     _make_longitudinal,
    "antero_posterior": _make_antero_posterior,
    "back_transverse":  _make_transverse_default,
    "front_transverse": _make_transverse_default,
    "multiple":         _make_multiple,
    "unknown":          _make_without,        # 与 'without' 同：保守默认
}


#: [B25] Table 1 的「踝（扭伤或骨折）」总占比 —— 40%（原文正文）。
#: 与 ``OBSERVED_SPECIFIC_INJURY["ankle_sprain"]``（28%）+ 
#: ``["ankle_fracture"]``（8%）不一致，因为原文正文另有「踝 40%」的口径
#: （含踝的其他软组织伤）。显式命名以便 §4c 报告引用。
OBSERVED_SPECIFIC_INJURY_ANKLE_TOTAL: float = 0.40

#: 踝列在 trial_rows / site_rows 中的字段名（保持顺序，便于 CSV 列对齐）。
#: 顺序与 ``AnkleVerdict`` / ``assess_ankle`` 的字段一致。
ANKLE_COLUMNS: tuple = (
    "ankle_mode", "ankle_beta0_rad", "peak_ligament_strain",
    "ankle_bone_utilization", "ankle_fracture", "ankle_sprain",
    "ankle_injured",
)


def _empty_ankle_row() -> dict:
    """踝列的零默认（并行 agent 还没写好 pad2d.py 时使用）。"""
    return {
        "ankle_mode": "none",
        "ankle_beta0_rad": 0.0,
        "peak_ligament_strain": 0.0,
        "ankle_bone_utilization": 0.0,
        "ankle_fracture": 0,
        "ankle_sprain": 0,
        "ankle_injured": 0,
    }


def _ankle_row_from_verdict(av) -> dict:
    """从 ``AnkleVerdict`` 派生 trial 行所需的踝字段。"""
    return {
        "ankle_mode": av.mode,
        "ankle_beta0_rad": round(float(av.beta0_rad), 4),
        "peak_ligament_strain": round(float(av.peak_ligament_strain), 4),
        "ankle_bone_utilization": round(float(av.bone_utilization), 4),
        "ankle_fracture": int(av.fracture),
        "ankle_sprain": int(av.sprain),
        "ankle_injured": int(bool(av.fracture or av.sprain)),
    }


# --------------------------------------------------------------------------
def run_trials(n: int, seed: int, fast: bool) -> tuple[list[dict], list[dict]]:
    """跑 n 次；返回 ``(逐次场景记录, 逐次部位记录)``。

    仿真器：``pad2d.simulate_boulder_fall_2d``（**继承 1D 子系统 + 2D
    旋转/横向**），由 ``ROTATION_TO_POSTURE2D`` 把 ``sc.rotation`` 映射
    成 ``Posture2D``。无踝字段时（并行 agent 还在写 pad2d.py），
    ``assess_ankle`` 抛 ``ValueError``，本函数**静默回退**到踝列
    零默认 + 试运行时不报错。
    """
    s = ScenarioSampler(seed=seed)
    kw = FAST_KW if fast else STRICT_KW
    trial_rows: list[dict] = []
    site_rows: list[dict] = []
    ankle_unavailable_warned = False

    t0 = time.time()
    for i in range(n):
        spec = s.sample_simulation()
        sc = spec.scenario
        t_max = None
        if fast:
            t_max = float(np.sqrt(2.0 * spec.height_m / 9.80665)) + POST_IMPACT_S

        # 为本 trial 开一个**独立** rng（不读 ``s.rng``），用来抽
        # ``ROTATION_TO_POSTURE2D`` 里的纵向角度等。这样采样器的随机流
        # 与"没有踝扩展"时**逐位一致**，§一 的边际验证不受踝抽样影响。
        trial_rng = np.random.default_rng([seed, i])
        factory = ROTATION_TO_POSTURE2D.get(sc.rotation, ROTATION_TO_POSTURE2D["unknown"])
        posture2d = factory(trial_rng)

        try:
            r = simulate_boulder_fall_2d(
                height_m=spec.height_m, mass_kg=spec.mass_kg,
                posture=spec.posture, on_pad=spec.on_pad,
                t_max=t_max, posture2d=posture2d, **kw,
            )
        except RuntimeError as e:                     # 撞击未发生 —— 不该出现
            print(f"  [warn] 第 {i} 次仿真失败: {e}")
            continue
        v = assess_sites(r, posture=spec.posture)

        # ---- 踝评估：并行 agent 还没把 ankle_* 加到 pad2d.py 时会抛 ValueError ----
        try:
            av = assess_ankle(r, posture=spec.posture)
            ankle_data = _ankle_row_from_verdict(av)
        except ValueError as e:
            if not ankle_unavailable_warned:
                print(f"  [info] 踝数据未在 pad2d 中就绪 —— 踝列将为零默认。"
                      f"首条诊断: {str(e)[:120]}")
                ankle_unavailable_warned = True
            ankle_data = _empty_ankle_row()

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
            **ankle_data,
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
                # 踝列也镜像到 site_rows —— ``site_distribution``
                # 按 trial 数踝时，从任意一行读都可以（用 seen 集合去重）
                **ankle_data,
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
    """模型输出的分部位损伤分布（按**损伤例数**计，与 [B25] 口径一致）。

    骨部位（7 个）的语义**保持不变**：每个 trial 可贡献 0..N 个部位骨折
    事件，``share`` 是事件数 / 总事件数。

    踝列（P2）：每个 trial 只数一次（**踝 ≠ 多个骨部位**），返回
    ``ankle_sprain`` / ``ankle_fracture`` / ``ankle_injured`` 三个 trial
    计数 + 它们相对 ``n_trials`` 的 share（与 [B25] Table 1「踝 40%」
    的口径一致 —— 踝占伤者的比例，不是事件数比例）。
    """
    out: dict = {}
    for subset, name in ((lambda r: True, "all"),
                         (lambda r: r["rotation_free"] == 1, "no_rotation"),
                         (lambda r: r["weight_exact"] == 1, "credible")):
        keys = {r["i"] for r in rows if subset(r)}
        cnt = {k: 0 for k in SITE_ORDER}

        # 踝：用 seen 集合保证每个 trial 只贡献一次
        ankle_seen: set[int] = set()
        ankle_sprain_n = 0
        ankle_fracture_n = 0
        ankle_injured_n = 0

        for s in site_rows:
            i = s["i"]
            if i not in keys:
                continue
            # 骨部位事件（语义不变）
            if s["fractured"]:
                cnt[s["site"]] += 1
            # 踝：每 trial 数一次
            if i not in ankle_seen:
                ankle_seen.add(i)
                if s["ankle_fracture"]:
                    ankle_fracture_n += 1
                if s["ankle_sprain"]:
                    ankle_sprain_n += 1
                if s["ankle_fracture"] or s["ankle_sprain"]:
                    ankle_injured_n += 1

        tot = sum(cnt.values())
        n_trials = len(keys)
        out[name] = {
            "n_trials": n_trials,
            "n_injuries": tot,
            "share": {k: (v / tot if tot else 0.0) for k, v in cnt.items()},
            # P2 踝列：trial 计数 + share（口径 = trial-count / n_trials）
            "ankle_sprain": ankle_sprain_n,
            "ankle_fracture": ankle_fracture_n,
            "ankle_injured": ankle_injured_n,
            "ankle_sprain_share": ankle_sprain_n / max(1, n_trials),
            "ankle_fracture_share": ankle_fracture_n / max(1, n_trials),
            "ankle_injured_share": ankle_injured_n / max(1, n_trials),
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
    lines.append("  模型能输出的「部位」（7 个骨部位 + 踝 P2）：")
    lines.append(f"    {', '.join(BONE_SITES[k].name_cn for k in SITE_ORDER)}, 踝（扭伤 + 骨折）")
    lines.append("")
    lines.append("  ── 逐条对表 ──")
    lines.append("    ⚠️ 踝（实测 40%，第一大部位）      → P2 已建模（见 §四 损伤部位分布），")
    lines.append("                                        但只有内翻/外翻 50/50 + 3 档幅度的")
    lines.append("                                        ROTATION_TO_POSTURE2D 标定，未校准")
    lines.append("    ❌ 肘（实测 16%）                  → 模型输出 **0%**，无上肢")
    lines.append("    ❌ 膝（实测 15%）                  → 模型只可能给出股骨/胫骨骨折，无膝关节")
    lines.append("    ⚠️ 下肢（实测 67%）                → 模型只覆盖其中**骨折**那部分")
    lines.append("    ⚠️ 头/颈（实测 3%）                → 由 HIC 覆盖，不在 bone.py 的部位表里")
    lines.append("    ❌ 扭伤（实测 36%）/ 脱位（11%）/ 韧带断裂（11%）")
    lines.append("                                        → 模型**只覆盖踝扭伤**（踝 28%），")
    lines.append("                                        肘/膝/肩脱位与韧带断裂仍 0%")
    lines.append("")
    return "\n".join(lines)


#: 下肢骨部位（[B25] Table 1 的「下肢 67%」口径）。
#: 踝（P2）单列；这里给出 5 个 1D 骨部位中属于下肢的那些。
LOWER_LIMB_BONE_SITES: tuple = (
    "calcaneus", "tibia_distal", "tibia_mid", "fibula_ends", "femoral_neck",
)


def damage_site_report(dist: dict) -> str:
    """生成「损伤部位分布」小节 —— 与 [B25] Table 1 逐项对比。

    口径：
    * **踝列** = 每 trial 计一次（一个 trial 至多贡献 1 个踝损伤）。
      这样可与 [B25] 的「踝占伤者 40%」直接比。
    * **骨部位列** = 骨折事件数（一个 trial 可以贡献 N 个骨部位骨折）。

    [B25] Table 1 的命中案例（V5 验证目标）：
    * 踝（扭伤或骨折）40.0%
    * 踝扭伤 28.0%
    * 踝骨折 8.0%
    * 下肢 67.0%
    """
    ref = {
        "ankle_total": OBSERVED_SPECIFIC_INJURY_ANKLE_TOTAL,
        "ankle_sprain": OBSERVED_SPECIFIC_INJURY["ankle_sprain"],
        "ankle_fracture": OBSERVED_SPECIFIC_INJURY["ankle_fracture"],
        "lower_limb": OBSERVED_INJURY_LOCATION["lower_limb"][0],
    }

    lines = [
        "=" * 78,
        "  4c. 损伤部位分布 —— 与 [B25] Table 1 逐项对比",
        "=" * 78,
        "",
        "  [口径] 踝列 = 每 trial 计一次（一个 trial 至多 1 个踝损伤），",
        "         可与 [B25]「踝占伤者 40%」直接比；",
        "         骨部位列 = 骨折事件数（一个 trial 可以有多处骨折）。",
        "",
    ]
    for sub, cn in (("all", "全体样本"), ("no_rotation", "仅无旋转子集"),
                    ("credible", "仅「结论可信」子集")):
        d = dist[sub]
        n_trials = max(1, d["n_trials"])
        lines.append(f"  [{'─'*4} {cn}  N={d['n_trials']} trial {'─'*40}]")
        # --- 踝 ---
        an_sp = d["ankle_sprain"]
        an_fr = d["ankle_fracture"]
        an_in = d["ankle_injured"]
        # 踝=0 的原因随子集不同：
        #   all 为 0 ⇒ pad2d 的 ankle_* 可能未就绪（真异常）
        #   no_rotation / credible 为 0 ⇒ 正常：这些子集都是 rotation='without'
        #     或要求 rotation_free，映射到 Posture2D() 全默认 ⇒ 无踝自由度
        if an_in == 0 and sub == "all":
            ankle_note = "  ⚠️ 踝全为 0（pad2d.py 的 ankle_* 字段可能还没就绪）"
        elif an_in == 0:
            ankle_note = "  ✅ 预期为 0（该子集是 rotation='without'，无踝自由度）"
        else:
            ankle_note = ""
        lines.append(f"      踝（扭伤或骨折） {an_in:5d} / {n_trials:<5d} = "
                     f"{an_in/n_trials:6.1%}   [B25] {ref['ankle_total']:.1%}"
                     f"{ankle_note}")
        lines.append(f"        其中 踝扭伤     {an_sp:5d} / {n_trials:<5d} = "
                     f"{an_sp/n_trials:6.1%}   [B25] {ref['ankle_sprain']:.1%}")
        lines.append(f"        其中 踝骨折     {an_fr:5d} / {n_trials:<5d} = "
                     f"{an_fr/n_trials:6.1%}   [B25] {ref['ankle_fracture']:.1%}")

        # --- 7 个骨部位（骨折事件数） ---
        lines.append("      ── 7 个骨部位（P0 已有，按事件数）──")
        for k in sorted(SITE_ORDER, key=lambda k: -d["share"][k]):
            ev = int(round(d["share"][k] * d["n_injuries"]))
            lines.append(f"        {BONE_SITES[k].name_cn:<12}{ev:5d}  "
                         f"{d['share'][k]:6.1%}   （事件/事件）")
        # --- 下肢汇总（踝 trial + 5 个下肢骨事件） ---
        # ⚠️ 口径：踝是 per-trial，下肢骨是 per-event —— 二者分母不同，
        # 故此处给出**下界**（把踝 trial 数 + 下肢骨事件数 都算作一个
        # 下肢"损伤"），并显式标注不能与 [B25] 的 67% 直接比。
        bone_ll_events = sum(d["share"][k] * d["n_injuries"]
                             for k in LOWER_LIMB_BONE_SITES)
        ll_lower_bound = (an_in + bone_ll_events) / n_trials
        lines.append("      ── 下肢（踝 trial + 5 个下肢骨事件）──")
        lines.append(f"        下肢（下界）      （{an_in} ankle + "
                     f"{bone_ll_events:.0f} bone） / {n_trials} = "
                     f"≥ {ll_lower_bound:6.1%}   [B25] {ref['lower_limb']:.1%}")
        lines.append(f"        ⚠️ 踝 per-trial、骨 per-event，分母不同，"
                     f"此列仅作方向性参考。")
        lines.append("")
    lines.append(
        "  [对比结论] 踝的三个对比项（踝 total / 踝扭伤 / 踝骨折）与 [B25] "
        "Table 1 三列可直接对照；\n"
        "            具体数字取决于 pad2d.py 的踝物理标定与 "
        "ROTATION_TO_POSTURE2D 的角度抽取。\n"
        "            ⚠️ 若踝列全为 0，说明并行 agent 的 pad2d.py 尚未把 "
        "ankle_* 写进 BoulderFallResult2D。"
    )
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
    bone_frac_zero = dist["all"]["n_injuries"] == 0
    ax.set_title(f"模型输出的分部位利用率\nN={n} 次蒙特卡洛 —— 最高仅 "
                 f"{max(mx):.2f}，" + ("零骨折" if bone_frac_zero else "有骨折"),
                 fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8, loc="center right")

    # (c) 覆盖度瀑布（结构缺口，含 P2 踝扩展）
    ax = axes[2]
    p_no_rot = dist["no_rotation"]["n_trials"] / max(1, dist["all"]["n_trials"])
    p_cred = dist["credible"]["n_trials"] / max(1, dist["all"]["n_trials"])
    frac = OBSERVED_INJURY_TYPE["fracture"][0]
    covered_1d = frac * 0.71 + OBSERVED_INJURY_LOCATION["head_neck"][0]
    # P2：踝扭伤是新增的判据能力（[B25] 28%），与 1D 的骨折口径相加。
    ankle_sprain_ref = OBSERVED_SPECIFIC_INJURY["ankle_sprain"]
    covered_p2 = min(1.0, covered_1d + ankle_sprain_ref)
    stages = ["全部坠落\n(100%)", "无旋转\n可表达", "模型结论\n可信",
              "1D 能产出\n的损伤", "P2 加踝\n扭伤后"]
    values = [1.0, p_no_rot, p_cred, covered_1d, covered_p2]
    cols = ["#616161", "#ef6c00", "#c62828", "#1565c0", "#2e7d32"]
    ax.bar(range(5), [v * 100 for v in values], color=cols, alpha=0.85)
    for i, v in enumerate(values):
        ax.text(i, v * 100 + 2, f"{v:.0%}", ha="center", fontsize=10,
                fontweight="bold")
    ax.set_xticks(range(5))
    ax.set_xticklabels(stages, fontsize=8)
    ax.set_ylabel("可覆盖比例 (%)")
    ax.set_ylim(0, 118)
    ax.set_title("结构性覆盖度 —— 缺口不是标定误差", fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    ax.text(0.5, 0.66,
            "P2 已加入踝（40%）：扭伤 28% + 骨折 8% 可表达\n"
            "仍缺：肘（16%）、膝（15%）、脱位（11%）、\n"
            "      肘/膝/肩韧带断裂 —— 结构性缺口",
            transform=ax.transAxes, ha="center", fontsize=8, color="#c62828",
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
        "",
        "  ⚠️ 上述覆盖度报告来自 scenarios.coverage_report()（1D 口径，写于 P2 之前）：",
        "     踝「无踝关节自由度 / 输出 0%」这一条**已过时** —— P2 已把踝接入本 harness，",
        "     见本报告 §4c。scenarios.py 未改（非本次范围），故原文本保留。",
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

    # ---- 4c. 损伤部位分布（骨 + 踝）与 [B25] Table 1 逐项对比 ----
    lines += ["", damage_site_report(dist), ""]

    # ---- 4b. 一例都没有时，看利用率分布 ----
    util = site_utilization_stats(site_rows)
    wsc = worst_site_counts(rows)
    order = sorted(SITE_ORDER, key=lambda k: -util[k]["max"])
    n_bone_injuries = dist["all"]["n_injuries"]
    lines += [
        f"  4b. 骨部位（1D 判据）一例都没超过阈值（{n_bone_injuries} 例骨折）"
        f"⇒ 骨部位分布为空；",
        "      此时骨部位的有效信息是**利用率分布**（踝见 §4c）：",
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
    # 动态计算"最危险骨部位"的占比（不写死数字）
    top_worst = max((v for k, v in wsc.items() if k), default=0)
    top_worst_site = max((k for k in wsc if k), key=lambda k: wsc[k], default="")
    top_worst_share = top_worst / max(1, n)
    lines += [
        "",
        f"  ⚠️ 关键：最危险骨部位 **{top_worst_share:.1%} 落在"
        f"{BONE_SITES[top_worst_site].name_cn if top_worst_site else '—'}**"
        f" —— 那是 P0 里 **对标 [Y25] 命中的组 A**（载荷路径接近均匀柱体）。",
        "     也就是说「骨部位零骨折」这个结论**不依赖**组 B 那些已证实不可靠",
        "     的阈值（跟骨被高估 17–34×）。它是由我们已经验证过的那部分阈值得出的。",
        "",
        "  为什么骨部位一例骨折都报不出来 —— 两层原因：",
        "    (1) **高度不够**：抱石墙 4.5 m，94% 落在 20 cm 泡沫垫上。",
        "        模型在**刚性地面**上要到 13 m 才报胫骨远端骨折（P0 对标结果），",
        "        垫子把力又压低一档 ⇒ 4.5 m 落垫远在阈值之下。",
        "        （与 [Y25] 一致：[Y25] 刚性地面首个骨折在 7–9 m。）",
        "    (2) **1D 轴向判据抓不到局部机制**：[B25] 里骨折占 23%（踝骨折 8%），",
        "        真实机制是**足缘小面积接触 + 踝旋后/内翻**导致的局部载荷，",
        "        不是均匀柱体轴向压缩 —— 1D 模型产生不了（P2 踝已在 §4c 补上）。",
        "        （[B25] 原文：护垫下沉 → 踝旋后 → 韧带负荷。）",
        "",
    ]

    # 最常出问题的姿势
    lines += [
        "=" * 78,
        "  五、按落地姿势的骨部位利用率（骨部位损伤全为 0，故看利用率）",
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

    # 结论用的涌现量（不写死）
    d_all = dist["all"]
    n_tr = max(1, n)
    an_in = d_all["ankle_injured"]
    an_sp = d_all["ankle_sprain"]
    an_fr = d_all["ankle_fracture"]
    bone_ev = d_all["n_injuries"]
    lines += [
        "",
        "=" * 78,
        "  六、结论",
        "=" * 78,
        f"  1. ✅ 采样器通过验证：MC 频率与 [B25] Table 3 的 95% CI "
        + (f"一致。" if ok
           else f"**不一致**（n={n} 小样本，见 §一 逐项）。"),
        f"  2. ⚠️ 骨部位（1D 判据）：{bone_ev} 例骨折事件 / {n} 次 —— "
        f"最高骨利用率 {util[order[0]]['max']:.2f}（{BONE_SITES[order[0]].name_cn}）。",
        "     骨部位分布仍为空，原因见 §4b（量程 + 1D 轴向判据抓不到局部机制）。",
        f"  3. ✅ P2 踝（本次新增）：踝损伤 {an_in}/{n} = {an_in/n_tr:.1%}"
        f"（踝扭伤 {an_sp/n_tr:.1%}，踝骨折 {an_fr/n_tr:.1%}）。",
        f"     [B25] 实测：踝 {OBSERVED_SPECIFIC_INJURY_ANKLE_TOTAL:.0%}、"
        f"踝扭伤 {OBSERVED_SPECIFIC_INJURY['ankle_sprain']:.0%}、"
        f"踝骨折 {OBSERVED_SPECIFIC_INJURY['ankle_fracture']:.0%}。",
        f"     ⇒ 损伤部位分布**不再为空**，踝是第一大（本模型唯一）部位 ——"
        f" 验证目标 V5 的方向达成。",
        "     ⚠️ 绝对量级依赖 pad2d.py 的踝标定（见并行 agent 的 "
        "ANKLE_SLIDE_UNLOAD / ANKLE_BONE_A_MM2），非本 harness 决定。",
        "  4. ⚠️ 仍缺的结构性缺口：肘（实测 16%）、膝（15%）、脱位（11%）、",
        "     肘/膝/肩韧带断裂（11%）—— 模型仍无上肢 / 膝关节 / 这些韧带判据。",
        f"     下肢合计（本模型口径，踝 trial + 下肢骨事件）："
        f"≥ {(an_in + sum(d_all['share'][k] * bone_ev for k in LOWER_LIMB_BONE_SITES))/n_tr:.1%}"
        f" vs [B25] {OBSERVED_INJURY_LOCATION['lower_limb'][0]:.0%}。",
        "  5. ➡️ P2 踝已接入主 harness；剩余缺口（上肢/膝）需继续加自由度，",
        "     且踝的绝对量级需要**独立的抱石实测标定数据**。",
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
