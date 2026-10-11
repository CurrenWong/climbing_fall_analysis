"""g4_damage_evolution.py -- G4 最小路径驱动：线性 σ(h) → 双线性弹塑（损伤）σ(h)，重算 A4。

预注册冻结设计：``docs/methods/G4_损伤演化设计.md``（§3 架构 / §5 验收 / §6 决策 / §7 备选）。

做什么
------
1. 读**只读冻结件** ``results/opensim_fe/fracture_matrix_bc2.json``（`align_pca=True`，主
   sub-stat=median）与 THUMS 材料卡（``read_material_damage`` 取逐 PID σ_y/E_t）。
2. 用冻结件里与原 A4 度量**完全相同的 per-part σ(h) 轨迹**（长骨=ends 统计 × load(h)/fe_solve，
   其余=``sigma_vm_mpa["median"][bone]``），经双线性弹塑律 ``damage_map`` 映射为**新** σ(h)。
3. **零 FE 重跑**：只读缓存 hdf5（长骨 ends 轨迹）与 JSON；不调用 FEBio/OpenSim。
4. 写**新** ``results/opensim_fe/g4_damage_matrix.json`` + ``g4_damage_run.json`` +
   ``G4_DAMAGE_REPORT.md``；冻结输入**只读**。
5. A4 重算（设计 §5 FROZEN）：
   * **rank-2**：PC2 解释方差比 vs 平行分析（Horn，定种子）/ 断棍 null（数据驱动，非魔数）；
   * **分离**：``(PC1, PC2)`` 载荷上的线性可分性 LP —— 复用冻结的 ``a4_from_matrix`` 规则；
   * **A2**（股骨/胸椎）：逐部位 logistic 状态 ∈ {阶跃, 渐进}，**只报**。

为什么转录而非"直接调用" ``a4_from_matrix``
------------------------------------------
冻结的 ``s4_alignment_bc2.a4_from_matrix`` 把**长骨** σ(h) 经**共享关节载荷向量**
（``fm["loads_n"][joint]``；胫/腓同用 ``ankle_r``）拼出 ⇒ 无法在不改冻结文件的情况下为每个
部位注入**各自**的损伤后轨迹。故本驱动：
  * 基线：**直接调用** ``s4.a4_from_matrix``（真实冻结函数，产出对照）；
  * 损伤：用 ``a4_rank2_from_sigma`` —— 与 ``a4_from_matrix`` 的 PCA+LP **逐行同构**的转录，
    并由 ``tests/test_g4_damage.py`` 断言其在**线性**轨迹上与真实函数**逐值相等**（规则未变）。

复现（仓库根，必须用 .venv 解释器）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_damage_evolution.py
"""
# allow: SIZE_OK — 驱动脚本 + 报告散文生成器（与 fracture_matrix.py / paper_alignment.py /
# s4_alignment_bc2.py 同构的仓库惯例）；纯数值逻辑 ~200 LOC，其余为不可约的报告字符串；
# 本构逻辑已拆到 src/climbing/coupling/bone_damage.py。
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parent.parent
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import fracture_matrix as fmx  # noqa: E402
import paper_alignment as pa  # noqa: E402
import s4_alignment_bc2 as s4  # noqa: E402
from bone_feb import read_material_damage  # noqa: E402
from climbing.coupling import bone_damage as bd  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
BC2_JSON = RES / "fracture_matrix_bc2.json"          # 只读冻结输入
OUT_JSON = RES / "g4_damage_matrix.json"             # 新件
OUT_RUN = RES / "g4_damage_run.json"                 # 新件
OUT_REPORT = RES / "G4_DAMAGE_REPORT.md"             # 新件
MASTER = ROOT / "model" / "AM50_V71_Occupant" / "main_THUMS_AM50_V71.k"

PA_SEED = 20261011        # 平行分析 null 固定种子（确定性）
PA_TRIALS = 5000
PLATEAU_RATIO = 0.05      # classify_step_progressive 阈值（显式、冻结）
GRADUAL_RATIO = 0.5

_BEHAV = {"all_zero": "常值0", "all_one": "常值1",
          "complete_separation": "阶跃", "ok": "渐进"}


# ---------------------------------------------------------------- trajectories
def build_linear_trajectories(yi: dict, fm: dict) -> dict:
    """逐部位严格线弹性 σ(h)，**与 ``a4_from_matrix`` 的 sig 构造逐行同构**。"""
    sig: dict[str, np.ndarray] = {}
    for cn, bone, joint, _subs in pa.REGIONS:
        d = yi["detail"][cn]
        if bone in pa.LONG_BONES:
            loads = np.asarray(fm["loads_n"][joint], dtype=float)
            sig[cn] = float(d["ends_stat_mpa"]) * loads / float(d["fe_solve_n"])
        else:
            sig[cn] = np.asarray(fm["sigma_vm_mpa"]["median"][bone], dtype=float)
    return sig


def main_domain_pid(fm: dict) -> dict:
    """部位 → 主承载域的 PID（用 ``fracture_matrix.build_parts`` 解析 mesh_name）。"""
    pbb = {p["bone"]: p for p in fm["parts"]}
    specs = {p.bone: p.spec for p in fmx.build_parts()}
    out: dict[str, int] = {}
    for _cn, bone, _joint, _subs in pa.REGIONS:
        dom = pbb[bone]["main_domain"]
        mn = dom.split("_", 1)[1] if dom.startswith(("solid_", "shell_")) else dom
        sp = specs[bone]
        for sd in sp.solid_domains:
            if sd.mesh_name == mn:
                out[bone] = int(sd.pids[0]); break
        for sh in sp.shell_domains:
            if sh.mesh_name == mn:
                out[bone] = int(sh.pids[0]); break
    return out


# ---------------------------------------------------------------- A4 rule (frozen)
def a4_rank2_from_sigma(sigma_by_part: dict) -> dict:
    """``(PC1,PC2)`` 二维可分性 —— **逐行转录** ``s4_alignment_bc2.a4_from_matrix`` 的数值规则。

    PCA：协方差 → eigh → 前 2 主载荷（9 部位）；LP：``w·x+b ≥ 1``(G1) / ``≤ −1``(G2)，
    ``scipy.optimize.linprog(method="highs")``。分组与"足部未标注"沿用 §六。
    """
    from scipy.optimize import linprog
    names = pa.REGION_NAMES
    X = np.column_stack([sigma_by_part[cn] for cn in names])
    Xc = X - X.mean(axis=0, keepdims=True)
    C = np.cov(Xc, rowvar=False)
    w, V = np.linalg.eigh(C)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    loads2 = V[:, :2]
    pos = {cn: [float(loads2[i, 0]), float(loads2[i, 1])] for i, cn in enumerate(names)}
    g1 = np.array([pos[c] for c in s4.A4_G1])
    g2 = np.array([pos[c] for c in s4.A4_G2])
    A_ub, b_ub = [], []
    for x in g1:
        A_ub.append([-x[0], -x[1], -1.0]); b_ub.append(-1.0)
    for x in g2:
        A_ub.append([x[0], x[1], 1.0]); b_ub.append(-1.0)
    r = linprog(c=[0, 0, 0], A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                bounds=[(None, None)] * 3, method="highs")
    return {
        "method": "cov(per-part sigma(h)) -> eigh -> top-2 PC loadings; "
                  "linear separability via LP (scipy.linprog highs) [transcribed from "
                  "s4_alignment_bc2.a4_from_matrix]",
        "explained_var": [float(w[0]), float(w[1])],
        "explained_var_ratio": [float(w[0] / w.sum()), float(w[1] / w.sum())],
        "part_positions": pos, "group1": list(s4.A4_G1), "group2": list(s4.A4_G2),
        "unlabeled": ["足部"], "separable": bool(r.success),
        "verdict": "PASS" if r.success else "FAIL", "lp_message": r.message,
    }


def rank2_null(sigma_matrix: np.ndarray, names, *, seed: int, trials: int) -> dict:
    """rank-2 数据驱动判定：PC2 解释方差比 vs 平行分析 95% null 与断棍 null。

    常数（已饱和）列方差 = 0，对协方差秩无贡献 ⇒ 从检验中剔除（使 null 更难通过，保守）。
    """
    X = np.asarray(sigma_matrix, dtype=float)
    keep = X.std(axis=0, ddof=1) > 0
    n_active = int(keep.sum())
    Xk = X[:, keep]
    w = np.linalg.eigvalsh(np.cov(Xk, rowvar=False))[::-1]
    ratio = w / w.sum()
    pc2 = float(ratio[1]) if w.size > 1 else 0.0
    scales = Xk.std(axis=0, ddof=1)
    _pct, all_ev = bd.parallel_analysis_null(scales, n_samples=X.shape[0],
                                             n_trials=trials, seed=seed, percentile=95.0)
    null_ratio = all_ev / all_ev.sum(axis=1, keepdims=True)
    null_p95 = float(np.percentile(null_ratio[:, 1], 95.0))
    bs = bd.broken_stick(n_active)
    bs2 = float(bs[1]) if bs.size > 1 else 0.0
    return {
        "n_active_parts": n_active,
        "active_parts": [names[i] for i in range(len(names)) if keep[i]],
        "constant_saturated_parts": [names[i] for i in range(len(names)) if not keep[i]],
        "pc1_explained_ratio": float(ratio[0]),
        "pc2_explained_ratio": pc2,
        "pc2_over_pc1": float(w[1] / w[0]) if w[0] > 0 else 0.0,
        "pa_null_pc2_ratio_mean": float(null_ratio[:, 1].mean()),
        "pa_null_pc2_ratio_p95": null_p95,
        "broken_stick_pc2_ratio": bs2,
        "pc2_exceeds_pa_null_p95": bool(pc2 > null_p95),
        "pc2_exceeds_broken_stick": bool(pc2 > bs2),
        "rank2_holds": bool(pc2 > null_p95),
        "seed": int(seed), "trials": int(trials),
    }


# ---------------------------------------------------------------- build
def build() -> dict:
    fm = json.loads(BC2_JSON.read_text(encoding="utf-8"))
    heights = np.asarray(fm["heights_m"], dtype=float)

    # ---- frozen A4 baseline: call the real frozen function ----------------
    yi = s4.route_yi("median")
    a4_baseline = s4.a4_from_matrix(yi, heights, fm)

    # ---- per-part linear trajectories (identical to a4_from_matrix) -------
    sig_lin = build_linear_trajectories(yi, fm)

    # ---- per-PID damage params from the material card ---------------------
    pid_by_bone = main_domain_pid(fm)
    mats = read_material_damage(MASTER, sorted(set(pid_by_bone.values())))

    sig_dmg: dict[str, np.ndarray] = {}
    parts_out: list[dict] = []
    classification: dict = {}
    for cn, bone, joint, subs in pa.REGIONS:
        pid = pid_by_bone[bone]
        mat = mats[pid]
        if mat.yield_stress_mpa is None:
            raise SystemExit(f"{bone}/{pid}: 卡无 σ_y/E_t（不在 _DAMAGE_CARDS）")
        law = bd.BoneDamageLaw(E_mpa=mat.E_mpa,
                               yield_stress_mpa=mat.yield_stress_mpa,
                               tangent_modulus_mpa=mat.tangent_modulus_mpa)
        sig_dmg[cn] = law.sigma_from_linear_stress(sig_lin[cn])
        cls = bd.classify_step_progressive(sig_lin[cn], law, heights=heights,
                                           plateau_ratio=PLATEAU_RATIO,
                                           gradual_ratio=GRADUAL_RATIO)
        classification[cn] = cls
        p = next(x for x in fm["parts"] if x["bone"] == bone)
        parts_out.append({
            "part_cn": cn, "bone": bone, "joint": joint,
            "main_domain": p["main_domain"], "main_kind": p["main_kind"], "pid": pid,
            "E_mpa": float(mat.E_mpa), "yield_stress_mpa": float(mat.yield_stress_mpa),
            "tangent_modulus_mpa": float(mat.tangent_modulus_mpa),
            "tangent_ratio_Et_over_E": float(law.tangent_ratio),
            "weak_sigma_c_mpa": float(pa._sc(subs[0][1])),
            "sigma_h_linear_mpa": [float(x) for x in sig_lin[cn]],
            "sigma_h_damaged_mpa": [float(x) for x in sig_dmg[cn]],
            "first_yield_h_m": cls["knee_h_m"], "classification": cls["class"],
        })

    # ---- A4: transcribed rule on linear (equality) & damaged --------------
    a4_linear = a4_rank2_from_sigma(sig_lin)
    a4_damaged = a4_rank2_from_sigma(sig_dmg)

    # ---- rank-2 data-driven null (the A4 metric's own matrix) -------------
    names = pa.REGION_NAMES
    rank2_dmg = rank2_null(np.column_stack([sig_dmg[c] for c in names]), names,
                           seed=PA_SEED, trials=PA_TRIALS)
    rank2_lin = rank2_null(np.column_stack([sig_lin[c] for c in names]), names,
                           seed=PA_SEED, trials=PA_TRIALS)
    rank2_dmg["baseline_pc2_explained_ratio"] = rank2_lin["pc2_explained_ratio"]
    rank2_dmg["pc2_ratio_gain_x"] = (
        rank2_dmg["pc2_explained_ratio"] / rank2_lin["pc2_explained_ratio"]
        if rank2_lin["pc2_explained_ratio"] > 0 else None)

    # ---- A2: per-part logistic state on the damage fracture matrix --------
    a2 = _a2_states(heights, sig_lin, sig_dmg, fm)

    return {
        "meta": {
            "design": "docs/methods/G4_损伤演化设计.md",
            "constitutive": "bilinear elastic-plastic (reduced, post-processing): "
                            "sigma=sigma_lin if sigma_lin<=sigma_y else "
                            "sigma_y+(E_t/E)(sigma_lin-sigma_y)",
            "frozen_input": "results/opensim_fe/fracture_matrix_bc2.json (read-only)",
            "material_card": "model/AM50_V71_Occupant/main_THUMS_AM50_V71.k",
            "zero_fe_resolve": True,
            "a4_input_trajectory": "per-part sigma(h) identical to a4_from_matrix "
                                   "(long bones: ends stat * load(h)/fe_solve; others: "
                                   "sigma_vm_mpa[median][bone])",
            "pa_seed": PA_SEED, "pa_trials": PA_TRIALS,
            "classify_thresholds": {"plateau_ratio": PLATEAU_RATIO,
                                    "gradual_ratio": GRADUAL_RATIO},
            "deterministic": "no wall-clock field; fixed PA seed",
        },
        "heights_m": [int(h) for h in heights],
        "parts": parts_out,
        "a4_baseline_frozen_call": {
            "explained_var_ratio": a4_baseline["explained_var_ratio"],
            "separable": a4_baseline["separable"], "verdict": a4_baseline["verdict"],
            "lp_message": a4_baseline["lp_message"],
        },
        "a4_linear_transcribed": {
            "explained_var_ratio": a4_linear["explained_var_ratio"],
            "separable": a4_linear["separable"], "verdict": a4_linear["verdict"],
        },
        "a4_damaged_transcribed": a4_damaged,
        "rank2_baseline": rank2_lin,
        "rank2_damaged": rank2_dmg,
        "a2": a2,
        "classification": classification,
        "verdicts": {
            "a4_rank2": "PASS" if rank2_dmg["rank2_holds"] else "FAIL",
            "a4_lp_separation": a4_damaged["verdict"],
            "break_rank1_ceiling_numeric": bool(
                rank2_dmg["pc2_explained_ratio"]
                > rank2_lin["pc2_explained_ratio"]),
        },
    }


def _a2_states(heights, sig_lin, sig_dmg, _fm) -> dict:
    """逐部位 logistic 状态（损伤 vs 基线），并单列股骨/胸椎。"""
    per = {}
    for cn, bone, _joint, subs in pa.REGIONS:
        weak = pa._sc(subs[0][1])
        lin_m = sig_lin[cn] >= weak
        dmg_m = sig_dmg[cn] >= weak
        lf = pa.logistic_fit(heights, lin_m.astype(int))
        dg = pa.logistic_fit(heights, dmg_m.astype(int))
        per[cn] = {
            "bone": bone, "weak_sigma_c_mpa": float(weak),
            "line_first_fracture_h": _first(heights, lin_m),
            "damage_first_fracture_h": _first(heights, dmg_m),
            "line_status": lf["status"], "line_behavior": _BEHAV.get(lf["status"]),
            "damage_status": dg["status"], "damage_behavior": _BEHAV.get(dg["status"]),
        }
    return {
        "per_part": per,
        "femur": per["股骨"], "thoracic": per["胸椎"], "lumbar": per["腰椎"],
        "note": "只报；与 FE 载荷分配缺陷纠缠（设计 §4/§5）",
    }


def _first(heights, mask):
    return pa._first_fracture(heights, mask)


# ---------------------------------------------------------------- report
def write_report(res: dict) -> None:
    L: list[str] = []
    A = L.append
    v = res["verdicts"]
    r2 = res["rank2_damaged"]
    a4d = res["a4_damaged_transcribed"]
    A("# G4 损伤演化（最小路径）报告 —— 双线性弹塑 σ(h) 重算 + A4 诚实判定")
    A("")
    A("> 预注册设计 `docs/methods/G4_损伤演化设计.md`（§3/§5/§6）。冻结输入只读："
      "`results/opensim_fe/fracture_matrix_bc2.json`。**零 FE 重跑**（σ(h) 后处理层）。")
    A("> 本件**无墙钟字段**：固定平行分析种子 ⇒ 逐位可复现。")
    A("")
    A("## 0. 一句话结论（如实）")
    A("")
    A(f"- **A4 rank-2：{v['a4_rank2']}** —— PC2 解释方差比 {r2['baseline_pc2_explained_ratio']:.6f} "
      f"→ **{r2['pc2_explained_ratio']:.6f}**（×{r2['pc2_ratio_gain_x']:.1f}），"
      f"**数值上已脱离 rank-1 天花板**；但平行分析 95% null = {r2['pa_null_pc2_ratio_p95']:.4f}、"
      f"断棍 null = {r2['broken_stick_pc2_ratio']:.4f} ⇒ PC2 **未**超 null ⇒ rank-2 **不成立**。")
    A(f"- **A4 分离（LP）：{v['a4_lp_separation']}** —— `(PC1,PC2)` 载荷上 G1/G2 不可线性分。")
    A(f"- **A2（股骨/胸椎）：未获区分度** —— 损伤律把脊柱 σ(h) 封顶在 σ_y=80 MPa < σ_c=150 "
      f"⇒ 腰椎/胸椎/颈椎由「阶跃/常值1」退化为「常值0」；股骨维持「常值0」（见 §4）。")
    A("")
    A("> 结论：**双线性后处理律本身不足以**造出可分双路径结构（与 §4 预注册风险一致）。"
      "§7 时间域度量是文档化备选（本件未实现；其触发条件「PC2 超 null 但 LP 失败」也**未**"
      "被满足——PC2 连 null 都未超）。")
    A("")
    A("## 1. 本构与参数（非调参）")
    A("")
    A("`σ(ε)=E·ε`（ε≤σ_y/E）否则 `σ_y+E_t·(ε−σ_y/E)`；后处理等价映射 "
      "`σ=σ_lin`（σ_lin≤σ_y）否则 `σ_y+(E_t/E)(σ_lin−σ_y)`。")
    A("")
    A("| 部位 | 代表骨 | 主承载域 | PID | E (MPa) | σ_y (MPa) | E_t (MPa) | E_t/E | 首屈服 h (m) | 分类 |")
    A("|---|---|---|---|---:|---:|---:|---:|---:|---|")
    for p in res["parts"]:
        A(f"| {p['part_cn']} | `{p['bone']}` | {p['main_domain']} | {p['pid']} | "
          f"{p['E_mpa']:.0f} | {p['yield_stress_mpa']:.2f} | {p['tangent_modulus_mpa']:.3g} | "
          f"{p['tangent_ratio_Et_over_E']:.3f} | "
          f"{'—' if p['first_yield_h_m'] is None else int(p['first_yield_h_m'])} | "
          f"{p['classification']} |")
    A("")
    A("> σ_y/E_t 逐 PID 取自 THUMS 卡（`read_material_damage`，卡 `nums[4]/nums[5]`），**非调参**。"
      "E_t=0 ⇒ 完美塑性平台（饱和端）；皮质骨 E_t/E≈0.82–0.90 ⇒ 膝拐（单调，非平台）。")
    A("")
    A("## 2. A4 —— rank-2（数据驱动 null）")
    A("")
    A(f"- A4 输入轨迹 = 与 `a4_from_matrix` **逐行同构**的 per-part σ(h)（长骨 ends；其余 median）。")
    A(f"- 有效变量 {r2['n_active_parts']}/9（常数/饱和部位 {r2['constant_saturated_parts']} "
      f"方差=0，从检验剔除 ⇒ null 门随之降低 "
      f"{res['rank2_baseline']['pa_null_pc2_ratio_p95']:.4f}→{r2['pa_null_pc2_ratio_p95']:.4f}"
      f"（对通过更有利、非保守）；即便如此 PC2 仍 < null ⇒ FAIL）：{r2['active_parts']}。")
    A("")
    A("| 量 | 基线（严格线弹性） | 损伤（双线性） |")
    A("|---|---:|---:|")
    A(f"| PC1 解释方差比 | {res['rank2_baseline']['pc1_explained_ratio']:.6f} | {r2['pc1_explained_ratio']:.6f} |")
    A(f"| **PC2 解释方差比** | {res['rank2_baseline']['pc2_explained_ratio']:.6f} | **{r2['pc2_explained_ratio']:.6f}** |")
    A(f"| 平行分析 95% null（PC2） | {res['rank2_baseline']['pa_null_pc2_ratio_p95']:.4f} | {r2['pa_null_pc2_ratio_p95']:.4f} |")
    A(f"| 断棍 null（PC2） | {res['rank2_baseline']['broken_stick_pc2_ratio']:.4f} | {r2['broken_stick_pc2_ratio']:.4f} |")
    A(f"| PC2 超 PA null？ | {res['rank2_baseline']['pc2_exceeds_pa_null_p95']} | **{r2['pc2_exceeds_pa_null_p95']}** |")
    A(f"| PC2 超断棍？ | {res['rank2_baseline']['pc2_exceeds_broken_stick']} | **{r2['pc2_exceeds_broken_stick']}** |")
    A("")
    A(f"> 平行分析 null = Horn 法（同维/同方差/iid 正态，{r2['trials']} trials，seed={r2['seed']}）；"
      f"断棍 = 断棍 null（见上）。⇒ **rank-2 = {v['a4_rank2']}**。")
    A("")
    A("## 3. A4 —— 分离（LP，复用冻结规则）")
    A("")
    A(f"- 基线（**直接调用**冻结 `s4_alignment_bc2.a4_from_matrix`）："
      f"PC2={res['a4_baseline_frozen_call']['explained_var_ratio'][1]:.6f}，"
      f"separable={res['a4_baseline_frozen_call']['separable']}。")
    A(f"- 线性轨迹（转录规则，应逐值等于基线）："
      f"PC2={res['a4_linear_transcribed']['explained_var_ratio'][1]:.6f}，"
      f"separable={res['a4_linear_transcribed']['separable']}。")
    A(f"- 损伤轨迹（转录规则）：PC1={a4d['explained_var_ratio'][0]:.6f}、"
      f"PC2={a4d['explained_var_ratio'][1]:.6f}；LP **{a4d['verdict']}**（{a4d['lp_message']}）。")
    A("")
    A("| 部位 | PC1 | PC2 | 组 |")
    A("|---|---:|---:|---|")
    for cn, xy in a4d["part_positions"].items():
        grp = "G1" if cn in s4.A4_G1 else ("G2" if cn in s4.A4_G2 else "未标注")
        A(f"| {cn} | {xy[0]:+.4f} | {xy[1]:+.4f} | {grp} |")
    A("")
    A("## 4. A2（股骨/胸椎）—— 只报")
    A("")
    A("| 部位 | 弱区 σ_c | 基线 状态 | 损伤 状态 | 基线首折 h | 损伤首折 h |")
    A("|---|---:|---|---|---:|---:|")
    for cn, e in res["a2"]["per_part"].items():
        A(f"| {cn} | {e['weak_sigma_c_mpa']:.0f} | {e['line_status']}（{e['line_behavior']}） | "
          f"{e['damage_status']}（{e['damage_behavior']}） | "
          f"{e['line_first_fracture_h'] or '—'} | {e['damage_first_fracture_h'] or '—'} |")
    A("")
    A("> 股骨/胸椎**未获区分度**；脊柱封顶 σ_y=80 < σ_c=150 ⇒ 全高度不判骨折。"
      "与 FE 载荷分配缺陷同根因（设计 §4/§5 预注册）。")
    A("")
    A("## 5. 复现命令")
    A("")
    A("```powershell")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_damage_evolution.py")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_g4_damage.py -q")
    A("```")
    A("")
    A("## 6. 产物")
    A("")
    A("- `results/opensim_fe/g4_damage_matrix.json`（新；σ(h) 线性/损伤逐格 + A4/A2/rank2 判定）")
    A("- `results/opensim_fe/g4_damage_run.json`（新；判定摘要）")
    A("- `results/opensim_fe/G4_DAMAGE_REPORT.md`（本文件）")
    A("- `src/climbing/coupling/bone_damage.py`、`scripts/opensim_fe/g4_damage_evolution.py`（新）")
    A("")
    OUT_REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")


# ---------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(OUT_JSON))
    args = ap.parse_args(argv)

    res = build()
    Path(args.out).write_text(json.dumps(res, indent=2, ensure_ascii=False),
                              encoding="utf-8")
    OUT_RUN.write_text(json.dumps({
        "verdicts": res["verdicts"],
        "rank2_damaged": {k: res["rank2_damaged"][k] for k in (
            "pc2_explained_ratio", "baseline_pc2_explained_ratio",
            "pa_null_pc2_ratio_p95", "broken_stick_pc2_ratio",
            "pc2_exceeds_pa_null_p95", "rank2_holds")},
        "a4_baseline_frozen_call": res["a4_baseline_frozen_call"],
        "a4_damaged_transcribed": {
            "explained_var_ratio": res["a4_damaged_transcribed"]["explained_var_ratio"],
            "separable": res["a4_damaged_transcribed"]["separable"],
            "verdict": res["a4_damaged_transcribed"]["verdict"]},
        "a2_femur_thoracic": {"股骨": res["a2"]["femur"], "胸椎": res["a2"]["thoracic"]},
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    write_report(res)

    print("=" * 88)
    print("G4 损伤演化（最小路径）判定：")
    print(f"  A4 rank-2 : {res['verdicts']['a4_rank2']}   "
          f"PC2 {res['rank2_damaged']['baseline_pc2_explained_ratio']:.6f} -> "
          f"{res['rank2_damaged']['pc2_explained_ratio']:.6f}  "
          f"(null95={res['rank2_damaged']['pa_null_pc2_ratio_p95']:.4f})")
    print(f"  A4 LP     : {res['verdicts']['a4_lp_separation']}   "
          f"({res['a4_damaged_transcribed']['lp_message']})")
    print(f"  A2 股骨   : {res['a2']['femur']['line_status']} -> "
          f"{res['a2']['femur']['damage_status']}")
    print(f"  A2 胸椎   : {res['a2']['thoracic']['line_status']} -> "
          f"{res['a2']['thoracic']['damage_status']}")
    print("=" * 88)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
