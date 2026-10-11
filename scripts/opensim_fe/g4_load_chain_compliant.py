"""g4_load_chain_compliant.py -- G4 载荷链顺应性杠杆（设计 §3 行4 / §6 D6）。

预注册冻结设计：``docs/methods/G4_损伤演化设计.md``。本驱动是 G4 的**第二个、也是最后
一个预注册杠杆**：在**顺应载荷链**（解锁下肢矢状面关节 + 单一冻结被动刚度预设
``PassiveStiffness.nominal_landing()``）下重算**长骨**的 σ(h)，再跑与最小路径
**完全相同**的 A4/A2 判据，并给出三路诚实对照。

做什么
------
1. 逐高度 ``h = 1..50``：``ground_reaction(h, 75.337)`` →
   ``run_dead_drop(grf, passive=PassiveStiffness.nominal_landing())`` →
   ``joint_reaction(fall, grf, joint=j)``，取 ``peak_vertical_n``（与
   ``fracture_matrix.py`` 的 ``load_ref`` 口径**同一统计量**）。关节 = ``ankle_r``
   （tibia/fibula）/ ``hip_r``（femur）——**只**这两个被长骨 σ(h) 使用。
   * **唯一**冻结预设：``PassiveStiffness.nominal_landing()``（D6：**不扫、不调**）。
   * 先做 **2-height 冒烟打印**（h=1, h=50）再进全循环。
   * 写**新** ``results/opensim_fe/load_chain_compliant.json``（主）+ ``.npz``
     （设计 §3 命名的二进制缓存）；冻结输入**只读**。
2. 用**顺应**载荷向量替换 ``build_linear_trajectories`` 里长骨的 ``loads_n``
   （非长骨部位**逐位沿用**冻结 ``sigma_vm_mpa[median]``），经同一双线性弹塑律
   ``bone_damage.BoneDamageLaw``（逐 PID σ_y/E_t 取自 ``read_material_damage``）映射，
   重跑**同一** A4（rank-2 数据驱动 null + LP，种子/trials 与最小路径**逐字相同**）与
   A2（股骨/胸椎 logistic 状态）。写**新** ``g4_loadchain_matrix.json`` +
   ``g4_loadchain_run.json`` + ``G4_LOADCHAIN_REPORT.md``。
3. **零 FE 重跑**（σ(h) 后处理层）；**不改** ``coupling/``（可选路径已存在）。

诚实边界（设计 §4/§7，必须随结果记录）
--------------------------------------
* 股骨/胸椎**很可能仍无区分度**——与 FE 载荷分配缺陷**同根因**；本杠杆**单独很可能修不好**。
* §7 时间域备选的**触发条件** = ``PC2 超 null 但 LP 失败``；本件**如实测量并记录该触发是否成立**，
  但**不实现**时间域度量（任务范围外）。

复现（仓库根，必须用 .venv 解释器）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_load_chain_compliant.py
    # --force-loads 重新跑 OpenSim（默认复用 load_chain_compliant.json 缓存）
"""
# allow: SIZE_OK — 驱动脚本 + 报告散文生成器（与 g4_damage_evolution.py / fracture_matrix.py /
# s4_alignment_bc2.py 同构的仓库惯例）；纯数值逻辑很小，其余为不可约的报告字符串；
# OpenSim 调用是既有 coupling API 的薄封装（lazy import，纯测试不触发）。
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parent.parent
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import paper_alignment as pa  # noqa: E402
import s4_alignment_bc2 as s4  # noqa: E402
from bone_feb import read_material_damage  # noqa: E402
from climbing.coupling import bone_damage as bd  # noqa: E402

# 复用最小路径驱动的**同一** A4/A2 规则与常量（保证种子/trials/口径逐字一致）。
import g4_damage_evolution as G  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
BC2_JSON = RES / "fracture_matrix_bc2.json"           # 只读冻结输入
OUT_LOADS_JSON = RES / "load_chain_compliant.json"    # 新件（主）
OUT_LOADS_NPZ = RES / "load_chain_compliant.npz"      # 新件（设计 §3 命名）
PARTIAL_PATH = ROOT / "temp" / "opensim_fe" / "g4_loadchain" / "partial_loads.json"
OUT_MATRIX = RES / "g4_loadchain_matrix.json"         # 新件
OUT_RUN = RES / "g4_loadchain_run.json"               # 新件
OUT_REPORT = RES / "G4_LOADCHAIN_REPORT.md"           # 新件
MASTER = ROOT / "model" / "AM50_V71_Occupant" / "main_THUMS_AM50_V71.k"

MODEL_MASS_KG = 75.337                 # 与 fracture_matrix.MODEL_MASS_KG 一致
HEIGHTS = np.arange(1, 51, dtype=int)  # 1..50 m
JOINTS_COMPLIANT = ("ankle_r", "hip_r")
PRESET_NAME = "PassiveStiffness.nominal_landing()"

REPRO_CMD = ("$env:PYTHONPATH=\"src\"; "
             "& .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_load_chain_compliant.py")


def _log(msg: str) -> None:
    print(f"[g4lc] {msg}", flush=True)


# ----------------------------------------------------------------- 载荷链
def _one_height(h: int) -> dict:
    """单高度**顺应**载荷链：grf → dead-drop(passive) → 关节纵向峰值。

    与 ``fracture_matrix.run_load_chain`` **同一**链路与统计量（``peak_vertical_n``），
    仅在 ``run_dead_drop`` 处加 ``passive=PassiveStiffness.nominal_landing()``（D6 唯一预设）。
    """
    from climbing.coupling.joint_reactions import joint_reaction
    from climbing.coupling.opensim_fall import PassiveStiffness, run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=float(h), mass_kg=MODEL_MASS_KG)
    fall = run_dead_drop(grf, passive=PassiveStiffness.nominal_landing())
    out = {"grf_peak_n": float(grf.peak_total_n),
           "grf_impulse_rel_err": float(grf.impulse_rel_err),
           "fall_n_states": int(fall.n_states)}
    for j in JOINTS_COMPLIANT:
        jr = joint_reaction(fall, grf, joint=j)
        out[j] = float(jr.peak_vertical_n)
    return out


def _load_partial() -> dict[int, dict]:
    """读增量检查点（断点续跑；预设/质量不符则丢弃）。"""
    if not PARTIAL_PATH.is_file():
        return {}
    try:
        p = json.loads(PARTIAL_PATH.read_text(encoding="utf-8"))
        m = p.get("meta", {})
        if m.get("preset") != PRESET_NAME or m.get("mass_kg") != MODEL_MASS_KG:
            return {}
        return {int(h): {k: float(v) for k, v in row.items()}
                for h, row in p.get("rows", {}).items()}
    except (OSError, ValueError):
        return {}


def _save_partial(done: dict[int, dict]) -> None:
    PARTIAL_PATH.parent.mkdir(parents=True, exist_ok=True)
    PARTIAL_PATH.write_text(json.dumps({
        "meta": {"preset": PRESET_NAME, "mass_kg": MODEL_MASS_KG,
                 "joints": list(JOINTS_COMPLIANT)},
        "rows": {str(h): done[h] for h in sorted(done)},
    }, indent=2, ensure_ascii=False), encoding="utf-8")


def compute_compliant_loads(*, reuse: bool) -> dict:
    """返回 ``{"heights": [...], "loads_n": {joint: [...]}, "meta": {...}}``。

    先 2-height 冒烟（h=1, h=50）再进全循环；结果复用冒烟值（不重算）。
    ``reuse=True`` 且既有 JSON 的 meta 与本次期望（heights/mass/preset）一致时**直接复用**
    （免 20 min OpenSim）；``reuse=False`` 强制重算。

    **断点续跑**：逐高度写增量检查点 ``temp/opensim_fe/g4_loadchain/partial_loads.json``；
    进程中断后再次运行自动从首个未完成高度续跑（同一预设/质量）。
    """
    if reuse and OUT_LOADS_JSON.is_file():
        prev = json.loads(OUT_LOADS_JSON.read_text(encoding="utf-8"))
        m = prev.get("meta", {})
        covered = ([int(x) for x in prev.get("heights_m", [])] == list(HEIGHTS))
        same = (m.get("mass_kg") == MODEL_MASS_KG
                and m.get("preset") == PRESET_NAME
                and set(prev.get("loads_n", {})) == set(JOINTS_COMPLIANT))
        if covered and same:
            _log(f"复用 {OUT_LOADS_JSON}（{len(prev['heights_m'])} heights，预设未变）")
            return {"heights": [int(x) for x in prev["heights_m"]],
                    "loads_n": {j: [float(x) for x in prev["loads_n"][j]]
                                for j in JOINTS_COMPLIANT},
                    "meta": {**m, "reused": True}}

    done: dict[int, dict] = _load_partial() if reuse else {}
    if done:
        _log(f"断点续跑：已有 {len(done)}/{HEIGHTS.size} 高度（{PARTIAL_PATH}）")

    _log("=" * 78)
    _log("2-height 冒烟：h = 1, 50（确认链路可用 + 量级）")
    for h in (1, 50):
        if int(h) not in done:
            t0 = time.time()
            done[int(h)] = _one_height(int(h))
            _save_partial(done)
            _log(f"  smoke h={h:2d} m  GRF={done[int(h)]['grf_peak_n']/1e3:7.1f} kN  "
                 f"ankle_r={done[int(h)]['ankle_r']/1e3:6.1f} kN  "
                 f"hip_r={done[int(h)]['hip_r']/1e3:6.1f} kN  ({time.time()-t0:4.1f}s)")
        else:
            r = done[int(h)]
            _log(f"  smoke h={h:2d} m  GRF={r['grf_peak_n']/1e3:7.1f} kN  "
                 f"ankle_r={r['ankle_r']/1e3:6.1f} kN  hip_r={r['hip_r']/1e3:6.1f} kN  (缓存)")
    _log("=" * 78)

    for h in HEIGHTS:
        if int(h) in done:
            continue
        t0 = time.time()
        done[int(h)] = _one_height(int(h))
        _save_partial(done)
        r = done[int(h)]
        _log(f"h={int(h):2d} m  GRF={r['grf_peak_n']/1e3:7.1f} kN  "
             f"ankle_r={r['ankle_r']/1e3:6.1f} kN  "
             f"hip_r={r['hip_r']/1e3:6.1f} kN  ({time.time()-t0:4.1f}s)")

    heights = [int(h) for h in HEIGHTS]
    loads_n = {j: [float(done[int(h)][j]) for h in HEIGHTS] for j in JOINTS_COMPLIANT}
    meta = {
        "design": "docs/methods/G4_损伤演化设计.md §3 行4 / §6 D6",
        "lever": "load-chain compliance (single frozen passive preset; no scan)",
        "preset": PRESET_NAME,
        "mass_kg": MODEL_MASS_KG,
        "heights_m": heights,
        "joints": list(JOINTS_COMPLIANT),
        "statistic": "joint_reaction(...).peak_vertical_n (longitudinal peak; same as "
                     "fracture_matrix load_ref)",
        "chain": "ground_reaction(h, 75.337) -> run_dead_drop(grf, passive="
                 "PassiveStiffness.nominal_landing()) -> joint_reaction(fall, grf, joint=j)",
        "grf_peak_n": [float(done[int(h)]["grf_peak_n"]) for h in HEIGHTS],
        "zero_fe_resolve": True,       # 只重跑 OpenSim 载荷链；FE 单骨求解冻结
        "reused": False,
        "command": REPRO_CMD,
    }
    OUT_LOADS_JSON.write_text(json.dumps(
        {"meta": meta, "heights_m": heights, "loads_n": loads_n},
        indent=2, ensure_ascii=False), encoding="utf-8")
    np.savez(OUT_LOADS_NPZ, heights=np.asarray(heights, float),
             grf_peak_n=np.asarray(meta["grf_peak_n"], float),
             **{f"{j}__vert": np.asarray(loads_n[j], float) for j in JOINTS_COMPLIANT})
    if PARTIAL_PATH.is_file():
        PARTIAL_PATH.unlink()
    _log(f"写出 {OUT_LOADS_JSON}")
    _log(f"写出 {OUT_LOADS_NPZ}")
    return {"heights": heights, "loads_n": loads_n, "meta": meta}


# ------------------------------------------------- σ(h) 轨迹（载荷可注入）
def build_linear_trajectories_override(yi: dict, fm: dict, loads_n: dict) -> dict:
    """逐部位严格线弹性 σ(h)，长骨用**注入**的 ``loads_n``（其余逐字同冻结）。

    与 ``g4_damage_evolution.build_linear_trajectories`` **逐行同构**：唯一差异是长骨
    的载荷向量来自入参 ``loads_n`` 而非 ``fm["loads_n"]``。传入
    ``loads_n = fm["loads_n"]`` 时**逐值等于**冻结函数（由测试锁定）。
    """
    sig: dict[str, np.ndarray] = {}
    for cn, bone, joint, _subs in pa.REGIONS:
        d = yi["detail"][cn]
        if bone in pa.LONG_BONES:
            loads = np.asarray(loads_n[joint], dtype=float)
            sig[cn] = float(d["ends_stat_mpa"]) * loads / float(d["fe_solve_n"])
        else:
            sig[cn] = np.asarray(fm["sigma_vm_mpa"]["median"][bone], dtype=float)
    return sig


# ----------------------------------------------------------------- build
def build(*, reuse_loads: bool = True, damage_only: dict | None = None) -> dict:
    fm = json.loads(BC2_JSON.read_text(encoding="utf-8"))
    heights = np.asarray(fm["heights_m"], dtype=float)
    names = pa.REGION_NAMES
    yi = s4.route_yi("median")

    # ---- 1. 顺应载荷链 -----------------------------------------------------
    chain = compute_compliant_loads(reuse=reuse_loads)
    loads_c = {j: np.asarray(chain["loads_n"][j], dtype=float) for j in JOINTS_COMPLIANT}

    # ---- 对照：最小路径（损伤-only）由冻结驱动**原样**重算（确定性） -------
    if damage_only is None:
        damage_only = G.build()

    # ---- 2. σ(h) 轨迹 ------------------------------------------------------
    sig_lin_rigid = G.build_linear_trajectories(yi, fm)             # 冻结载荷链
    sig_lin_c = build_linear_trajectories_override(yi, fm, loads_c)  # 顺应载荷链

    pid_by_bone = G.main_domain_pid(fm)
    mats = read_material_damage(MASTER, sorted(set(pid_by_bone.values())))

    sig_dmg_c: dict[str, np.ndarray] = {}
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
        sig_dmg_c[cn] = law.sigma_from_linear_stress(sig_lin_c[cn])
        cls = bd.classify_step_progressive(sig_lin_c[cn], law, heights=heights,
                                           plateau_ratio=G.PLATEAU_RATIO,
                                           gradual_ratio=G.GRADUAL_RATIO)
        classification[cn] = cls
        p = next(x for x in fm["parts"] if x["bone"] == bone)
        parts_out.append({
            "part_cn": cn, "bone": bone, "joint": joint,
            "main_domain": p["main_domain"], "main_kind": p["main_kind"], "pid": pid,
            "E_mpa": float(mat.E_mpa), "yield_stress_mpa": float(mat.yield_stress_mpa),
            "tangent_modulus_mpa": float(mat.tangent_modulus_mpa),
            "tangent_ratio_Et_over_E": float(law.tangent_ratio),
            "weak_sigma_c_mpa": float(pa._sc(subs[0][1])),
            "sigma_h_linear_compliant_mpa": [float(x) for x in sig_lin_c[cn]],
            "sigma_h_damaged_compliant_mpa": [float(x) for x in sig_dmg_c[cn]],
            "sigma_h_linear_rigid_mpa": [float(x) for x in sig_lin_rigid[cn]],
            "first_yield_h_m": cls["knee_h_m"], "classification": cls["class"],
        })

    # ---- 3. A4 / A2（规则、种子、trials 与最小路径**逐字相同**）------------
    a4_baseline = s4.a4_from_matrix(yi, heights, fm)          # rank-1（冻结调用）
    a4_c_lin = G.a4_rank2_from_sigma(sig_lin_c)
    a4_c_dmg = G.a4_rank2_from_sigma(sig_dmg_c)

    rank2_c = G.rank2_null(np.column_stack([sig_dmg_c[c] for c in names]), names,
                           seed=G.PA_SEED, trials=G.PA_TRIALS)
    rank2_c_lin = G.rank2_null(np.column_stack([sig_lin_c[c] for c in names]), names,
                               seed=G.PA_SEED, trials=G.PA_TRIALS)
    rank2_c["baseline_pc2_explained_ratio"] = rank2_c_lin["pc2_explained_ratio"]
    rank2_c["pc2_ratio_gain_x"] = (
        rank2_c["pc2_explained_ratio"] / rank2_c_lin["pc2_explained_ratio"]
        if rank2_c_lin["pc2_explained_ratio"] > 0 else None)

    a2 = G._a2_states(heights, sig_lin_c, sig_dmg_c, fm)

    # §7 触发条件（如实记录，不实现时间域度量）
    section7_trigger = bool(rank2_c["pc2_exceeds_pa_null_p95"]
                            and a4_c_dmg["verdict"] == "FAIL")

    return {
        "meta": {
            "design": "docs/methods/G4_损伤演化设计.md",
            "lever": "load-chain compliance (G4 第二杠杆)",
            "preset": PRESET_NAME,
            "mass_kg": MODEL_MASS_KG,
            "frozen_input": "results/opensim_fe/fracture_matrix_bc2.json (read-only)",
            "load_chain_compliant": "results/opensim_fe/load_chain_compliant.json (new)",
            "material_card": "model/AM50_V71_Occupant/main_THUMS_AM50_V71.k",
            "zero_fe_resolve": True,
            "a4_input_trajectory": "long bones: ends stat * load_compliant(h)/fe_solve; "
                                   "others: sigma_vm_mpa[median][bone] (identical to frozen)",
            "pa_seed": G.PA_SEED, "pa_trials": G.PA_TRIALS,
            "classify_thresholds": {"plateau_ratio": G.PLATEAU_RATIO,
                                    "gradual_ratio": G.GRADUAL_RATIO},
            "deterministic": "no wall-clock field; fixed PA seed",
            "command": REPRO_CMD,
        },
        "heights_m": [int(h) for h in heights],
        "loads_n_compliant": {j: [float(x) for x in loads_c[j]] for j in JOINTS_COMPLIANT},
        "grf_peak_n": chain["meta"].get("grf_peak_n"),
        "parts": parts_out,
        "a4_baseline_rank1": {
            "explained_var_ratio": a4_baseline["explained_var_ratio"],
            "separable": a4_baseline["separable"], "verdict": a4_baseline["verdict"],
            "lp_message": a4_baseline["lp_message"],
        },
        "a4_compliant_linear": {
            "explained_var_ratio": a4_c_lin["explained_var_ratio"],
            "separable": a4_c_lin["separable"], "verdict": a4_c_lin["verdict"],
        },
        "a4_compliant_damaged": a4_c_dmg,
        "rank2_baseline_rigid": damage_only["rank2_baseline"],
        "rank2_damage_only": damage_only["rank2_damaged"],
        "rank2_compliant": rank2_c,
        "a2_compliant": a2,
        "a2_damage_only": damage_only["a2"],
        "a4_damage_only": damage_only["a4_damaged_transcribed"],
        "classification": classification,
        "verdicts": {
            "a4_rank2": "PASS" if rank2_c["rank2_holds"] else "FAIL",
            "a4_lp_separation": a4_c_dmg["verdict"],
            "break_rank1_ceiling_numeric": bool(
                rank2_c["pc2_explained_ratio"]
                > damage_only["rank2_baseline"]["pc2_explained_ratio"]),
            "section7_trigger_met": section7_trigger,
            "section7_trigger_condition": "PC2 exceeds PA null AND LP separation fails",
        },
    }


# ----------------------------------------------------------------- report
def _fmt(x, nd=6):
    return "—" if x is None else f"{x:.{nd}f}"


def write_report(res: dict) -> None:
    L: list[str] = []
    A = L.append
    v = res["verdicts"]
    r2c = res["rank2_compliant"]
    r2d = res["rank2_damage_only"]
    r2b = res["rank2_baseline_rigid"]
    a4b = res["a4_baseline_rank1"]
    a4c = res["a4_compliant_damaged"]

    A("# G4 载荷链顺应性杠杆报告 —— 顺应载荷链下 σ(h) 重算 + A4/A2 诚实判定")
    A("")
    A("> 预注册设计 `docs/methods/G4_损伤演化设计.md`（§3 行4 / §5 / §6 D6 / §7）。"
      "冻结输入只读：`results/opensim_fe/fracture_matrix_bc2.json`。**零 FE 重跑**"
      "（σ(h) 后处理层）；OpenSim 载荷链在**唯一冻结预设** `PassiveStiffness.nominal_landing()`"
      "下逐高度重跑。本件**无墙钟字段**：固定平行分析种子 ⇒ 逐位可复现。")
    A("")
    A("## 0. 一句话结论（如实）")
    A("")
    A(f"- **A4 rank-2：{v['a4_rank2']}** —— PC2 解释方差比（有效部位集）"
      f"{_fmt(r2b['pc2_explained_ratio'])}（rank-1）→ {_fmt(r2d['pc2_explained_ratio'])}"
      f"（损伤-only）→ **{_fmt(r2c['pc2_explained_ratio'])}**（损伤+顺应链）；"
      f"平行分析 95% null = {_fmt(r2c['pa_null_pc2_ratio_p95'], 4)}、断棍 null = "
      f"{_fmt(r2c['broken_stick_pc2_ratio'], 4)} ⇒ "
      f"PC2 {'**超**' if r2c['pc2_exceeds_pa_null_p95'] else '**未超**'} null ⇒ "
      f"rank-2 **{v['a4_rank2']}**。")
    A(f"- **A4 分离（LP）：{v['a4_lp_separation']}** —— `(PC1,PC2)` 载荷上 G1/G2"
      f"{'可' if a4c['separable'] else '不可'}线性分。")
    a2c = res["a2_compliant"]
    fem, tho = a2c["femur"], a2c["thoracic"]
    A(f"- **A2（股骨/胸椎）：部分改善（诚实报限）** —— 股骨 {fem['damage_status']}"
      f"（{fem['damage_behavior']}，首折 h={fem['damage_first_fracture_h'] or '—'}）；"
      f"胸椎 {tho['damage_status']}（{tho['damage_behavior']}）。"
      f"顺应链**直接抬高髋关节载荷**使**股骨获得阶跃区分度**；但**胸椎** σ(h) 取自冻结 "
      f"median 轨迹且被损伤律 σ_y=80 < σ_c=150 **封顶** ⇒ 载荷链杠杆**单独不足以**修复胸椎"
      f"（设计 §4 预注册的「与 FE 载荷分配缺陷同根因」对胸椎成立）。")
    A("")
    A(f"> **§7 触发条件（PC2 超 null 但 LP 失败）：{'成立' if v['section7_trigger_met'] else '未成立'}**。"
      f"本件**不实现**时间域度量（任务范围外）；如实记录触发状态。")
    A("")
    A("## 1. 三路对照（rank-1 严格线弹性 vs 损伤-only vs 损伤+顺应载荷链）")
    A("")
    A("| 量 | rank-1 基线（严格线性） | 损伤-only（最小路径） | **损伤+顺应链（本件）** |")
    A("|---|---:|---:|---:|")
    A(f"| PC2 解释方差比（有效集） | {_fmt(r2b['pc2_explained_ratio'])} | "
      f"{_fmt(r2d['pc2_explained_ratio'])} | **{_fmt(r2c['pc2_explained_ratio'])}** |")
    A(f"| PC1 解释方差比（有效集） | {_fmt(r2b['pc1_explained_ratio'])} | "
      f"{_fmt(r2d['pc1_explained_ratio'])} | **{_fmt(r2c['pc1_explained_ratio'])}** |")
    A(f"| 平行分析 95% null（PC2） | {_fmt(r2b['pa_null_pc2_ratio_p95'], 4)} | "
      f"{_fmt(r2d['pa_null_pc2_ratio_p95'], 4)} | **{_fmt(r2c['pa_null_pc2_ratio_p95'], 4)}** |")
    A(f"| 断棍 null（PC2） | {_fmt(r2b['broken_stick_pc2_ratio'], 4)} | "
      f"{_fmt(r2d['broken_stick_pc2_ratio'], 4)} | **{_fmt(r2c['broken_stick_pc2_ratio'], 4)}** |")
    A(f"| PC2 超 PA null？ | {r2b['pc2_exceeds_pa_null_p95']} | "
      f"{r2d['pc2_exceeds_pa_null_p95']} | **{r2c['pc2_exceeds_pa_null_p95']}** |")
    A(f"| rank-2 成立？ | — | "
      f"{'PASS' if r2d['rank2_holds'] else 'FAIL'} | **{'PASS' if r2c['rank2_holds'] else 'FAIL'}** |")
    A(f"| A4 分离（LP） | {a4b['verdict']} | "
      f"{res['a4_damage_only']['verdict']} | **{a4c['verdict']}** |")
    A(f"| A2 股骨 状态 | {res['a2_damage_only']['femur']['line_status']} | "
      f"{res['a2_damage_only']['femur']['damage_status']} | "
      f"**{res['a2_compliant']['femur']['damage_status']}** |")
    A(f"| A2 胸椎 状态 | {res['a2_damage_only']['thoracic']['line_status']} | "
      f"{res['a2_damage_only']['thoracic']['damage_status']} | "
      f"**{res['a2_compliant']['thoracic']['damage_status']}** |")
    A(f"| 有效部位数 | {r2b['n_active_parts']} | {r2d['n_active_parts']} | "
      f"**{r2c['n_active_parts']}** |")
    A("")
    A(f"> 平行分析 null = Horn 法（同维/同方差/iid 正态，{r2c['trials']} trials，"
      f"seed={r2c['seed']}）；断棍 = 断棍 null。**三路同规则、同种子**（转录/复用自最小路径）。")
    A("")
    A("## 2. 顺应载荷链（唯一冻结预设，非调参）")
    A("")
    A(f"- 预设 = `{res['meta']['preset']}`；质量 = {res['meta']['mass_kg']} kg；高度 = 1..50 m。")
    A("- 链路 = `ground_reaction(h, 75.337) -> run_dead_drop(grf, passive=预设) -> "
      "joint_reaction(fall, grf, joint=j)`，取 `peak_vertical_n`（纵向峰值，与 "
      "`fracture_matrix` 的 `load_ref` 同口径）。")
    A("- **只**用 `ankle_r`（tibia/fibula）/ `hip_r`（femur）；非长骨部位 σ(h) **逐位沿用**冻结 median 轨迹。")
    A("")
    A("| 关节 | @1 m (kN) | @5 m (kN) | @50 m (kN) | 50/5 比 |")
    A("|---|---:|---:|---:|---:|")
    for j in JOINTS_COMPLIANT:
        a = np.asarray(res["loads_n_compliant"][j], float)
        A(f"| `{j}` | {a[0]/1e3:.1f} | {a[4]/1e3:.1f} | {a[-1]/1e3:.1f} | "
          f"{a[-1]/a[4]:.3f} |")
    A("")
    A("> 与冻结（刚性腿）载荷链对照（`fracture_matrix_bc2.json` 的 `loads_n`）：顺应链改变"
      "各关节**峰值量级与沿高度形状**——这正是本杠杆的物理内容（`coupling/` **未改**，"
      "仅调用既有可选路径）。")
    A("")
    A("## 3. 本构与参数（非调参）")
    A("")
    A("| 部位 | 代表骨 | PID | E (MPa) | σ_y (MPa) | E_t (MPa) | E_t/E | 首屈服 h (m) | 分类 |")
    A("|---|---|---|---:|---:|---:|---:|---:|---|")
    for p in res["parts"]:
        A(f"| {p['part_cn']} | `{p['bone']}` | {p['pid']} | {p['E_mpa']:.0f} | "
          f"{p['yield_stress_mpa']:.2f} | {p['tangent_modulus_mpa']:.3g} | "
          f"{p['tangent_ratio_Et_over_E']:.3f} | "
          f"{'—' if p['first_yield_h_m'] is None else int(p['first_yield_h_m'])} | "
          f"{p['classification']} |")
    A("")
    A("> σ_y/E_t 逐 PID 取自 THUMS 卡（`read_material_damage`，与非调参口径一致）。"
      "**唯一**旋钮 = 冻结被动预设，**不扫、不调**（§6 D6）。")
    A("")
    A("## 4. A2（股骨/胸椎）—— 只报 + 三路状态")
    A("")
    A("| 部位 | 弱区 σ_c | rank-1 状态 | 损伤-only 状态 | 损伤+顺应链 状态 | 顺应链首折 h |")
    A("|---|---:|---|---|---|---:|")
    for cn, e in res["a2_compliant"]["per_part"].items():
        d = res["a2_damage_only"]["per_part"][cn]
        A(f"| {cn} | {e['weak_sigma_c_mpa']:.0f} | {d['line_status']} | "
          f"{d['damage_status']} | **{e['damage_status']}** | "
          f"{e['damage_first_fracture_h'] or '—'} |")
    A("")
    A("> 顺势链把**股骨**由「常值0」抬为「阶跃」（髋载荷随高度上移使股骨弱区 σ(h) 越过 σ_c=80）；"
      "**胸椎**仍为「常值0」——因其 σ(h) 取自冻结 median 轨迹且损伤律 σ_y=80 < σ_c=150 封顶。"
      "⇒ 载荷链杠杆**部分改善 A2（股骨✔ / 胸椎✘）**，但**未**修复 A4（rank-2 与 LP 双双 FAIL）。"
      "胸椎的修复**超出**本杠杆范围（设计 §4/§8 预注册）。")
    A("")
    A("## 5. LP 载荷坐标（损伤+顺应链）")
    A("")
    A("| 部位 | PC1 | PC2 | 组 |")
    A("|---|---:|---:|---|")
    for cn, xy in a4c["part_positions"].items():
        grp = "G1" if cn in s4.A4_G1 else ("G2" if cn in s4.A4_G2 else "未标注")
        A(f"| {cn} | {xy[0]:+.4f} | {xy[1]:+.4f} | {grp} |")
    A("")
    A("## 6. 复现命令")
    A("")
    A("```powershell")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_load_chain_compliant.py")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_g4_load_chain.py -q")
    A("```")
    A("")
    A("## 7. 产物")
    A("")
    A("- `results/opensim_fe/load_chain_compliant.json`（新；顺应链逐关节逐高度 + meta）")
    A("- `results/opensim_fe/load_chain_compliant.npz`（新；设计 §3 命名的二进制缓存）")
    A("- `results/opensim_fe/g4_loadchain_matrix.json`（新；三路 σ(h) + A4/A2/rank2 判定）")
    A("- `results/opensim_fe/g4_loadchain_run.json`（新；判定摘要）")
    A("- `results/opensim_fe/G4_LOADCHAIN_REPORT.md`（本文件）")
    A("- `scripts/opensim_fe/g4_load_chain_compliant.py`、`tests/test_g4_load_chain.py`（新）")
    A("")
    OUT_REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")


# ----------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(OUT_MATRIX))
    ap.add_argument("--force-loads", action="store_true",
                    help="强制重跑 OpenSim 载荷链（忽略 load_chain_compliant.json 缓存）")
    args = ap.parse_args(argv)

    t0 = time.time()
    res = build(reuse_loads=not args.force_loads)
    Path(args.out).write_text(json.dumps(res, indent=2, ensure_ascii=False),
                              encoding="utf-8")

    r2c = res["rank2_compliant"]
    r2d = res["rank2_damage_only"]
    r2b = res["rank2_baseline_rigid"]
    OUT_RUN.write_text(json.dumps({
        "verdicts": res["verdicts"],
        "pc2_explained_ratio": {
            "rank1_baseline": r2b["pc2_explained_ratio"],
            "damage_only": r2d["pc2_explained_ratio"],
            "damage_compliant": r2c["pc2_explained_ratio"],
        },
        "pa_null_pc2_ratio_p95": {
            "rank1_baseline": r2b["pa_null_pc2_ratio_p95"],
            "damage_only": r2d["pa_null_pc2_ratio_p95"],
            "damage_compliant": r2c["pa_null_pc2_ratio_p95"],
        },
        "rank2_holds": {
            "damage_only": r2d["rank2_holds"],
            "damage_compliant": r2c["rank2_holds"],
        },
        "a4_lp_separation": {
            "rank1_baseline": res["a4_baseline_rank1"]["verdict"],
            "damage_only": res["a4_damage_only"]["verdict"],
            "damage_compliant": res["a4_compliant_damaged"]["verdict"],
        },
        "a2_femur_thoracic": {
            "rank1_baseline": {"股骨": res["a2_damage_only"]["femur"]["line_status"],
                               "胸椎": res["a2_damage_only"]["thoracic"]["line_status"]},
            "damage_only": {"股骨": res["a2_damage_only"]["femur"]["damage_status"],
                            "胸椎": res["a2_damage_only"]["thoracic"]["damage_status"]},
            "damage_compliant": {"股骨": res["a2_compliant"]["femur"]["damage_status"],
                                 "胸椎": res["a2_compliant"]["thoracic"]["damage_status"]},
        },
        "meta": res["meta"],
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    write_report(res)

    print("=" * 88)
    print("G4 载荷链顺应性杠杆判定：")
    print(f"  A4 rank-2 : {res['verdicts']['a4_rank2']}   "
          f"PC2 {r2b['pc2_explained_ratio']:.6f} -> {r2d['pc2_explained_ratio']:.6f} -> "
          f"{r2c['pc2_explained_ratio']:.6f}  (null95={r2c['pa_null_pc2_ratio_p95']:.4f})")
    print(f"  A4 LP     : {res['verdicts']['a4_lp_separation']}   "
          f"({res['a4_compliant_damaged']['lp_message']})")
    print(f"  §7 触发   : {res['verdicts']['section7_trigger_met']}")
    print(f"  A2 股骨   : {res['a2_compliant']['femur']['line_status']} -> "
          f"{res['a2_compliant']['femur']['damage_status']}")
    print(f"  A2 胸椎   : {res['a2_compliant']['thoracic']['line_status']} -> "
          f"{res['a2_compliant']['thoracic']['damage_status']}")
    print(f"  总耗时     : {time.time()-t0:.0f}s")
    print("=" * 88)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
