# -*- coding: utf-8 -*-
"""P2 · 端部 BC（P5 式解剖关节面接触）+ 中轴载荷分配 —— 能力门探针。

预注册：`docs/methods/P2_端部BC与载荷分配_预注册.md`（2026-10-15 冻结，跑前）。
侦察：`results/opensim_fe/P2_BC_SCOPE_INVENTORY.md`（提交 e603d7e）。

**为什么是新模块而不是改 `p5_head_contact.py`**：后者是 P5.0 的已收口记录，
硬编码 `femur_r`。本模块**只复用、零改动**，并把它泛化到任意 `bone`。
泛化的结构依据：`p4_pb_probe.head_joint_nodes` 的约定（+Z 端 = cap 接触承载端，
−Z 端 = 分布载荷端，力向 +Z 压缩）**与解剖无关** —— 股骨头球冠和椎体终板
是同一个拓扑角色，故公式逐字成立。

**冻结配方**（承预注册 §2，逐字继承 P5.0）：
  GAP_MM −0.3 / TOP_FRAC 0.30 / CAP_AXDIV_MIN 0.5 / fric 0.0
  lateral_support [(0.05,0.95)] / 板锁 RxRy / cap 端 3 节点横向锚点
  板压入 PLATE_DISP_MM 0.5

**冻结 penalty 公式**（预注册 §4）：penalty = 50 × F_ref / 57941.3
  锚点 = P5.0 实测 PASS（penalty=50 @ 股骨 22923−57941 N）。
  一次给出全部部位，**零自由度可调**；不收敛 = 结果，不是待调参数。

**红线**：不得调参凑 A1/A2；冻结输入只读；产出新件。
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import bone_feb as bf  # noqa: E402
import g4_fe_p1b as p1b  # noqa: E402
import p4_pb_probe  # noqa: E402
import p5_head_contact as p5  # noqa: E402
import paper_alignment as pa  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
BASE = ROOT / "temp" / "opensim_fe" / "p2_head_contact"

# ---- 冻结常量（承预注册 §2） ----------------------------------------------
PLATE_DISP_MM = p5.PLATE_DISP_MM      # 0.5
GAP_MM = p5.GAP_MM                    # -0.3
TOP_FRAC = p5.TOP_FRAC                # 0.30
LOAD_FRAC = 0.20                      # 载荷端带份额（承股骨 HEAD_FRAC）
FRIC = 0.0
LATERAL = [(0.05, 0.95)]
N_ANCHOR = 3
TIME_STEPS = 10

# ---- 冻结 penalty 公式（预注册 §4） ----------------------------------------
PENALTY_REF_LOAD_N = 57941.3          # P5.0 实测 PASS 的最大载荷
PENALTY_REF_VALUE = 50.0              # P5.0 实测 PASS 的 penalty


def penalty_for(f_ref: float) -> float:
    """预注册 §4 公式。唯一自变量是载荷（罚刚度 ∝ penalty，接触承载 ⟺ 罚够刚）。"""
    return PENALTY_REF_VALUE * float(f_ref) / PENALTY_REF_LOAD_N


# ---- 缺陷 B 载荷分配（预注册 §3；库内 risk_1d_loadshare 逐椎质量） -----------
LOAD_SHARE = {
    "L3": 1.0000, "T6": 0.6104, "C5": 0.1840, "parietal_r": 0.1442,
    "tibia_r": 0.7759, "fibula_r": 0.2241,
    "femur_r": 1.0, "calcaneus_r": 1.0, "R_HIPBONE": 1.0,
}

# ---- 预注册默认子集（#3 颅骨排除 / #4 腰 L3 排除，均可翻转） ----------------
SUBSET = ["T6", "C5", "femur_r", "tibia_r", "fibula_r"]

HEIGHTS = [10, 20, 30, 50]

# 能力门阈值（预注册 §5.1，跑前冻结）
BEAR_LO, BEAR_HI = 0.80, 1.20
R2_MIN = 0.20


# ---------------------------------------------------------------------------
# 泛化的两端上下文（结构同 p4_pb_probe.head_joint_nodes，与解剖无关）
# ---------------------------------------------------------------------------
def end_nodes(mesh: dict, *, load_frac: float = LOAD_FRAC) -> dict:
    """+Z 端 = cap（接触承载）端；−Z 端 = 分布载荷端；力向 +Z（压缩）。

    与 `head_joint_nodes` 的差别仅是：载荷端带宽由参数给出（股骨 0.20），
    且**不返回股骨专属的 band_eccentricity**（那是 P4-a′ 诊断量）。
    """
    nodes_rot = bf._align_pca_rotate(mesh["nodes"])
    axis_idx, a_min, a_max, top_local, _bot = bf._select_top_bottom(
        nodes_rot, axis=None, top_frac=TOP_FRAC, bot_frac=TOP_FRAC)
    span = float(a_max - a_min)
    coord = nodes_rot[:, axis_idx]
    load_thr = a_min + load_frac * span
    load_sel = (np.flatnonzero(coord <= load_thr) + 1).astype(np.int64)  # 1-based
    e_axis = np.zeros(3, dtype=float)
    e_axis[axis_idx] = 1.0
    return {
        "axis_idx": int(axis_idx),
        "a_min_mm": float(a_min), "a_max_mm": float(a_max), "span_mm": span,
        "load_nodes_1based": load_sel.tolist(),
        "n_load_nodes": int(len(load_sel)),
        "cap_nodes_1based": [int(x) for x in top_local],
        "n_cap_nodes": int(len(top_local)),
        # bone_feb 语义：load_vec = -load_n × load_dir_vec ⇒ 传 -e_axis 得 +e_axis 力
        "load_dir_vec": (-e_axis).tolist(),
    }


def anchor_nodes(cap_nodes, n: int = N_ANCHOR) -> np.ndarray:
    arr = np.asarray(cap_nodes, dtype=np.int64)
    idx = np.unique(np.linspace(0, len(arr) - 1, n).round().astype(int))
    return arr[idx]


# ---------------------------------------------------------------------------
# cap 几何 —— **壳扩展版**（P5.0 的 `cap_geom` 不能直接用于椎体，见下）
# ---------------------------------------------------------------------------
def _combined_elems_all(mesh: dict, spec) -> list[list[int]]:
    """solid **+ shell** 全部单元（1-based node-id）。

    与 `p5_head_contact._combined_elems` 的唯一差别：**并入壳单元**。

    为什么必须并（2026-10-15 T6 实测 4/4 负 Jacobian 根因）：
      - 长骨（股/胫/腓）的皮质 `*_CORT` 是 **solid hex**，本就在 solid 域内
        ⇒ P5.0 的口径恰好落在**外表面**，故股骨能跑通（接触承载 100.0001%）。
      - 中轴骨（椎/骨盆/颅骨）的 `*_CORT` 是 **`*ELEMENT_SHELL`**（solid 数 = 0）
        ⇒ 只取 solid 会把 cap 落进**松质骨内部**，压板压的是内部、
        皮质壳被从外侧压溃 ⇒ 负 Jacobian。
    ⇒ 这是 P2 推广中继「9 部位非同构」之后的**第二个非同性点**。
    """
    elems: list[list[int]] = []
    for sd in spec.solid_domains:
        b = mesh["by_solid"][sd.mesh_name]
        for row in b["tets"]:
            elems.append([int(x) for x in row])
        for row in b["hexes"]:
            elems.append([int(x) for x in row])
    for sh in spec.shell_domains:
        b = mesh["by_shell"][sh.mesh_name]
        for row in b["q4"]:
            elems.append([int(x) for x in row])
        for row in b["t3"]:
            elems.append([int(x) for x in row])
    return elems


def cap_geom(mesh: dict, spec, top_frac: float = TOP_FRAC,
             include_shell: bool = True) -> dict:
    """`p5_head_contact.cap_geom` 的壳扩展版 —— **筛选口径逐字相同**。

    口径不变：`band_facets` 取端带 → `axdiv = normal·axis` →
    只保留 `axdiv > CAP_AXDIV_MIN`（10-14 根因修复，只取真端面）。
    仅把**候选单元集**从「只有 solid」扩为「solid + shell」。

    `include_shell=False` 时退回 `p5` 原口径（solid-only），用于
    **单变量对照**（2026-10-15 T6 诊断）。
    """
    nodes_rot = bf._align_pca_rotate(mesh["nodes"])
    axis_idx, a_min, a_max, top_local, _bot = bf._select_top_bottom(
        nodes_rot, axis=None, top_frac=top_frac, bot_frac=top_frac)
    if include_shell:
        elems = _combined_elems_all(mesh, spec)
    else:
        elems = p5._combined_elems(mesh, spec)
    facets = p1b.oriented_facets(nodes_rot, elems)
    band = p1b.band_facets(facets, top_local)
    axis = np.zeros(3)
    axis[axis_idx] = 1.0
    for f in band:
        f["axdiv"] = float(np.dot(f["normal"], axis))
    cap = [f for f in band if f["axdiv"] > p5.CAP_AXDIV_MIN]
    return {
        "axis_idx": int(axis_idx), "axis": axis.tolist(),
        "a_min": float(a_min), "a_max": float(a_max),
        "n_top_nodes": len(top_local), "top_local": top_local,
        "band": band, "cap": cap,
        "n_cap_facets": int(len(cap)),
        "a_cap_mm2": float(sum(f["area"] for f in cap)),
        "n_cap_tri": int(sum(1 for f in cap if len(f["poly"]) == 3)),
        "n_cap_quad": int(sum(1 for f in cap if len(f["poly"]) == 4)),
    }


def official_domain(bone: str) -> tuple[str, str]:
    """该部位的**官方 A2 判定域**（照 BC2 的 `main_domain`，不用 solid 域）。

    胸椎 T6 例：BC2 判 `shell_t6_cort`（皮质壳 712 单元），
    而 `load_ctx` 给的 `solid_dom` 是 `t6_spon`（松质骨）—— **两者不是一回事**。
    判皮质强度必须用皮质域。
    """
    fmj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    row = next((r for r in fmj["parts"] if r["bone"] == bone), None)
    if not row:
        return "", ""
    return str(row.get("main_domain", "")), str(row.get("main_kind", ""))


# ---------------------------------------------------------------------------
# 后处理（泛化：p5 的版本硬编码 displacement/femur_r_nodes）
# ---------------------------------------------------------------------------
def load_bear_stats(h5: Path, load_n: float, bone: str) -> dict:
    """接触承载验证 + 横向滑移哨兵。

    口径与 P5.0 逐字一致：读 FEBio 原生 ``surface_data/contact force`` 的
    **单面**（`bone_top`）—— 两面为作用力/反作用力，相加恒 0。
    """
    import h5py
    with h5py.File(h5, "r") as h:
        sid = sorted(h["states"].keys(), key=int)[-1]
        sd = h[f"states/{sid}"]
        dis = np.asarray(sd[f"node_data/displacement/{bone}_nodes"], dtype=float)
        if "surface_data" not in sd or "contact force" not in sd["surface_data"]:
            raise KeyError(
                f"h5 无原生 contact force 输出（states/{sid}/surface_data/contact force）"
                " —— .feb 的 <Output> 段必须含 <var type=\"contact force\"/>")
        cf = sd["surface_data/contact force"]
        f_bone = np.asarray(cf["bone_top"], dtype=float).sum(axis=0)
        f_plate = np.asarray(cf["plate_bottom"], dtype=float).sum(axis=0)
    r_z = float(f_bone[2])
    return {
        "R_z_contact_n": r_z,
        "contact_f_bone_n": [float(x) for x in f_bone],
        "contact_f_plate_n": [float(x) for x in f_plate],
        "contact_newton3_resid_n": float(np.abs(f_bone + f_plate).max()),
        "R_axis_check": abs(abs(r_z) - load_n) / load_n,
        "bear_frac": abs(r_z) / load_n if load_n else float("nan"),
        "dis_xy_max_mm": float(np.sqrt(dis[:, 0] ** 2 + dis[:, 1] ** 2).max()),
        "dis_z_max_mm": float(np.abs(dis[:, 2]).max()),
    }


def domain_stat(h5: Path, dom: str) -> dict:
    """全域统计（单子区域部位用；口径 = median 为官方 Route 乙，另报 p95/max）。

    弯曲指标 `r2_linear` = σ_zz 对单元横向坐标 (x,y) 的线性回归 r²
    （沿用 P4-a′ 口径：基线 0.291 vs P4-a′ 0.054）。
    """
    import h5py
    with h5py.File(h5, "r") as h:
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(np.float64)
        conn = h[f"meshes/0/domains/{dom}"][:]
        sids = sorted(h["states"].keys(), key=int)
        stress = h[f"states/{sids[-1]}/element_data/stress/{dom}"][:].astype(np.float64)
    idx = conn[:, 1:] - 1
    cent = coords[idx].mean(axis=1)
    vm = pa.vm_voigt(stress)
    szz = stress[:, 2]
    # 横向坐标（相对质心，轴向分量清零 —— 承 P4 实测 e=365mm 虚高的修正）
    c = cent - cent.mean(axis=0)
    ax = int(np.argmax(coords.max(0) - coords.min(0)))
    perp = c.copy()
    perp[:, ax] = 0.0
    # 单变量投影后回归（与 P4-a′ 同式：对 perp 的两个主方向做最小二乘）
    r2 = 0.0
    try:
        A = np.column_stack([perp[:, 0], perp[:, 1], np.ones(len(perp))])
        coef, *_ = np.linalg.lstsq(A, szz, rcond=None)
        pred = A @ coef
        ss_res = float(((szz - pred) ** 2).sum())
        ss_tot = float(((szz - szz.mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    except Exception:  # noqa: BLE001
        r2 = 0.0
    return {
        "n_elem": int(len(vm)),
        "median_mpa": float(np.median(vm)),
        "p95_mpa": float(np.percentile(vm, 95)),
        "max_mpa": float(vm.max()),
        "szz_median_mpa": float(np.median(szz)),
        "szz_min_mpa": float(szz.min()),
        "szz_max_mpa": float(szz.max()),
        "r2_linear": float(r2),
    }


# ---------------------------------------------------------------------------
# 构建 / 运行 / 判定
# ---------------------------------------------------------------------------
def build_one(bone: str, feb: Path, load: float, ctx, en, geom, pg, penalty: float,
              with_shell: bool = True) -> Path:
    part, mesh, solid_names, shell_specs, dom, joint, loads, heights = ctx
    bf.build_bone_feb(
        part.spec, feb, load_n=load, top_frac=TOP_FRAC, bot_frac=TOP_FRAC,
        axis=None, with_shell=with_shell, time_steps=TIME_STEPS, mesh=mesh,
        align_pca=True,
        lateral_support=LATERAL,
        # ⚠️ 必须 False（逐字承 P5.0 `p5_pen`）：先生成**全固定** bottom_fix，
        #    再由 `_anchor_xy_patch` 只把 z_dof 置 0 ⇒ 锚点只约束 x/y、轴向走接触。
        #    误传 True 会先生成 x/y 固定、z 自由，随后 patch 匹配到 z 已是 0 的块
        #    而把 x/y 一并清零 ⇒ bottom_fix 三个 dof 全自由
        #    ⇒ FEBio 报 "No degrees of freedom are constrained"（2026-10-15 T6 4/4 FAIL 实测）。
        #    FEBio 约定：x_dof=1 **固定** / 0 自由（由 lat_support_0 的 x/y=1,z=0 反证）。
        bottom_fix_axial_only=False,
        load_nodes_override=np.asarray(en["load_nodes_1based"], dtype=np.int64),
        fix_nodes_override=anchor_nodes(en["cap_nodes_1based"]),
        load_dir_vec=np.asarray(en["load_dir_vec"], dtype=float),
    )
    text = feb.read_text(encoding="ISO-8859-1")
    text = p5._inject_plate_contact(
        text, pg, geom, fric=FRIC, disp_mm=PLATE_DISP_MM,
        first_extra_id=int(mesh["nodes"].shape[0]) + 1, penalty=penalty)
    text = p5._anchor_xy_patch(text)   # z_dof 置 0 ⇒ 轴向全走接触
    feb.write_text(text, encoding="ISO-8859-1")
    return feb


def judge_case(rec: dict) -> dict:
    """预注册 §5.1 能力门（跑前写死）。不通过 = 该部位不可交付，不回调 penalty。"""
    reasons = []
    if not rec.get("converged"):
        reasons.append("不收敛")
    bear = rec.get("bear_frac")
    if bear is None or not (BEAR_LO <= bear <= BEAR_HI):
        reasons.append(f"接触承载 {bear:.3f} ∉ [{BEAR_LO},{BEAR_HI}]"
                       if bear is not None else "无接触承载读数")
    bend = rec.get("r2_linear", 0.0) is not None and (
        rec.get("r2_linear", 0.0) >= R2_MIN or rec.get("szz_max_mpa", 0.0) > 0)
    if not bend:
        reasons.append(f"弯曲路径缺失 r2={rec.get('r2_linear', float('nan')):.3f} "
                       f"szz_max={rec.get('szz_max_mpa', float('nan')):.1f}")
    return {"pass": not reasons, "reasons": reasons,
            "bear_ok": bear is not None and BEAR_LO <= bear <= BEAR_HI,
            "bend_ok": bool(bend)}


def run_bone(bone: str, heights: list[int], force_loads: list[float] | None = None,
             with_shell: bool = True) -> dict:
    tag = "" if with_shell else "   <<< no-shell 对照（与股骨同构）"
    print(f"\n{'='*92}\n[part] {bone}{tag}\n{'='*92}", flush=True)
    ctx = p4_pb_probe.load_ctx(bone)
    part, mesh, solid_names, shell_specs, dom, joint, joint_loads, h_arr = ctx
    share = LOAD_SHARE[bone]
    print(f"[{bone}] mesh={mesh['n_nodes']} nodes dom={dom} joint={joint} "
          f"share={share:.4f} with_shell={with_shell}", flush=True)

    if force_loads:
        # 诊断模式：显式载荷（不按高度取）。仅用于**判别 BC 形式 vs 载荷/网格上限**，
        # 不作为交付结果、不得用于凑判据（预注册 §4 红线）。
        plan = [(0, float(h_arr[0]), L * share) for L in force_loads]
    else:
        plan = [(hno, float(h_arr[hno - 1]), float(joint_loads[hno - 1]) * share)
                for hno in heights]

    en = end_nodes(mesh)
    geom = cap_geom(mesh, part.spec, top_frac=TOP_FRAC, include_shell=with_shell)
    pg = p5.plate_geom_from_mesh(bf._align_pca_rotate(mesh["nodes"]), geom)
    dom_official, dom_kind = official_domain(bone)
    # no-shell 对照里 shell 域不存在，判定域必须退回 solid 域
    stat_dom = (dom_official or dom) if with_shell else dom
    shells_for_post = shell_specs if with_shell else []
    f_max = float(joint_loads[max(heights) - 1]) * share
    pen = penalty_for(f_max)
    print(f"[{bone}] cap={geom['n_cap_facets']}f "
          f"({geom['n_cap_tri']}tri/{geom['n_cap_quad']}quad) A={geom['a_cap_mm2']:.1f}mm²"
          f" | 载荷端 n={en['n_load_nodes']} cap端 n={en['n_cap_nodes']} "
          f"| F_max={f_max:.1f}N penalty={pen:.2f} "
          f"| 判定域={stat_dom}({dom_kind}) solid域={dom}", flush=True)

    rows = []
    run_tag = "" if with_shell else "_noshell"
    for hno, h_m, load in plan:
        d = BASE / bone / (f"diag{int(load)}N{run_tag}" if force_loads
                           else f"h{hno:02d}{run_tag}")
        d.mkdir(parents=True, exist_ok=True)
        feb = d / f"{bone}.feb"
        h5 = feb.with_suffix(".hdf5")
        rec = {"bone": bone, "h_no": hno, "h_m": float(h_m),
               "load_n": load, "load_share": share, "penalty": penalty_for(load),
               "domain": stat_dom, "solid_domain": dom}
        if not h5.is_file():
            t = time.time()
            try:
                build_one(bone, feb, load, ctx, en, geom, pg, rec["penalty"],
                          with_shell=with_shell)
                bf.run_bone_feb(feb, workdir=feb.parent, timeout=1800)
                if feb.with_suffix(".xplt").is_file():
                    bf.post_process(feb.with_suffix(".xplt"), h5,
                                    solid_names=solid_names,
                                    shell_specs=shells_for_post, label=bone)
                rec["converged"] = True
            except Exception as exc:  # noqa: BLE001
                rec["converged"] = False
                rec["error"] = f"{type(exc).__name__}: {exc}"
            rec["seconds"] = round(time.time() - t, 1)
        else:
            rec["converged"] = True
            rec["seconds"] = 0.0
        if h5.is_file():
            try:
                rec.update(domain_stat(h5, stat_dom))
                rec.update(load_bear_stats(h5, load, bone))
            except Exception as exc:  # noqa: BLE001
                rec["post_error"] = f"{type(exc).__name__}: {exc}"
        rec.update(judge_case(rec))
        rows.append(rec)
        print(f"[{bone}] h={hno:2d} F={load:9.1f}N pen={rec['penalty']:6.2f} "
              f"{'OK ' if rec.get('converged') else 'FAIL'} "
              + (f"med={rec['median_mpa']:8.2f} p95={rec['p95_mpa']:8.2f} "
                 f"szzmax={rec['szz_max_mpa']:8.1f} r2={rec['r2_linear']:6.3f} "
                 f"bear={rec['bear_frac']:6.4f} [{rec['seconds']}s] "
                 f"pass={rec['pass']}" if rec.get("median_mpa") is not None
                 else rec.get("error", rec.get("post_error", "?"))[:90]),
              flush=True)
    return {"bone": bone, "joint": joint, "domain": dom, "load_share": share,
            "penalty_at_fmax": pen, "cap_facets": geom["n_cap_facets"],
            "cap_area_mm2": geom["a_cap_mm2"], "rows": rows}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bones", nargs="*", default=SUBSET)
    ap.add_argument("--loads", nargs="*", type=float, default=None,
                    help="诊断模式：显式关节载荷 N（乘 share 后求解）。"
                         "仅用于判别 BC 形式 vs 载荷/网格上限，不得用于凑判据。")
    ap.add_argument("--no-shell", action="store_true",
                    help="单变量对照：把皮质壳整个从 .feb 去掉（with_shell=False），"
                         "cap 退回 solid-only 口径 ⇒ T6 变成与股骨同构的单一 solid 域。")
    args = ap.parse_args()
    t0 = time.time()
    out = []
    for b in args.bones:
        out.append(run_bone(b, HEIGHTS, force_loads=args.loads,
                            with_shell=not args.no_shell))
    doc = {"probe": "P2 端部 BC + 载荷分配 · 能力门",
           "prereg": "docs/methods/P2_端部BC与载荷分配_预注册.md",
           "frozen": {"gap_mm": GAP_MM, "top_frac": TOP_FRAC,
                      "cap_axdiv_min": p5.CAP_AXDIV_MIN, "fric": FRIC,
                      "lateral_support": LATERAL, "plate_disp_mm": PLATE_DISP_MM,
                      "load_frac": LOAD_FRAC,
                      "penalty_formula": f"{PENALTY_REF_VALUE} * F_ref / {PENALTY_REF_LOAD_N}",
                      "bear_gate": [BEAR_LO, BEAR_HI], "r2_min": R2_MIN},
           "load_share": LOAD_SHARE, "heights_m": HEIGHTS, "parts": out,
           "elapsed_s": round(time.time() - t0, 1)}
    (RES / "P2_HEAD_CONTACT.json").write_text(
        json.dumps(doc, indent=2, ensure_ascii=False), encoding="utf-8")
    print("\n" + "=" * 92)
    for p in out:
        rs = p["rows"]
        npass = sum(1 for r in rs if r.get("pass"))
        print(f"  {p['bone']:<12} 能力门 {npass}/{len(rs)}"
              f"  penalty={p['penalty_at_fmax']:.2f}  cap={p['cap_facets']}f/"
              f"{p['cap_area_mm2']:.1f}mm²")
        for r in rs:
            if not r.get("pass") and r.get("converged"):
                print(f"      h={r['h_no']:2d} 未过：{'; '.join(r['reasons'])}")
    print(f"\n落盘 {RES / 'P2_HEAD_CONTACT.json'}  总耗时 {doc['elapsed_s']}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
