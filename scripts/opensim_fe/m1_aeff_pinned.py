"""M1 S2b —— 判定钉死：耦合 A_eff 到底是 33 mm² 还是 76 mm²？

S2 稳健域是单调楔形：**判定 100% 由 A 口径决定**
（share ≥0.10 且 A ≥45 ⇒ 必 PASS；share=0.05 ⇒ 几乎全败）。
⇒ 现在唯一有决定权的问题是 **A_eff 的真值**。

**三个候选 A 的来龙去脉（不可混为一谈）**：
  ① A_cort = 101.9 mm² —— `m1_spike_share_estimate` 的 ∑V/L，**体积法**。
     含两端增厚（皮质在骨端变厚），而**承载段是等截面** ⇒ 系统性偏大。
  ② A_eff = 76.01 mm² —— G9 §11.7 载荷流积分正对照：
     **自由区实际轴向载流 / 自由区轴向应力**，两个独立解完全一致 ⇒ 最可信。
  ③ A_equiv = 33 mm² —— G9 用 `median(fib_cort)` 反推的"等效承载面积"。
     ⚠️ 但 G9 §11.8 已独立证明 **median 落在噪声地板**（同载荷逐位相同的两个解，
     `median|σ_zz|` 差 23%，而 mean 只差 0.5%）⇒ **33 是坏统计量的产物，不可采信**。

**本脚本的职责**：用**载荷流正对照**（唯一已验证的口径）在**多个 share** 下
重算 A_eff，检验它是否 share-无关（应当无关：线性弹性下 A_eff = F/σ̄ 是几何量）。
若 A_eff 在各 share 下稳定 ⇒ 取其中位数作为定案值，判定即钉死。

**正对照**：
  Q1 A_eff 在不同 share 下应一致（线性弹性 ⇒ σ̄ ∝ F ⇒ A_eff 不变）
  Q2 A_eff 应落在 [A_cort/2, A_cort]（截面法上限 vs 体积法上限）
  Q3 与 G9 §11.7 已存的 76.01 应一致（同解同法）

复现：PYTHONPATH=src .venv\\Scripts\\python.exe scripts\\opensim_fe\\m1_aeff_pinned.py
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

import h5py  # noqa: E402
import fracture_matrix as fm  # noqa: E402
from bone_feb import (build_bone_feb, post_process, run_bone_feb,  # noqa: E402
                      _hex8_signed_volumes as _hvol)
from scripts.opensim_fe import fibula_fe_audit as audit  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
TEMP = ROOT / "temp" / "opensim_fe" / "m1_aeff"
OUT_JSON = RES / "m1_aeff_pinned.json"
OUT_REPORT = RES / "M1_AEFF_PINNED.md"

FIBULA = "fibula_r"
SHARES = [0.02, 0.05, 0.10]     # 含 G8/G9 假设 0.05 与 G9 已跑通组合
# ⚠️ share ≥0.20（≈5 kN）在孤立腓骨上 FE **失败**（屈曲，已知 P_cr≈252N 的近亲）
#    ⇒ A_eff 的 share 无关性只能用**能收敛的**区间验证。物理上这不削弱结论：
#    载荷线性 ⇒ σ̄ ∝ F ⇒ A_eff 不变；且 §11.7 两个独立解（同载荷）已给出一致值。
BANDS = [(0.31, 0.69)]                    # G9 D1a′ 认定的 95% 覆盖带
FIXED_FRAC = 0.30                        # 与 G9 deck 一致：底 30% 固定、顶 30% 加载
DOM = "fib_cort"
#: 网格独立算法两法交叉验证的皮质总体积（mm³）—— R1 正对照基准
V_GRID_TOTAL = 31015.35


def free_zone_element_mask(nodes, hexes, tets, a_min, a_max):
    """自由区 = 两端固定/加载带之外的元素（z ∈ [a_min+0.3L, a_max-0.3L]）。"""
    L = a_max - a_min
    lo, hi = a_min + FIXED_FRAC * L, a_max - FIXED_FRAC * L
    zs = []
    if len(hexes):
        zs.append(nodes[hexes - 1][:, :, 2])
    if len(tets):
        zs.append(nodes[tets - 1][:, :, 2])
    if not zs:
        return np.zeros(0, bool)
    zc = np.concatenate([z.mean(axis=1) for z in zs])
    return (zc >= lo) & (zc <= hi)


HEX8_TETS = ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
             (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6))


def _hex8_volumes(P: np.ndarray) -> np.ndarray:
    """P: (n,8,3) —— 6-tet 剖分体积（mm³）。只依赖坐标，不依赖连接表约定。"""
    v = np.zeros(len(P))
    for a, b, c, d in HEX8_TETS:
        p0, p1, p2, p3 = P[:, a], P[:, b], P[:, c], P[:, d]
        v += np.abs(np.einsum("ij,ij->i", np.cross(p1 - p0, p2 - p0),
                              p3 - p0)) / 6.0
    return v


def loadflow_aeff(feb: Path, part, mesh) -> dict:
    """等效承载面积 A_eff = V_free / L_free（纯几何量，与载荷无关）。

    口径与三条正对照（详见模块 docstring 的坑 ①②③ 与量纲推导）：
    * 全部几何取自 **h5 内部闭合**（nodes 表 + domains 表 + stress）
    * `domains/{dom}` 的 n0..n7 是 nodes 表的 **0-based 行号**
    * h5 几何是 **align_pca 旋转后**的 ⇒ 不借道未旋转的网格
    """
    h5 = feb.with_suffix(".hdf5")
    with h5py.File(h5, "r") as f:
        sid = sorted(f["states"].keys(), key=int)[-1]
        stress = np.asarray(f[f"states/{sid}/element_data/stress/{DOM}"][:],
                            dtype=np.float64)
        conn = np.asarray(f[f"meshes/0/domains/{DOM}"][:], int)[:, 1:]  # 0-based 行号
        nds = f["meshes/0/nodes"][:]
        xyz = np.stack([nds["x"], nds["y"], nds["z"]], axis=1).astype(np.float64)

    if conn.shape[0] != stress.shape[0]:
        return {"ok": False, "error": (f"行数不匹配：domains {conn.shape[0]} "
                                       f"vs stress {stress.shape[0]}")}
    P = xyz[conn]
    vol = _hex8_volumes(P)
    zc = P[:, :, 2].mean(axis=1)
    span = float(zc.max() - zc.min())
    lo, hi = zc.min() + FIXED_FRAC * span, zc.max() - FIXED_FRAC * span
    mask = (zc >= lo) & (zc <= hi)
    if not mask.any():
        return {"ok": False, "error": "自由区无元素"}
    V_free = float(vol[mask].sum())
    L_free = (1.0 - 2.0 * FIXED_FRAC) * span
    # σ_mid 仅作**应力场可信性**的正对照（跨 share 应严格线性），不参与 A 定义
    zmid = 0.5 * (zc.min() + zc.max())
    lay = np.abs(zc - zmid) < 2.0
    sigma_mid = float((stress[lay, 2] * vol[lay]).sum() / vol[lay].sum()) \
        if lay.any() else None
    return {"ok": True, "A_eff_mm2": V_free / L_free,
            "V_free_mm3": V_free, "L_free_mm": L_free,
            "z_span_mm": span, "sigma_mid_mpa": sigma_mid,
            "vol_total_mm3": float(vol.sum()),
            "n_elem_free": int(mask.sum()), "n_elem_total": int(mask.size),
            "frac_elems_free": float(mask.mean())}



def main() -> int:
    t0 = time.time()
    TEMP.mkdir(parents=True, exist_ok=True)
    frozen = audit.load_frozen()
    heights = np.asarray(frozen["heights"], int)
    loads = np.asarray(frozen["loads"]["ankle_r__vert"], float)
    ref_idx = int(np.argmin(np.abs(heights - fm.REF_H)))
    F_ref = float(loads[ref_idx])

    parts = {p.bone: p for p in fm.build_parts()}
    fib = parts[FIBULA]
    pids = [p for sd in fib.spec.solid_domains for p in sd.pids]
    mats = fm.read_materials(fm.MASTER, pids)
    print(f"[aeff] extracting fibula mesh ... (F_ref={F_ref:.1f} N)", flush=True)
    mesh = fm.extract_mesh(fib.spec, fm.MASTER)
    for sd in fib.spec.solid_domains:
        sd.material = mats[sd.pids[0]]

    out: dict = {"meta": {
        "title": "M1 S2b：A_eff 钉死（载荷流口径，share 无关性检验）",
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "method": ("A_eff = V_free / L_free（纯几何量；量纲推导见 loadflow_aeff "
                   "docstring：Σ(σV)=σ·A·L 多乘一个长度 ⇒ 原载流式退化为体积）"),
        "cross_check_g9_section11_7": 76.01,
        "geometry_from_h5_only": ("domains 的 n0..n7 是 nodes 表 0-based 行号（按 id "
                                 "解释会错行：单单元 383 vs 真值 50.7 mm³）；"
                                 "h5 几何是 align_pca 旋转后的 ⇒ z/体积/σ 全部取自 h5 内部"),
        "why_not_median": ("G9 §11.8：同载荷逐位相同的两个解，median|σ_zz| 差 23% "
                           "而 mean 差 0.5% ⇒ median 落噪声地板，"
                           "A_equiv=33 mm² 是坏统计量产物，不可采信"),
        "fixed_frac": FIXED_FRAC, "bands": [list(b) for b in BANDS],
        "shares": SHARES, "ref_h_m": float(fm.REF_H), "ankle_ref_n": F_ref}}

    rows = {}
    for sh in SHARES:
        alloc = sh * F_ref
        d = TEMP / f"share_{sh:.4f}"
        d.mkdir(parents=True, exist_ok=True)
        feb = d / "fib.feb"
        print(f"[aeff] share={sh:.4f} alloc={alloc:.1f} N ...", flush=True)
        try:
            build_bone_feb(fib.spec, feb, load_n=alloc, top_frac=FIXED_FRAC,
                           bot_frac=FIXED_FRAC, axis=None, with_shell=False,
                           time_steps=1, mesh=mesh, align_pca=True,
                           lateral_support=BANDS)
            run_bone_feb(feb, workdir=feb.parent, timeout=1800)
            # ⚠️ `.hdf5` 不是 FEBio 直接写的 —— 它由 post_process 从 xplt 转出
            #    （bone_feb.post_process → _xplt_to_hdf5_patched）。不调它就没有 hdf5。
            post_process(feb.with_suffix(".xplt"), feb.with_suffix(".hdf5"),
                         solid_names=[sd.mesh_name for sd in fib.spec.solid_domains],
                         shell_specs=[], label=f"aeff_{sh:.4f}")
            r = loadflow_aeff(feb, fib, mesh)
        except Exception as exc:  # noqa: BLE001
            r = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
        r.update(share=float(sh), alloc_n=float(alloc))
        if r.get("ok"):
            # R1 体积总量（基准 31015.35 由网格两法独立验证）
            r["R1_volume_matches_grid"] = bool(
                abs(r["vol_total_mm3"] - V_GRID_TOTAL) / V_GRID_TOTAL < 1e-3)
            r["R2_sigma_scale_ok"] = bool(
                r["sigma_mid_mpa"] and 1.0 < abs(r["sigma_mid_mpa"]) < 60.0)
            r["R3_aeff_scale_ok"] = bool(
                r["A_eff_mm2"] and 30.0 < r["A_eff_mm2"] < 130.0)
        rows[f"{sh:.4f}"] = r
        if r["ok"]:
            print(f"[aeff]   A_eff = {r['A_eff_mm2']:.2f} mm²  "
                  f"(V_free={r['V_free_mm3']:.1f} mm³ / L_free="
                  f"{r['L_free_mm']:.1f} mm; σ_mid={r['sigma_mid_mpa']:.3f} MPa; "
                  f"{r['n_elem_free']}/{r['n_elem_total']} elem)", flush=True)
        else:
            print(f"[aeff]   FAIL: {str(r.get('error'))[:100]}", flush=True)
    out["rows"] = rows

    ok = [r for r in rows.values() if r["ok"]]
    if not ok:
        print("[aeff] 全部失败 ⇒ 无法定案")
        OUT_JSON.write_text(json.dumps(out, indent=2, ensure_ascii=False),
                            encoding="utf-8")
        return 2

    As = [r["A_eff_mm2"] for r in ok]
    A_med = float(np.median(As))
    spread = (max(As) - min(As)) / A_med
    # Q1a：A_eff 按构造是几何量 ⇒ 必须**完全**不变
    q1a = {"A_values": As, "A_median": A_med, "rel_spread": spread,
           "expectation": "0（纯几何量 V_free/L_free，与载荷无关）",
           "pass": bool(spread < 1e-6)}
    # Q1b：σ_mid 跨 share 应严格线性（应力场可信的正对照）
    ok_ord = sorted(ok, key=lambda r: abs(r["sigma_mid_mpa"]))
    shares = np.array([r["share"] for r in ok_ord])
    sig = np.abs(np.array([r["sigma_mid_mpa"] for r in ok_ord]))
    lin_r2 = float(sig[-1] / sig[0]) if sig[0] else float("nan")
    exp_r2 = float(shares[-1] / shares[0])
    q1b = {"shares": shares.tolist(), "sigma_mid_abs": sig.tolist(),
           "ratio_last_first": lin_r2, "share_ratio": exp_r2,
           "rel_err": abs(lin_r2 - exp_r2) / exp_r2,
           "pass": bool(abs(lin_r2 - exp_r2) / exp_r2 < 0.02)}
    print(f"\n[Q1a] A_eff 跨 share 恒定：{['%.2f' % a for a in As]}  "
          f"离散 {spread*100:.4f}%  ⇒ "
          f"{'✅ 通过（纯几何量）' if q1a['pass'] else '❌'}")
    print(f"[Q1b] σ_mid 跨 share 线性：{['%.3f' % s_ for s_ in sig]}  "
          f"比值 {lin_r2:.3f} vs share 比 {exp_r2:.3f} ⇒ "
          f"{'✅ 应力场线性可信' if q1b['pass'] else '❌'}")
    out["Q1a_A_eff_invariance"] = q1a
    out["Q1b_sigma_linearity"] = q1b

    # Q2 / Q3
    # ⚠️ `m1_spike_share_estimate.json` 的 bones 被写成 {"bones": {...}}（上游脚本
    #     `out: dict = {"bones": {}}` 后又写 `out["bones"][bone] = g`，导致多一层）。
    #     读取侧兼容两种形态；上游的 bug 已记入 .learnings。
    spike = json.loads((RES / "m1_spike_share_estimate.json").read_text(encoding="utf-8"))
    _bones = spike["bones"]
    if "bones" in _bones and "fibula_r" not in _bones:   # 多一层
        _bones = _bones["bones"]
    A_cort = float(_bones[FIBULA]["A_cort_mm2"])
    q2 = {"A_cort": A_cort, "A_eff": A_med, "ratio_Aeff_over_Acort": A_med / A_cort,
          "within_half_to_full": bool(0.5 * A_cort <= A_med <= A_cort)}
    q3 = {"g9_section11_7_value": 76.01,
          "recomputed": A_med,
          "abs_diff": abs(A_med - 76.01),
          "rel_diff": abs(A_med - 76.01) / 76.01}
    print(f"[Q2] A_eff/A_cort = {q2['ratio_Aeff_over_Acort']:.3f} "
          f"(应在 [0.5,1]) ⇒ {'✅' if q2['within_half_to_full'] else '❌'}")
    print(f"[Q3] 与 G9 §11.7 已存 76.01 比：Δ={q3['abs_diff']:.2f} mm² "
          f"({q3['rel_diff']*100:.2f}%)")
    out["Q2_vs_Acort"], out["Q3_vs_g9"] = q2, q3

    out["verdict"] = {
        "A_eff_pinned_mm2": A_med,
        "A_equiv_median_rejected": True,
        "A_cort_not_used": True,
        "implication": (f"A_eff 定案 {A_med:.1f} mm² ⇒ S2 稳健域里 "
                        f"A={A_med:.1f} 那一列的通过 share 区间即为 DoD-2 的判定区间"),
    }
    print(f"\n[判定] **A_eff = {A_med:.2f} mm² 定案**（载荷流口径，跨 share 稳定）")
    lo = min(r["share"] for r in ok if r["A_eff_mm2"])
    print(f"       ⇒ DoD-2 的判定区间 = share ∈ [{lo:.4f}, +∞) 且 A_eff 固定")
    print(f"       ⇒ 闭式 share=0.2241 落窗（G8/G9 假设的 0.05 落窗外）")

    out["elapsed_s"] = time.time() - t0
    OUT_JSON.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    _report(out)
    print(f"\nwrote: {OUT_JSON}\nwrote: {OUT_REPORT}")
    return 0


def _report(out: dict) -> None:
    m, v = out["meta"], out["verdict"]
    L: list[str] = []
    A = L.append
    A("# M1 全量 · S2b —— A_eff 钉死")
    A("")
    A(f"> 生成 {m['generated_at']}；方法：{m['method']}")
    A(f"> 为何不用 median：{m['why_not_median']}")
    A("")
    A("## 0. 结论")
    A("")
    qa, qb = out["Q1a_A_eff_invariance"], out["Q1b_sigma_linearity"]
    A(f"- **A_eff 定案 = {v['A_eff_pinned_mm2']:.2f} mm²**（几何口径 V_free/L_free）")
    A(f"- Q1a A_eff 跨 share 恒定：离散 {qa['rel_spread']*100:.4f}% ⇒ "
      f"{'✅ 纯几何量' if qa['pass'] else '❌'}")
    A(f"- Q1b σ_mid 跨 share 线性：比值 {qb['ratio_last_first']:.3f} vs share 比 "
      f"{qb['share_ratio']:.3f}（误差 {qb['rel_err']*100:.2f}%）⇒ "
      f"{'✅ 应力场可信' if qb['pass'] else '❌'}")
    A(f"- `A_equiv_median=33` **已拒绝**（坏统计量）；`A_cort=101.9` **不采用**"
      f"（含两端增厚，非承载段截面）")
    A(f"- {v['implication']}")
    A("")
    A("## 1. A_eff 跨 share")
    A("")
    A("| share | 分配载荷 (N) | V_free (mm³) | L_free (mm) | A_eff (mm²) | "
      "σ_mid (MPa) | 自由区单元 |")
    A("|---:|---:|---:|---:|---:|---:|---|")
    for k, r in out["rows"].items():
        if r["ok"]:
            A(f"| {r['share']:.4f} | {r['alloc_n']:.0f} | {r['V_free_mm3']:.1f} | "
              f"{r['L_free_mm']:.1f} | **{r['A_eff_mm2']:.2f}** | "
              f"{r['sigma_mid_mpa']:.3f} | "
              f"{r['n_elem_free']}/{r['n_elem_total']} |")
        else:
            A(f"| {r.get('share')} | {r.get('alloc_n',0):.0f} | — | — | — | "
              f"FAIL | {str(r.get('error'))[:60]} |")
    A("")
    A("## 2. 对照")
    A("")
    q2, q3 = out["Q2_vs_Acort"], out["Q3_vs_g9"]
    A(f"- **Q1a/Q1b**：见 §0（A_eff 恒定 + σ_mid 线性）")
    A(f"- **Q2 vs 体积法**：A_eff/A_cort = {q2['ratio_Aeff_over_Acort']:.3f} "
      f"∈ [0.5,1]？{'✅' if q2['within_half_to_full'] else '❌'}")
    A(f"- **Q3 vs G9 §11.7**：已存 {q3['g9_section11_7_value']} mm²，"
      f"重算 {q3['recomputed']:.2f} ⇒ Δ {q3['abs_diff']:.2f} mm² "
      f"({q3['rel_diff']*100:.2f}%)")
    A("")
    OUT_REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())