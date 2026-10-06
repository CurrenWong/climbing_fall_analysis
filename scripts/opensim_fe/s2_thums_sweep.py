"""S2 · THUMS AM50 calcaneus rigid-ground HEIGHT sweep, with Achilles activation.

Option C follow-up: the crude whole-foot geometry (``r_foot.vtp``, the lumped
IITD foot mesh) is replaced by the registered **AM50 THUMS calcaneus** two-domain
model (``thums_feb.py``); the OpenSim load pipeline is reused **verbatim**.

Load side (identical chain to ``s2_height_sweep.py`` / ``s3_achilles_sweep.py``)::

    h  -> climbing.coupling.opensim_grf.ground_reaction(h)
       -> climbing.coupling.opensim_fall.run_dead_drop(grf)
       -> climbing.coupling.joint_loads.subtalar_reaction(fall, grf, side="r")
       -> peak_vertical_n  = F_subt(h)          (真跑正动力学，不靠 √h 外推)
    transfer(F_subt, achilles_activation=a)      -> F_subt dir, F_ach = a·F_MAX

FE side: the mesh is built **once** (``thums_feb.load_thums_mesh``) and re-solved
per ``(h, a)`` via ``thums_feb.build_thums_feb(..., use_rigid=False)`` — i.e.
the same ``PressureLoad p = F/A`` convention the old pipeline uses, so the two
are comparable; Achilles (``a>0``) is applied on the (now non-empty)
``achilles`` surface as the reference ``TractionLoad``.

Metric (fracture): ``gauge_max`` on the **cortical** domain (``"calcaneus"``,
hex8) using the mixed-mesh gauge helper ``thums_feb.gauge_von_mises_mixed``
(process-zone regularised peak, ℓ = 4 mm).  The first height where
``gauge_max > 150 MPa`` is the **first-fracture height** for that activation.
This is a **post-processing** criterion only — no element deletion, no fracture
as a solver action.

Outputs (``results/opensim_fe/``)::

    s2_thums_height_sweep.json / .csv / .png   and   S2_THUMS_REPORT.md

Usage::

    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\s2_thums_sweep.py
    ... --heights 1 5 --acts 0 1
    ... --gauge-mm 4 --outdir results\\opensim_fe

Units: mm-N-MPa-s.
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
# scripts/opensim_fe/s2_thums_sweep.py -> scripts/opensim_fe -> scripts -> repo root
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import thums_feb as T  # noqa: E402  (new Option-C module, imported verbatim)

#: default height grid (mirrors ``s2_height_sweep.py``'s DEFAULT_GRID).
DEFAULT_GRID = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50]
#: triceps-surae activation levels (a=0 no Achilles, a=1 max isometric).
DEFAULT_ACTS = [0.0, 1.0]
#: [Y25]-tabled cortical compression strength -> fracture threshold (MPa).
THRESHOLD_MPA = 150.0
#: process-zone regularisation length ℓ (mm) == gauge radius default.
REG_LEN_MM = 4.0
#: Rajagopal2015 triceps-surae F_max sum (N), from ``load_transfer``.
F_MAX_TRICEPS_N = 10885.0
#: [Y25] foot-fracture cluster (m) our first-fracture height is compared with.
Y25_LO_M, Y25_HI_M = 7.0, 9.0
#: OLD whole-foot pipeline first-fracture height (smooth mesh, gauge ℓ=4 mm).
OLD_FIRST_FRACTURE_M = 26.5

CSV_FIELDS = [
    "h", "act", "subtalar_kN", "achilles_kN", "t1_rel_err",
    "max", "p99", "p95", "mean", "gauge_max",
]


# ---------------------------------------------------------------------------
# load side (verbatim mirror of s2_height_sweep._load_at)
# ---------------------------------------------------------------------------
def _load_at(height: float) -> tuple[float, float]:
    """Run the multibody half-drop -> (subtalar peak vertical force N, T1 err)."""
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=height)
    fall = run_dead_drop(grf)
    jr = subtalar_reaction(fall, grf, side="r")
    return float(jr.peak_vertical_n), float(grf.impulse_rel_err)


# ---------------------------------------------------------------------------
# FE post side (mirrors thums_feb_run._read_domain_vm)
# ---------------------------------------------------------------------------
def _read_domain_vm(xplt: Path, domain: str, hdf5: Path) -> np.ndarray:
    """Last-state per-element von Mises (MPa) for ``domain``; writes HDF5 if absent."""
    import h5py

    from climbing.coupling import fe_post
    from pyfebio import xplt as xplt_mod

    if not Path(hdf5).is_file():
        xplt_mod.to_hdf5(str(xplt), str(hdf5))
    with h5py.File(hdf5, "r") as h5:
        stress, used = fe_post._read_last_stress(h5, domain)
    return fe_post.von_mises_from_voigt(stress)


def _solve_case(
    mesh: dict,
    cort_centroids: np.ndarray,
    cort_vol: np.ndarray,
    h: float,
    a: float,
    f_sub: float,
    t1: float,
    gauge_mm: float,
    workdir: Path,
) -> dict:
    """Build -> solve -> post one ``(h, a)`` case; returns the CSV/JSON row.

    ``f_sub``/``t1`` come from ``_load_at(h)`` (once per height, shared by all
    activations) — exactly the loop shape of ``s3_achilles_sweep.py``.
    """
    from climbing.coupling import fe_post
    from climbing.coupling.febio_run import run_febio
    from climbing.coupling.load_transfer import transfer

    spec = transfer(f_sub, achilles_activation=a)

    feb = (workdir / f"s2_thums_h{h:g}_a{a:g}.feb").resolve()
    xplt = feb.with_suffix(".xplt")
    hdf5 = feb.with_suffix(".hdf5")
    for stale in (xplt, hdf5):
        if stale.exists():
            stale.unlink()

    T.build_thums_feb(
        mesh,
        feb,
        load_n=spec.subtalar_n,
        load_dir=spec.subtalar_dir,
        achilles_n=spec.achilles_n,
        achilles_dir=spec.achilles_dir,
        use_rigid=False,          # PressureLoad p=F/A: same convention as old pipeline
        time_steps=1,
    )
    rc = run_febio(feb, workdir=feb.parent, timeout=1800)  # raises on rc != 0
    if rc != 0:  # defensive: run_febio already raises, but be explicit
        raise RuntimeError(f"FEBio rc={rc} for h={h} a={a}")

    # Headline cortical peak via the shared fe_post entry point (writes HDF5);
    # then read the full per-element field for the mixed-mesh gauge.
    peak_public = fe_post.peak_von_mises(xplt, hdf5_path=hdf5, domain="calcaneus")
    vm = _read_domain_vm(xplt, "calcaneus", hdf5)
    if len(vm) != len(cort_centroids):
        raise RuntimeError(
            f"CORT 单元数不匹配：vm={len(vm)} vs hex={len(cort_centroids)} (h={h} a={a})"
        )
    if abs(float(vm.max()) - float(peak_public)) > 1e-6 * max(1.0, float(peak_public)):
        raise RuntimeError(
            f"fe_post.peak_von_mises={peak_public} 与逐单元峰值 {float(vm.max())} 不一致"
        )

    gauge = T.gauge_von_mises_mixed(vm, cort_centroids, cort_vol, gauge_mm)
    row = {
        "h": float(h),
        "act": float(a),
        "subtalar_kN": spec.subtalar_n / 1e3,
        "achilles_kN": spec.achilles_n / 1e3,
        "t1_rel_err": t1,
        "max": float(vm.max()),
        "p99": float(np.percentile(vm, 99)),
        "p95": float(np.percentile(vm, 95)),
        "mean": float(vm.mean()),
        "gauge_max": float(gauge.max()),
    }
    return row


# ---------------------------------------------------------------------------
# reporting helpers
# ---------------------------------------------------------------------------
def _first_fracture(rows: list[dict], acts: list[float]) -> dict[float, dict]:
    """Per activation: first height where ``gauge_max`` exceeds the threshold."""
    out: dict[float, dict] = {}
    for a in acts:
        sub = sorted((r for r in rows if r["act"] == a), key=lambda r: r["h"])
        triggered = [r for r in sub if r["gauge_max"] > THRESHOLD_MPA]
        if not triggered:
            out[a] = {
                "act": a,
                "achilles_kN": a * F_MAX_TRICEPS_N / 1e3,
                "h_first_m": None,
                "h_lo_m": None,
                "h_interp_m": None,
                "gauge_at_first": None,
                "triggered": False,
                "max_h_m": sub[-1]["h"] if sub else None,
            }
            continue
        first = triggered[0]
        # bracket [previous sample, first] for a (rough) linear crossing estimate
        prev = [r for r in sub if r["h"] < first["h"]]
        h_interp = None
        h_lo = None
        if prev:
            p = prev[-1]
            h_lo = p["h"]
            g0, g1 = p["gauge_max"], first["gauge_max"]
            if g1 > g0:
                h_interp = p["h"] + (first["h"] - p["h"]) * (THRESHOLD_MPA - g0) / (g1 - g0)
        out[a] = {
            "act": a,
            "achilles_kN": a * F_MAX_TRICEPS_N / 1e3,
            "h_first_m": first["h"],
            "h_lo_m": h_lo,
            "h_interp_m": h_interp,
            "gauge_at_first": first["gauge_max"],
            "triggered": True,
            "max_h_m": sub[-1]["h"] if sub else None,
        }
    return out


def _write_plot(rows: list[dict], acts: list[float], out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(acts), 1)))
    for i, a in enumerate(acts):
        sub = sorted((r for r in rows if r["act"] == a), key=lambda r: r["h"])
        h = [r["h"] for r in sub]
        c = colors[i]
        ax.plot(h, [r["gauge_max"] for r in sub], "-", color=c, marker="o", ms=4,
                label=f"gauge_max  a={a:g}  (F_ach={a*F_MAX_TRICEPS_N/1e3:.1f} kN)")
        ax.plot(h, [r["max"] for r in sub], ":", color=c, lw=1.0, alpha=0.55,
                marker="^", ms=3, label=f"raw max  a={a:g}")

    ax.axhline(THRESHOLD_MPA, color="r", lw=1.3, ls="--",
               label=f"皮质骨折阈值 {THRESHOLD_MPA:.0f} MPa ([Y25])")
    ax.axvspan(Y25_LO_M, Y25_HI_M, color="green", alpha=0.10,
               label=f"[Y25] 足部骨折簇 {Y25_LO_M:.0f}–{Y25_HI_M:.0f} m")
    ax.set_xlabel("drop height h (m)")
    ax.set_ylabel("von Mises (MPa)  ·  域 `calcaneus` (CORT)")
    ax.set_title("S2-THUMS · AM50 跟骨 σ_vm(h) · 刚性地面 · 距下关节压力载荷 + 跟腱")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=1, loc="upper left")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def _write_report(
    rows: list[dict],
    acts: list[float],
    ff: dict[float, dict],
    mesh: dict,
    gauge_mm: float,
    outdir: Path,
    n_solves: int,
    *,
    tag: str = "",
    source_npz: str | None = None,
    joint_area_mm2: float | None = None,
    old_joint_area_mm2: float = 222.0,
) -> None:
    def _ff_txt(s: dict) -> str:
        if not s["triggered"]:
            return f"未触发（≤ {s['max_h_m']:.0f} m）"
        t = f"{s['h_first_m']:.0f} m"
        if s["h_interp_m"] is not None:
            t += f"（线性内插 ≈ {s['h_interp_m']:.2f} m）"
        return t

    src = source_npz or "temp/opensim_fe/thums_calcaneus/calcaneus_r_calcnframe.npz"
    tag_txt = tag.strip("_")
    title = "S2-THUMS · AM50 跟骨刚性地面全高度扫描（Option C）"
    if tag_txt:
        title += f" · {tag_txt}"
    lines: list[str] = []
    lines += [
        f"# {title}",
        "",
        f"**模型**：`{src}`（AM50 THUMS 右跟骨，mm，OpenSim `calcn_r` 帧；"
        "原点 = 距下关节面质心）。",
        "双材料域：皮质 `calcaneus`（hex8，E=15000 MPa，ν=0.30）+ 松质 `trabecular`"
        "（tet4，E=73.4 MPa，ν=0.45）。",
        "",
        f"- 单元：CORT hex8 = {mesh['n_elements_cort']}，SPON tet4 = {mesh['n_elements_spon']}；"
        f"节点 {mesh['n_nodes']}。",
        f"- 具名面：subtalar_joint {mesh['n_joint_faces']} 面 / plantar {mesh['n_plantar_faces']} 面 / "
        f"achilles {mesh['n_achilles_faces']} 面。",
    ]
    if joint_area_mm2 is not None:
        lines += [
            f"- **距下关节载荷面**：dist<25 mm 且 n·(+Y)>0.2 → {joint_area_mm2:.1f} mm²"
            f"（旧 foot-frame 配准 = {old_joint_area_mm2:.1f} mm²，14 面）。",
        ]
    lines += [
        f"- 求解：FEBio 4.13 STATIC × {n_solves} 次，**全部退出码 0**；距下关节载荷用 "
        f"`PressureLoad p=F/A`（`use_rigid=False`，与旧管线同口径），跟腱用 `TractionLoad`。",
        f"- 判据：皮质域 `gauge_max`（体平均正则化，半径 R={gauge_mm:g} mm = ℓ）> "
        f"**{THRESHOLD_MPA:.0f} MPa**（[Y25] 皮质压缩强度）。**后处理判据，无单元删除。**",
        "",
        "## 载荷链（复用 OpenSim 管线原样）",
        "",
        "```",
        "h -> ground_reaction(h) -> run_dead_drop(grf) -> subtalar_reaction(fall, grf, side='r')",
        "  -> peak_vertical_n = F_subt(h)      # 真跑正动力学",
        "  -> load_transfer.transfer(F_subt, achilles_activation=a)",
        "       subtalar_dir=(0,-1,0), achilles_dir=(0.15,0.985,0), F_MAX=10885 N",
        "```",
        "",
        "## 首次骨折高度（gauge_max > 150 MPa）",
        "",
        "| 激活 a | 跟腱力 (kN) | 首次骨折高度 | 触发点 gauge_max (MPa) |",
        "|---|---|---|---|",
    ]
    for a in acts:
        s = ff[a]
        g = "—" if s["gauge_at_first"] is None else f"{s['gauge_at_first']:.1f}"
        lines.append(f"| {a:g} | {s['achilles_kN']:.2f} | {_ff_txt(s)} | {g} |")

    lines += [
        "",
        f"**与 [Y25] 足部骨折簇（{Y25_LO_M:.0f}–{Y25_HI_M:.0f} m）对比**：",
        "",
    ]
    for a in acts:
        s = ff[a]
        if not s["triggered"]:
            rel = f"在本扫描上界（{s['max_h_m']:.0f} m）内未触发，晚于 [Y25] 簇"
        else:
            hh = s["h_interp_m"] if s["h_interp_m"] is not None else s["h_first_m"]
            if hh < Y25_LO_M:
                rel = f"≈{hh:.2f} m，**早于** [Y25] {Y25_LO_M:.0f}–{Y25_HI_M:.0f} m"
            elif hh <= Y25_HI_M:
                rel = f"≈{hh:.2f} m，**落在** [Y25] {Y25_LO_M:.0f}–{Y25_HI_M:.0f} m 内"
            else:
                rel = f"≈{hh:.2f} m，**晚于** [Y25] {Y25_LO_M:.0f}–{Y25_HI_M:.0f} m"
        lines.append(f"- **a={a:g}**：首次骨折高度 {_ff_txt(s)} → {rel}。")
    lines += [
        "",
        f"- 旧整体足网格（均匀 E=15000，平滑网格，gauge ℓ=4 mm）首次骨折高度 ≈ "
        f"**{OLD_FIRST_FRACTURE_M:.1f} m**（文档记录，`docs/OpenSim_FE交接.md`）。"
        "新模型改用真实 AM50 跟骨 + CORT/SPON 双域，绝对高度**不可直接移植**，只作量级/方向核对。",
        "",
        "## 逐点结果",
        "",
        "| h (m) | a | F_subt (kN) | F_ach (kN) | max | p99 | p95 | mean | gauge_max |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in sorted(rows, key=lambda x: (x["act"], x["h"])):
        lines.append(
            f"| {r['h']:g} | {r['act']:g} | {r['subtalar_kN']:.2f} | {r['achilles_kN']:.2f} | "
            f"{r['max']:.2f} | {r['p99']:.2f} | {r['p95']:.2f} | {r['mean']:.3f} | "
            f"{r['gauge_max']:.2f} |"
        )

    lines += [
        "",
        "## 诚实说明 / 既有 caveat（必读）",
        "",
        f"- **目标几何替换**：旧 S1–S3 的 `\"calcaneus\"` 域实际是 `r_foot.vtp` 的**整足**"
        "网格（121.7 cm³），不是单一跟骨；本模型的 `\"calcaneus\"` 是真实 AM50 跟骨"
        f"（CORT {mesh['volume_cort_mm3']:.0f} + SPON {mesh['volume_spon_mm3']:.0f} mm³）。"
        "两者数值**不可直接比较**，只能看趋势/方向。",
        "- **骨位置**：AM50 跟骨按其**解剖标志**置于 `calcn_r` 帧 —— 距下关节面质心映射到"
        "原点，跖面朝 −Y，后结节朝 −X；并非随意移动。",
        "- **距下关节载荷块**：本配准把**距下关节面质心**放到原点，故 `subtalar_joint` "
        "面（dist<25 mm 且 n·(+Y)>0.2）覆盖整个距下关节区，形成正常的载荷块"
        + (f"（面积 {joint_area_mm2:.1f} mm²，旧 foot-frame 仅 {old_joint_area_mm2:.1f} mm²）。"
           if joint_area_mm2 is not None else "。"),
        "- **峰值受边界奇异支配**：载荷面 / 跖面固定棱边处 σ_vm 随网格发散，故报告以 "
        "`gauge_max`（过程区 ℓ=4 mm 正则化）为强度判据，裸 `max` 仅作诊断。",
        "- 跟腱 `a=1` 施加 10.885 kN（Rajagopal2015 三头肌 F_max 之和的等长上界），"
        "属**参数化锚定**，非具体激活时程。",
        "- 本扫描是**筛选性**研究：绝对应力不是临床预测；σ_vm≥150 MPa 判据仅在后处理中应用，"
        "不做单元删除 / 不把断裂作为求解动作。",
        "",
        "## 产物",
        "",
        f"- `s2_thums_height_sweep{tag}.json` / `.csv`",
        f"- `s2_thums_height_sweep{tag}.png`",
        f"- `S2_THUMS_REPORT{tag.upper() if tag else ''}.md`（本文件）",
    ]
    rep_name = f"S2_THUMS_REPORT{tag.upper() if tag else ''}.md"
    (outdir / rep_name).write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="S2-THUMS 跟骨全高度扫描（刚性地面，含跟腱激活）")
    ap.add_argument("--heights", type=float, nargs="+", default=None, help="高度列表 (m)")
    ap.add_argument("--acts", type=float, nargs="+", default=None, help="三头肌激活 a∈[0,1]")
    ap.add_argument("--gauge-mm", type=float, default=REG_LEN_MM,
                    help=f"过程区正则化半径 (mm)；默认 = ℓ = {REG_LEN_MM:g} mm")
    ap.add_argument("--outdir", type=Path, default=ROOT / "results" / "opensim_fe")
    ap.add_argument("--mesh", type=Path, default=None,
                    help="替代配准网格 .npz（默认 temp/.../calcaneus_r_calcnframe.npz）；"
                         "解剖配准用 calcaneus_r_anatframe.npz")
    ap.add_argument("--tag", type=str, default="",
                    help="输出文件名后缀（如 _anat），保持旧输出不变；默认空 = 旧行为")
    args = ap.parse_args(argv)

    heights = [float(h) for h in (args.heights if args.heights else DEFAULT_GRID)]
    acts = [float(a) for a in (args.acts if args.acts else DEFAULT_ACTS)]
    for a in acts:
        if not 0.0 <= a <= 1.0:
            raise SystemExit(f"activation 必须在 [0,1]，得到 {a!r}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    workdir = ROOT / "temp" / "opensim_fe" / f"thums_sweep{args.tag}"
    workdir.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("S2-THUMS · AM50 跟骨刚性地面全高度扫描（CORT hex8 + SPON tet4）")
    print(f"  heights ({len(heights)}): {heights}")
    print(f"  acts ({len(acts)}): {acts}  (F_MAX={F_MAX_TRICEPS_N:.0f} N)")
    print(f"  gauge R={args.gauge_mm:g} mm   threshold={THRESHOLD_MPA:.0f} MPa")
    print("=" * 72, flush=True)

    t0 = time.time()
    mesh = T.load_thums_mesh(args.mesh) if args.mesh else T.load_thums_mesh()
    print(f"[mesh] source={mesh['source_npz']}")
    print(f"[mesh] 节点 {mesh['n_nodes']}  CORT hex8={mesh['n_elements_cort']}  "
          f"SPON tet4={mesh['n_elements_spon']}  "
          f"joint={mesh['n_joint_faces']}f  plantar={mesh['n_plantar_faces']}f  "
          f"achilles={mesh['n_achilles_faces']}f", flush=True)
    cort_centroids, cort_vol = T.hex_centroids_volumes(mesh["nodes"], mesh["elements_cort"])

    rows: list[dict] = []
    n_solves = 0
    for h in heights:
        f_sub, t1 = _load_at(h)
        for a in acts:
            row = _solve_case(
                mesh, cort_centroids, cort_vol, h, a, f_sub, t1, args.gauge_mm, workdir
            )
            rows.append(row)
            n_solves += 1
            print(f"h={h:5.1f} m  a={a:g}  F_sub={row['subtalar_kN']:6.2f} kN  "
                  f"F_ach={row['achilles_kN']:6.2f} kN  max={row['max']:7.2f}  "
                  f"p99={row['p99']:6.2f}  p95={row['p95']:6.2f}  "
                  f"gauge(R{args.gauge_mm:g})={row['gauge_max']:7.2f} MPa  "
                  f"[rc=0, T1={row['t1_rel_err']*100:+.1f}%, {time.time()-t0:.0f}s]",
                  flush=True)

    ff = _first_fracture(rows, acts)

    tag = args.tag
    json_out = args.outdir / f"s2_thums_height_sweep{tag}.json"
    csv_out = args.outdir / f"s2_thums_height_sweep{tag}.csv"
    png_out = args.outdir / f"s2_thums_height_sweep{tag}.png"
    rep_out = args.outdir / f"S2_THUMS_REPORT{tag.upper() if tag else ''}.md"

    rows_sorted = sorted(rows, key=lambda r: (r["act"], r["h"]))
    json_out.write_text(
        json.dumps(rows_sorted, indent=2, ensure_ascii=False, default=float),
        encoding="utf-8",
    )
    with csv_out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in rows_sorted:
            w.writerow({k: r[k] for k in CSV_FIELDS})

    # joint-patch area of the actual loaded mesh (exact FE-builder rule)
    from climbing.coupling.febio_model import joint_area_mm2
    joint_area = float(joint_area_mm2(mesh, "subtalar_joint"))

    _write_plot(rows, acts, png_out)
    _write_report(rows, acts, ff, mesh, args.gauge_mm, args.outdir, n_solves,
                  tag=tag, source_npz=mesh["source_npz"],
                  joint_area_mm2=joint_area, old_joint_area_mm2=222.0)

    print("\n" + "=" * 72)
    print("S2-THUMS 报告：首次骨折高度（gauge_max > 150 MPa，域 `calcaneus`/CORT）")
    for a in acts:
        s = ff[a]
        if s["triggered"]:
            extra = "" if s["h_interp_m"] is None else f"  (内插≈{s['h_interp_m']:.2f} m)"
            print(f"  a={a:g} (F_ach={s['achilles_kN']:.2f} kN): {s['h_first_m']:.0f} m{extra}"
                  f"   gauge@first={s['gauge_at_first']:.1f} MPa")
        else:
            print(f"  a={a:g} (F_ach={s['achilles_kN']:.2f} kN): 未触发（≤ {s['max_h_m']:.0f} m）")
    print(f"[Y25] 足部骨折簇 = {Y25_LO_M:.0f}–{Y25_HI_M:.0f} m；"
          f"旧整体足网格首次骨折 ≈ {OLD_FIRST_FRACTURE_M:.1f} m（文档值，不可直接移植）。")
    print(f"[meta] FEBio 求解次数 = {n_solves}（全部 rc=0）  "
          f"gauge R={args.gauge_mm:g} mm  threshold={THRESHOLD_MPA:.0f} MPa")
    print(f"[out] {json_out}")
    print(f"[out] {csv_out}")
    print(f"[out] {png_out}")
    print(f"[out] {rep_out}")
    print(f"[done] {time.time()-t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
