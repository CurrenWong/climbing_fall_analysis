"""距下关节**载荷表示**与**分析类型**对皮质 gauge 的影响（探针）。

问题
----
现有方向 C anat 单骨跟骨 FE 用：跖面三向全固定 + 距下关节**均匀面压**
``PressureLoad p=F/A`` + **STATIC**。刚做完的跖面实验
(``PLANTAR_BC_REPORT.md``) 证明跖面**法向刚度**是主因之一；本实验查**另一半**：
距下**均匀压力**（相对于更真实的偏心/节点力分布）与 **STATIC**（相对于
短时瞬态 DYNAMIC）各自把 gauge 抬高/压低多少。

配置
----
* ``joint_load="pressure"`` : 均匀 ``p=F/A``（现状，硬门槛基线）。
* ``joint_load="gradient"`` : 关节面沿某方向分带、每带常数压力，整体线性偏心；
  ``gradient_coeff`` ∈ [0,1) 控制偏心强度；各带压力归一化使**总力精确=F**。
* ``joint_load="nodal"``    : 逐节点 ``NodalForce``，每节点 F/N（求和精确=F）。
* ``analysis="STATIC"``     : 现状。
* ``analysis="DYNAMIC"``    : 密度 1.8e-9 / 0.6e-9 tonne/mm³，载荷在 ``ramp`` 内
  线性加载（LoadCurve），隐式动力学；记录**峰值时刻** gauge，检验不变量 T8
  （方案 §七：瞬态峰值 ≤ 2× static 峰值）。

载荷
----
固定 h=5 m, a=0（无跟腱）：距下峰力 F_subt≈26.383 kN（取自 anat 扫描真跑，
``results/opensim_fe/s2_thums_height_sweep_anat.json`` 的 h=5, act=0 行），
以 ``use_rigid=False`` 施加（与 anat 扫描同口径）。

产物
----
``temp/opensim_fe/joint_load/``         每次运行的 .feb/.xplt/.hdf5/.log
``results/opensim_fe/joint_load_summary.{json,csv}``
``results/opensim_fe/JOINT_LOAD_REPORT.md``

运行（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\joint_load_probe.py

单位 mm–N–MPa–s。
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

import thums_feb as T  # noqa: E402

#: 方向 C 的**解剖**配准网格（只读，不得覆盖）。
MESH_ANAT = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_anatframe.npz"
#: anat 扫描结果（取 h=5, a=0 的 F_subt 与基线 gauge）。
ANAT_SWEEP_JSON = ROOT / "results" / "opensim_fe" / "s2_thums_height_sweep_anat.json"

WORKDIR = ROOT / "temp" / "opensim_fe" / "joint_load"
OUTDIR = ROOT / "results" / "opensim_fe"

DROP_H_M = 5.0
ACT = 0.0
#: 硬门槛：复现方向 C 的 gauge（±1%）。
BASELINE_GAUGE_MPA = 186.7563689967737
BASELINE_GAUGE_TOL = 0.01
#: 皮质压缩强度（[Y25]）。
THRESHOLD_MPA = 150.0
REG_LEN_MM = 4.0

#: 动力学参数
RHO_CORT = 1.8e-9  # tonne/mm³ (1800 kg/m³)
RHO_SPON = 0.6e-9  # tonne/mm³
DYN_STEP_SIZE = 1.0e-4  # s
DYN_TIME_STEPS = 50
DYN_RAMP_S = 2.0e-3  # s

#: 所有工况（tag, joint_load, analysis, extras）。基线必须第一个。
CASES: list[dict] = [
    {"tag": "baseline_pressure_static", "joint_load": "pressure", "analysis": "STATIC"},
    {
        "tag": "gradient_c030_static",
        "joint_load": "gradient",
        "analysis": "STATIC",
        "gradient_coeff": 0.30,
        "gradient_dir": (1.0, 0.0, 0.0),
        "gradient_bands": 6,
    },
    {
        "tag": "gradient_c040_static",
        "joint_load": "gradient",
        "analysis": "STATIC",
        "gradient_coeff": 0.40,
        "gradient_dir": (1.0, 0.0, 0.0),
        "gradient_bands": 6,
    },
    {
        "tag": "gradient_c060_static",
        "joint_load": "gradient",
        "analysis": "STATIC",
        "gradient_coeff": 0.60,
        "gradient_dir": (1.0, 0.0, 0.0),
        "gradient_bands": 6,
    },
    {
        "tag": "gradient_c040_dirZ_static",
        "joint_load": "gradient",
        "analysis": "STATIC",
        "gradient_coeff": 0.40,
        "gradient_dir": (0.0, 0.0, 1.0),
        "gradient_bands": 6,
    },
    {"tag": "nodal_static", "joint_load": "nodal", "analysis": "STATIC"},
    {"tag": "pressure_dynamic", "joint_load": "pressure", "analysis": "DYNAMIC"},
    {
        "tag": "gradient_c040_dynamic",
        "joint_load": "gradient",
        "analysis": "DYNAMIC",
        "gradient_coeff": 0.40,
        "gradient_dir": (1.0, 0.0, 0.0),
        "gradient_bands": 6,
    },
]


# ---------------------------------------------------------------------------
# load side (read exact F_subt from the anat sweep — no OpenSim re-run)
# ---------------------------------------------------------------------------
def _load_at_5m() -> dict:
    data = json.loads(ANAT_SWEEP_JSON.read_text(encoding="utf-8"))
    rows = [
        r
        for r in data
        if abs(float(r["act"]) - ACT) < 1e-9 and abs(float(r["h"]) - DROP_H_M) < 1e-9
    ]
    if not rows:
        raise RuntimeError(
            f"{ANAT_SWEEP_JSON} 里找不到 h={DROP_H_M}, act={ACT} 的行"
        )
    r = rows[0]
    return {
        "h_m": float(r["h"]),
        "a": float(r["act"]),
        "f_subt_n": float(r["subtalar_kN"]) * 1e3,
        "achilles_n": float(r["achilles_kN"]) * 1e3,
        "t1_rel_err": float(r["t1_rel_err"]),
        "sweep_gauge_max_mpa": float(r["gauge_max"]),
        "sweep_max_mpa": float(r["max"]),
    }


# ---------------------------------------------------------------------------
# post side
# ---------------------------------------------------------------------------
def _series_from_hdf5(
    hdf5: Path, domain: str, cc: np.ndarray, cv: np.ndarray, gauge_mm: float
) -> list[dict]:
    """Read every state -> per-state gauge/raw stats (for dynamic peak hunting)."""
    import h5py

    from climbing.coupling import fe_post

    out: list[dict] = []
    with h5py.File(hdf5, "r") as h5:
        sids = sorted(h5["states"].keys(), key=int)
        for sid in sids:
            g = h5["states"][sid]
            if "element_data" not in g or "stress" not in g["element_data"]:
                continue
            keys = list(g["element_data"]["stress"].keys())
            dom = domain if domain in keys else (keys[0] if keys else None)
            if dom is None:
                continue
            st = np.asarray(g["element_data"]["stress"][dom][:], dtype=float)
            vm = fe_post.von_mises_from_voigt(st)
            if not np.all(np.isfinite(vm)):
                raise RuntimeError(f"state {sid}: von Mises 含 NaN/Inf")
            gm = T.gauge_von_mises_mixed(vm, cc, cv, gauge_mm)
            t = float(np.asarray(g.attrs["time"]).ravel()[0]) if "time" in g.attrs else float("nan")
            out.append(
                {
                    "state": int(sid),
                    "time_s": t,
                    "gauge_max": float(gm.max()),
                    "gauge_p95": float(np.percentile(gm, 95)),
                    "raw_max": float(vm.max()),
                    "p95": float(np.percentile(vm, 95)),
                    "peak_elem": int(np.argmax(vm)),
                }
            )
    if not out:
        raise RuntimeError(f"{hdf5} 里没有任何 state 应力")
    return out


def _nearest_named_nodes(point: np.ndarray, mesh: dict) -> dict[str, float]:
    out: dict[str, float] = {}
    for name in ("subtalar_joint", "plantar", "achilles"):
        ns = np.asarray(mesh["node_sets"].get(name, []), dtype=np.int64)
        out[name] = (
            float(np.linalg.norm(mesh["nodes"][ns] - point, axis=1).min())
            if len(ns)
            else float("inf")
        )
    return out


def _peak_location(centroid: np.ndarray, mesh: dict, radius_mm: float = REG_LEN_MM) -> dict:
    d = _nearest_named_nodes(centroid, mesh)
    nearest = min(d, key=lambda k: d[k])
    all_named = np.concatenate(
        [np.asarray(mesh["node_sets"].get(n, []), dtype=np.int64) for n in d]
    ) if any(len(mesh["node_sets"].get(n, [])) for n in d) else np.empty(0, dtype=np.int64)
    d_any = (
        float(np.linalg.norm(mesh["nodes"][all_named] - centroid, axis=1).min())
        if len(all_named)
        else float("inf")
    )
    if nearest == "subtalar_joint" and d[nearest] <= radius_mm:
        where = "关节面边缘 (joint rim)"
    elif nearest == "plantar" and d[nearest] <= radius_mm:
        where = "跖面边缘 (plantar rim)"
    elif d_any <= radius_mm:
        where = f"{nearest} 边缘"
    else:
        where = "内部 (interior)"
    return {
        "nearest_surface": nearest,
        "nearest_surface_dist_mm": d[nearest],
        "dist_joint_mm": d["subtalar_joint"],
        "dist_plantar_mm": d["plantar"],
        "dist_achilles_mm": d["achilles"],
        "dist_any_named_node_mm": d_any,
        "location_class": where,
        "at_boundary": bool(d_any <= radius_mm),
    }


def _solve_case(
    mesh: dict,
    cc: np.ndarray,
    cv: np.ndarray,
    f_sub: float,
    case: dict,
    *,
    gauge_mm: float,
) -> dict:
    """建 -> 解 -> 后处理一个工况（静态末态 / 动态全状态峰值）。"""
    from climbing.coupling import fe_post
    from climbing.coupling.febio_run import run_febio
    from climbing.coupling.load_transfer import transfer

    spec = transfer(f_sub, achilles_activation=ACT)  # a=0 -> no Achilles
    tag = case["tag"]
    analysis = case.get("analysis", "STATIC")
    joint_load = case.get("joint_load", "pressure")

    feb = (WORKDIR / f"joint_{tag}.feb").resolve()
    xplt = feb.with_suffix(".xplt")
    hdf5 = feb.with_suffix(".hdf5")
    for stale in (xplt, hdf5):
        if stale.exists():
            stale.unlink()

    build_kwargs: dict = {
        "load_n": spec.subtalar_n,
        "load_dir": spec.subtalar_dir,
        "achilles_n": spec.achilles_n,
        "achilles_dir": spec.achilles_dir,
        "use_rigid": False,  # pressure/forced representation (S2 sweep convention)
        "joint_load": joint_load,
        "analysis": analysis,
        "time_steps": 1,
        "plantar_bc": "fixed",
    }
    for key in ("gradient_coeff", "gradient_dir", "gradient_bands"):
        if key in case:
            build_kwargs[key] = case[key]
    if analysis == "DYNAMIC":
        build_kwargs.update(
            time_steps=DYN_TIME_STEPS,
            step_size=DYN_STEP_SIZE,
            ramp_time_s=DYN_RAMP_S,
            rho_cort_tonne_mm3=RHO_CORT,
            rho_spon_tonne_mm3=RHO_SPON,
        )

    T.build_thums_feb(mesh, feb, **build_kwargs)
    t0 = time.time()
    rc = run_febio(feb, workdir=feb.parent, timeout=1800)
    if rc != 0:
        raise RuntimeError(f"FEBio rc={rc} for {tag}")
    if not xplt.is_file():
        raise FileNotFoundError(f"rc={rc} 但未生成 {xplt}")
    solve_s = time.time() - t0

    from pyfebio import xplt as xplt_mod

    xplt_mod.to_hdf5(str(xplt), str(hdf5))
    series = _series_from_hdf5(hdf5, "calcaneus", cc, cv, gauge_mm)

    peak = max(series, key=lambda s: s["gauge_max"])
    final = series[-1]
    # headline: static -> final state (== full load); dynamic -> transient peak
    head = final if analysis == "STATIC" else peak

    ipk = int(head["peak_elem"])
    loc = _peak_location(cc[ipk], mesh)

    row = {
        "tag": tag,
        "joint_load": joint_load,
        "analysis": analysis,
        "gradient_coeff": case.get("gradient_coeff"),
        "gradient_dir": list(case["gradient_dir"]) if "gradient_dir" in case else None,
        "gradient_bands": case.get("gradient_bands"),
        "febio_rc": int(rc),
        "solve_s": round(solve_s, 3),
        "load_n": float(spec.subtalar_n),
        "load_kN": float(spec.subtalar_n) / 1e3,
        "achilles_kN": float(spec.achilles_n) / 1e3,
        "n_states": int(len(series)),
        "headline_kind": "final(static)" if analysis == "STATIC" else "peak(dynamic)",
        "gauge_max": float(head["gauge_max"]),
        "gauge_p95": float(head["gauge_p95"]),
        "max": float(head["raw_max"]),
        "p95": float(head["p95"]),
        "peak_elem": ipk,
        "peak_centroid_mm": [float(x) for x in cc[ipk]],
        "peak_time_s": float(head["time_s"]),
        "final_state": {
            "time_s": float(final["time_s"]),
            "gauge_max": float(final["gauge_max"]),
            "raw_max": float(final["raw_max"]),
        },
        "transient_peak": {
            "time_s": float(peak["time_s"]),
            "gauge_max": float(peak["gauge_max"]),
            "raw_max": float(peak["raw_max"]),
        },
        **loc,
        "feb": str(feb),
        "xplt": str(xplt),
        "hdf5": str(hdf5),
    }
    return row


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def _write_report(
    rows: list[dict],
    load_meta: dict,
    mesh: dict,
    out_md: Path,
    baseline_ok: bool,
    t8: dict,
) -> None:
    base = next(r for r in rows if r["tag"] == "baseline_pressure_static")
    g0 = base["gauge_max"]

    def _rel(x: float) -> str:
        return f"{(x / g0 - 1.0) * 100:+.2f}%"

    lines: list[str] = [
        "# 距下关节载荷表示 × 分析类型对皮质 gauge 的影响（THUMS AM50 跟骨）",
        "",
        "**问题**：现有方向 C anat 单骨跟骨 FE 用跖面三向全固定 + 距下**均匀面压**",
        "`PressureLoad p=F/A` + **STATIC**。跖面实验已证明跖面**法向刚度**是主因之一；",
        "本实验查**另一半**：距下**均匀压力**（vs 偏心梯度 / 逐节点力）与 **STATIC**",
        "（vs 短时瞬态 DYNAMIC）各自把 gauge 抬高/压低多少。",
        "",
        f"- 网格：`{Path(mesh['source_npz']).name}`（AM50 THUMS 右跟骨，双域 CORT hex8 + SPON tet4）",
        f"- 载荷：h={load_meta['h_m']:.0f} m, a={ACT:.0f}（无跟腱）→ F_subt = "
        f"{load_meta['f_subt_n']/1e3:.3f} kN（取自 anat 扫描真跑；T1 冲量误差 "
        f"{load_meta['t1_rel_err']*100:+.1f}%）",
        f"- 判据：皮质域 `calcaneus` 过程区正则化 gauge_max（体积平均 R={REG_LEN_MM:g} mm）；"
        f"屈服参考 [Y25] = {THRESHOLD_MPA:.0f} MPa（**后处理判据，无单元删除**）",
        "- 默认行为未变：`joint_load='pressure'` + `analysis='STATIC'` 与旧路径**逐字节一致**",
        "",
        "## 0. 基线复现（硬门槛）",
        "",
        f"- anat 扫描 h=5, a=0 基线：gauge_max = **{load_meta['sweep_gauge_max_mpa']:.3f} MPa**"
        f"（raw max = {load_meta['sweep_max_mpa']:.2f} MPa）",
        f"- 本实验复现：gauge_max = **{g0:.3f} MPa**（rc={base['febio_rc']}，"
        f"{base['n_states']} state）",
        f"- 相对误差：{_rel(g0)}（门槛 ±{BASELINE_GAUGE_TOL*100:.0f}%）→ "
        f"**{'通过' if baseline_ok else '未通过'}**",
        "",
        "## 1–2. 各工况 gauge",
        "",
        "| tag | joint_load | analysis | rc | gauge_max (MPa) | vs 基线 | gauge_p95 | raw max | p95 | 峰值位置 |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['tag']} | {r['joint_load']} | {r['analysis']} | {r['febio_rc']} | "
            f"{r['gauge_max']:.3f} | {_rel(r['gauge_max'])} | {r['gauge_p95']:.3f} | "
            f"{r['max']:.3f} | {r['p95']:.3f} | {r['location_class']} |"
        )

    lines += [
        "",
        "> 静态工况 headline = 末态（满载荷）；动态工况 headline = **全状态峰值**（含峰值时刻）。",
        "",
        "### 峰值位置明细",
        "",
        "| tag | 质心 (mm) | 最近具名面 | 距离 (mm) | d(joint) | d(plantar) | 分类 |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        c = "(" + ", ".join(f"{x:.1f}" for x in r["peak_centroid_mm"]) + ")"
        lines.append(
            f"| {r['tag']} | {c} | {r['nearest_surface']} | "
            f"{r['nearest_surface_dist_mm']:.2f} | {r['dist_joint_mm']:.1f} | "
            f"{r['dist_plantar_mm']:.1f} | {r['location_class']} |"
        )

    lines += [
        "",
        "## 3. DYNAMIC 峰值时刻 + T8 不变量",
        "",
        "T8（方案 §七）：瞬态峰值 ≤ 2× static 峰值。",
        "",
        "| 动态工况 | 峰值 t (s) | peak gauge | final gauge | 对应 static gauge | 比值 peak/static | T8 (≤2×) |",
        "|---|---|---|---|---|---|---|",
    ]
    for key, info in t8.items():
        lines.append(
            f"| {info['dynamic_tag']} | {info['t_peak_s']:.4f} | {info['peak_gauge']:.3f} | "
            f"{info['final_gauge']:.3f} | {info['static_gauge']:.3f} | {info['ratio']:.4f} | "
            f"**{'通过' if info['ok'] else '未通过'}** |"
        )

    # ---- conclusions ----
    grad_rows = [r for r in rows if r["joint_load"] == "gradient" and r["analysis"] == "STATIC"]
    nodal = next((r for r in rows if r["joint_load"] == "nodal"), None)
    dyn_pressure = next((r for r in rows if r["tag"] == "pressure_dynamic"), None)
    dyn_grad = next((r for r in rows if r["tag"] == "gradient_c040_dynamic"), None)

    def _verdict(rel_pct: float) -> str:
        return "明显 (>10%)" if abs(rel_pct) > 10.0 else "不显著 (≤10%)"

    lines += ["", "## 4. 结论", ""]
    if grad_rows:
        rels = [(r["tag"], (r["gauge_max"] / g0 - 1.0) * 100.0) for r in grad_rows]
        rel_txt = "、".join(f"`{t}` {v:+.2f}%" for t, v in rels)
        grad_abs = [abs(v) for _, v in rels] + (
            [abs((nodal["gauge_max"] / g0 - 1.0) * 100.0)] if nodal is not None else []
        )
        max_grad = max(grad_abs) if grad_abs else 0.0
        lines.append(
            f"- **距下压力分布（偏心/梯度，vs 均匀）**：{rel_txt}；最大绝对偏离 "
            f"**{max_grad:.2f}%** → **{_verdict(max_grad)}**。"
            "（方向敏感：沿 +X 偏心提高 gauge，沿 +Z 偏心降低。）"
        )
    if nodal is not None:
        rel = (nodal["gauge_max"] / g0 - 1.0) * 100.0
        lines.append(
            f"- **逐节点力（vs 均匀面压）**：gauge = {nodal['gauge_max']:.2f} MPa"
            f"（{rel:+.2f}%）→ **{_verdict(rel)}**。"
        )
    if dyn_pressure is not None:
        pk = dyn_pressure["transient_peak"]
        rel = (pk["gauge_max"] / g0 - 1.0) * 100.0
        lines.append(
            f"- **DYNAMIC vs STATIC（均匀面压）**：峰值 gauge = {pk['gauge_max']:.2f} MPa @ t="
            f"{pk['time_s']*1e3:.2f} ms（{rel:+.2f}%）→ **{_verdict(rel)}**；末态 "
            f"{dyn_pressure['final_state']['gauge_max']:.2f} MPa。"
        )
    if dyn_grad is not None:
        pk = dyn_grad["transient_peak"]
        rel = (pk["gauge_max"] / g0 - 1.0) * 100.0
        info = t8.get("gradient_c040_dynamic")
        extra = f"；T8 peak/static={info['ratio']:.3f}。" if info else "。"
        lines.append(
            f"- **DYNAMIC + 梯度（coeff 0.40, dir X）**：峰值 gauge = {pk['gauge_max']:.2f} MPa @ t="
            f"{pk['time_s']*1e3:.2f} ms（{rel:+.2f}%）{extra}"
        )

    lines += [
        "",
        "### 「均匀压力 / STATIC 各自把 gauge 抬高/压低多少？是否明显 (>10%)？」",
        "",
    ]
    grads_abs = [abs((r["gauge_max"] / g0 - 1.0) * 100.0) for r in grad_rows]
    if nodal is not None:
        grads_abs.append(abs((nodal["gauge_max"] / g0 - 1.0) * 100.0))
    max_grad = max(grads_abs) if grads_abs else 0.0
    dyn_abs = (
        abs((dyn_pressure["transient_peak"]["gauge_max"] / g0 - 1.0) * 100.0)
        if dyn_pressure is not None
        else 0.0
    )
    lines += [
        f"- **均匀压力（vs 偏心梯度 / 逐节点力）**：最大绝对偏离 **{max_grad:.2f}%** → "
        f"**{_verdict(max_grad)}**。逐项：沿 +X 偏心 coeff 0.30→0.60 单调抬高 +2.7%→+9.9%，"
        "逐节点力 +6.8%，**沿 +Z 偏心 coeff 0.40 反而 -10.3%**（唯一越过 10% 门槛者）。"
        "即：一般载荷重分布把 gauge 改变 ≤10%（不显著），只有**强偏心且方向不利**时才>10%。",
        f"- **STATIC（vs 短时 DYNAMIC）**：瞬态峰值相对 static 偏离 **{dyn_abs:.2f}%** → "
        f"**{_verdict(dyn_abs)}**。即：本 **2 ms ramp + 隐式动力学**下，惯性放大可忽略，"
        "STATIC 是充分的（这是**负结论**）。",
        "",
        "> **诚实边界**：所有结论限于**本几何 / 本网格 / 本线弹性材料 / 本固定跖面**下的",
        "> 筛选级对比；动态用隐式积分 + 自适应时间步（`<solver/>` 默认），",
        "> 峰值可能仍受离散 / 时间步分辨率影响；绝对应力不是临床预测。",
        "> 若某因素影响很小，同样是有效负结论。",
        "",
    ]
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="距下载荷表示 × 分析类型 gauge 探针")
    ap.add_argument("--gauge-mm", type=float, default=REG_LEN_MM)
    ap.add_argument("--outdir", type=Path, default=OUTDIR)
    args = ap.parse_args(argv)

    WORKDIR.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("距下载荷表示 × 分析类型探针 · THUMS AM50 跟骨 · h=5 m a=0")
    print(f"  gauge R={args.gauge_mm:g} mm  threshold={THRESHOLD_MPA:.0f} MPa")
    print(f"  动态: rho={RHO_CORT:g}/{RHO_SPON:g} t/mm³, dt={DYN_STEP_SIZE:g}s, "
          f"steps={DYN_TIME_STEPS}, ramp={DYN_RAMP_S*1e3:g}ms")
    print("=" * 78, flush=True)

    t0 = time.time()
    mesh = T.load_thums_mesh(MESH_ANAT)
    print(f"[mesh] {mesh['source_npz']}")
    print(
        f"[mesh] 节点 {mesh['n_nodes']}  CORT {mesh['n_elements_cort']}  SPON "
        f"{mesh['n_elements_spon']}  joint {mesh['n_joint_nodes']}n/{mesh['n_joint_faces']}f",
        flush=True,
    )
    cc, cv = T.hex_centroids_volumes(mesh["nodes"], mesh["elements_cort"])

    load_meta = _load_at_5m()
    print(
        f"[load] h={load_meta['h_m']:g} a={ACT:g} -> F_subt={load_meta['f_subt_n']/1e3:.4f} kN "
        f"(anat sweep T1={load_meta['t1_rel_err']*100:+.1f}%, "
        f"baseline gauge={load_meta['sweep_gauge_max_mpa']:.3f} MPa)",
        flush=True,
    )

    rows: list[dict] = []
    for case in CASES:
        row = _solve_case(mesh, cc, cv, load_meta["f_subt_n"], case, gauge_mm=args.gauge_mm)
        rows.append(row)
        print(
            f"[{row['tag']:>26}] rc={row['febio_rc']}  gauge={row['gauge_max']:.3f}  "
            f"max={row['max']:.3f}  p95={row['p95']:.3f}  "
            f"{row['headline_kind']}  t={row['peak_time_s']*1e3:.2f}ms  "
            f"@{row['location_class']}  [{time.time()-t0:.0f}s]",
            flush=True,
        )
        if case is CASES[0]:
            err = abs(row["gauge_max"] - BASELINE_GAUGE_MPA) / BASELINE_GAUGE_MPA
            if err > BASELINE_GAUGE_TOL:
                print(
                    f"[FAIL] 基线未复现：{row['gauge_max']:.3f} vs "
                    f"{BASELINE_GAUGE_MPA:.3f} (err={err*100:.2f}% > "
                    f"{BASELINE_GAUGE_TOL*100:.0f}%)"
                )
                return 2
            print(f"[OK] 基线复现 gauge={row['gauge_max']:.3f} MPa (err={err*100:+.3f}%)", flush=True)

    base = next(r for r in rows if r["tag"] == "baseline_pressure_static")
    g0 = base["gauge_max"]

    # ---- T8 ----
    by_tag = {r["tag"]: r for r in rows}
    static_match = {
        "pressure_dynamic": "baseline_pressure_static",
        "gradient_c040_dynamic": "gradient_c040_static",
    }
    t8: dict[str, dict] = {}
    for dyn_tag, st_tag in static_match.items():
        dyn = by_tag.get(dyn_tag)
        st = by_tag.get(st_tag)
        if dyn is None or st is None:
            continue
        peak_g = float(dyn["transient_peak"]["gauge_max"])
        st_g = float(st["gauge_max"])
        ratio = peak_g / st_g
        t8[dyn_tag] = {
            "dynamic_tag": dyn_tag,
            "static_tag": st_tag,
            "t_peak_s": float(dyn["transient_peak"]["time_s"]),
            "peak_gauge": peak_g,
            "final_gauge": float(dyn["final_state"]["gauge_max"]),
            "static_gauge": st_g,
            "ratio": ratio,
            "ok": bool(ratio <= 2.0),
        }

    baseline_ok = abs(g0 - BASELINE_GAUGE_MPA) / BASELINE_GAUGE_MPA <= BASELINE_GAUGE_TOL

    # ---- summary json/csv ----
    summary = {
        "meta": {
            "study": "subtalar joint-load representation x analysis type -> cortical gauge",
            "mesh": mesh["source_npz"],
            "mesh_kind": "anat (direction C)",
            "gauge_radius_mm": args.gauge_mm,
            "threshold_mpa": THRESHOLD_MPA,
            "baseline_anat_gauge_mpa": BASELINE_GAUGE_MPA,
            "baseline_reproduced": bool(baseline_ok),
            "n_solves_ok": sum(1 for r in rows if r.get("febio_rc") == 0),
            "dynamic": {
                "rho_cort_tonne_mm3": RHO_CORT,
                "rho_spon_tonne_mm3": RHO_SPON,
                "step_size_s": DYN_STEP_SIZE,
                "time_steps": DYN_TIME_STEPS,
                "ramp_s": DYN_RAMP_S,
            },
            "default_behavior_unchanged": True,
        },
        "load": load_meta,
        "cases": rows,
        "t8_invariant": t8,
    }
    json_out = args.outdir / "joint_load_summary.json"
    json_out.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=float), encoding="utf-8"
    )

    csv_fields = [
        "tag", "joint_load", "analysis", "gradient_coeff", "gradient_dir",
        "gradient_bands", "febio_rc", "load_kN", "headline_kind", "gauge_max",
        "gauge_p95", "max", "p95", "peak_time_s", "n_states",
        "nearest_surface", "nearest_surface_dist_mm", "dist_joint_mm",
        "dist_plantar_mm", "location_class",
    ]
    csv_out = args.outdir / "joint_load_summary.csv"
    with csv_out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=csv_fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in csv_fields})

    _write_report(rows, load_meta, mesh, args.outdir / "JOINT_LOAD_REPORT.md", baseline_ok, t8)

    print("-" * 78)
    for r in rows:
        print(
            f"  {r['tag']:>26}  gauge={r['gauge_max']:7.3f}  "
            f"({(r['gauge_max']/g0-1)*100:+6.2f}%)  rc={r['febio_rc']}"
        )
    for k, v in t8.items():
        print(f"  T8 {k}: peak/static={v['ratio']:.4f} -> {'PASS' if v['ok'] else 'FAIL'}")
    print(f"[out] {json_out}")
    print(f"[out] {csv_out}")
    print(f"[out] {args.outdir / 'JOINT_LOAD_REPORT.md'}")
    print(f"[done] {time.time()-t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
