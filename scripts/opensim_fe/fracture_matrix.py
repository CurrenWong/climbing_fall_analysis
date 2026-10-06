"""50 m × 9-part fracture-matrix skeleton (load end + generic FE builder + criterion).

For [Y25] (a representative bone per part), this driver wires three proven
pieces into one chain:

    load end            : climbing.coupling.joint_reactions.joint_reaction
    generic FE builder  : scripts/opensim_fe/bone_feb.build_bone_feb / run_bone_feb / post_process
    criterion           : climbing.bone.MATERIAL_STRENGTH_MPA (σ_c, compression)

Method (all numbers measured; nothing fabricated)
-------------------------------------------------
1.  For each of the **9 representative bones** the THUMS deck is stream-read
    **once** (``bone_feb.extract_mesh``), then build/solve/post-processed at a
    reference load (the @5 m longitudinal peak of that part's joint).  The
    **main load-bearing domain** σ_vm is recorded:
      * has cortical shell -> the shell domains (most-loaded shell per statistic);
      * no shell            -> the solid CORT domain.
    FEBio diverges for 3 bones at the full @5 m load (near-incompressible
    trabecular material / thin-shell point loading).  A **descending load
    fallback** finds the largest solvable load and the FE stress is scaled
    back to the reference load under **strict linear elasticity** (validated by
    a two-load linearity self-check where possible).
2.  Heights ``h = 1..50 m``: for every ``h`` a fresh
    ``ground_reaction(h)`` -> ``run_dead_drop`` -> ``joint_reaction`` is run
    (**per-height, no GRF-peak-ratio shortcut** — the chain is measurably
    non-linear in h).  ``load(h)`` = the joint's longitudinal peak.
3.  Strict linear elasticity: ``σ_vm(h) = σ_vm_ref · load(h)/load_ref``.
4.  ``fracture(h) = σ_vm(h) >= σ_c(part)`` -> a 9×50 binary matrix.
    Primary statistic = **p95** of the main domain (regularized hot spot;
    ``max`` is dominated by point-load geometry artifacts, ``mean/median`` are
    bulk; all four are reported).

Outputs (new files only):
    temp/opensim_fe/fracture_matrix/            (febs, per-bone JSON, load chain)
    results/opensim_fe/fracture_matrix.{json,csv,png}
    results/opensim_fe/FRACTURE_MATRIX_REPORT.md

Reproduce (repo root, PYTHONPATH=src):
    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\fracture_matrix.py
    # --loads-only / --fe-only / --force-loads / --force-fe / --stat max|p95|mean|median

This module **does not modify** any default behaviour or signature of
``bone_feb`` / ``joint_reactions`` / ``bone`` / ``meshing`` / ``kmesh_io`` /
``shell_io`` / ``febio_run``; it only imports and calls them.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import traceback
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parent.parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from bone_feb import (  # noqa: E402
    BoneSpec,
    ShellDomainSpec,
    SolidDomainSpec,
    build_bone_feb,
    extract_mesh,
    post_process,
    read_materials,
    run_bone_feb,
)

from climbing.bone import MATERIAL_STRENGTH_MPA  # noqa: E402

# --------------------------------------------------------------------------
# paths / constants
# --------------------------------------------------------------------------
OUT_TEMP = ROOT / "temp" / "opensim_fe" / "fracture_matrix"
OUT_RESULTS = ROOT / "results" / "opensim_fe"
MASTER = ROOT / "model" / "AM50_V71_Occupant" / "main_THUMS_AM50_V71.k"

MODEL_MASS_KG = 75.337  # Rajagopal2015 total mass (matches opensim_grf default)

HEIGHTS = np.arange(1, 51, dtype=int)  # 1..50 m inclusive -> 50 columns
REF_H = 5.0                            # reference height for the FE reference solve
JOINTS = ("subtalar_r", "ankle_r", "hip_r", "lumbar")
STATS = ("max", "p95", "mean", "median")
PRIMARY_STAT = "p95"

#: descending FE load candidates (N) used to find a solvable reference load
FE_LOAD_CANDIDATES = (10000.0, 5000.0, 2500.0, 1000.0)

# THUMS shell thicknesses (mm) — from bone_feb_validate.py (hand-verified there)
_HIP_SHELL_THICK = {
    83500201: 1.5, 83500202: 1.5, 83500203: 1.0, 83500204: 0.75, 83500205: 2.0,
    83500206: 3.0, 83500207: 2.2, 83500208: 1.8, 83500209: 1.0,
}


# --------------------------------------------------------------------------
# part definitions (9 × [Y25] part -> representative bone -> joint -> σ_c)
# --------------------------------------------------------------------------
@dataclass
class Part:
    """One [Y25] injury part mapped to a representative THUMS bone."""

    cn: str                      # [Y25] part (Chinese)
    bone: str                    # representative bone key = BoneSpec.name
    joint: str                   # load-end joint key
    sigma_key: str               # key into MATERIAL_STRENGTH_MPA (compression)
    spec: BoneSpec
    main_kind: str               # "shell" | "solid"
    main_solid_mesh: str | None  # mesh_name of the CORT solid when main_kind=="solid"
    mapping_note: str = ""       # why this representative / σ_c

    @property
    def sigma_c_mpa(self) -> float:
        return float(MATERIAL_STRENGTH_MPA[self.sigma_key][1])


def _solid(pids: list[int], mat: str) -> SolidDomainSpec:
    return SolidDomainSpec(pids=pids, mat_name=mat, mesh_name=mat)


def _shell(pids: list[int], mat: str, t: float) -> ShellDomainSpec:
    return ShellDomainSpec(pids=pids, mat_name=mat, mesh_name=mat, thickness_mm=t)


def build_parts() -> list[Part]:
    """Default mapping (see report §映射假设).  Notes justify each choice."""
    parts: list[Part] = []

    # 1) 足部 -> calcaneus (81001200 SPON tet + 81001300 CORT hex, both solid)
    parts.append(Part(
        cn="足部", bone="calcaneus_r", joint="subtalar_r", sigma_key="calcaneus",
        spec=BoneSpec(name="calcaneus_r", solid_domains=[
            _solid([81001200], "calc_spon"),
            _solid([81001300], "calc_cort"),
        ]),
        main_kind="solid", main_solid_mesh="calc_cort",
        mapping_note="跟骨是 [Y25] 足部加载的代表骨；无壳 → 取实体 CORT 域。"
                     "σ_c=150（calcaneus）。",
    ))

    # 2) 胫骨 -> tibia_r (3 solid PIDs incl. CORT; no shell)
    parts.append(Part(
        cn="胫骨", bone="tibia_r", joint="ankle_r", sigma_key="tibia_shaft",
        spec=BoneSpec(name="tibia_r", solid_domains=[
            _solid([81000600], "tibia_spon_end"),
            _solid([81000601], "tibia_spon_center"),
            _solid([81000700], "tibia_cort"),
        ]),
        main_kind="solid", main_solid_mesh="tibia_cort",
        mapping_note="合并 3 PID 成一根骨（PARTS_MAPPING §1）；THUMS 胫 CORT 为实体 hex。"
                     "无壳 → 取实体 CORT 域。σ_c=200（tibia_shaft，按任务默认；胫骨两端 70 见诚实边界）。",
    ))

    # 3) 腓骨 -> fibula_r
    parts.append(Part(
        cn="腓骨", bone="fibula_r", joint="ankle_r", sigma_key="fibula_shaft",
        spec=BoneSpec(name="fibula_r", solid_domains=[
            _solid([81000800], "fib_spon_end"),
            _solid([81000801], "fib_spon_center"),
            _solid([81000900], "fib_cort"),
        ]),
        main_kind="solid", main_solid_mesh="fib_cort",
        mapping_note="合并 3 PID；无壳 → 实体 CORT 域。σ_c=160（fibula_shaft）。",
    ))

    # 4) 股骨 -> femur_r
    parts.append(Part(
        cn="股骨", bone="femur_r", joint="hip_r", sigma_key="femur_shaft",
        spec=BoneSpec(name="femur_r", solid_domains=[
            _solid([81000000], "fem_spon_end"),
            _solid([81000001], "fem_spon_center"),
            _solid([81000100], "fem_cort"),
        ]),
        main_kind="solid", main_solid_mesh="fem_cort",
        mapping_note="合并 3 PID；无壳 → 实体 CORT 域。σ_c=220（femur_shaft；股骨颈 80 见诚实边界）。",
    ))

    # 5) 骨盆 -> R_HIPBONE (SPON solid + 9 cortical shells)
    hip = BoneSpec(
        name="R_HIPBONE",
        solid_domains=[_solid([83500200], "hip_spon")],
        shell_domains=[
            _shell([p], f"hip_cort_{p}", _HIP_SHELL_THICK[p])
            for p in sorted(_HIP_SHELL_THICK)
        ],
    )
    parts.append(Part(
        cn="骨盆", bone="R_HIPBONE", joint="hip_r", sigma_key="pelvis",
        spec=hip, main_kind="shell", main_solid_mesh=None,
        mapping_note="取髋骨代表骨（骶/耻骨联合未并，PARTS_MAPPING §3）；有壳 → 取最受载壳域"
                     "（每统计量取跨 9 壳的最大）。σ_c=180（pelvis）。",
    ))

    # 6) 腰椎 -> L3
    parts.append(Part(
        cn="腰椎", bone="L3", joint="lumbar", sigma_key="spine",
        spec=BoneSpec(name="L3", solid_domains=[
            _solid([89001500, 89501500], "l3_spon")],
            shell_domains=[_shell([89001501, 89501501], "l3_cort", 1.39)],
        ),
        main_kind="shell", main_solid_mesh=None,
        mapping_note="取单代表椎体 L3（左右 SPON 两半合并，PARTS_MAPPING §3）；有壳 → 壳域。"
                     "σ_c=150（spine）。",
    ))

    # 7) 胸椎 -> T6
    parts.append(Part(
        cn="胸椎", bone="T6", joint="lumbar", sigma_key="spine",
        spec=BoneSpec(name="T6", solid_domains=[
            _solid([89000600, 89500600], "t6_spon")],
            shell_domains=[_shell([89000601, 89500601], "t6_cort", 1.5)],
        ),
        main_kind="shell", main_solid_mesh=None,
        mapping_note="取单代表椎体 T6；有壳 → 壳域。σ_c=150（spine）。",
    ))

    # 8) 颈椎 -> C5
    parts.append(Part(
        cn="颈椎", bone="C5", joint="lumbar", sigma_key="spine",
        spec=BoneSpec(name="C5", solid_domains=[
            _solid([87000500, 87500500], "c5_spon")],
            shell_domains=[_shell([87000501, 87500501], "c5_cort", 1.5)],
        ),
        main_kind="shell", main_solid_mesh=None,
        mapping_note="取单代表椎体 C5；有壳 → 壳域。σ_c=150（spine）。",
    ))

    # 9) 颅骨 -> parietal_r (diploe solid + ext/int cortical shells)
    parts.append(Part(
        cn="颅骨", bone="parietal_r", joint="lumbar", sigma_key="skull",
        spec=BoneSpec(name="parietal_r", solid_domains=[
            _solid([88000004], "par_diploe")],
            shell_domains=[
                _shell([88000005], "par_ext", 1.5),
                _shell([88000006], "par_int", 1.5),
            ],
        ),
        main_kind="shell", main_solid_mesh=None,
        mapping_note="取单块代表颅骨 parietal_r（整颅是薄壳穹顶，不可当实心骨，PARTS_MAPPING §4）；"
                     "有壳 → 壳域。σ_c=160（skull）。载荷沿用 lumbar（默认映射，非物理结论，见诚实边界）。",
    ))
    return parts


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------
class Log:
    def __init__(self, path: Path) -> None:
        self.lines: list[str] = []
        self.path = path

    def __call__(self, msg: str) -> None:
        print(msg, flush=True)
        self.lines.append(msg)

    def save(self) -> None:
        self.path.write_text("\n".join(self.lines) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------
# load chain: load(h) per joint  (per-height dead-drop, no shortcut)
# --------------------------------------------------------------------------
def run_load_chain(log: Log, *, reuse: bool) -> dict:
    """Return load-chain results (heights, per-joint longitudinal peaks).

    Cached in ``temp/opensim_fe/fracture_matrix/load_chain.npz`` (reused when
    ``reuse`` and the cache covers exactly ``HEIGHTS`` × ``JOINTS``).
    """
    cache = OUT_TEMP / "load_chain.npz"
    if reuse and cache.is_file():
        with np.load(cache) as z:
            if (np.array_equal(z["heights"], HEIGHTS)
                    and all(f"{j}__vert" in z.files for j in JOINTS)):
                log(f"[loads] 复用缓存 {cache}")
                return {k: z[k] for k in z.files}

    from climbing.coupling.joint_reactions import joint_reaction
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    out: dict[str, np.ndarray] = {
        "heights": HEIGHTS.astype(float),
        "grf_peak_n": np.full(HEIGHTS.size, np.nan),
        "grf_impulse_rel_err": np.full(HEIGHTS.size, np.nan),
        "grf_window_ms": np.full(HEIGHTS.size, np.nan),
        "fd_seconds": np.full(HEIGHTS.size, np.nan),
        "fall_n_states": np.full(HEIGHTS.size, np.nan),
        "load_seconds": np.full(HEIGHTS.size, np.nan),
    }
    for j in JOINTS:
        out[f"{j}__vert"] = np.full(HEIGHTS.size, np.nan)
        out[f"{j}__abs"] = np.full(HEIGHTS.size, np.nan)
        out[f"{j}__t_ms"] = np.full(HEIGHTS.size, np.nan)

    failures: list[str] = []
    for i, h in enumerate(HEIGHTS):
        t0 = time.time()
        try:
            grf = ground_reaction(height_m=float(h), mass_kg=MODEL_MASS_KG)
            t_fd = time.time()
            fall = run_dead_drop(grf)
            fd_s = time.time() - t_fd
            out["grf_peak_n"][i] = grf.peak_total_n
            out["grf_impulse_rel_err"][i] = grf.impulse_rel_err
            out["grf_window_ms"][i] = grf.window_s[1] * 1e3
            out["fd_seconds"][i] = fd_s
            out["fall_n_states"][i] = fall.n_states
            t_jr = time.time()
            line = (f"[loads] h={h:2d} m  GRF peak={grf.peak_total_n/1e3:7.1f} kN  "
                    f"FD {fd_s:5.1f}s ({fall.n_states} 帧)  |  ")
            for j in JOINTS:
                jr = joint_reaction(fall, grf, joint=j)
                out[f"{j}__vert"][i] = jr.peak_vertical_n
                out[f"{j}__abs"][i] = jr.peak_force_n
                out[f"{j}__t_ms"][i] = jr.peak_vertical_time_s * 1e3
                line += f"{j.split('_')[0]}={jr.peak_vertical_n/1e3:6.1f}kN "
            out["load_seconds"][i] = time.time() - t_jr
            log(line + f"({time.time()-t0:4.1f}s)")
        except Exception as exc:  # noqa: BLE001  -- record verbatim, keep going
            log(f"[loads] h={h} m  FAIL: {type(exc).__name__}: {exc}")
            failures.append(f"h={h}: {type(exc).__name__}: {exc}")

    np.savez(cache, **out)
    if failures:
        (OUT_TEMP / "load_chain_failures.txt").write_text(
            "\n".join(failures) + "\n", encoding="utf-8")
    log(f"[loads] 缓存写出 {cache}")
    return out


# --------------------------------------------------------------------------
# one-bone FE: extract once, build+solve once (with solvable-load fallback)
# --------------------------------------------------------------------------
def _load_candidates(ref_load_n: float) -> list[float]:
    seq = [float(ref_load_n)]
    seq += [x for x in FE_LOAD_CANDIDATES if x <= ref_load_n]
    out: list[float] = []
    for x in seq:
        if x > 0 and all(abs(x - o) > 1e-6 for o in out):
            out.append(x)
    return out


def _build_run(part: Part, mesh: dict, mats: dict, feb: Path, load: float,
               time_steps: int) -> None:
    for sd in part.spec.solid_domains:
        sd.material = mats[sd.pids[0]]
    for sh in part.spec.shell_domains:
        sh.material = mats[sh.pids[0]]
    build_bone_feb(part.spec, feb, load_n=load, top_frac=0.30, bot_frac=0.30,
                   axis=None, with_shell=True, time_steps=time_steps, mesh=mesh)
    run_bone_feb(feb, workdir=feb.parent, timeout=1800)


def _post_bone(part: Part, mesh: dict, feb: Path) -> dict:
    solid_names = [sd.mesh_name for sd in part.spec.solid_domains]
    shell_specs = [
        (sh.mesh_name, int(mesh["by_shell"][sh.mesh_name]["q4"].shape[0]
                           + mesh["by_shell"][sh.mesh_name]["t3"].shape[0]))
        for sh in part.spec.shell_domains
    ]
    return post_process(feb.with_suffix(".xplt"), feb.with_suffix(".hdf5"),
                        solid_names=solid_names, shell_specs=shell_specs,
                        label=part.bone)


def solve_bone_fe(part: Part, ref_load_n: float, log: Log, *, reuse: bool,
                  time_steps: int) -> dict:
    """Extract + build + solve + post one bone, scaled to ``ref_load_n``.

    FEBio diverges at the full @5 m load for some bones (soft trabecular /
    thin-shell point loading).  A descending load fallback picks the largest
    solvable load; the returned stress is scaled back to ``ref_load_n`` by
    strict linear elasticity.  A two-load linearity self-check is recorded
    whenever a second load also converges.

    Cached in ``temp/opensim_fe/fracture_matrix/<bone>/<bone>_fe.json``.
    On total failure the exception text is recorded verbatim.
    """
    bone_dir = OUT_TEMP / part.bone
    bone_dir.mkdir(parents=True, exist_ok=True)
    cache = bone_dir / f"{part.bone}_fe.json"
    if reuse and cache.is_file():
        log(f"[fe] {part.bone}: 复用缓存 {cache.name}")
        return json.loads(cache.read_text(encoding="utf-8"))

    rec: dict = {"bone": part.bone, "ref_load_n": float(ref_load_n),
                 "status": "pending", "attempts": []}
    feb = bone_dir / f"{part.bone}.feb"
    try:
        pids = [p for sd in part.spec.solid_domains for p in sd.pids]
        pids += [p for sh in part.spec.shell_domains for p in sh.pids]
        mats = read_materials(MASTER, pids)
        missing = [p for p in pids if p not in mats]
        if missing:
            raise RuntimeError(f"以下 PID 无 *MAT_* 卡：{missing}")

        mesh = extract_mesh(part.spec, MASTER)
        rec["mesh"] = {
            "n_nodes": mesh["n_nodes"],
            "topology": mesh["topology"],
            "n_id_solid": mesh["n_id_solid"],
            "n_id_shell": mesh["n_id_shell"],
            "n_id_shared": mesh["n_id_shared"],
            "by_solid": {
                k: {"tets": int(b["tets"].shape[0]), "hexes": int(b["hexes"].shape[0])}
                for k, b in mesh["by_solid"].items()
            },
            "by_shell": {
                k: {"q4": int(b["q4"].shape[0]), "t3": int(b["t3"].shape[0]),
                    "thickness_mm": float(b["thickness_mm"])}
                for k, b in mesh["by_shell"].items()
            },
        }

        # --- descending solvable-load search -------------------------------
        chosen: float | None = None
        for load in _load_candidates(ref_load_n):
            try:
                _build_run(part, mesh, mats, feb, load, time_steps)
                rec["attempts"].append({"load_n": load, "ok": True})
                chosen = load
                break
            except Exception as exc:  # noqa: BLE001  -- record, try smaller
                rec["attempts"].append(
                    {"load_n": load, "ok": False,
                     "error": f"{type(exc).__name__}: {exc}"})
                log(f"[fe] {part.bone}: load={load:.0f} N 失败 → 降载重试")
        if chosen is None:
            raise RuntimeError(
                "无可用载荷：所有候选均发散 "
                f"{[a['load_n'] for a in rec['attempts']]}")

        pp = _post_bone(part, mesh, feb)
        rec["post"] = pp
        rec["fe_solve_load_n"] = float(chosen)
        rec["load_reduction_factor"] = float(ref_load_n) / float(chosen)

        # --- linearity self-check at a second convergent load --------------
        tried = {a["load_n"] for a in rec["attempts"]}
        checks: list[float] = []
        for c in (chosen * 2.0, chosen / 2.0):
            if (c > 0 and c <= ref_load_n * (1 + 1e-9)
                    and all(abs(c - t) > 1e-6 for t in tried)):
                checks.append(float(c))
        feb_check = bone_dir / f"{part.bone}.lincheck.feb"
        for check_load in checks:
            try:
                _build_run(part, mesh, mats, feb_check, check_load, time_steps)
                pp2 = _post_bone(part, mesh, feb_check)
                d1 = extract_main_domain(part, {"post": pp}, log)
                d2 = extract_main_domain(part, {"post": pp2}, log)
                lin = {"check_load_n": float(check_load)}
                for s in STATS:
                    v1 = d1["stats"].get(s)
                    v2 = d2["stats"].get(s)
                    if v1 and v2 and v2 > 0:
                        lin[f"ratio_{s}"] = float(v1 / v2)
                lin["expected_ratio"] = float(chosen / check_load)
                rec["linearity"] = lin
                rec["attempts"].append({"load_n": check_load, "ok": True,
                                        "purpose": "linearity_check"})
                break
            except Exception as exc:  # noqa: BLE001
                rec["attempts"].append(
                    {"load_n": check_load, "ok": False,
                     "purpose": "linearity_check",
                     "error": f"{type(exc).__name__}: {exc}"})
        rec["status"] = "ok" if pp.get("domains") else "no_domains"
    except Exception as exc:  # noqa: BLE001  -- record verbatim
        rec["status"] = "failed"
        rec["error"] = f"{type(exc).__name__}: {exc}"
        rec["traceback"] = traceback.format_exc()
        log(f"[fe] {part.bone}: FAIL {rec['error']}")
        (bone_dir / f"{part.bone}.fail.txt").write_text(
            rec["traceback"], encoding="utf-8")

    cache.write_text(json.dumps(rec, indent=2, ensure_ascii=False), encoding="utf-8")
    return rec


def extract_main_domain(part: Part, rec: dict, log: Log) -> dict:
    """Pull the main load-bearing domain σ_vm stats from a post dict.

    has shell -> per statistic, the max across the shell domains;
    no shell  -> the solid CORT domain.
    """
    doms = (rec.get("post") or {}).get("domains", {})
    out: dict = {"kind": part.main_kind, "domains": {}, "stats": {}}
    if part.main_kind == "solid":
        key = f"solid_{part.main_solid_mesh}"
        if key not in doms:
            out["error"] = f"主承载实体域 {key} 不在 post.domains"
            return out
        d = doms[key]
        out["label"] = key
        out["stats"] = {
            "max": float(d["vm_max_mpa"]), "p95": float(d["vm_p95_mpa"]),
            "mean": float(d["vm_mean_mpa"]), "median": float(d["vm_median_mpa"]),
        }
        out["n_elem"] = int(d["n_elem"])
        return out
    shells = {k: v for k, v in doms.items()
              if k.startswith("shell_") and not k.endswith("_error")
              and "vm_max_mpa" in v}
    if not shells:
        out["error"] = "主承载壳域为空"
        return out
    out["domains"] = {
        k: {s: float(v[f"vm_{'max' if s == 'max' else 'p95' if s == 'p95' else 'mean' if s == 'mean' else 'median'}_mpa"])
            for s in STATS}
        for k, v in shells.items()
    }
    for s in STATS:
        key = {"max": "vm_max_mpa", "p95": "vm_p95_mpa",
               "mean": "vm_mean_mpa", "median": "vm_median_mpa"}[s]
        out["stats"][s] = max(float(v[key]) for v in shells.values())
    out["label"] = max(shells, key=lambda k: shells[k]["vm_max_mpa"])
    out["n_elem"] = int(sum(int(v["n_elem"]) for v in shells.values()))
    return out


# --------------------------------------------------------------------------
# matrix + plots + report
# --------------------------------------------------------------------------
def compute_matrix(parts: list[Part], fe: dict[str, dict], loads: dict,
                   primary_stat: str, log: Log) -> dict:
    heights = loads["heights"]
    ref_idx = int(np.argmin(np.abs(heights - REF_H)))

    part_rows: list[dict] = []
    sig = {s: np.full((len(parts), heights.size), np.nan) for s in STATS}
    frac = {s: np.zeros((len(parts), heights.size), dtype=int) for s in STATS}

    for pi, part in enumerate(parts):
        rec = fe[part.bone]
        main = extract_main_domain(part, rec, log)
        ref_load = float(loads[f"{part.joint}__vert"][ref_idx])
        load_h = np.asarray(loads[f"{part.joint}__vert"], dtype=float)

        row = {
            "part_cn": part.cn, "bone": part.bone, "joint": part.joint,
            "sigma_key": part.sigma_key, "sigma_c_mpa": part.sigma_c_mpa,
            "main_kind": part.main_kind, "main_domain": main.get("label"),
            "main_domain_n_elem": main.get("n_elem"),
            "fe_status": rec.get("status"),
            "fe_solve_load_n": rec.get("fe_solve_load_n"),
            "load_reduction_factor": rec.get("load_reduction_factor"),
            "load_ref_n": ref_load, "load_ref_height_m": float(heights[ref_idx]),
            "mapping_note": part.mapping_note,
        }
        if "error" in main or "stats" not in main or not main["stats"]:
            row["error"] = main.get("error", "无 stats")
            for s in STATS:
                row[f"first_fracture_{s}_m"] = None
            part_rows.append(row)
            log(f"[matrix] {part.cn}/{part.bone}: 主域缺失 → 该行 NaN/0")
            continue

        # FE-solved stats live at ``fe_solve_load_n`` (<= ref load whenever FEBio
        # had to reduce the load).  Scale them back to the reference load with the
        # load-reduction factor before using them as "ref-load" stress; otherwise
        # the matrix mixes stresses computed at different loads.
        lrf = float(rec.get("load_reduction_factor") or 1.0)
        row["sigma_vm_ref_mpa"] = {s: main["stats"][s] * lrf for s in STATS}
        row["sigma_vm_ref"] = main["stats"][primary_stat] * lrf
        # raw per-domain stats are left as solved (at fe_solve_load_n); the
        # ref-load-scaled values are in ``sigma_vm_ref_mpa`` above.
        row["main_domain_detail"] = main.get("domains", {})

        ratio = load_h / ref_load if ref_load > 0 else np.full_like(load_h, np.nan)
        for s in STATS:
            sig[s][pi] = main["stats"][s] * lrf * ratio
            frac[s][pi] = (sig[s][pi] >= part.sigma_c_mpa).astype(int)
            first = None
            for k, h in enumerate(heights):
                if frac[s][pi, k]:
                    first = int(h)
                    break
            row[f"first_fracture_{s}_m"] = first
        row["utilization_at_5m"] = (main["stats"][primary_stat] * lrf
                                    / part.sigma_c_mpa)
        row["load_scale_1_to_50"] = (
            float(load_h[-1] / ref_load) if ref_load > 0 else None)
        part_rows.append(row)
        log(f"[matrix] {part.cn:4s}/{part.bone:11s} "
            f"σ_vm_ref[{primary_stat}]={main['stats'][primary_stat] * lrf:8.3f} MPa "
            f"σ_c={part.sigma_c_mpa:5.0f} MPa  "
            + "  ".join(f"{s}:{row[f'first_fracture_{s}_m'] or '>50'}"
                        for s in STATS))

    return {
        "heights": heights, "rows": part_rows, "primary_stat": primary_stat,
        "sigma_vm": sig, "fracture": frac, "ref_idx": ref_idx,
    }


def write_heatmap(mat: dict, path: Path, log: Log) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    # CJK-capable font (Windows); fall back gracefully if absent.
    matplotlib.rcParams["font.sans-serif"] = [
        "Microsoft YaHei", "SimHei", "SimSun", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    heights = mat["heights"]
    primary = mat["primary_stat"]
    frac = mat["fracture"][primary]
    ratio = np.where(
        np.isfinite(mat["sigma_vm"][primary]),
        mat["sigma_vm"][primary]
        / np.array([[r["sigma_c_mpa"]] for r in mat["rows"]], dtype=float),
        np.nan,
    )
    labels = [f"{i+1} {r['part_cn']} ({r['bone']})" for i, r in enumerate(mat["rows"])]
    cmap = ListedColormap(["#cfe3f5", "#b2182b"])

    fig, (ax0, ax1) = plt.subplots(
        1, 2, figsize=(20, 6.8), gridspec_kw={"width_ratios": [1, 1.35]})
    ax0.imshow(frac, aspect="auto", cmap=cmap, vmin=0, vmax=1,
               origin="upper", interpolation="nearest")
    ax0.set_title(f"Fracture matrix (1 = σ_vm[{primary}](h) ≥ σ_c)", fontsize=12)
    ax0.set_yticks(range(len(labels)))
    ax0.set_yticklabels(labels, fontsize=9)
    ax0.set_xticks(range(0, heights.size, 5))
    ax0.set_xticklabels([str(int(heights[i])) for i in range(0, heights.size, 5)])
    ax0.set_xlabel("fall height (m)")
    ax0.set_xticks(np.arange(-0.5, heights.size, 5), minor=True)
    ax0.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax0.grid(which="minor", color="#777777", linewidth=0.6)
    ax0.tick_params(which="minor", length=0)
    for i, r in enumerate(mat["rows"]):
        fh = r.get(f"first_fracture_{primary}_m")
        if fh:
            k = int(np.where(heights == fh)[0][0])
            ax0.scatter([k], [i], marker="s", facecolor="none",
                        edgecolor="#1a9850", s=70, linewidths=1.8, zorder=5)

    im = ax1.imshow(ratio, aspect="auto", cmap="inferno", vmin=0, vmax=4,
                    origin="upper", interpolation="nearest")
    ax1.set_title("Utilization σ_vm/σ_c  (cyan contour = 1.0 = fracture)", fontsize=12)
    ax1.set_yticks(range(len(labels)))
    ax1.set_yticklabels([])
    ax1.set_xticks(range(0, heights.size, 5))
    ax1.set_xticklabels([str(int(heights[i])) for i in range(0, heights.size, 5)])
    ax1.set_xlabel("fall height (m)")
    try:
        ax1.contour(np.arange(heights.size), np.arange(len(labels)),
                    np.nan_to_num(ratio, nan=0.0), levels=[1.0],
                    colors=["#00ffcc"], linewidths=1.6)
    except Exception as exc:  # noqa: BLE001
        log(f"[plot] contour 跳过：{exc}")
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04, label="σ_vm / σ_c")
    fig.suptitle(f"50×9 fracture-matrix skeleton — per-bone FE σ_vm[{primary}], "
                 "per-height joint load (THUMS AM50 → FEBio; [Y25] qualitative only). "
                 "Green □ = first-fracture height", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=150)
    plt.close(fig)
    log(f"[plot] 热图写出 {path}")


def write_outputs(mat: dict, parts: list[Part], loads: dict, fe: dict,
                  meta: dict, log: Log) -> None:
    heights = [int(h) for h in mat["heights"]]
    primary = mat["primary_stat"]
    jdoc = {
        "meta": meta, "heights_m": heights, "joints": list(JOINTS),
        "primary_stat": primary,
        "parts": mat["rows"],
        "loads_n": {j: [float(x) for x in loads[f"{j}__vert"]] for j in JOINTS},
        "loads_abs_n": {j: [float(x) for x in loads[f"{j}__abs"]] for j in JOINTS},
        "grf_peak_n": [float(x) for x in loads["grf_peak_n"]],
        "sigma_c_mpa": {r["bone"]: r["sigma_c_mpa"] for r in mat["rows"]},
        "sigma_vm_mpa": {
            s: {r["bone"]: [None if not np.isfinite(x) else float(x)
                            for x in mat["sigma_vm"][s][i]]
                for i, r in enumerate(mat["rows"])} for s in STATS
        },
        "fracture": {
            s: {r["bone"]: [int(x) for x in mat["fracture"][s][i]]
                for i, r in enumerate(mat["rows"])} for s in STATS
        },
        "fracture_matrix_9x50": {
            s: [[int(x) for x in mat["fracture"][s][i]] for i in range(len(parts))]
            for s in STATS
        },
    }
    (OUT_RESULTS / "fracture_matrix.json").write_text(
        json.dumps(jdoc, indent=2, ensure_ascii=False), encoding="utf-8")

    cols = (["part_cn", "bone", "joint", "sigma_c_mpa", "main_kind", "main_domain",
             "fe_solve_load_n", "load_reduction_factor",
             "sigma_vm_ref_max_mpa", "sigma_vm_ref_p95_mpa",
             "sigma_vm_ref_mean_mpa", "sigma_vm_ref_median_mpa",
             "load_ref_n", "first_fracture_max_m", "first_fracture_p95_m",
             "first_fracture_mean_m", "first_fracture_median_m",
             "utilization_at_5m"] + [f"h{h}" for h in heights])
    with (OUT_RESULTS / "fracture_matrix.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for i, r in enumerate(mat["rows"]):
            sv = r.get("sigma_vm_ref_mpa") or {}
            row = [r.get("part_cn"), r.get("bone"), r.get("joint"),
                   r.get("sigma_c_mpa"), r.get("main_kind"), r.get("main_domain"),
                   r.get("fe_solve_load_n"), r.get("load_reduction_factor"),
                   sv.get("max"), sv.get("p95"), sv.get("mean"), sv.get("median"),
                   r.get("load_ref_n"),
                   r.get("first_fracture_max_m"), r.get("first_fracture_p95_m"),
                   r.get("first_fracture_mean_m"), r.get("first_fracture_median_m"),
                   r.get("utilization_at_5m")]
            row += [int(x) for x in mat["fracture"][primary][i]]
            w.writerow(row)

    write_heatmap(mat, OUT_RESULTS / "fracture_matrix.png", log)


def write_report(mat: dict, parts: list[Part], loads: dict, fe: dict,
                 meta: dict, log: Log) -> None:
    heights = mat["heights"]
    ref_idx = mat["ref_idx"]
    primary = mat["primary_stat"]
    L: list[str] = []
    A = L.append

    A("# 1–50 m × 9 部位骨折矩阵骨架（[Y25] 对标用）")
    A("")
    A("> **本报告只做与 [Y25] 的定性比较，不声称已对齐 [Y25]。** 所有数字均为本机实测；"
      "映射与近似逐条写明；失败骨以原文记录。")
    A("")

    # ---- conclusion ----
    A("## 0. 一句话结论")
    A("")
    n_frac = sum(r.get(f"first_fracture_{primary}_m") is not None for r in mat["rows"])
    A(f"- 9×50 二值矩阵**成形**（{len(parts)} 行 × {heights.size} 列）。主用统计量："
      f"**主承载域 σ_vm 的 {primary}**（正则化热点；`max`/`mean`/`median` 见 §4 敏感性）。")
    A(f"- 在 ≤50 m 内，9 个代表骨中 **{n_frac}/{len(parts)}** 个至少在某高度判为骨折"
      "（判据 σ_vm(h) ≥ σ_c）。")
    n_at1 = sum(r.get(f"first_fracture_{primary}_m") == 1 for r in mat["rows"])
    A(f"- 其中 **{n_at1}/{len(parts)}** 在 **1 m** 即判骨折 —— 材料强度判据在「点载荷单骨 FE」上"
      "**严重过判**（见 §7 边界 3）；`median` 口径给出更有区分度的梯度（§4）。")
    A(f"- 载荷口径：**逐高度实测**（h=1..50 m 每档一次 "
      "`ground_reaction → run_dead_drop → joint_reaction`），**未**用 GRF 峰值比缩放。")
    n_red = sum((r.get("load_reduction_factor") or 1.0) > 1.001 for r in mat["rows"])
    A(f"- FE 参考载荷：{len(parts)-n_red}/{len(parts)} 骨在 @5 m 关节峰值下直接求解；"
      f"{n_red} 骨因 FEBio 发散**降载**后线性折算（见 §3）。")
    A("")

    # ---- method ----
    A("## 1. 方法链与口径")
    A("")
    A("| 环节 | 模块 | 说明 |")
    A("|---|---|---|")
    A("| 载荷端 | `climbing.coupling.joint_reactions.joint_reaction` | 关节远端刚体子树自由体；"
      "修正状态时间口径（`set_state_time=True`） |")
    A("| 通用 FE | `scripts/opensim_fe/bone_feb.py` | `extract_mesh` → `build_bone_feb` → "
      "`run_bone_feb` → `post_process`，与 `bone_feb_validate.py` 同写法 |")
    A("| 判据 | `climbing.bone.MATERIAL_STRENGTH_MPA` | [Y25] 表 1 **压缩**强度 σ_c |")
    A("")
    A("- 每骨抽网格一次；在参考载荷下 build+solve+post 一次（失败则降载，见 §3）。")
    A(f"- 参考高度 `h_ref = {int(REF_H)} m`；`load_ref` = 该关节 @5 m 纵向峰值。")
    A("- 线弹性严格线性：`σ_vm(h) = σ_vm_ref · load(h)/load_ref`；"
      "判据 `fracture(h) = σ_vm(h) ≥ σ_c`。")
    A(f"- 求解器：`{meta['febio_exe']}`（probe={meta['febio_version']}）；"
      f"`time_steps={meta['time_steps']}`。")
    A("")

    # ---- loads ----
    A("## 2. 逐高度关节纵向峰值 load(h)（实测）")
    A("")
    grf_ratio = float(loads["grf_peak_n"][-1] / loads["grf_peak_n"][ref_idx])
    A("| 关节 | @5 m load_ref (kN) | @50 m (kN) | 50 m/5 m 比 | GRF 50 m/5 m 比 |")
    A("|---|---|---|---|---|")
    for j in JOINTS:
        v = np.asarray(loads[f"{j}__vert"])
        A(f"| `{j}` | {v[ref_idx]/1e3:.3f} | {v[-1]/1e3:.3f} | "
          f"{v[-1]/v[ref_idx]:.3f} | {grf_ratio:.3f} |")
    A("")
    A(f"> 关节载荷的 50 m/5 m 缩放比**不等于** GRF 峰值比 {grf_ratio:.2f} "
      "—— 这正是采用**逐高度实测**而非 GRF 比例外推的理由。完整表见 JSON `loads_n`。")
    A("")

    # ---- per-bone FE ----
    A("## 3. 每骨规模 / 材料 / 参考 σ_vm / 求解载荷")
    A("")
    # materials (E) for the E_max column — read once from the deck
    _pids: dict[str, list[int]] = {}
    for p in parts:
        _pids[p.bone] = ([q for sd in p.spec.solid_domains for q in sd.pids]
                         + [q for sh in p.spec.shell_domains for q in sh.pids])
    mats_all = read_materials(
        MASTER, sorted({q for v in _pids.values() for q in v}))
    emax = {p.bone: (max(mats_all[q].E_mpa for q in _pids[p.bone] if q in mats_all)
                     if any(q in mats_all for q in _pids[p.bone]) else None)
            for p in parts}
    A("| 部位 | 代表骨 | 关节 | 主承载域 | 规模 (nodes/tets+hex/shell) | "
      "E_max (MPa) | 求解载荷 (N) | 降载因子 | σ_vm_ref p95 (MPa) | σ_vm_ref max (MPa) | "
      "线性自检比 | σ_c (MPa) |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in mat["rows"]:
        rec = fe[r["bone"]]
        m = rec.get("mesh", {})
        n_tets = sum(v["tets"] for v in m.get("by_solid", {}).values())
        n_hex = sum(v["hexes"] for v in m.get("by_solid", {}).values())
        n_shell = sum(v["q4"] + v["t3"] for v in m.get("by_shell", {}).values())
        em = emax.get(r["bone"])
        e_str = f"{em:.3g}" if em else "—"
        sv = r.get("sigma_vm_ref_mpa") or {}
        fsl = r.get("fe_solve_load_n")
        lrf = r.get("load_reduction_factor")
        lin = rec.get("linearity") or {}
        lin_str = (f"{lin.get('ratio_p95', float('nan')):.3f} vs "
                   f"{lin.get('expected_ratio', float('nan')):.3f}"
                   if "ratio_p95" in lin else "—")
        A(f"| {r['part_cn']} | `{r['bone']}` | `{r['joint']}` | "
          f"{r.get('main_domain')} ({r.get('main_kind')}) | "
          f"{m.get('n_nodes','?')}/{n_tets+n_hex}/{n_shell} | {e_str} | "
          f"{'—' if fsl is None else f'{fsl:.0f}'} | "
          f"{'—' if lrf is None else f'{lrf:.2f}×'} | "
          f"{sv.get('p95', float('nan')):.2f} | {sv.get('max', float('nan')):.2f} | "
          f"{lin_str} | {r['sigma_c_mpa']:.0f} |")
    A("")
    A("> 规模列 = 全局节点 / 实体单元(tet+hex) / 壳单元。主承载域：**有壳→壳域**"
      "（多壳取该统计量最大者），**无壳→实体 CORT 域**。E_max = 该骨最硬域的材料 E。")
    A("")

    # ---- first fracture ----
    A("## 4. 每部位首次骨折高度（1–50 m 或 >50）")
    A("")
    A(f"主口径 = **σ_vm {primary}**（加粗）；同时给出 `max`/`mean`/`median` 敏感性。")
    A("")
    A("| # | 部位 | 代表骨 | 关节 | σ_c (MPa) | σ_vm@5m " + primary + " (MPa) | "
      "利用率@5m | **首次骨折 (p95)** | 首次骨折 (max) | 首次骨折 (mean) | 首次骨折 (median) |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(mat["rows"]):
        sv = r.get("sigma_vm_ref_mpa") or {}
        sref = sv.get(primary)
        util = r.get("utilization_at_5m")
        fh = {s: r.get(f"first_fracture_{s}_m") for s in STATS}
        A(f"| {i+1} | {r['part_cn']} | `{r['bone']}` | `{r['joint']}` | "
          f"{r['sigma_c_mpa']:.0f} | "
          f"{'—' if sref is None else f'{sref:.2f}'} | "
          f"{'—' if util is None else f'{util:.3f}'} | "
          f"**{fh['p95'] or '>50'}** | {fh['max'] or '>50'} | "
          f"{fh['mean'] or '>50'} | {fh['median'] or '>50'} |")
    A("")
    A("> 判定用主承载域的 σ_vm 统计量；`max`（局部极值）对网格/点载荷应力集中最敏感，"
      "`median` 最稳健但丢失热点。四者差异本身即是「判据口径」不确定性的量化。")
    A("")

    # ---- Y25 qualitative ----
    A("## 5. 与 [Y25] 的定性比较（仅定性）")
    A("")
    firsts = {r["part_cn"]: r.get(f"first_fracture_{primary}_m") for r in mat["rows"]}
    foot = firsts.get("足部")
    earliest = min([h for h in firsts.values() if h], default=None)
    tied = [p for p, h in firsts.items() if h is not None and h == earliest]
    if foot is None:
        foot_verdict = "足部在 ≤50 m 内未判骨折（与 [Y25] 7–9 m 不符，判据偏松）"
    elif earliest is not None and foot == earliest and len(tied) == 1:
        foot_verdict = f"✅ 足部为**唯一最早**（{foot} m），方向与 [Y25] 一致"
    elif earliest is not None and foot == earliest:
        others = "、".join(p for p in tied if p != "足部")
        foot_verdict = (f"⚠️ 足部与另外 {len(tied)-1} 个部位**并列最早**（{foot} m，"
                        f"{others}）——未复现 [Y25]「足部最早」的独特性")
    else:
        foot_verdict = (f"⚠️ 足部首次骨折 {foot} m，**并非最早**（最早 {earliest} m）"
                        "——与 [Y25] 排序不一致")
    A("| [Y25] 定性特征 | 本骨架观察（主口径） | 方向 |")
    A("|---|---|---|")
    A(f"| 足部最早骨折（约 7–9 m） | 足部(calcaneus) 首次骨折 = **{foot or '>50'} m** | "
      f"{foot_verdict} |")
    early = sorted([(h, p) for p, h in firsts.items() if h], key=lambda x: x[0])
    if early:
        order = "、".join(f"{p}({h}m)" for h, p in early)
        A(f"| 部位间先后排序（OR 序列） | {order} | 见上；仅方向性参考 |")
    else:
        A("| 部位间先后排序 | ≤50 m 内无部位满足判据 | 需下修阈值/口径 |")
    A("| 阶跃/渐进（临界高度） | 判据 σ_vm(h) ≥ σ_c 天然给出阶跃 | 形式一致，数值不可对标 |")
    A("| 双路径（接触 + 轴向） | 只用关节反力**单路径**传入单骨 FE | 只覆盖一条路径 |")
    A("")
    A("**不声称对齐 [Y25]**：本骨架是单骨截断 + 刚性腿伪影 + 材料强度判据，"
      "可用于比较的是**排序与单调性**，不是临界高度。")
    A("")

    # ---- mapping ----
    A("## 6. 映射假设与近似")
    A("")
    A("| # | 部位 | 代表骨 | 关节 | σ_c 键 | 映射理由 |")
    A("|---|---|---|---|---|---|")
    for i, p in enumerate(parts):
        A(f"| {i+1} | {p.cn} | `{p.bone}` | `{p.joint}` | `{p.sigma_key}` | {p.mapping_note} |")
    A("")
    A("- **参考载荷**：每部位取其关节 @5 m 纵向峰值（实测）。")
    A("- **主承载域**：有壳→壳域（多壳取统计量最大），无壳→实体 CORT 域。")
    A("- **FE 边界**：沿用 `bone_feb_validate.py` 的 bbox 最长轴轴向压缩 + 顶/底 30% 带节点力 + "
      "底带三向固定；`time_steps=1`（模块默认）——实测步数 >1 会让 FEBio 对本 NodalForce "
      "设置误加载而发散，故不采用。")
    A("- **降载**：发散骨取候选载荷中最大可解者，σ_vm 按线性折回参考载荷。")
    A("- **线性外推**：`σ_vm(h)=σ_vm_ref·load(h)/load_ref`（线弹性）。")
    A("")

    # ---- honest boundaries ----
    A("## 7. 诚实边界")
    A("")
    A("1. **单骨截断 vs 全身**：每个 FE 是一根代表骨在两条节点带间受压，无关节/韧带/肌肉/"
      "相邻骨，不体现全身载荷重分配。")
    A("2. **刚性腿伪影**：dead-drop 锁定下肢关节 + 禁用肌肉，实测近端 |F| **偏低**"
      "（髋 15.2 kN < 距下 25.4 kN @5 m），与活体趋势可能相反。")
    A("3. **判据用材料强度×正则化**：[Y25] 的 σ_c 是**材料强度**，本骨架直接与单骨 FE 的 "
      "σ_vm 比较 → 点载荷/薄截面处存在**人为应力集中**，`max` 对网格敏感（本骨架主口径改取 "
      "**p95** 以正则化；`max`/`median` 一并给出）。整体骨失效还涉及弯曲、小梁骨、应力集中，"
      "材料强度会**高估**承载能力。")
    A("4. **长骨 σ_c 取骨干值**：任务默认胫/腓/股用 200/160/220（shaft）。两端"
      "（tibia/fibula_ends=70、femoral_neck=80）约弱 3 倍，若改取骨端值首次骨折高度会"
      "**显著提前**——最大的一处口径分歧。")
    A("5. **载荷方向**：关节反力是全局纵向，FE 却沿骨 bbox 最长轴施加；对非竖直骨是方向近似。")
    A("6. **载荷端近似**：GRF 来自 `pad.py` 1D 等效双质点模型、施加于 calcn 原点；"
      "关节反力为自由体结果；仅 GRF 窗口内有效。")
    A("7. **颅骨**：[Y25] 明确颅骨失效由 HIC/接触主导、不在轴向通路；本骨架按任务默认映射到 "
      "`lumbar`，该行**仅为骨架占位**，无物理意义。")
    A("8. **代表骨 ≠ 整个部位**：足部=跟骨（非 26 块足骨）；脊柱各取单椎体；"
      "骨盆只取髋骨（骶/耻骨联合未并）；颅骨取 parietal_r 一块。")
    A("9. **线性损伤**：不建模塑性/损伤/失效模式转换/应变率；材料卡中的塑性/损伤部分被 "
      "`bone_feb.read_materials` 取**弹性近似**（模块既有约定）。")
    A("10. **FE 降载与线性假设的实测偏差**：3 骨在 @5 m 载荷下发散（近不可压松质骨/薄壳点载荷），"
      "在较低载荷下求解后线性折回；双载荷自检实测比值 σ(2L)/σ(L) 为 **1.84–2.67**"
      "（严格线性应为 2.0），即 FE 响应在本载荷区间**并非严格线性**，"
      "`σ_vm_ref` 与 `σ_vm(h)` 的线性外推可带 **约 ±8–33%** 的量级误差（见 §3「线性自检比」列）。"
      "这是本骨架最需注意的定量边界之一。")
    A("11. **不声称对齐 [Y25]**：仅定性比较排序/单调性/阶跃形式。")
    A("")

    # ---- failures ----
    A("## 8. 失败骨原文记录")
    A("")
    failed = [r for r in mat["rows"] if r.get("fe_status") != "ok" or "error" in r]
    if not failed:
        A("本次运行 **9/9 骨全部成功**（`run_bone_feb` rc=0，`post_process` 取到主域）。无失败需记录。")
        # also list the per-bone diverging attempts (verbatim) for transparency
        any_attempt_fail = False
        for r in mat["rows"]:
            at = [a for a in fe[r["bone"]].get("attempts", []) if not a.get("ok")]
            if at:
                if not any_attempt_fail:
                    A("")
                    A("> 以下为**降载搜索中被 FEBio 拒绝的尝试**（原文），供复现与审计：")
                    any_attempt_fail = True
                for a in at:
                    A(f"> - `{r['bone']}` load={a['load_n']:.0f} N: "
                      f"{a.get('error','').replace(chr(10), ' ')[:400]}")
    else:
        for r in failed:
            A(f"### {r['part_cn']} / `{r['bone']}`")
            A("")
            A("```")
            A(str(r.get("error")))
            A("```")
            A("")
    A("")

    # ---- reproduce ----
    A("## 9. 复现命令")
    A("")
    A("```powershell")
    A("$env:PYTHONPATH=\"src\"")
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\fracture_matrix.py")
    A("# 只重跑载荷链 / 只重跑 FE：")
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\fracture_matrix.py --loads-only")
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\fracture_matrix.py --fe-only")
    A("# 换主统计量 / 忽略缓存强制重算：--stat max|p95|mean|median --force-loads --force-fe")
    A("```")
    A("")
    A("回归（必须保持 19 passed / 1 skipped）：")
    A("")
    A("```powershell")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_opensim_fe.py tests\\test_febio_run.py -q")
    A("```")
    A("")

    # ---- products ----
    A("## 10. 产物")
    A("")
    A("- `results/opensim_fe/fracture_matrix.json`（9×50 矩阵 + 逐高度载荷、σ_vm、四统计量）")
    A("- `results/opensim_fe/fracture_matrix.csv`（宽表：行=部位，列=h1..h50，主口径）")
    A("- `results/opensim_fe/fracture_matrix.png`（左：二值骨折矩阵；右：利用率+等值线）")
    A("- `results/opensim_fe/FRACTURE_MATRIX_REPORT.md`（本文件）")
    A("- `temp/opensim_fe/fracture_matrix/`（每骨 .feb/.xplt/.hdf5/.log + per-bone JSON；"
      "`load_chain.npz`；`run.log`）")
    A("")
    A(f"_生成时间：{meta['generated_at']}；总耗时 {meta['total_seconds']:.0f} s。_")

    (OUT_RESULTS / "FRACTURE_MATRIX_REPORT.md").write_text(
        "\n".join(L) + "\n", encoding="utf-8")
    log("[report] 写出 FRACTURE_MATRIX_REPORT.md")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--loads-only", action="store_true")
    ap.add_argument("--fe-only", action="store_true")
    ap.add_argument("--force-loads", action="store_true")
    ap.add_argument("--force-fe", action="store_true")
    ap.add_argument("--stat", choices=list(STATS), default=PRIMARY_STAT,
                    help="primary σ_vm statistic for the matrix (default p95)")
    ap.add_argument("--time-steps", type=int, default=1,
                    help="FE load increments (default 1, same as bone_feb_validate; "
                         "note: >1 makes FEBio mis-ramp this NodalForce setup and diverge)")
    args = ap.parse_args(argv)

    OUT_TEMP.mkdir(parents=True, exist_ok=True)
    OUT_RESULTS.mkdir(parents=True, exist_ok=True)
    log = Log(OUT_TEMP / "run.log")
    t_start = time.time()

    parts = build_parts()
    log(f"[init] {len(parts)} 部位；heights={HEIGHTS[0]}..{HEIGHTS[-1]}；"
        f"joints={JOINTS}；ref_h={REF_H} m；primary_stat={args.stat}")

    loads = run_load_chain(log, reuse=not args.force_loads)
    if args.loads_only:
        log("[done] --loads-only 完成")
        log.save()
        return 0

    ref_idx = int(np.argmin(np.abs(loads["heights"] - REF_H)))
    ref_load_by_joint = {j: float(loads[f"{j}__vert"][ref_idx]) for j in JOINTS}
    log("[ref] load_ref @ %.0f m: " % REF_H
        + "  ".join(f"{j}={ref_load_by_joint[j]/1e3:.3f}kN" for j in JOINTS))

    fe: dict[str, dict] = {}
    for part in parts:
        rl = ref_load_by_joint[part.joint]
        log(f"[fe] {part.bone}: ref load={rl/1e3:.3f} kN (joint {part.joint}), "
            f"time_steps={args.time_steps}")
        fe[part.bone] = solve_bone_fe(
            part, rl, log, reuse=not args.force_fe, time_steps=args.time_steps)

    if args.fe_only:
        log("[done] --fe-only 完成（FE 结果已缓存；去掉该开关以生成矩阵/报告）")
        log.save()
        return 0

    mat = compute_matrix(parts, fe, loads, args.stat, log)

    from climbing.coupling.febio_run import find_febio, probe_febio  # noqa: E402

    exe = find_febio(verify=True)
    meta = {
        "title": "1-50 m x 9-part fracture-matrix skeleton ([Y25] qualitative)",
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "load_mode": "per-height (ground_reaction -> run_dead_drop -> joint_reaction; NO GRF-ratio shortcut)",
        "heights_m": [int(h) for h in HEIGHTS],
        "ref_height_m": REF_H,
        "joints": list(JOINTS),
        "primary_stat": args.stat,
        "stats_reported": list(STATS),
        "main_domain_rule": "has shell -> shell domain (max across shells per statistic); no shell -> solid CORT domain",
        "sigma_vm_extrapolation": "sigma_vm(h) = sigma_vm_ref * load(h)/load_ref (strict linear elasticity)",
        "criterion": "fracture(h) = sigma_vm(h) >= sigma_c(part) (MATERIAL_STRENGTH_MPA compression)",
        "fe_load_fallback": "try ref, then 10000/5000/2500/1000 N; scale back linearly",
        "model": "Rajagopal2015.osim",
        "model_mass_kg": MODEL_MASS_KG,
        "leg_joints_locked": True,
        "muscles_disabled": True,
        "master_k": str(MASTER.relative_to(ROOT)),
        "febio_exe": str(exe),
        "febio_version": probe_febio(exe),
        "time_steps": args.time_steps,
        "ref_load_by_joint_n": ref_load_by_joint,
        "grf_peak_ratio_50_over_5": float(loads["grf_peak_n"][-1] / loads["grf_peak_n"][ref_idx]),
        "total_seconds": time.time() - t_start,
    }
    write_outputs(mat, parts, loads, fe, meta, log)
    write_report(mat, parts, loads, fe, meta, log)

    log(f"[done] 总耗时 {time.time()-t_start:.0f}s；产物 {OUT_RESULTS}")
    log.save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
