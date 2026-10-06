"""Validation driver for bone_feb.py on 4 THUMS bones.

For each bone we (re-)build the ``.feb`` with the new generic builder, run
FEBio, and write a per-bone JSON + a top-level summary.  Numbers go to
``temp/opensim_fe/bone_feb/<bone>/...`` and the report goes to
``results/opensim_fe/BONE_FEB_BUILDER_REPORT.md`` (built separately by hand
or by ``main()``).

Bones
-----
* L3            : 1 solid SPON + 1 shell CORT (regression vs l3_shell_feb.py)
* tibia_r       : 2 solids (SPON_end + SPON_center) + 0 shells (solid-only)
* parietal_r    : 1 solid diploë + 2 shells (external + internal)
* R_HIPBONE     : 1 solid SPON + 9 variant-thickness shells (the hardest case)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(_HERE))

from bone_feb import (  # noqa: E402
    BoneSpec,
    Material,
    ShellDomainSpec,
    SolidDomainSpec,
    build_bone_feb,
    extract_mesh,
    post_process,
    read_materials,
    run_bone_feb,
)

OUT_DIR = ROOT / "temp" / "opensim_fe" / "bone_feb"
MASTER = ROOT / "model" / "AM50_V71_Occupant" / "main_THUMS_AM50_V71.k"

LOAD_N_DEFAULT = 200.0  # N; same as l3_shell_feb.py for L3 regression


# ---------------------------------------------------------------------------
# bone definitions
# ---------------------------------------------------------------------------
def l3_spec() -> BoneSpec:
    """L3 vertebra: L_L3 (89001500) + R_L3 (89501500) SPON (solid) and
    L_L3_CORT (89001501) + R_L3_CORT (89501501) cortical shell."""
    return BoneSpec(
        name="L3",
        solid_domains=[
            SolidDomainSpec(
                pids=[89001500, 89501500],
                mat_name="bone_spon",
                mesh_name="L3_spon",
            )
        ],
        shell_domains=[
            ShellDomainSpec(
                pids=[89001501, 89501501],
                mat_name="bone_cort",
                mesh_name="L3_cort",
                thickness_mm=1.39,
            )
        ],
    )


def tibia_r_spec() -> BoneSpec:
    """Right tibia: SPON (end + center) and CORT solid (no shell)."""
    return BoneSpec(
        name="tibia_r",
        solid_domains=[
            SolidDomainSpec(
                pids=[81000600],
                mat_name="bone_spon_end",
                mesh_name="tibia_spon_end",
            ),
            SolidDomainSpec(
                pids=[81000601],
                mat_name="bone_spon_center",
                mesh_name="tibia_spon_center",
            ),
            SolidDomainSpec(
                pids=[81000700],
                mat_name="bone_cort",
                mesh_name="tibia_cort",
            ),
        ],
        shell_domains=[],
    )


def parietal_r_spec() -> BoneSpec:
    """Right parietal: diploë solid + external/internal cortical shells."""
    return BoneSpec(
        name="parietal_r",
        solid_domains=[
            SolidDomainSpec(
                pids=[88000004],
                mat_name="bone_diploe",
                mesh_name="parietal_diploe",
            )
        ],
        shell_domains=[
            ShellDomainSpec(
                pids=[88000005],
                mat_name="bone_cort",
                mesh_name="parietal_ext",
                thickness_mm=1.5,
            ),
            ShellDomainSpec(
                pids=[88000006],
                mat_name="bone_cort",
                mesh_name="parietal_int",
                thickness_mm=1.5,
            ),
        ],
    )


def hipbone_spec() -> BoneSpec:
    """Right hipbone: 1 SPON + 9 variant-thickness cortical shells.

    Thickness from PARTS_MAPPING.md §1 and _all_names.txt inspection:
        83500201  CORT_1.5                  t=1.5
        83500202  CORT_ACETABLUM1.5         t=1.5
        83500203  CORT_ACETABLUM1.0         t=1.0
        83500204  CORT_0.75                 t=0.75
        83500205  CORT_2                    t=2.0
        83500206  CORT_3                    t=3.0
        83500207  CORT_2.2                  t=2.2
        83500208  CORT_1.8                  t=1.8
        83500209  CORT_1                    t=1.0
    """
    return BoneSpec(
        name="R_HIPBONE",
        solid_domains=[
            SolidDomainSpec(
                pids=[83500200],
                mat_name="bone_spon",
                mesh_name="hip_spon",
            )
        ],
        shell_domains=[
            ShellDomainSpec(
                pids=[83500201], mat_name="bone_cort",
                mesh_name="hip_cort_1p5",   thickness_mm=1.5,
            ),
            ShellDomainSpec(
                pids=[83500202], mat_name="bone_cort",
                mesh_name="hip_cort_ac1p5", thickness_mm=1.5,
            ),
            ShellDomainSpec(
                pids=[83500203], mat_name="bone_cort",
                mesh_name="hip_cort_ac1p0", thickness_mm=1.0,
            ),
            ShellDomainSpec(
                pids=[83500204], mat_name="bone_cort",
                mesh_name="hip_cort_0p75",  thickness_mm=0.75,
            ),
            ShellDomainSpec(
                pids=[83500205], mat_name="bone_cort",
                mesh_name="hip_cort_2",     thickness_mm=2.0,
            ),
            ShellDomainSpec(
                pids=[83500206], mat_name="bone_cort",
                mesh_name="hip_cort_3",     thickness_mm=3.0,
            ),
            ShellDomainSpec(
                pids=[83500207], mat_name="bone_cort",
                mesh_name="hip_cort_2p2",   thickness_mm=2.2,
            ),
            ShellDomainSpec(
                pids=[83500208], mat_name="bone_cort",
                mesh_name="hip_cort_1p8",   thickness_mm=1.8,
            ),
            ShellDomainSpec(
                pids=[83500209], mat_name="bone_cort",
                mesh_name="hip_cort_1",     thickness_mm=1.0,
            ),
        ],
    )


def all_pids(spec: BoneSpec) -> list[int]:
    out: list[int] = []
    for sd in spec.solid_domains:
        out.extend(sd.pids)
    for sh in spec.shell_domains:
        out.extend(sh.pids)
    return out


def _materialize(spec: BoneSpec, mats: dict[int, Material]) -> None:
    for sd in spec.solid_domains:
        # all PIDs in one SolidDomainSpec are assumed to share one *MAT_*
        m = mats.get(int(sd.pids[0]))
        if m is None:
            raise RuntimeError(
                f"{spec.name}: solid {sd.mat_name!r} PID {sd.pids[0]} "
                f"没有 *MAT_* 卡（甲板中可能存在引用但未明确定义）"
            )
        sd.material = m
    for sh in spec.shell_domains:
        m = mats.get(int(sh.pids[0]))
        if m is None:
            raise RuntimeError(
                f"{spec.name}: shell {sh.mat_name!r} PID {sh.pids[0]} "
                f"没有 *MAT_* 卡"
            )
        sh.material = m


def run_one(name: str, spec: BoneSpec, *, load_n: float = LOAD_N_DEFAULT,
            top_frac: float = 0.30, bot_frac: float = 0.30,
            with_shell: bool | None = None,
            axis: int | None = None,
            time_steps: int = 1) -> dict:
    """Build + solve + post-process one bone, return the result dict."""
    bone_dir = OUT_DIR / name
    bone_dir.mkdir(parents=True, exist_ok=True)
    feb = bone_dir / f"{name}.feb"

    # materials
    mats = read_materials(MASTER, all_pids(spec))
    _materialize(spec, mats)

    # mesh (extracted once per bone; .feb writes downstream of this)
    mesh = extract_mesh(spec, MASTER)
    mesh_summary = {
        "n_nodes": mesh["n_nodes"],
        "n_id_solid": mesh["n_id_solid"],
        "n_id_shell": mesh["n_id_shell"],
        "n_id_shared": mesh["n_id_shared"],
        "n_id_shell_only": mesh["n_id_shell_only"],
        "topology": mesh["topology"],
        "by_solid_keys": sorted(mesh["by_solid"].keys()),
        "by_shell_keys": sorted(mesh["by_shell"].keys()),
        "by_solid": {
            mn: {"n_solid": b["n_solid"], "tets": int(b["tets"].shape[0]),
                 "hexes": int(b["hexes"].shape[0]),
                 "n_flipped_tets": int(b.get("n_flipped_tets", 0)),
                 "n_flipped_hexes": int(b.get("n_flipped_hexes", 0))}
            for mn, b in mesh["by_solid"].items()
        },
        "by_shell": {
            mn: {"n_shell": b["n_shell"], "q4": int(b["q4"].shape[0]),
                 "t3": int(b["t3"].shape[0]), "thickness_mm": b["thickness_mm"]}
            for mn, b in mesh["by_shell"].items()
        },
    }

    use_shell = (len(spec.shell_domains) > 0) if with_shell is None else bool(with_shell)

    print("=" * 68, flush=True)
    print(f"BONE: {name}   load={load_n} N  with_shell={use_shell}", flush=True)
    print("=" * 68, flush=True)
    for pid in all_pids(spec):
        m = mats.get(pid)
        if m:
            print(f"  mat PID {pid}: E={m.E_mpa:>10.3f}  ν={m.nu:.3f}  "
                  f"({m.src_keyword}; {m.note})", flush=True)

    meta = build_bone_feb(
        spec, feb,
        load_n=load_n,
        top_frac=top_frac,
        bot_frac=bot_frac,
        axis=axis,
        with_shell=use_shell,
        time_steps=time_steps,
        mesh=mesh,
    )

    # run FEBio
    run_info = run_bone_feb(feb, workdir=feb.parent, timeout=900)
    rc = run_info["rc"]

    # post-process
    solid_names = [sd.mesh_name for sd in spec.solid_domains]
    shell_specs = []
    if use_shell:
        for sh in spec.shell_domains:
            b = mesh["by_shell"][sh.mesh_name]
            shell_specs.append((sh.mesh_name, int(b["q4"].shape[0] + b["t3"].shape[0])))

    pp = {"label": name, "domains": {}, "rc": rc}
    if feb.with_suffix(".xplt").is_file() and rc == 0:
        pp = post_process(
            feb.with_suffix(".xplt"),
            feb.with_suffix(".hdf5"),
            solid_names=solid_names,
            shell_specs=shell_specs,
            label=name,
        )
        pp["rc"] = rc
    else:
        pp["rc"] = rc
        pp["note"] = "FEBio did not produce .xplt"

    result = {
        "bone": name,
        "feb": str(feb),
        "febio_rc": rc,
        "load_n": load_n,
        "with_shell": use_shell,
        "mesh": mesh_summary,
        "materials": {
            mn: {"pid": m.pid, "E_mpa": m.E_mpa, "nu": m.nu,
                 "src_keyword": m.src_keyword, "note": m.note}
            for mn, m in mats.items()
        },
        "feb_meta": meta,
        "run": run_info,
        "post": pp,
    }
    out_json = bone_dir / f"{name}_result.json"
    out_json.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[out] {out_json}", flush=True)
    return result


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary: dict[str, dict] = {}
    cases = [
        ("L3",         l3_spec(),        {"load_n": 200.0}),
        ("tibia_r",    tibia_r_spec(),   {"load_n": 200.0}),
        ("parietal_r", parietal_r_spec(),{"load_n": 200.0}),
        ("R_HIPBONE",  hipbone_spec(),   {"load_n": 200.0}),
    ]
    t_start = time.time()
    for name, spec, kw in cases:
        try:
            r = run_one(name, spec, **kw)
        except Exception as exc:  # noqa: BLE001
            print(f"[FAIL:{name}] {type(exc).__name__}: {exc}", flush=True)
            r = {"bone": name, "error": repr(exc)}
        summary[name] = {
            "rc": r.get("febio_rc"),
            "load_n": r.get("load_n"),
            "with_shell": r.get("with_shell"),
            "error": r.get("error"),
        }
    summary["_meta"] = {
        "total_elapsed_s": time.time() - t_start,
        "master_k": str(MASTER),
    }
    (OUT_DIR / "_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print("\n" + "=" * 68, flush=True)
    print("BONE SUMMARY", flush=True)
    print("=" * 68, flush=True)
    for name, r in summary.items():
        if name == "_meta":
            continue
        print(f"  {name:12s}  rc={r.get('rc')}  load={r.get('load_n')}  "
              f"shell={r.get('with_shell')}  err={r.get('error')}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())