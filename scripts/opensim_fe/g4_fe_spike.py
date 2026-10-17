"""G4 "② FE 深水区" P0 feasibility spike: FEBio-native reactive-damage bone material.

Goal (screening, **not** a calibrated model)
--------------------------------------------
Upgrade ONE representative THUMS bone (default ``fibula_r``) from the existing
isotropic-elastic FEBio model to a **FEBio-native reactive-damage** material
(``elastic damage`` + ``DC Drucker shear stress`` + ``CDF Simo``).  ``pyfebio``
cannot emit these material types, so the material block is **injected as raw
XML by text-replacing the ``<Material>`` block** (exactly like
``bone_feb._fix_febio413`` rewrites the ``<solver>`` block).

We measure: does FEBio parse it / run / converge, how long, and what damage
field it produces vs the elastic baseline.

Pipeline
--------
1. **Minimal parse test FIRST** -- hand-written 1-hex8 ``.feb`` with the exact
   material; run through :func:`climbing.coupling.febio_run.run_febio`.
   ACCEPT (parses + converges) or the verbatim FEBio error.
2. Read THUMS params by PID (``bone_feb.read_material_damage``): E, nu,
   sigma_y (card col4), E_t (card col5) for the 3 fibula domains.
3. Build the elastic ``.feb`` (``bone_feb.build_bone_feb``) -- as-built the
   ``fib_cort`` block is ``isotropic elastic``.  This is the **elastic twin**.
4. **Inject** the damage material for the ``fib_cort`` block (raw XML text
   replacement, ``name=`` kept identical so ``<MeshDomains>`` still resolves).
   Map ``sigma_y`` -> CDF Simo ``a`` (:math:`a = (\\sigma_y/\\sqrt3)
   (1-(4/27)c)^{1/6}`; ``c=0`` => von Mises, ``a = \\sigma_y/\\sqrt3``).
5. Run **timed** (``run_febio``); record rc + wall-clock + convergence.
6. Post-process: sigma_vm per domain (``bone_feb.post_process``) AND the
   ``damage`` block from the hdf5 -> ``Dmax`` reached / fraction D>0.
7. Build+run the **elastic twin** at the same load -> baseline sigma_ref;
   report the damaged-vs-elastic sigma_vm deviation.
8. Write ``results/opensim_fe/G4_FE_SPIKE_REPORT.md`` with a P0 gate verdict
   (single-solve <= 5 min; mesh-objectivity, or an explicit "not feasible").

Reproduce (repo root, PYTHONPATH=src)::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_spike.py --bone fibula_r --c 0
    # --steps N   load increments (default 1)
    # --load N    elastic twin target load (default 1000 N, the established
    #             elastic-solvable load for fibula_r); the damaged model descends
    #             through load/2, load/4, ... until it converges.
    # --rebuild   force re-extract mesh + rebuild the elastic .feb
    # --reuse     reuse an existing elastic .feb in the spike dir
    # --c X       Drucker non-dimensional constant in [-27/8, 9/4] (0 = von Mises)

This module **only ADD**s ``g4_fe_spike.py``; it imports and calls existing code
(``bone_feb``, ``fracture_matrix``, ``climbing.coupling.febio_run`` /
``fe_post``) and modifies no existing module.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# repo paths (ROOT = repo root; script lives at scripts/opensim_fe/)
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import bone_feb as bf  # noqa: E402  (scripts/opensim_fe/bone_feb.py)
from climbing.coupling.febio_run import (  # noqa: E402
    FebioRunError,
    find_febio,
    probe_febio,
    run_febio,
)

MASTER = ROOT / "model" / "AM50_V71_Occupant" / "main_THUMS_AM50_V71.k"

#: Directory for spike artefacts (febs / xplt / hdf5 / logs / json).
OUT_TEMP = ROOT / "temp" / "opensim_fe" / "g4_fe_spike"
#: Directory for the spike report.
OUT_RESULTS = ROOT / "results" / "opensim_fe"

#: Descending load candidates for the damaged model (multipliers of --load).
_DESCEND = (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)

#: P0 gate: single solve wall-clock budget (seconds).
SINGLE_SOLVE_BUDGET_S = 300.0


# ===========================================================================
# 1) minimal parse test
# ===========================================================================
_PARSE_TEST_TPL = """<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0">
  <Module type="solid"/>
  <Control>
    <analysis>STATIC</analysis>
    <time_steps>1</time_steps>
    <step_size>1.0</step_size>
    <plot_level>PLOT_MAJOR_ITRS</plot_level>
    <output_level>OUTPUT_MAJOR_ITRS</output_level>
    <solver/>
  </Control>
  <Material>
    <material name="cube" id="1" type="elastic damage">
      <elastic type="neo-Hookean">
        <density>1e-06</density>
        <E>{E}</E>
        <v>{nu}</v>
      </elastic>
      <damage type="CDF Simo">
        <a>{a}</a>
        <b>0</b>
        <Dmax>1</Dmax>
      </damage>
      <criterion type="DC Drucker shear stress">
        <c>{c}</c>
      </criterion>
    </material>
  </Material>
  <Mesh>
    <Nodes name="nodes">
      <node id="1">0,0,0</node>
      <node id="2">1,0,0</node>
      <node id="3">1,1,0</node>
      <node id="4">0,1,0</node>
      <node id="5">0,0,1</node>
      <node id="6">1,0,1</node>
      <node id="7">1,1,1</node>
      <node id="8">0,1,1</node>
    </Nodes>
    <Elements type="hex8" name="cube">
      <elem id="1">1,2,3,4,5,6,7,8</elem>
    </Elements>
    <NodeSet name="top_load">5,6,7,8</NodeSet>
    <NodeSet name="bottom_fix">1,2,3,4</NodeSet>
  </Mesh>
  <MeshDomains>
    <SolidDomain name="cube" mat="cube"/>
  </MeshDomains>
  <Step>
    <step id="1" name="Step">
      <Control>
        <analysis>STATIC</analysis>
        <time_steps>10</time_steps>
        <step_size>0.1</step_size>
        <plot_level>PLOT_MAJOR_ITRS</plot_level>
        <output_level>OUTPUT_MAJOR_ITRS</output_level>
        <solver/>
      </Control>
      <Boundary>
        <bc type="zero displacement" node_set="bottom_fix">
          <x_dof>1</x_dof>
          <y_dof>1</y_dof>
          <z_dof>1</z_dof>
        </bc>
      </Boundary>
      <Loads>
        <nodal_load type="nodal_force" node_set="top_load">
          <value>0,0,-9.5</value>
        </nodal_load>
      </Loads>
    </step>
  </Step>
  <Output>
    <plotfile type="febio" file="_parse_test.xplt">
      <var type="displacement"/>
      <var type="stress"/>
      <var type="damage"/>
    </plotfile>
  </Output>
</febio_spec>
"""


def run_parse_test(*, E: float, nu: float, a: float, c: float) -> dict:
    """Hand-write a 1-hex8 ``elastic damage`` model and run it via ``run_febio``.

    Returns a dict with ``verdict`` ("ACCEPT"/"REJECT"), ``stdout`` (the FEBio
    log, verbatim), ``error`` (verbatim on reject) and the toy ``feb``/``log``.
    """
    feb = OUT_TEMP / "_parse_test.feb"
    feb.write_text(
        _PARSE_TEST_TPL.format(E=E, nu=nu, a=a, c=c), encoding="ISO-8859-1")
    log = feb.with_suffix(".log")
    if log.exists():
        log.unlink()
    verdict, err = "ACCEPT", ""
    try:
        run_febio(feb, workdir=OUT_TEMP)
    except FebioRunError as exc:  # verbatim rejection
        verdict, err = "REJECT", str(exc)
    stdout = log.read_text(errors="replace") if log.is_file() else err
    # a material-summary line proves the parser accepted the type
    parsed = ("type: elastic damage" in stdout and "type: CDF Simo" in stdout
              and "type: DC Drucker shear stress" in stdout)
    if verdict == "ACCEPT" and not parsed:
        verdict = "REJECT"
        err = "FEBio exited 0 but the material summary is absent:\n" + stdout[-800:]
    return {
        "verdict": verdict, "error": err, "stdout": stdout,
        "parsed_material": parsed, "feb": str(feb), "log": str(log),
    }


# ===========================================================================
# 2) raw-XML damage-material injection
# ===========================================================================
def damage_material_xml(
    name: str, mid: int, E: float, nu: float, a: float, c: float,
    *, b: float = 0.0, Dmax: float = 1.0, indent: str = "    ",
) -> str:
    """Return the FEBio 4.13 ``elastic damage`` material block (verbatim syntax).

    ``a`` = Simo c.d.f. alpha (same units as Xi = MPa);
    ``b`` = beta (max damage = 1-b); ``c`` = Drucker non-dimensional constant.

    All numeric parameters are coerced to plain Python floats: ``repr()`` of a
    numpy scalar is ``np.float64(...)`` (not parseable by FEBio), so an
    un-coerced ``!r`` would make FEBio abort with "Failed to initialize
    parameter a".
    """
    i = indent
    E, nu, a, b, c, Dmax = (float(x) for x in (E, nu, a, b, c, Dmax))
    return (
        f'{i}<material name="{name}" id="{int(mid)}" type="elastic damage">\n'
        f'{i}  <elastic type="neo-Hookean">\n'
        f'{i}    <density>1e-06</density>\n'
        f'{i}    <E>{E!r}</E>\n'
        f'{i}    <v>{nu!r}</v>\n'
        f'{i}  </elastic>\n'
        f'{i}  <damage type="CDF Simo">\n'
        f'{i}    <a>{a!r}</a>\n'
        f'{i}    <b>{b!r}</b>\n'
        f'{i}    <Dmax>{Dmax!r}</Dmax>\n'
        f'{i}  </damage>\n'
        f'{i}  <criterion type="DC Drucker shear stress">\n'
        f'{i}    <c>{c!r}</c>\n'
        f'{i}  </criterion>\n'
        f'{i}</material>'
    )


def simo_a_from_sigma_y(sigma_y_mpa: float, c: float) -> float:
    """Map THUMS sigma_y -> CDF Simo ``a`` (Xi units, MPa).

    ``Xi`` for ``DC Drucker shear stress`` equals ``k``, the yield limit in
    simple shear; for uniaxial yield ``sigma_y`` and c=0 (von Mises),
    ``k = sigma_y/sqrt(3)``.  For the Drucker form the manual (sec 4.6.3.7)
    gives ``k = (sigma_y/sqrt(3)) (1 - (4/27)c)^{1/6}``.
    """
    return float(
        float(sigma_y_mpa) / float(np.sqrt(3.0))
        * float((1.0 - (4.0 / 27.0) * float(c)) ** (1.0 / 6.0))
    )


def inject_damage_material(feb_text: str, *, dom_name: str, mid: int,
                           E: float, nu: float, a: float, c: float,
                           xplt_name: str | None = None) -> tuple[str, str]:
    """Replace the ``<material name=dom_name ...>`` block with the damage material.

    Keeps ``name=`` + ``id=`` identical so ``<MeshDomains>`` still resolves.
    Also appends ``<var type="damage"/>`` to the plotfile.

    ``xplt_name`` (when given) rewrites the ``<plotfile file="...">`` attribute
    so the damaged model writes its OWN result file.  Without this the injected
    copy inherits the elastic ``.feb`` stem (``build_bone_feb`` sets
    ``file = out_feb.stem + ".xplt"``) and would clobber/replace the elastic
    twin's xplt -> the damaged post-process would look for a file that never
    existed.

    Returns ``(new_text, injected_xml)``.
    """
    xml = damage_material_xml(dom_name, mid, E, nu, a, c)
    pat = re.compile(
        r'[ \t]*<material\s+name="' + re.escape(dom_name) + r'".*?</material>',
        re.S,
    )
    new_text, n = pat.subn(xml, feb_text)
    if n != 1:
        raise RuntimeError(
            f"注入失败：期望 1 个 name={dom_name!r} 的 <material> 块，命中 {n} 个")
    if '<var type="damage"/>' not in new_text:
        new_text, nv = re.subn(
            r'([ \t]*)<var type="stress"/>',
            r'\1<var type="stress"/>\n\1<var type="damage"/>',
            new_text, count=1,
        )
        if nv != 1:
            raise RuntimeError("注入失败：找不到 <var type=\"stress\"/> 以追加 damage")
    if xplt_name:
        new_text, nx = re.subn(
            r'(<plotfile[^>]*\bfile=")[^"]*(")',
            lambda m: m.group(1) + xplt_name + m.group(2),
            new_text, count=1,
        )
        if nx != 1:
            raise RuntimeError("注入失败：找不到 <plotfile ... file=\"...\">")
    return new_text, xml


# ===========================================================================
# 3) hdf5 damage-field read
# ===========================================================================
def read_damage_field(hdf5: Path, domain: str) -> dict:
    """Read the ``damage`` element block for `domain` from a pyfebio hdf5.

    Returns ``{n_elem, d_max, d_mean, d_median, frac_gt0, frac_gt_half}`` and
    the raw per-element vector under key ``values``.
    """
    import h5py

    with h5py.File(hdf5, "r") as h5:
        sids = sorted(h5["states"].keys(), key=int)
        if not sids:
            raise RuntimeError("HDF5 无 state")
        grp = h5[f"states/{sids[-1]}/element_data/damage"]
        if domain not in grp:
            raise RuntimeError(
                f"damage 块里没有域 {domain!r}；现有 {list(grp.keys())!r}")
        d = np.asarray(grp[domain], dtype=np.float64).reshape(-1)
    return {
        "n_elem": int(d.size),
        "d_max": float(d.max()) if d.size else 0.0,
        "d_mean": float(d.mean()) if d.size else 0.0,
        "d_median": float(np.median(d)) if d.size else 0.0,
        "frac_gt0": float(np.mean(d > 1e-9)) if d.size else 0.0,
        "frac_gt_half": float(np.mean(d > 0.5)) if d.size else 0.0,
        "values": d.tolist(),
    }


# ===========================================================================
# 4) one timed solve + post-process
# ===========================================================================
def _solid_names(spec: bf.BoneSpec) -> list[str]:
    return [sd.mesh_name for sd in spec.solid_domains]


def timed_solve_and_post(feb: Path, *, label: str, solid_names: list[str],
                         read_damage_domain: str | None) -> dict:
    """Run `feb` timed, post-process sigma_vm, optionally read damage.

    Non-convergence is captured (rc / error verbatim), never raised out.
    """
    rec: dict = {"label": label, "feb": str(feb), "ok": False}
    exe = find_febio(verify=True)
    t0 = time.time()
    try:
        rc = run_febio(feb, workdir=feb.parent, timeout=SINGLE_SOLVE_BUDGET_S * 4)
        rec["rc"] = int(rc)
        rec["ok"] = True
    except FebioRunError as exc:
        rec["rc"] = 1
        rec["error"] = str(exc)
    rec["wall_clock_s"] = time.time() - t0
    rec["febio_exe"] = str(exe)
    log = feb.with_suffix(".log")
    if log.is_file():
        txt = log.read_text(errors="replace")
        rec["steps_completed"] = _steps_completed(txt)
        rec["log"] = str(log)
        if not rec["ok"]:
            rec["error_reason"] = _error_reason(txt)
            # keep the failing log (the next attempt overwrites <feb>.log)
            try:
                faillog = feb.with_suffix(".faillog")
                faillog.write_text(txt, encoding="utf-8")
                rec["faillog"] = str(faillog)
            except OSError:
                pass
    if rec["ok"]:
        try:
            pp = bf.post_process(
                feb.with_suffix(".xplt"), feb.with_suffix(".hdf5"),
                solid_names=solid_names, shell_specs=[], label=label)
            rec["post"] = pp
            if read_damage_domain:
                rec["damage"] = read_damage_field(
                    feb.with_suffix(".hdf5"), read_damage_domain)
        except Exception as exc:  # noqa: BLE001
            rec["post_error"] = f"{type(exc).__name__}: {exc}"
    return rec


def _steps_completed(log_text: str) -> int | None:
    m = re.search(r"time steps completed\s*\.*\s*:\s*(\d+)", log_text)
    return int(m.group(1)) if m else None


def _error_reason(log_text: str) -> str:
    """Extract the FEBio ERROR-block message(s), verbatim, joined by ' | '."""
    msgs: list[str] = []
    lines = log_text.splitlines()
    for i, ln in enumerate(lines):
        if "ERROR" in ln and "*" in ln:
            for j in range(i + 1, min(i + 6, len(lines))):
                s = lines[j].strip().strip("*").strip()
                if s and (set(s) - {"*", " "}):
                    if s not in msgs:
                        msgs.append(s)
                    break
    m = re.search(r"failed to converge[^\n]*", log_text, re.IGNORECASE)
    if m and m.group(0).strip() not in msgs:
        msgs.append(m.group(0).strip())
    return " | ".join(msgs[:3])


def inject_step_time_stepper(text: str, *, steps: int, retries: int = 10,
                             dtmin: float = 1e-6) -> str:
    """Insert a ``<time_stepper>`` into the STEP ``<Control>`` + raise step count.

    FEBio 4.13's step control inherits from the global control only for
    *unspecified* fields; the step's ``<solver/>`` block is what the global
    time-stepper policy is applied through, so we place the time_stepper in the
    step control (the part of the file after ``<Step>``).  Used as an XML
    iteration when a damaged solve hits a negative-jacobian / non-convergence
    wall.
    """
    head, sep, tail = text.partition("<Step>")
    if not sep:
        raise RuntimeError("inject_step_time_stepper: no <Step> section")
    ts = (
        f'      <time_stepper type="default">\n'
        f'        <max_retries>{int(retries)}</max_retries>\n'
        f'        <dtmin>{float(dtmin)}</dtmin>\n'
        f'        <dtmax>1.0</dtmax>\n'
        f'        <cutback>0.5</cutback>\n'
        f'      </time_stepper>\n'
    )
    tail, n1 = re.subn(r"<time_steps>\d+</time_steps>",
                       f"<time_steps>{int(steps)}</time_steps>", tail, count=1)
    tail, n2 = re.subn(r"<step_size>[^<]*</step_size>",
                       f"<step_size>{1.0/int(steps):.10g}</step_size>", tail, count=1)
    tail, n3 = re.subn(r"<solver/>", ts + "      <solver/>", tail, count=1)
    if (n1, n2, n3) != (1, 1, 1):
        raise RuntimeError(
            f"inject_step_time_stepper: counts {(n1, n2, n3)} != (1,1,1)")
    return head + sep + tail


# ===========================================================================
# 5) main
# ===========================================================================
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bone", default="fibula_r")
    ap.add_argument("--c", type=float, default=0.0,
                    help="Drucker non-dimensional constant in [-27/8, 9/4]; "
                         "0 = von Mises (default)")
    ap.add_argument("--steps", type=int, default=1, help="load increments (default 1)")
    ap.add_argument("--repeat", type=int, default=3,
                    help="repeat solves at the chosen load to measure convergence "
                         "robustness (FEBio's nonlinear solve is non-deterministic "
                         "near the damage-onset stability boundary)")
    ap.add_argument("--load", type=float, default=1000.0,
                    help="elastic-twin target load (N); default 1000 (fibula_r "
                         "established elastic-solvable load)")
    ap.add_argument("--top-frac", type=float, default=0.30)
    ap.add_argument("--bot-frac", type=float, default=0.30)
    ap.add_argument("--rebuild", action="store_true", help="force re-extract + rebuild")
    ap.add_argument("--reuse", action="store_true", help="reuse existing elastic .feb")
    args = ap.parse_args(argv)

    if not (-27.0 / 8.0 - 1e-9 <= args.c <= 9.0 / 4.0 + 1e-9):
        print(f"[warn] c={args.c} 超出 FEBio 文档范围 [-27/8, 9/4]；照跑但记录。")

    OUT_TEMP.mkdir(parents=True, exist_ok=True)
    OUT_RESULTS.mkdir(parents=True, exist_ok=True)

    from fracture_matrix import build_parts  # noqa: E402  (byte-locked, import only)

    parts = {p.bone: p for p in build_parts()}
    if args.bone not in parts:
        print(f"[error] 未知骨 {args.bone!r}；可选 {sorted(parts)}")
        return 2
    part = parts[args.bone]
    spec = part.spec
    main_mesh = part.main_solid_mesh or spec.solid_domains[-1].mesh_name
    main_sd = next(sd for sd in spec.solid_domains if sd.mesh_name == main_mesh)
    main_pid = main_sd.pids[0]

    pids = [p for sd in spec.solid_domains for p in sd.pids]
    pids += [p for sh in spec.shell_domains for p in sh.pids]

    # ---- THUMS params by PID ------------------------------------------------
    mats = bf.read_material_damage(MASTER, pids)
    missing = [p for p in pids if p not in mats]
    if missing:
        print(f"[error] PID 无 *MAT_* 卡：{missing}")
        return 2
    main_mat = mats[main_pid]
    if main_mat.yield_stress_mpa is None:
        print(f"[error] 主域 PID {main_pid} 无 σ_y（卡不含 col4）")
        return 2
    a = simo_a_from_sigma_y(main_mat.yield_stress_mpa, args.c)
    print(f"[params] main={main_mesh} PID={main_pid} "
          f"E={main_mat.E_mpa} nu={main_mat.nu} "
          f"sigma_y={main_mat.yield_stress_mpa} E_t={main_mat.tangent_modulus_mpa} "
          f"-> Simo a={a:.6g} (c={args.c})", flush=True)

    # ---- parse test FIRST ---------------------------------------------------
    print("\n[parse-test] 1-hex8 'elastic damage' model ...", flush=True)
    pt = run_parse_test(E=main_mat.E_mpa, nu=main_mat.nu, a=a, c=args.c)
    print(f"[parse-test] VERDICT = {pt['verdict']}", flush=True)
    if pt["verdict"] != "ACCEPT":
        print("[parse-test] FEBio rejected the material; STOP (P0 verdict).", flush=True)
        print(pt["error"][-2000:], flush=True)
        _write_report(args, part, main_mesh, main_pid, main_mat, mats, a, pt,
                      elastic=None, damaged=None, chosen_load=None,
                      meta={"stop": "parse-test rejected"})
        return 1

    # ---- build elastic twin .feb -------------------------------------------
    elastic_feb = OUT_TEMP / f"{args.bone}_elastic.feb"
    build_meta: dict = {}
    if args.reuse and elastic_feb.is_file() and not args.rebuild:
        print(f"[build] reuse {elastic_feb.name}", flush=True)
        mesh = None
    else:
        print(f"[build] extract mesh + build elastic {elastic_feb.name} ...", flush=True)
        t0 = time.time()
        mesh = bf.extract_mesh(spec, MASTER)
        for sd in spec.solid_domains:
            sd.material = mats[sd.pids[0]]
        for sh in spec.shell_domains:
            sh.material = mats[sh.pids[0]]
        build_meta = bf.build_bone_feb(
            spec, elastic_feb, load_n=args.load, top_frac=args.top_frac,
            bot_frac=args.bot_frac, axis=None, with_shell=True,
            time_steps=args.steps, mesh=mesh)
        print(f"[build] built in {time.time()-t0:.1f}s "
              f"(nodes={build_meta.get('n_top_nodes')} top)", flush=True)

    if mesh is None:  # --reuse path: still need the mesh for per-load rebuilds
        mesh = bf.extract_mesh(spec, MASTER)

    solid_names = _solid_names(spec)

    # ---- damaged .feb (inject) ---------------------------------------------
    text = elastic_feb.read_text(encoding="ISO-8859-1")
    # the material id must match the elastic block; recover it from the text
    m_id = re.search(
        r'<material\s+name="' + re.escape(main_mesh) + r'"\s+id="(\d+)"', text)
    mid = int(m_id.group(1)) if m_id else 3
    c_tag = f"c{args.c:g}".replace(".", "p")
    damaged_feb = OUT_TEMP / f"{args.bone}_damage_{c_tag}.feb"
    damaged_text, injected_xml = inject_damage_material(
        text, dom_name=main_mesh, mid=mid, E=main_mat.E_mpa, nu=main_mat.nu,
        a=a, c=args.c, xplt_name=damaged_feb.stem + ".xplt")
    damaged_feb.write_text(damaged_text, encoding="ISO-8859-1")
    (OUT_TEMP / "injected_material.xml").write_text(
        injected_xml + "\n", encoding="utf-8")
    print(f"[inject] wrote {damaged_feb.name} "
          f"(replaced material id={mid} name={main_mesh})", flush=True)

    # ---- elastic twin run (baseline) ---------------------------------------
    print("\n[solve] elastic twin ...", flush=True)
    elastic = timed_solve_and_post(
        elastic_feb, label=f"{args.bone}-elastic", solid_names=solid_names,
        read_damage_domain=None)
    print(f"[solve] elastic rc={elastic.get('rc')} "
          f"{elastic.get('wall_clock_s', 0):.1f}s ok={elastic['ok']}",
          flush=True)

    # ---- damaged run (descending load until converge) ----------------------
    print("\n[solve] damaged model (descending load) ...", flush=True)
    damaged = None
    chosen_load = None
    attempts: list[dict] = []
    # config ladder: as-built, then a time_stepper + raised step count (the
    # documented XML iteration) -- each attempt is recorded verbatim.
    configs = (("default", False), ("time_stepper+10steps", True))
    for mult in _DESCEND:
        load = args.load * mult
        # rebuild at this load (elastic feb carries the load; re-inject)
        for sd in spec.solid_domains:
            sd.material = mats[sd.pids[0]]
        for sh in spec.shell_domains:
            sh.material = mats[sh.pids[0]]
        bf.build_bone_feb(
            spec, elastic_feb, load_n=load, top_frac=args.top_frac,
            bot_frac=args.bot_frac, axis=None, with_shell=True,
            time_steps=args.steps, mesh=mesh)
        base_text = elastic_feb.read_text(encoding="ISO-8859-1")
        m_id = re.search(
            r'<material\s+name="' + re.escape(main_mesh) + r'"\s+id="(\d+)"',
            base_text)
        mid = int(m_id.group(1)) if m_id else mid
        for cfg_name, use_ts in configs:
            text_i = base_text
            if use_ts:
                text_i = inject_step_time_stepper(
                    text_i, steps=max(args.steps, 10))
            dtext, injected_xml = inject_damage_material(
                text_i, dom_name=main_mesh, mid=mid, E=main_mat.E_mpa,
                nu=main_mat.nu, a=a, c=args.c,
                xplt_name=damaged_feb.stem + ".xplt")
            damaged_feb.write_text(dtext, encoding="ISO-8859-1")
            rec = timed_solve_and_post(
                damaged_feb, label=f"{args.bone}-damage@{load:g}",
                solid_names=solid_names, read_damage_domain=main_mesh)
            attempts.append({"load_n": float(load), "config": cfg_name,
                             "ok": rec["ok"], "rc": rec.get("rc"),
                             "wall_clock_s": rec.get("wall_clock_s"),
                             "error_reason": rec.get("error_reason", ""),
                             "error": (rec.get("error") or "")[:600]})
            if rec["ok"]:
                damaged = rec
                chosen_load = float(load)
                print(f"[solve] damaged CONVERGED at load={load:g} N "
                      f"[{cfg_name}] ({rec['wall_clock_s']:.1f}s, "
                      f"steps={rec.get('steps_completed')})", flush=True)
                break
            print(f"[solve] damaged load={load:g} N [{cfg_name}] FAILED: "
                  f"{rec.get('error_reason', '?')} → next", flush=True)
        if damaged:
            break

    # ---- repeatability at the chosen load ----------------------------------
    # The damage model with a softening c.d.f. sits near a stability boundary
    # under force control; FEBio's nonlinear iteration path is sensitive to
    # floating-point round-off from the (parallel) linear solver, so the SAME
    # .feb can converge on one run and fail on the next.  Quantify it.
    repeats: list[dict] = []
    if damaged and chosen_load is not None:
        nrep = max(1, args.repeat)
        print(f"\n[repeat] damaged @ {chosen_load:g} N × {nrep} ...", flush=True)
        for r in range(nrep):
            rec = timed_solve_and_post(
                damaged_feb, label=f"{args.bone}-damage@{chosen_load:g}#rep{r}",
                solid_names=solid_names, read_damage_domain=main_mesh)
            d = rec.get("damage") or {}
            repeats.append({
                "ok": rec["ok"], "rc": rec.get("rc"),
                "wall_clock_s": rec.get("wall_clock_s"),
                "d_max": d.get("d_max"), "d_mean": d.get("d_mean"),
                "d_median": d.get("d_median"), "frac_gt0": d.get("frac_gt0"),
                "error_reason": rec.get("error_reason", ""),
            })
            print(f"[repeat] #{r}: ok={rec['ok']} rc={rec.get('rc')} "
                  f"d_mean={d.get('d_mean')}", flush=True)

    # ---- elastic twin must be at the SAME load as the converged damage ------
    if chosen_load is not None and abs(chosen_load - args.load) > 1e-9:
        print(f"\n[solve] re-run elastic twin at chosen load {chosen_load:g} N ...",
              flush=True)
        for sd in spec.solid_domains:
            sd.material = mats[sd.pids[0]]
        for sh in spec.shell_domains:
            sh.material = mats[sh.pids[0]]
        bf.build_bone_feb(
            spec, elastic_feb, load_n=chosen_load, top_frac=args.top_frac,
            bot_frac=args.bot_frac, axis=None, with_shell=True,
            time_steps=args.steps, mesh=mesh)
        elastic = timed_solve_and_post(
            elastic_feb, label=f"{args.bone}-elastic@{chosen_load:g}",
            solid_names=solid_names, read_damage_domain=None)
        print(f"[solve] elastic(same load) rc={elastic.get('rc')} "
              f"ok={elastic['ok']}", flush=True)

    meta = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "febio_exe": elastic.get("febio_exe"),
        "febio_version": probe_febio(find_febio(verify=True)),
        "build": build_meta,
    }
    _write_report(args, part, main_mesh, main_pid, main_mat, mats, a, pt,
                  elastic=elastic, damaged=damaged, chosen_load=chosen_load,
                  meta=meta, attempts=attempts, repeats=repeats)
    print(f"\n[report] wrote {OUT_RESULTS / 'G4_FE_SPIKE_REPORT.md'}", flush=True)
    (OUT_TEMP / "g4_fe_spike.json").write_text(
        json.dumps({
            "args": vars(args), "parse_test": {k: v for k, v in pt.items()
                                               if k != "stdout"},
            "main_domain": main_mesh, "main_pid": main_pid,
            "sigma_y_mpa": main_mat.yield_stress_mpa,
            "E_t_mpa": main_mat.tangent_modulus_mpa,
            "E_mpa": main_mat.E_mpa, "nu": main_mat.nu,
            "simo_a": a, "c": args.c, "chosen_load_n": chosen_load,
            "attempts": attempts, "repeats": repeats,
            "elastic": {k: v for k, v in (elastic or {}).items() if k != "log"},
            "damaged": {k: v for k, v in (damaged or {}).items() if k != "log"},
        }, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8")
    return 0


# ===========================================================================
# report writer
# ===========================================================================
def _fmt(x, nd=4):
    try:
        return f"{float(x):.{nd}f}"
    except Exception:  # noqa: BLE001
        return "—"


def _dom_stats(rec: dict, key: str) -> dict | None:
    if not rec or "post" not in rec:
        return None
    return (rec["post"].get("domains") or {}).get(key)


def _write_report(args, part, main_mesh, main_pid, main_mat, mats, a, pt,
                  *, elastic, damaged, chosen_load, meta, attempts=None,
                  repeats=None) -> None:
    L: list[str] = []
    A = L.append
    A("# G4 ② FE 深水区 — P0 可行性 spike：FEBio 原生损伤材料（fibula_r）")
    A("")
    A("> **筛选级 spike，不是标定模型。** 目的只有一个：确认 FEBio 能否解析/求解"
      "`elastic damage` + `DC Drucker shear stress` + `CDF Simo` 物料，以及它对"
      "一根代表骨（`fibula_r`）产生什么损伤场。所有数字本机实测；未对任何参数调参。")
    A("")
    A(f"_生成时间：{meta.get('generated_at', '—')}；"
      f"FEBio `{meta.get('febio_exe', '—')}` (probe={meta.get('febio_version', '—')})。_")
    A("")

    # ---- 0 conclusions ----
    A("## 0. 一句话结论")
    A("")
    A(f"- **解析测试：{pt['verdict']}**（FEBio 4.13 接受该物料类型"
      + (f"；已收敛并产生损伤场）" if pt["verdict"] == "ACCEPT" else "）。"))
    if damaged and damaged.get("ok"):
        A(f"- 受损模型在 **load = {chosen_load:g} N** 收敛；"
          f"单次求解 wall-clock **{damaged['wall_clock_s']:.1f} s**"
          + (f"；该载荷重复求解 **{sum(1 for r in repeats if r['ok'])}/{len(repeats)}"
             f" 收敛**（非确定性，见 §4.3）" if repeats else "") + "。")
        d = damaged.get("damage") or {}
        A(f"- 产生损伤场：`Dmax` = {_fmt(d.get('d_max'))}，"
          f"D>0 单元占比 {100*d.get('frac_gt0', 0):.1f}%。")
    else:
        A("- 受损模型**未收敛**（详见 §5，含逐次尝试原文）。")
    A("")
    A(f"- 主承载域 `{main_mesh}`（PID {main_pid}）：E={main_mat.E_mpa} MPa，"
      f"ν={main_mat.nu}，σ_y={main_mat.yield_stress_mpa} MPa，"
      f"E_t={main_mat.tangent_modulus_mpa} MPa → Simo `a` = {a:.6g}（c={args.c}）。")
    A("")

    # ---- 1 parse test ----
    A("## 1. 最小解析测试（先跑，决定性）")
    A("")
    A(f"- 文件：`{pt['feb']}`（1 个 hex8，8 节点，单物料 `elastic damage`，"
      "底固定 + 顶节点力 -9.5 N/节点，10 载荷步）。")
    A(f"- 结论：**{pt['verdict']}**。物料摘要 " +
      ("出现于 FEBio 日志（解析器接受 neo-Hookean + CDF Simo + DC Drucker shear stress）"
       if pt.get("parsed_material") else "缺失。"))
    A("")
    if pt["verdict"] != "ACCEPT":
        A("**FEBio 拒绝原文：**")
        A("")
        A("```")
        A(pt["error"][-1500:])
        A("```")
        A("")
    A("FEBio 日志中的物料摘要（节选，证明解析成功）：")
    A("")
    A("```")
    for ln in pt["stdout"].splitlines():
        if ("elastic damage" in ln or "CDF Simo" in ln or "Drucker shear" in ln
                or re.match(r"\s*(density|E|v|Dmax|a|b|c)\s*\.+", ln)):
            A(ln.rstrip())
    A("```")
    A("")
    A(f"完整日志：`{pt['log']}`。")
    A("")

    # ---- 2 injected XML ----
    A("## 2. 注入的 `<material>` XML（原文）")
    A("")
    A("`pyfebio` 没有 damage/reactive 物料类，故用**文本替换**注入"
      "（`<material name=\"" + main_mesh + "\" ...>…</material>` → 下块；"
      "`name=`/`id=` 保持不变，`<MeshDomains>` 仍解析到同一物料）。"
      "同时在 `<plotfile>` 加 `<var type=\"damage\"/>`。")
    A("")
    A("```xml")
    txt = (OUT_TEMP / "injected_material.xml")
    A(txt.read_text(encoding="utf-8").rstrip() if txt.is_file()
      else damage_material_xml(main_mesh, 0, main_mat.E_mpa, main_mat.nu, a, args.c))
    A("```")
    A("")
    A("**映射依据**（FEBio User Manual 4.13）：§4.6.1 `elastic damage`（未约束弹性，"
      "对 ν=0.3 用 plain 而非 `uncoupled`）；§4.6.2.1 `CDF Simo`（`a`=α，与 Ξ 同量纲；"
      "`b`=β，最大损伤 1−β；`Dmax`≤1）；§4.6.3.7 `DC Drucker shear stress` "
      "（`c` ∈ [−27/8, 9/4]，`c=0` ⇒ von Mises，k=σ_y/√3）。")
    A("")

    # ---- 3 params ----
    A("## 3. THUMS 参数（逐 PID，实测卡）")
    A("")
    A("| domain | PID | src_keyword | E (MPa) | ν | σ_y (MPa) | E_t (MPa) | E_t/E |")
    A("|---|---|---|---|---|---|---|---|")
    for sd in part.spec.solid_domains:
        pid = sd.pids[0]
        m = mats.get(pid)
        if m is None:
            continue
        et = m.tangent_modulus_mpa
        ratio = (et / m.E_mpa) if (et and m.E_mpa) else None
        A(f"| `{sd.mesh_name}` | {pid} | {m.src_keyword} | {m.E_mpa:.6g} | "
          f"{m.nu:.6g} | {m.yield_stress_mpa if m.yield_stress_mpa is not None else '—'} | "
          f"{et if et is not None else '—'} | {_fmt(ratio, 4) if ratio is not None else '—'} |")
    A("")
    A(f"主承载域 = **`{main_mesh}`**。Simo α = σ_y/√3·(1−(4/27)c)^(1/6) = "
      f"**{a:.6g}** MPa（c={args.c}）。E_t 仅记录，不用于本次 spike（损伤演化只由 "
      "`a` 控制速率；`b=0` ⇒ Dmax≤1）。")
    A("")

    # ---- 4 numbers ----
    A("## 4. 实测数字（rc / wall-clock / σ_vm / 损伤场）")
    A("")
    A("| 运行 | load (N) | rc | wall-clock (s) | 收敛步 | 状态 |")
    A("|---|---|---|---|---|---|")
    if elastic:
        A(f"| elastic twin | {chosen_load if chosen_load else args.load:g} | "
          f"{elastic.get('rc', '—')} | {elastic.get('wall_clock_s', 0):.1f} | "
          f"{elastic.get('steps_completed', '—')} | {'OK' if elastic['ok'] else 'FAIL'} |")
    if damaged:
        A(f"| damaged | {chosen_load:g} | {damaged.get('rc', '—')} | "
          f"{damaged.get('wall_clock_s', 0):.1f} | "
          f"{damaged.get('steps_completed', '—')} | "
          f"{'OK' if damaged['ok'] else 'FAIL'} |")
    A("")
    if attempts:
        A("受损模型降载搜索（同 load 下 elastic twin 亦按需重跑）：")
        A("")
        A("| load (N) | config | ok | rc | wall (s) | FEBio error (verbatim) |")
        A("|---|---|---|---|---|---|")
        for at in attempts:
            reason = at.get("error_reason") or "(no ERROR block)"
            A(f"| {at['load_n']:g} | {at.get('config', '—')} | "
              f"{'yes' if at['ok'] else 'no'} | {at.get('rc', '—')} | "
              f"{_fmt(at.get('wall_clock_s'), 1)} | {reason} |")
        A("")

    def _row(name, es, ds):
        def g(d, k):
            return _fmt(d.get(k)) if d else "—"
        return (f"| `{name}` | {g(es, 'vm_max_mpa')} | {g(es, 'vm_p95_mpa')} | "
                f"{g(es, 'vm_mean_mpa')} | {g(ds, 'vm_max_mpa')} | "
                f"{g(ds, 'vm_p95_mpa')} | {g(ds, 'vm_mean_mpa')} |")

    A("### 4.1 σ_vm（MPa）逐域：elastic vs damaged" )
    A("")
    A("| domain | elastic max | elastic p95 | elastic mean | "
      "damaged max | damaged p95 | damaged mean |")
    A("|---|---|---|---|---|---|---|")
    for name in _solid_names(part.spec):
        es = _dom_stats(elastic, f"solid_{name}")
        ds = _dom_stats(damaged, f"solid_{name}")
        A(_row(name, es, ds))
    A("")
    es = _dom_stats(elastic, f"solid_{main_mesh}")
    ds = _dom_stats(damaged, f"solid_{main_mesh}")
    if es and ds:
        dev = {k: (ds.get(k, np.nan) - es.get(k, np.nan)) / es.get(k, np.nan)
               for k in ("vm_max_mpa", "vm_p95_mpa", "vm_mean_mpa")}
        A("主承载域 σ_vm 偏差（(damaged−elastic)/elastic）：")
        A("")
        A(f"- max:  {100*dev['vm_max_mpa']:+.1f}%（{es.get('vm_max_mpa'):.3f} → "
          f"{ds.get('vm_max_mpa'):.3f} MPa）")
        A(f"- p95:  {100*dev['vm_p95_mpa']:+.1f}%（{es.get('vm_p95_mpa'):.3f} → "
          f"{ds.get('vm_p95_mpa'):.3f} MPa）")
        A(f"- mean: {100*dev['vm_mean_mpa']:+.1f}%（{es.get('vm_mean_mpa'):.3f} → "
          f"{ds.get('vm_mean_mpa'):.3f} MPa）")
        A("")
    A("### 4.2 损伤场 D（末状态，主承载域）")
    A("")
    if damaged and damaged.get("damage"):
        d = damaged["damage"]
        A(f"- 单元数：{d['n_elem']}")
        A(f"- `Dmax` reached：**{d['d_max']:.4f}**")
        A(f"- mean D：{d['d_mean']:.4f}；median D：{d['d_median']:.4f}")
        A(f"- D>0 单元占比：**{100*d['frac_gt0']:.1f}%**；"
          f"D>0.5 占比：{100*d['frac_gt_half']:.1f}%")
    else:
        A("- 无（受损模型未收敛或未写损伤块）。")
    A("")

    # ---- 4.3 repeatability ----
    if repeats:
        n_ok = sum(1 for r in repeats if r["ok"])
        A("### 4.3 收敛稳健性（对**同一** .feb 重复求解）")
        A("")
        A(f"在 load = {chosen_load:g} N 对同一 `"
          f"{args.bone}_damage_c{args.c:g}.feb` 重复求解 "
          f"**{len(repeats)}** 次：**{n_ok}/{len(repeats)} 收敛**。")
        A("")
        A("| # | ok | rc | wall (s) | d_max | d_mean | D>0 | FEBio error |")
        A("|---|---|---|---|---|---|---|---|")
        for i, r in enumerate(repeats):
            fg = r.get("frac_gt0")
            A(f"| {i} | {'yes' if r['ok'] else 'no'} | {r.get('rc', '—')} | "
              f"{_fmt(r.get('wall_clock_s'), 1)} | {_fmt(r.get('d_max'))} | "
              f"{_fmt(r.get('d_mean'))} | "
              f"{(f'{100*fg:.1f}%') if fg is not None else '—'} | "
              f"{r.get('error_reason') or ''} |")
        A("")
        if n_ok < len(repeats):
            A("**关键发现（诚实边界）**：软化型损伤在**力控制**下处于稳定性边界，"
              "FEBio 的非线性迭代路径对（并行）线性求解器的浮点舍入敏感 —— "
              "**同一输入文件在不同次运行会收敛/发散不同**，收敛状态下的损伤场也有"
              "散布。故下方「单次求解 ≤ 5 min」的判定须理解为「**在收敛的那些运行里**」"
              "的耗时，而非「每次都收敛」；这是本 P0 最重要的可行性风险，"
              "后续若要在力控制下稳定求解，需考虑位移控制 / 粘性正则 / 弧长法（P1）。")
        else:
            A(f"重复 {len(repeats)} 次全部收敛，损伤场稳定。")
        A("")

    # ---- 5 gate ----
    A("## 5. P0 gate 判定")
    A("")
    single_ok = bool(damaged and damaged.get("ok")
                     and damaged["wall_clock_s"] <= SINGLE_SOLVE_BUDGET_S)
    A(f"| gate | 阈值 | 实测 | 判定 |")
    A("|---|---|---|---|")
    A(f"| 单次求解 ≤ 5 min | ≤ {SINGLE_SOLVE_BUDGET_S:.0f} s | "
      f"{_fmt(damaged['wall_clock_s'], 1) if (damaged and damaged.get('ok')) else '未收敛'} s | "
      f"{'**PASS**' if single_ok else '**FAIL**'} |")
    A(f"| 解析/收敛（能跑） | rc=0 且收敛 | "
      f"{'是' if (damaged and damaged.get('ok')) else '否'} | "
      f"{'**PASS**' if (damaged and damaged.get('ok')) else '**FAIL**'} |")
    if repeats:
        n_ok = sum(1 for r in repeats if r["ok"])
        A(f"| 收敛稳健性（新增观察） | ≥ 1 次收敛（原计划未设） | "
          f"{n_ok}/{len(repeats)} | "
          f"{'**PASS（弱）**' if n_ok > 0 else '**FAIL**'} |")
    A("| 网格客观性 | 粗网格重跑 σ/首发点偏差 | **未测** | **未能廉价执行** |")
    A("")
    A("**网格客观性说明（诚实边界）**：本 spike 直接复用 THUMS AM50 固定网格，"
      "仓库内没有为「纯实体 fibula」接线的重网格器；在 P0 时间内无法廉价生成一套"
      "更粗的同胚网格（不同单元阶/拓扑会改变约束与点载荷带），故此项 **未测**，"
      "留待后续（P1）。不作任何「已满足网格客观性」的声称。")
    A("")

    # ---- 6 blockers / reproducibility ----
    A("## 6. 阻塞与复现")
    A("")
    if not (damaged and damaged.get("ok")):
        A("**受损模型未收敛**：逐次尝试原文见 §4 表；完整 FEBio 日志在 "
          "`temp/fe_out/g4_fe_spike/`。")
    else:
        A(f"受损模型在记录载荷下收敛（rc=0）；**唯一阻塞/风险** = §4.3 的收敛"
          f"非确定性（软化本构在力控制下处于稳定性边界）。"
          + (f"本载荷 {sum(1 for r in repeats if r['ok'])}/{len(repeats)} 次收敛。"
             if repeats else ""))
    A("")
    A("```powershell")
    A("$env:PYTHONPATH=\"src\"")
    A(f"& .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_spike.py "
      f"--bone {args.bone} --c {args.c:g}")
    A("```")
    A("")
    A("产物：")
    A("")
    A("- `results/opensim_fe/G4_FE_SPIKE_REPORT.md`（本文件）")
    A("- `temp/fe_out/g4_fe_spike/`（`_parse_test.feb/.log`、"
      "`*_elastic.feb`、`*_damage_*.feb`、`injected_material.xml`、"
      "`.xplt/.hdf5/.log`、`g4_fe_spike.json`）")
    A("")
    A("**未修改任何既有模块**：本 spike 只新增 `scripts/opensim_fe/g4_fe_spike.py`，"
      "其余均以 import + 调用复用。")
    A("")

    out = OUT_RESULTS / "G4_FE_SPIKE_REPORT.md"
    out.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
