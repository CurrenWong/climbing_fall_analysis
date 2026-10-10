"""E1 probe: does auto_penalty=1 explain the p>=2 negative-Jacobian failures?

Method (single-variable, artifact-verified):
  control arm = run build_feb.py as-is   -> deck has <auto_penalty>1</auto_penalty>
  testarm     = same deck, ONLY <auto_penalty>1</...> -> 0, then FEBio re-run in an
                isolated dir.  Diff assertion: exactly one line may differ.

Why: FEBio's auto penalty makes <penalty> a *scale factor* on an auto estimate
(primary surface stiffness + mesh size); the official webinar warns that mode can
generate a large enough penalty to invert elements.  Our project rule
(ankle_lig_feb.py:127-128) says penalty sweeps must disable it -- but the sweep
driver went through build_feb.py:696 which re-enables it.

Usage:
    .venv\\Scripts\\python.exe scripts/ankle_fe/penalty_auto0_probe.py
    .venv\\Scripts\\python.exe scripts/ankle_fe/penalty_auto0_probe.py --pen 2 4
Outputs: results/opensim_fe/ankle_penalty_auto0.json  (+ stdout table)
"""
from __future__ import annotations

import argparse
import difflib
import itertools
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
DECK = ROOT / "temp/thums/febio/ankle_tibiotalar.feb"
OUTDIR = ROOT / "temp/thums/febio/_e1"
JSONOUT = ROOT / "results/opensim_fe/ankle_penalty_auto0.json"

# same configuration as the original sweep (temp/pyfebio_demo/pen_converge.py)
BASE_ENV = {
    "TIE": "0",
    "LAUGON": "PENALTY",
    "TIME_STEPS": "20",
    "DRIVER": "disp",
    "DISP_MM": "0.1472",
}
RC_RE = re.compile(r"FEBio rc=(-?\d+)\s+用时\s+([\d.]+)s")


def signals(text: str) -> dict:
    flat = text.replace("  ", " ")
    m = RC_RE.search(text)
    tj = re.search(r"Elapsed time :\s*(\d+):(\d+):(\d+)", text)
    steps = re.search(r"Number of time steps completed \.+: (\d+)", text)
    ft = re.findall(r"failed to converge at time :\s*([\d.eE+-]+)", text)
    dtm = re.findall(r"Time step size\s*[:=]\s*([\d.eE+-]+)", text)
    neg = re.findall(r"(\d+) negative jacobians detected|Negative jacobian detected",
                     text)
    return {
        "rc": int(m.group(1)) if m else None,
        "wall_s": float(m.group(2)) if m else None,
        "normal_term": "N O R M A L" in flat,
        "neg_jac": bool(neg),
        "neg_jac_msg": (neg[0] if neg and isinstance(neg[0], str) else
                        (f"{neg[0]} negative jacobians" if neg else None)),
        "steps_done": int(steps.group(1)) if steps else None,
        "fail_t": max((float(x) for x in ft), default=None),
        "fail_t_last": sorted({float(x) for x in ft})[-4:] or None,
        "fe_elapsed": (f"{tj.group(1)}:{tj.group(2)}:{tj.group(3)}" if tj else None),
        "diverged": "Problem is diverging" in text,
    }


def run_control(pen: str, timeout: int) -> dict:
    env = dict(os.environ)
    env.update(BASE_ENV)
    env["PENALTY"] = pen
    t0 = time.time()
    r = subprocess.run([str(PY), "-u", "scripts/ankle_fe/build_feb.py"],
                       cwd=str(ROOT), env=env, capture_output=True,
                       text=True, timeout=timeout)
    out = (r.stdout or "") + (r.stderr or "")
    rec = {"pen": pen, "arm": "control(auto=1)", **signals(out),
           "driver_wall_s": round(time.time() - t0, 1)}
    if DECK.exists():
        txt = DECK.read_text(encoding="latin-1")
        rec["deck_auto_penalty"] = (re.search(r"<auto_penalty>([^<]*)</auto_penalty>",
                                              txt) or [None, "?"])[1]
    return rec


def patch_and_run(pen: str, timeout: int, diag_negj: bool = False) -> dict:
    """Copy deck, flip ONLY auto_penalty 1 -> 0 (assert exactly one line), re-run."""
    if not DECK.exists():
        raise FileNotFoundError(DECK)
    orig = DECK.read_text(encoding="latin-1")
    # guard: the deck must be the one just emitted for THIS penalty, otherwise we
    # would be testing the previous case's <penalty> with a flipped auto flag.
    mp = re.search(r"<penalty>([^<]*)</penalty>", orig)
    if not mp or abs(float(mp.group(1)) - float(pen)) > 1e-12:
        raise RuntimeError(
            f"deck <penalty>={mp.group(1) if mp else None} != requested {pen} "
            f"-- stale deck, run the control arm first (--arms both)")
    n = orig.count("<auto_penalty>1</auto_penalty>")
    if n != 1:
        raise RuntimeError(
            f"deck auto_penalty=1 count = {n} (expected 1) -- deck format changed?")
    patched = orig.replace("<auto_penalty>1</auto_penalty>",
                           "<auto_penalty>0</auto_penalty>")
    negj_elems: list[str] = []
    if diag_negj:
        # official switch: make FEBio print every element/gauss point that inverts
        flag = "<output_negative_jacobians>1</output_negative_jacobians>"
        if flag not in patched:
            if "</Output>" not in patched:
                raise RuntimeError("no </Output> to inject diagnostic flag into")
            patched = patched.replace("</Output>", f"  {flag}\n  </Output>", 1)
    # line-level guard: every changed hunk must be the auto_penalty line, the
    # </Output> line, or the injected diagnostic flag -- nothing else moves.
    o_lines, p_lines = orig.splitlines(), patched.splitlines()
    allowed = ("<auto_penalty>", "</Output>", "<output_negative_jacobians>")
    touched: list[str] = []
    for tag_, i1, i2, j1, j2 in difflib.SequenceMatcher(
            None, o_lines, p_lines).get_opcodes():
        if tag_ == "equal":
            continue
        touched += o_lines[i1:i2] + p_lines[j1:j2]
    bad = [ln for ln in touched if not any(k in ln for k in allowed)]
    assert not bad, f"patch leaked outside allowed tags: {bad[:3]}"
    assert any("<auto_penalty>0<" in ln for ln in touched), "auto flip missing"

    tag = f"pen{pen.replace('.', 'p')}_auto0_ts{BASE_ENV['TIME_STEPS']}" \
          + ("_negj" if diag_negj else "")
    cdir = OUTDIR / tag
    if cdir.exists():
        shutil.rmtree(cdir)
    cdir.mkdir(parents=True)
    deck2 = cdir / "ankle_tibiotalar.feb"
    deck2.write_text(patched, encoding="latin-1")

    t0 = time.time()
    r = subprocess.run([str(FEBIO), "-i", str(deck2.name)], cwd=str(cdir),
                       capture_output=True, text=True, timeout=timeout)
    out = (r.stdout or "") + (r.stderr or "")
    (cdir / "febio.log").write_text(out, encoding="utf-8", errors="replace")
    negj_elems = sorted(set(re.findall(
        r"Negative jacobian was detected at element (\d+) at gauss point (\d+)",
        out)), key=lambda t: (int(t[0]), int(t[1])))
    return {"pen": pen, "arm": "test(auto=0)", **signals(out),
            "driver_wall_s": round(time.time() - t0, 1),
            "patched_lines": len(touched) // 2, "case_dir": str(cdir),
            "negj_elem_gauss": [f"e{a}@gp{b}" for a, b in negj_elems][:40],
            "negj_elem_count": len(negj_elems)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pen", nargs="*", default=["0.001", "0.01", "0.1", "1",
                                                 "2", "4", "8"])
    ap.add_argument("--timeout", type=int, default=2400)
    ap.add_argument("--time-steps", default=None,
                    help="override TIME_STEPS (deck step_size = 1/this); E2 lever")
    ap.add_argument("--diag-negj", action="store_true",
                    help="also inject <output_negative_jacobians>1 (output-only, "
                         "lists every inverting element/gauss point)")
    ap.add_argument("--arms", default="both", choices=["both", "control", "test"])
    ap.add_argument("--out", default=str(JSONOUT))
    args = ap.parse_args()

    if not FEBIO.exists():
        raise SystemExit(f"FEBio missing: {FEBIO}")
    if args.time_steps:
        BASE_ENV["TIME_STEPS"] = str(args.time_steps)
    print(f"[E1] config={BASE_ENV}", flush=True)

    rows: list[dict] = []
    for pen in args.pen:
        if args.arms in ("both", "control"):
            rec = run_control(pen, args.timeout)
            rows.append(rec)
            print(f"  [control auto=1] PENALTY={pen:>5s} rc={rec['rc']} "
                  f"normal={rec['normal_term']} negJ={rec['neg_jac']} "
                  f"steps={rec['steps_done']} {rec['wall_s']}s", flush=True)
        if args.arms in ("both", "test"):
            try:
                rec = patch_and_run(pen, args.timeout, diag_negj=args.diag_negj)
            except subprocess.TimeoutExpired:
                rec = {"pen": pen, "arm": "test(auto=0)", "rc": None,
                       "normal_term": False, "neg_jac": None,
                       "timeout": True, "driver_wall_s": args.timeout}
            rows.append(rec)
            print(f"  [test    auto=0] PENALTY={pen:>5s} rc={rec['rc']} "
                  f"normal={rec['normal_term']} negJ={rec['neg_jac']} "
                  f"steps={rec['steps_done']} {rec.get('wall_s')}s", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"config": BASE_ENV, "rows": rows,
                               "generated": time.strftime("%Y-%m-%d %H:%M:%S")},
                              indent=2, ensure_ascii=False), encoding="utf-8")

    print("\n" + "=" * 84)
    print(f"{'PENALTY':>9} {'arm':>14} {'rc':>3} {'NORMAL':>7} {'negJ':>5} "
          f"{'steps':>6} {'time(s)':>9}")
    for r in rows:
        print(f"{r['pen']:>9} {r['arm']:>14} {str(r['rc']):>3} "
              f"{str(r['normal_term']):>7} {str(r['neg_jac']):>5} "
              f"{str(r['steps_done']):>6} {str(r.get('wall_s')):>9}")
    print(f"\njson -> {out}")


if __name__ == "__main__":
    main()