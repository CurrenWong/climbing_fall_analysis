"""R13: 翻案「AUGLAG 更差（穿透 0.46 → 0.65 mm）」。

原结论（A1+A2，2026-10-06）来自**苹果对橘子**的比较：
  SUBTALAR_LAUGON=PENALTY -> ankle_lig_feb.py:1828 写 auto_penalty=1（自动估计刚度）
  SUBTALAR_LAUGON=AUGLAG  -> 同一行写 auto_penalty=0 + 原始 <penalty>=1（极软）
⇒ 0.46 vs 0.65 很可能不是「AUGLAG 表述更差」，而是「penalty=1 的原始罚太软」。

本脚本单变量复验：tibiotalar 保持生产原样（PENALTY+auto=1），只改 deck 里
**subtalar 那一块** 的 <penalty>（2/8/16），逐档重跑 FEBio。

用法:
  .venv\\Scripts\\python.exe -u temp/pyfebio_demo/r13_subtalar_auglag.py
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
DECK = WORK / "ankle_lig_WK.feb"
PY = ROOT / ".venv/Scripts/python.exe"
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
OUT = ROOT / "temp/pyfebio_demo/r13"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
OUT.mkdir(parents=True, exist_ok=True)

DISP = os.environ.get("DISP_MM", "+1.0")
BASE = dict(DRIVER="dispx", DISP_MM=DISP, DET_MIN="0.0", TRUNC_MM="60",
            AXIS_SRC="full", ALIGN_AXIS="1", FIX_TALUS="transverse",
            CART_NEW="1", TIE="1", CART_DIR="temp/thums/cart_m1",
            POSE_IE="0", LOG_DATA="1", SUBTALAR="1", TALUS_SUPPORT="spring",
            TALUS_K="500", MAX_UPS="0", PYTHONIOENCODING="utf-8")
PATCH_PENS = [p for p in os.environ.get("PATCH_PENS", "2,8,16").split(",") if p]
PROD_TIMEOUT = int(os.environ.get("PROD_TIMEOUT", "2400"))
FEB_TIMEOUT = int(os.environ.get("FEB_TIMEOUT", "1500"))


def sum_top(p: Path) -> float:
    """两顶面 z 向反力合力（与 a1a2_matrix.py 同口径）。"""
    if not p.exists():
        return 0.0
    cur_t, buf, last = None, [], {}
    for l in p.read_text(encoding="utf-8", errors="replace").split("\n"):
        if l.startswith("*Time"):
            if buf:
                last[cur_t] = np.array(buf, float)[:, 1:4].sum(0)
            cur_t = float(l.split("=")[1]); buf = []
        elif l.strip() and not l.startswith("*"):
            try:
                buf.append([float(x) for x in l.split(",")])
            except ValueError:
                pass
    if buf:
        last[cur_t] = np.array(buf, float)[:, 1:4].sum(0)
    return float(last[max(last)][2]) if last else 0.0


def febio_metrics(d: Path) -> dict:
    log = d / "febio.log"
    txt = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    steps = re.findall(r"Number of time steps completed \.+: (\d+)", txt)
    ft = re.findall(r"failed to converge at time :\s*([0-9.eE+-]+)", txt)
    neg = re.findall(r"(\d+) negative jacobians detected", txt)
    el = re.findall(r"Elapsed time :\s*([\d:]+)", txt)
    return {"normal": "N O R M A L" in txt.replace("  ", " "),
            "steps_done": int(steps[-1]) if steps else None,
            "fail_t": max((float(x) for x in ft), default=None),
            "neg_jac": len(neg), "elapsed": el[-1] if el else None,
            "diverged": txt.count("Problem is diverging")}


def contact_blocks(txt: str) -> dict:
    out = {}
    for b in re.findall(r"<contact\b.*?</contact>", txt, re.S):
        nm = re.search(r'name="([^"]*)"', b)
        if not nm:
            continue
        g = lambda t: (re.search(rf"<{t}>([^<]*)</{t}>", b) or [None, "-"])[1]
        out[nm.group(1)] = {"laugon": g("laugon"), "penalty": g("penalty"),
                            "auto_penalty": g("auto_penalty")}
    return out


def run_production(tag: str, extra: dict) -> dict:
    env = dict(os.environ)
    env.update(BASE)
    env.update(extra)
    t0 = time.time()
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=PROD_TIMEOUT)
    out = (r.stdout or "") + (r.stderr or "")
    (OUT / f"{tag}.stdout.txt").write_text(out, encoding="utf-8")
    m = re.search(r"FEBio rc=(-?\d+)", out)
    rec = {"tag": tag, "env": extra, "rc": int(m.group(1)) if m else None,
           "wall_s": round(time.time() - t0, 1),
           "F_top_z": sum_top(WORK / "reac_tibia_WK.txt")
           + sum_top(WORK / "reac_fibula_WK.txt")}
    if DECK.exists():
        shutil.copy(DECK, OUT / f"{tag}.deck.feb")
        rec["contacts"] = contact_blocks(DECK.read_text(encoding="latin-1"))
    if (WORK / "febio_WK.log").exists():
        shutil.copy(WORK / "febio_WK.log", OUT / f"{tag}.febio.log")
        rec.update({k: v for k, v in febio_metrics(OUT).items() if False})
        txt = (OUT / f"{tag}.febio.log").read_text(encoding="utf-8", errors="replace")
        steps = re.findall(r"Number of time steps completed \.+: (\d+)", txt)
        rec["steps_done"] = int(steps[-1]) if steps else None
        rec["normal_log"] = "N O R M A L" in txt.replace("  ", " ")
    print(f"[{tag}] rc={rec['rc']} F_top_z={rec['F_top_z']:.4f} N "
          f"wall={rec['wall_s']}s contacts={rec.get('contacts')}", flush=True)
    return rec


def patch_subtalar_and_run(template: Path, pen: str) -> dict:
    txt = template.read_text(encoding="latin-1")
    blocks = list(re.finditer(r"<contact name=\"subtalar\".*?</contact>", txt, re.S))
    assert len(blocks) == 1, f"expected 1 subtalar contact, got {len(blocks)}"
    blk = blocks[0].group(0)
    cur = contact_blocks(blk)["subtalar"]
    assert cur["laugon"].upper() == "AUGLAG" and cur["auto_penalty"] == "0", \
        f"template subtalar contact not AUGLAG/auto0: {cur}"
    m = re.search(r"<penalty>[^<]*</penalty>", blk)
    assert m, "no <penalty> in subtalar block"
    new_blk = blk.replace(m.group(0), f"<penalty>{float(pen):g}</penalty>")
    patched = txt[:blocks[0].start()] + new_blk + txt[blocks[0].end():]
    diff = [x for x in re.findall(r"<penalty>[^<]*</penalty>",
                                  new_blk)]  # only one penalty inside block
    assert diff == [f"<penalty>{float(pen):g}</penalty>"], diff

    d = OUT / f"aug_p{pen.replace('.', 'p')}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    (d / "ankle_lig_WK.feb").write_text(patched, encoding="latin-1")
    t0 = time.time()
    r = subprocess.run([str(FEBIO), "-i", "ankle_lig_WK.feb"], cwd=str(d),
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=FEB_TIMEOUT)
    (d / "febio.log").write_text((r.stdout or "") + (r.stderr or ""),
                                encoding="utf-8")
    rec = {"tag": f"aug_p{pen}", "patched_penalty": float(pen),
           "subtalar_contact": contact_blocks(new_blk)["subtalar"],
           "wall_s": round(time.time() - t0, 1),
           "F_top_z": sum_top(d / "reac_tibia_WK.txt")
           + sum_top(d / "reac_fibula_WK.txt"), **febio_metrics(d)}
    print(f"[aug_p{pen}] normal={rec['normal']} F_top_z={rec['F_top_z']:.4f} N "
          f"steps={rec['steps_done']} fail_t={rec['fail_t']} negJ={rec['neg_jac']} "
          f"wall={rec['wall_s']}s", flush=True)
    return rec


def main() -> None:
    rows = []
    print("=== R13 step 1: production baseline (SUBTALAR_LAUGON=PENALTY, auto=1) "
          "— 生产现状", flush=True)
    rows.append(run_production("base_pen_auto1",
                               {"SUBTALAR_LAUGON": "PENALTY", "PENALTY": "1"}))
    print("=== R13 step 2: emit AUGLAG template (auto=0, penalty=1) ===",
          flush=True)
    rows.append(run_production("auglag_template",
                               {"SUBTALAR_LAUGON": "AUGLAG", "PENALTY": "1"}))
    template = OUT / "auglag_template.deck.feb"
    if not template.exists():
        raise SystemExit("template deck missing — cannot patch")
    for pen in PATCH_PENS:
        try:
            rows.append(patch_subtalar_and_run(template, pen))
        except subprocess.TimeoutExpired:
            print(f"[aug_p{pen}] TIMEOUT >{FEB_TIMEOUT}s", flush=True)
            rows.append({"tag": f"aug_p{pen}", "patched_penalty": float(pen),
                         "timeout": True})

    (OUT / "r13.json").write_text(json.dumps(
        {"disp_mm": DISP, "rows": rows,
         "generated": time.strftime("%Y-%m-%d %H:%M:%S")},
        indent=2, ensure_ascii=False), encoding="utf-8")
    print("\n" + "=" * 78)
    print(f"{'tag':>18} {'rc/normal':>10} {'F_top_z(N)':>12} "
          f"{'steps':>6} {'fail_t':>7} {'negJ':>5} {'wall(s)':>8}")
    for r in rows:
        nm = r.get("normal", r.get("normal_log"))
        print(f"{r['tag']:>18} {str(r.get('rc'))+'/'+str(nm):>10} "
              f"{r.get('F_top_z', float('nan')):>12.4f} "
              f"{str(r.get('steps_done')):>6} {str(r.get('fail_t')):>7} "
              f"{str(r.get('neg_jac')):>5} {str(r.get('wall_s')):>8}")
    print(f"\njson -> {OUT/'r13.json'}")


if __name__ == "__main__":
    main()