"""R13c: 距下接触 AUGLAG 停滞的真堵点 = 增广容差（augtol）卡死。

R13b 证据（V2_refs100 日志 t=0.755307 状态块）：
    augmentation # 1  sliding interface # 2
       D multiplier : 1.851805e-02   REQUIRED 2.000000e-02   ← 没到容差
       maximum gap  : 4.324446e-02   ***                     ← 0.043mm 间隙未闭合
`REQUIRED 2.0e-02` 恰等于 deck 里硬写的 <tolerance>0.02</tolerance>
（ankle_lig_feb.py:1830）⇒ 该 tag 就是增广容差。

本脚本只改 subtalar 接触块内的容差（与 max_aug_iterations），其余逐字不动。
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
TEMPLATE = ROOT / "temp/pyfebio_demo/r13/aug_p8/ankle_lig_WK.feb"
OUT = ROOT / "temp/pyfebio_demo/r13c"
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
TIMEOUT = 1800
VARIANTS = [
    ("W1_tol0.2", {"tolerance": 0.2}),
    ("W2_tol1.0", {"tolerance": 1.0}),
    ("W3_tol0.2_aug20", {"tolerance": 0.2, "max_aug_iterations": 20}),
    ("W4_tol0.2_reloc1", {"tolerance": 0.2, "node_reloc": 1}),
    ("W5_reloc1", {"node_reloc": 1}),
]
ONLY = [t for t in os.environ.get("ONLY", "").split(",") if t]


def sum_top(p: Path) -> float:
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


def metrics(d: Path) -> dict:
    log = d / "febio.log"
    txt = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    steps = re.findall(r"Number of time steps completed \.+: (\d+)", txt)
    ft = re.findall(r"failed to converge at time :\s*([0-9.eE+-]+)", txt)
    neg = re.findall(r"(\d+) negative jacobians detected", txt)
    dm = re.findall(r"D multiplier :\s*([\d.eE+-]+)\s*([\d.eE+-]+)", txt)
    gp = re.findall(r"maximum gap\s*:\s*([\d.eE+-]+)", txt)
    el = re.findall(r"Elapsed time :\s*([\d:]+)", txt)
    bad_tag = re.search(r'tag "([^"]+)"[^\n]*unrecognized tag', txt)
    return {"normal": "N O R M A L" in txt.replace("  ", " "),
            "steps_done": int(steps[-1]) if steps else None,
            "fail_t": max((float(x) for x in ft), default=None),
            "neg_jac": len(neg), "elapsed": el[-1] if el else None,
            "last_Dmult": dm[-1] if dm else None,
            "last_maxgap": float(gp[-1]) if gp else None,
            "bad_tag": bad_tag.group(1) if bad_tag else None,
            "size": log.stat().st_size if log.exists() else 0}


def run(tag: str, mods: dict) -> dict:
    txt = TEMPLATE.read_text(encoding="latin-1")
    m = list(re.finditer(r"<contact name=\"subtalar\".*?</contact>", txt, re.S))
    assert len(m) == 1
    blk = m[0].group(0)
    assert "<tolerance>0.02</tolerance>" in blk, blk[-400:]
    new = blk.replace("<tolerance>0.02</tolerance>",
                      f"<tolerance>{mods['tolerance']:g}</tolerance>")
    if "node_reloc" in mods:
        assert "<node_reloc>0</node_reloc>" in blk, blk[-400:]
        new = new.replace("<node_reloc>0</node_reloc>",
                          f"<node_reloc>{mods['node_reloc']}</node_reloc>")
    if "max_aug_iterations" in mods:
        new = new.replace("</contact>",
                          f"<max_aug_iterations>{mods['max_aug_iterations']}"
                          f"</max_aug_iterations>\n    </contact>")
    patched = txt[:m[0].start()] + new + txt[m[0].end():]

    d = OUT / tag
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    (d / "ankle_lig_WK.feb").write_text(patched, encoding="latin-1")
    t0 = time.time()
    try:
        r = subprocess.run([str(FEBIO), "-i", "ankle_lig_WK.feb"], cwd=str(d),
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=TIMEOUT)
        (d / "febio.log").write_text((r.stdout or "") + (r.stderr or ""),
                                    encoding="utf-8")
        mm = metrics(d)
    except subprocess.TimeoutExpired:
        mm = {"normal": False, "timeout": True}
    rec = {"tag": tag, "mods": mods, "wall_s": round(time.time() - t0, 1),
           "F_top_z": sum_top(d / "reac_tibia_WK.txt")
           + sum_top(d / "reac_fibula_WK.txt"), **mm}
    print(f"[{tag}] normal={rec['normal']} fail_t={rec.get('fail_t')} "
          f"negJ={rec.get('neg_jac')} Dmult={rec.get('last_Dmult')} "
          f"maxgap={rec.get('last_maxgap')} bad_tag={rec.get('bad_tag')} "
          f"F={rec['F_top_z']:.4f} wall={rec['wall_s']}s", flush=True)
    return rec


def main() -> None:
    assert TEMPLATE.exists(), TEMPLATE
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [run(t, mo) for t, mo in VARIANTS if not ONLY or t in ONLY]
    (OUT / "r13c.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False),
                                   encoding="utf-8")
    print("\n" + "=" * 78)
    print(f"{'tag':>18} {'normal':>7} {'fail_t':>9} {'negJ':>5} {'maxgap':>10} "
          f"{'F_top_z':>10} {'wall':>7}")
    for r in rows:
        print(f"{r['tag']:>18} {str(r['normal']):>7} {str(r.get('fail_t')):>9} "
              f"{str(r.get('neg_jac')):>5} {str(r.get('last_maxgap')):>10} "
              f"{r['F_top_z']:>10.4f} {r['wall_s']:>7}")


if __name__ == "__main__":
    main()