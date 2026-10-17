"""R13b: 距下接触 AUGLAG(penalty=8) 停滞于 t≈0.77 的处方试跑（roadmap §判读）。

停滞特征 = 残差卡住 + 零负 J（不是翻转）⇒ 处方是"迭代次数/时间步不够"：
  V1 只细分步长  time_steps 20 -> 40（step_size 同步 0.05 -> 0.025）
  V2 只加 reform  max_refs 50 -> 100
  V3 两者同时
基线 deck = R13 产出的 aug_p8（subtalar AUGLAG/auto0/penalty=8，其余与生产一致）。

用法: .venv\\Scripts\\python.exe -u temp/pyfebio_demo/r13b_stall_probe.py
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
OUT = ROOT / "temp/pyfebio_demo/r13b"
TEMPLATE = ROOT / "temp/pyfebio_demo/r13/aug_p8/ankle_lig_WK.feb"
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
TIMEOUT = int(1800)
VARIANTS = [
    ("V1_ts40", {"time_steps": 40}),
    ("V2_refs100", {"max_refs": 100}),
    ("V3_ts40_refs100", {"time_steps": 40, "max_refs": 100}),
]


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
    res = re.findall(r"residual\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+", txt)
    last_res = float(res[-1][1]) if res else None
    el = re.findall(r"Elapsed time :\s*([\d:]+)", txt)
    return {"normal": "N O R M A L" in txt.replace("  ", " "),
            "steps_done": int(steps[-1]) if steps else None,
            "fail_t": max((float(x) for x in ft), default=None),
            "neg_jac": len(neg), "last_residual": last_res,
            "diverged": txt.count("Problem is diverging"),
            "elapsed": el[-1] if el else None}


def run(tag: str, mods: dict) -> dict:
    txt = TEMPLATE.read_text(encoding="latin-1")
    if "time_steps" in mods:
        ts = mods["time_steps"]
        assert txt.count("<time_steps>20</time_steps>") == 1
        assert txt.count("<step_size>0.050000</step_size>") == 1
        txt = txt.replace("<time_steps>20</time_steps>",
                          f"<time_steps>{ts}</time_steps>")
        txt = txt.replace("<step_size>0.050000</step_size>",
                          f"<step_size>{1.0/ts:.6f}</step_size>")
    if "max_refs" in mods:
        assert txt.count("<max_refs>50</max_refs>") == 1
        txt = txt.replace("<max_refs>50</max_refs>",
                          f"<max_refs>{mods['max_refs']}</max_refs>")
    d = OUT / tag
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    (d / "ankle_lig_WK.feb").write_text(txt, encoding="latin-1")
    t0 = time.time()
    try:
        r = subprocess.run([str(FEBIO), "-i", "ankle_lig_WK.feb"], cwd=str(d),
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=TIMEOUT)
        (d / "febio.log").write_text((r.stdout or "") + (r.stderr or ""),
                                    encoding="utf-8")
        m = metrics(d)
    except subprocess.TimeoutExpired:
        m = {"normal": False, "timeout": True}
    rec = {"tag": tag, "mods": mods, "wall_s": round(time.time() - t0, 1),
           "F_top_z": sum_top(d / "reac_tibia_WK.txt")
           + sum_top(d / "reac_fibula_WK.txt"), **m}
    print(f"[{tag}] normal={rec['normal']} fail_t={rec.get('fail_t')} "
          f"steps={rec.get('steps_done')} negJ={rec.get('neg_jac')} "
          f"last_res={rec.get('last_residual')} F={rec['F_top_z']:.4f} "
          f"wall={rec['wall_s']}s", flush=True)
    return rec


def main() -> None:
    assert TEMPLATE.exists(), TEMPLATE
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [run(t, m) for t, m in VARIANTS]
    (OUT / "r13b.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False),
                                   encoding="utf-8")
    print("\n" + "=" * 70)
    print(f"{'tag':>20} {'normal':>7} {'fail_t':>8} {'steps':>6} {'negJ':>5} "
          f"{'last_res':>10} {'F_top_z':>10} {'wall':>7}")
    for r in rows:
        print(f"{r['tag']:>20} {str(r['normal']):>7} {str(r.get('fail_t')):>8} "
              f"{str(r.get('steps_done')):>6} {str(r.get('neg_jac')):>5} "
              f"{str(r.get('last_residual')):>10} {r['F_top_z']:>10.4f} "
              f"{r['wall_s']:>7}")


if __name__ == "__main__":
    main()