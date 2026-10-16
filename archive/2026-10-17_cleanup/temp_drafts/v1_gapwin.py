"""gap 窗口扫描（penalty=5000 固定）。

推理链：
  * 失败与**外加载荷无关**（×1e-6=0.026 N 仍失败）⇒ 不是载荷造成的
  * 失败与 **penalty 无关**（0.1→5000 同款 8 negJ）⇒ 不是单一接触刚度造成的
  * 但 `gap=-0.05`（初始过盈）下，**初始过盈力 = gap × penalty × 面积 与外加载荷无关**
    ⇒ penalty=5000 时 ≈ 0.05×5000×2700 ≈ 675 kN；penalty=0.1 时 ≈ 13.5 N
      —— 两者量级天差地别却**同款失败** ⇒ 说明**两个不同机制**都给出 ~8 个倒置
        （软 penalty：欠支撑刚体模态；硬 penalty：初始过盈压溃）
  * ⇒ 需要的是**中间窗口**：gap 接近 0 但不为负 ⇒ 无过盈力，且接触一加载即闭合

本跑：penalty=5000，gap ∈ {-0.05(对照), -0.005, 0.0, +0.001, +0.01, +0.05}
"""
import json
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT / "src"), str(ROOT / "scripts" / "opensim_fe")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import thums_feb as T  # noqa: E402
from s2_thums_sweep import _read_domain_vm  # noqa: E402
from climbing.coupling.febio_run import run_febio  # noqa: E402

MESH = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_anatframe.npz"
WORK = ROOT / "temp" / "opensim_fe" / "v1_gapwin"
F0, REF = 26383.133775894516, 186.7563689967737
GAPS = [-0.05, -0.005, 0.0, 0.001, 0.01, 0.05]


def run(mesh, gap, tag):
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    feb = (wd / f"{tag}.feb").resolve()
    for ext in (".xplt", ".hdf5", ".log"):
        p = feb.with_suffix(ext)
        if p.exists():
            p.unlink()
    T.build_thums_feb(mesh, feb, load_n=F0, use_rigid=False, time_steps=1,
                      plantar_bc="contact", contact_gap_mm=gap,
                      contact_hold_spring_k=None,
                      contact_pad_e_mpa=0.185, contact_pad_thickness_mm=15.4,
                      contact_penalty=5000.0, contact_auto_penalty=0)
    t0 = time.time()
    try:
        run_febio(feb, workdir=wd, timeout=1200)
        rc = 0
    except Exception:                              # noqa: BLE001
        rc = 1
    el = time.time() - t0
    lg = feb.with_suffix(".log")
    txt = lg.read_text(encoding="utf-8", errors="replace") if lg.is_file() else ""
    neg = [int(x) for x in re.findall(r'(\d+)\s+negative jacobians detected', txt)]
    steps = re.findall(r'beginning time step \d+ : ([\d.eE+-]+)', txt)
    d = {"gap": gap, "rc": rc, "s": el, "n_neg": len(neg),
         "neg_first": neg[0] if neg else None, "steps": len(steps),
         "last_t": steps[-1] if steps else None}
    if rc == 0:
        vm = _read_domain_vm(feb.with_suffix(".xplt"), "calcaneus",
                             feb.with_suffix(".hdf5"))
        d["sigma_max"] = float(vm.max())
        d["rel"] = float(vm.max()) / REF
    return d


if __name__ == "__main__":
    mesh = T.load_thums_mesh(MESH)
    print(f"h=5 m  penalty=5000  pad_E=0.185  hold=None  STATIC")
    print(f"参照 fixed = {REF:.4f} MPa\n")
    rows = []
    for g in GAPS:
        tag = f"gap{g:+g}"
        r = run(mesh, g, tag)
        rows.append(r)
        s = (f"✅ σ_max={r['sigma_max']:.4f} (fixed×{r['rel']:.4f})"
             if "sigma_max" in r else
             f"❌ 负J×{r['n_neg']} 首={r['neg_first']} 进入步={r['steps']} "
             f"最后t={r['last_t']}")
        print(f"  gap={g:+g}  {s}  ({r['s']:.0f}s)", flush=True)
    out = ROOT / "results" / "opensim_fe" / "v1_gap_window_scan.json"
    out.write_text(json.dumps({"h_m": 5.0, "penalty": 5000.0, "rows": rows},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print("\n->", out)
