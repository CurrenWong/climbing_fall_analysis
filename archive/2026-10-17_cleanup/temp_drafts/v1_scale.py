"""决定性判别：那 8 个负 Jacobian 与**载荷量级**有关吗？

penalty 0.1→5000 完全无影响 ⇒ 怀疑与载荷/接触无关，是 deck/BC 的伪影。
本跑：载荷缩到 1e-6 / 1e-3 / 1e-1 / 1 × F_sub（gap=-0.05 保证 t=0 就接触；无 hold 旋钮）。
若连 1e-6 倍载荷都报同样的 `8 negative jacobians` ⇒ **与载荷无关的伪影**。
"""
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT / "src"), str(ROOT / "scripts" / "opensim_fe")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import thums_feb as T  # noqa: E402
from climbing.coupling.febio_run import run_febio  # noqa: E402

MESH = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_anatframe.npz"
WORK = ROOT / "temp" / "opensim_fe" / "v1_scale"
F0 = 26383.133775894516
SCALES = [1e-6, 1e-3, 1e-1, 1.0]


def run(mesh, scale, tag):
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    feb = (wd / f"{tag}.feb").resolve()
    for ext in (".xplt", ".hdf5", ".log"):
        p = feb.with_suffix(ext)
        if p.exists():
            p.unlink()
    T.build_thums_feb(mesh, feb, load_n=F0 * scale, use_rigid=False, time_steps=1,
                      plantar_bc="contact", contact_gap_mm=-0.05,
                      contact_hold_spring_k=None,
                      contact_pad_e_mpa=0.185, contact_pad_thickness_mm=15.4,
                      contact_penalty=5000.0, contact_auto_penalty=0)
    t0 = time.time()
    try:
        run_febio(feb, workdir=wd, timeout=900)
        rc = 0
    except Exception:                              # noqa: BLE001
        rc = 1
    el = time.time() - t0
    lg = feb.with_suffix(".log")
    import re
    txt = lg.read_text(encoding="utf-8", errors="replace") if lg.is_file() else ""
    neg = [int(x) for x in re.findall(r'(\d+)\s+negative jacobians detected', txt)]
    steps = re.findall(r'beginning time step \d+ : ([\d.eE+-]+)', txt)
    return {"scale": scale, "rc": rc, "s": el, "n_neg": len(neg),
            "neg_first": neg[0] if neg else None,
            "steps": len(steps), "last_t": steps[-1] if steps else None}


if __name__ == "__main__":
    mesh = T.load_thums_mesh(MESH)
    print(f"h=5 m  gap=-0.05  pad_E=0.185  penalty=5000  hold=None  STATIC")
    rows = []
    for s in SCALES:
        tag = f"sc{s:g}"
        print(f"### 载荷 ×{s:g}  (F = {F0*s:.4g} N)", flush=True)
        r = run(mesh, s, tag)
        rows.append(r)
        print(f"    rc={r['rc']}  负J 次数={r['n_neg']}  首={r['neg_first']}  "
              f"进入步={r['steps']}  最后 t={r['last_t']}  {r['s']:.0f}s", flush=True)
    out = ROOT / "results" / "opensim_fe" / "v1_load_scale_scan.json"
    out.write_text(json.dumps({"note": "判别 8 negJ 是否与载荷无关", "rows": rows},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print("\n->", out)
