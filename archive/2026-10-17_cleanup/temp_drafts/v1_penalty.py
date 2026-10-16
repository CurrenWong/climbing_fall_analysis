"""★★ 接触罚因子量级检验 —— 手册 §3.14 comment 2：ε_i = E·A/V（auto 计算）。

历史配方 `penalty=0.1`；皮质 E=15000 MPa、单元 ~3 mm ⇒ E/h ≈ **5000**
⇒ 历史值比正解小 ~5 万倍 ⇒ 表面自由穿透 ⇒ 等价无支撑 ⇒ 刚体模态。

本跑（STATIC，gap=-0.05 保证 t=0 有初始接触，`hold_spring_k=None` = 无旋钮）：
  penalty ∈ {0.1(对照), 1, 100, 1e3, 5e3, 1e5} × auto_penalty ∈ {0, 1}
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
from s2_thums_sweep import _read_domain_vm  # noqa: E402
from climbing.coupling.febio_run import run_febio  # noqa: E402

MESH = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_anatframe.npz"
WORK = ROOT / "temp" / "opensim_fe" / "v1_pen"
H, F_SUB, REF, T_PAD = 5.0, 26383.133775894516, 186.7563689967737, 15.4
GAP = -0.05
GRID = [(0.1, 0), (1.0, 0), (100.0, 0), (1000.0, 0), (5000.0, 0),
        (0.1, 1), (1.0, 1), (10.0, 1)]


def run(mesh, pen, ap, tag):
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    feb = (wd / f"{tag}.feb").resolve()
    for ext in (".xplt", ".hdf5", ".log"):
        p = feb.with_suffix(ext)
        if p.exists():
            p.unlink()
    T.build_thums_feb(mesh, feb, load_n=F_SUB, use_rigid=False, time_steps=1,
                      plantar_bc="contact", contact_gap_mm=GAP,
                      contact_hold_spring_k=None,
                      contact_pad_e_mpa=0.185, contact_pad_thickness_mm=T_PAD,
                      contact_penalty=pen, contact_auto_penalty=ap)
    t0 = time.time()
    try:
        rc = run_febio(feb, workdir=wd, timeout=1800)
    except Exception as exc:                       # noqa: BLE001
        return {"penalty": pen, "auto": ap, "rc": "ERR", "s": time.time() - t0}
    el = time.time() - t0
    if rc != 0:
        return {"penalty": pen, "auto": ap, "rc": int(rc), "s": el}
    vm = _read_domain_vm(feb.with_suffix(".xplt"), "calcaneus", feb.with_suffix(".hdf5"))
    return {"penalty": pen, "auto": ap, "rc": 0, "s": el,
            "sigma_max": float(vm.max()), "rel": float(vm.max()) / REF}


if __name__ == "__main__":
    mesh = T.load_thums_mesh(MESH)
    print(f"h={H} m  F_sub={F_SUB:.0f} N  gap={GAP}  pad_E=0.185  hold=None  STATIC")
    print(f"参照 fixed = {REF:.4f} MPa；手册正解 ε=E·A/V ≈ E/h ≈ 5000\n")
    rows = []
    for pen, ap in GRID:
        tag = f"pen{pen:g}_ap{ap:d}"
        print(f"### penalty={pen:g}  auto_penalty={ap}", flush=True)
        r = run(mesh, pen, ap, tag)
        rows.append(r)
        if "sigma_max" in r:
            print(f"    → ✅ 收敛  σ_max={r['sigma_max']:.4f} (fixed×{r['rel']:.4f})"
                  f"  {r['s']:.0f}s", flush=True)
        else:
            print(f"    → ❌ {r.get('rc')}  {r['s']:.0f}s", flush=True)
    out = ROOT / "results" / "opensim_fe" / "v1_penalty_scan.json"
    out.write_text(json.dumps({"h_m": H, "f_sub_n": F_SUB, "gap_mm": GAP,
                               "ref_fixed": REF, "rows": rows},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print("\n->", out)
    ok = [r for r in rows if r.get("rc") == 0]
    print(f"\n收敛 {len(ok)}/{len(rows)}")
