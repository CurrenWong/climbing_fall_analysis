"""★ 真实跟垫（可变形板）替代刚性地面 —— 手册 §8.12.4 的物理实现。

诊断基础（上一轮）：`contact` + 刚性板（`pad_e_mpa=200000`）在 DYNAMIC 下真积分到 t≈0.44
后 **460~520 个单元倒置** ⇒ 怀疑是**刚性地面把跖面局部压溃**。
真实跟垫是 **E≈0.185 MPa / t≈15.4 mm** 的软层（G5 已标定，PMC8925218）——
它正是防这件事的。而 `contact_pad_e_mpa` 默认 200000（数值刚性），**这条路从未跑过**。

本跑：`hold_spring_k=None`（无旋钮）
  板模量 ∈ {200000（对照·刚性）, 5.0, 0.185（跟垫）}
  gap  ∈ {+0.5（原）, -0.05（预接触）}
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
WORK = ROOT / "temp" / "opensim_fe" / "v1_padmat"
H, F_SUB, REF = 5.0, 26383.133775894516, 186.7563689967737
T_PAD = 15.4            # 跟垫厚度（G5）
GRID = [(0.5, 200000.0), (0.5, 0.185), (-0.05, 200000.0), (-0.05, 0.185),
        (-0.05, 5.0)]


def run(mesh, gap, e_pad, tag):
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    feb = (wd / f"{tag}.feb").resolve()
    for ext in (".xplt", ".hdf5", ".log"):
        p = feb.with_suffix(ext)
        if p.exists():
            p.unlink()
    T.build_thums_feb(mesh, feb, load_n=F_SUB, use_rigid=False, time_steps=1,
                      plantar_bc="contact", contact_gap_mm=gap,
                      contact_hold_spring_k=None,
                      contact_pad_e_mpa=e_pad, contact_pad_thickness_mm=T_PAD)
    d = feb.read_bytes().decode("iso-8859-1", errors="replace")
    import re
    m = re.search(r'name="pad_floor_mat".*?<E>([^<]+)</E>', d, re.S)
    print(f"    deck pad_floor E={m.group(1) if m else '?'}  "
          f"gap={gap:+g}", flush=True)
    t0 = time.time()
    try:
        rc = run_febio(feb, workdir=wd, timeout=1800)
    except Exception as exc:                       # noqa: BLE001
        return {"pad_e": e_pad, "gap": gap, "rc": "ERR", "s": time.time() - t0,
                "msg": str(exc)[-300:]}
    el = time.time() - t0
    if rc != 0:
        return {"pad_e": e_pad, "gap": gap, "rc": int(rc), "s": el}
    vm = _read_domain_vm(feb.with_suffix(".xplt"), "calcaneus", feb.with_suffix(".hdf5"))
    return {"pad_e": e_pad, "gap": gap, "rc": 0, "s": el,
            "sigma_max": float(vm.max()), "rel": float(vm.max()) / REF}


if __name__ == "__main__":
    mesh = T.load_thums_mesh(MESH)
    print(f"h={H} m  F_sub={F_SUB:.0f} N  hold_spring_k=None  pad_t={T_PAD} mm")
    print(f"参照 fixed = {REF:.4f} MPa\n")
    rows = []
    for gap, e in GRID:
        tag = f"gap{gap:+g}_e{e:g}"
        print(f"### gap={gap:+g}  pad_E={e:g} MPa", flush=True)
        r = run(mesh, gap, e, tag)
        rows.append(r)
        if "sigma_max" in r:
            print(f"    → ✅ σ_max={r['sigma_max']:.4f} (fixed×{r['rel']:.4f})"
                  f"  {r['s']:.0f}s", flush=True)
        else:
            print(f"    → ❌ {r.get('rc')}  {r['s']:.0f}s", flush=True)
            if r.get("msg"):
                print("      " + r["msg"].replace("\n", " ")[-200:], flush=True)
    out = ROOT / "results" / "opensim_fe" / "v1_pad_material_scan.json"
    out.write_text(json.dumps({"h_m": H, "f_sub_n": F_SUB, "pad_t_mm": T_PAD,
                               "ref_fixed": REF, "rows": rows},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print("\n->", out)
