"""★ 组合拳：手册 §8.4.2（DYNAMIC + `<solver>` 阻尼）+ 真实跟垫（可变形板）。

上一轮已知：
  * STATIC + 刚性板 + 无 hold 弹簧 ⇒ **1 s 秒退（8 negJ，刚体模态）**
  * DYNAMIC + 阻尼 + 刚性板    ⇒ **真积分到 t≈0.44，然后 460~520 单元倒置**
  * 板模量在 STATIC 下**无影响**（失败模式是刚体模态，不是压溃）

本跑把两者合起来：DYNAMIC + 阻尼 + 板模量 ∈ {200000(刚性), 0.185(跟垫)}，gap ∈ {+0.5, -0.05}
若跟垫能消掉倒置 ⇒ **「跖面按真实跟垫承载」在修正路径后是可行的**（推翻上一轮的「不可行」裁定）。
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
WORK = ROOT / "temp" / "opensim_fe" / "v1_combo"
H, F_SUB, REF, T_PAD = 5.0, 26383.133775894516, 186.7563689967737, 15.4
SOLVER_X = "<rhoi>-2</rhoi>\n<gamma>1.5</gamma>\n<beta>1</beta>\n"
GRID = [(-0.05, 200000.0, 0), (-0.05, 0.185, 0),
        (-0.05, 200000.0, 1), (-0.05, 0.185, 1), (0.5, 0.185, 1)]


def patch(txt):
    txt = txt.replace("<analysis>STATIC</analysis>", "<analysis>DYNAMIC</analysis>")
    txt = txt.replace("<time_steps>1</time_steps>", "<time_steps>20</time_steps>")
    txt = txt.replace("<step_size>1.0</step_size>", "<step_size>0.05</step_size>")
    if "<solver/>" in txt:
        txt = txt.replace("<solver/>", "<solver>\n" + SOLVER_X + "</solver>")
    else:
        txt = txt.replace("<solver>", "<solver>\n" + SOLVER_X, 99)
    return txt


def run(mesh, gap, e_pad, tag, ap=0):
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
                      contact_pad_e_mpa=e_pad, contact_pad_thickness_mm=T_PAD,
                      contact_auto_penalty=int(ap))
    feb.write_bytes(patch(feb.read_bytes().decode("iso-8859-1")).encode("iso-8859-1"))
    t0 = time.time()
    try:
        rc = run_febio(feb, workdir=wd, timeout=2400)
    except Exception as exc:                       # noqa: BLE001
        return {"pad_e": e_pad, "gap": gap, "rc": "ERR", "s": time.time() - t0,
                "msg": str(exc)[-400:]}
    el = time.time() - t0
    if rc != 0:
        return {"pad_e": e_pad, "gap": gap, "rc": int(rc), "s": el}
    vm = _read_domain_vm(feb.with_suffix(".xplt"), "calcaneus", feb.with_suffix(".hdf5"))
    return {"pad_e": e_pad, "gap": gap, "rc": 0, "s": el,
            "sigma_max": float(vm.max()), "rel": float(vm.max()) / REF}


if __name__ == "__main__":
    mesh = T.load_thums_mesh(MESH)
    print(f"h={H} m  F_sub={F_SUB:.0f} N  hold=None  DYNAMIC+阻尼  pad_t={T_PAD} mm")
    print(f"参照 fixed = {REF:.4f} MPa\n")
    rows = []
    for gap, e, ap in GRID:
        tag = f"g{gap:+g}_e{e:g}_ap{ap:d}"
        print(f"### DYNAMIC+阻尼  gap={gap:+g}  pad_E={e:g} MPa  auto_penalty={ap}",
              flush=True)
        r = run(mesh, gap, e, tag, ap)
        rows.append(r)
        if "sigma_max" in r:
            print(f"    → ✅ σ_max={r['sigma_max']:.4f} (fixed×{r['rel']:.4f})  {r['s']:.0f}s",
                  flush=True)
        else:
            print(f"    → ❌ {r.get('rc')}  {r['s']:.0f}s", flush=True)
            if r.get("msg"):
                print("      " + r["msg"].replace("\n", " ")[-160:], flush=True)
    out = ROOT / "results" / "opensim_fe" / "v1_dyn_pad_combo_scan.json"
    out.write_text(json.dumps({"h_m": H, "f_sub_n": F_SUB, "pad_t_mm": T_PAD,
                               "ref_fixed": REF, "rows": rows},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print("\n->", out)
