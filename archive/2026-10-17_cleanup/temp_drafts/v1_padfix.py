"""★ 验证：把对偶板从「完全刚体块」改成「可变形层 + 刚性地基」（只固定底面）。

发现的缺陷：`_build_pad_floor` 原先把**上下两排全部固定**
（`sorted(top_ids)+sorted(bot_ids)`）⇒ 板是完全刚体 ⇒ `pad_e_mpa` 形同无效
（实测 0.185 MPa 跟垫与 200 GPa 刚性板给出逐位相同的失败签名）。

本跑（`pad_fix_bottom_only=True`，opt-in）：
  板模量 ∈ {200000（刚性层）, 0.185（跟垫）}, gap ∈ {+0.05, -0.005}, penalty ∈ {0.1, 5000}
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
WORK = ROOT / "temp" / "opensim_fe" / "v1_padfix"
F0, REF = 26383.133775894516, 186.7563689967737
GRID = [(0.05, 0.185, 0.1), (0.05, 0.185, 5000.0), (0.05, 200000.0, 5000.0),
        (-0.005, 0.185, 5000.0)]


def run(mesh, gap, e, pen, tag):
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
                      contact_pad_e_mpa=e, contact_pad_thickness_mm=15.4,
                      contact_penalty=pen, contact_auto_penalty=0,
                      contact_pad_fix_bottom_only=True)
    d = feb.read_bytes().decode("iso-8859-1", errors="replace")
    nfix = re.search(r'<NodeSet name="pad_floor_fixed">(.*?)</NodeSet>', d, re.S)
    nfix = len(nfix.group(1).split(',')) if nfix else -1
    print(f"    固定节点数 = {nfix}（原为 50）", flush=True)
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
    out = {"gap": gap, "pad_e": e, "penalty": pen, "rc": rc, "s": el,
           "nfixed": nfix, "n_neg": len(neg), "neg_first": neg[0] if neg else None,
           "steps": len(steps), "last_t": steps[-1] if steps else None}
    if rc == 0:
        vm = _read_domain_vm(feb.with_suffix(".xplt"), "calcaneus", feb.with_suffix(".hdf5"))
        out["sigma_max"] = float(vm.max())
        out["rel"] = float(vm.max()) / REF
    return out


if __name__ == "__main__":
    mesh = T.load_thums_mesh(MESH)
    print(f"h=5 m  pad_fix_bottom_only=True  hold=None  STATIC")
    print(f"参照 fixed = {REF:.4f} MPa\n")
    rows = []
    for gap, e, pen in GRID:
        tag = f"g{gap:+g}_e{e:g}_p{pen:g}"
        print(f"### gap={gap:+g} pad_E={e:g} penalty={pen:g}", flush=True)
        r = run(mesh, gap, e, pen, tag)
        rows.append(r)
        s = (f"✅ σ_max={r['sigma_max']:.4f} (fixed×{r['rel']:.4f})"
             if "sigma_max" in r else
             f"❌ 负J×{r['n_neg']} 首={r['neg_first']} 步={r['steps']} t={r['last_t']}")
        print(f"    {s}  ({r['s']:.0f}s)", flush=True)
    out = ROOT / "results" / "opensim_fe" / "v1_pad_fix_bottom_only_scan.json"
    out.write_text(json.dumps({"h_m": 5.0, "rows": rows}, ensure_ascii=False,
                              indent=2), encoding="utf-8")
    print("\n->", out)
