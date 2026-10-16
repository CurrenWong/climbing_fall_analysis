"""★★ 双接触、双初始接触、**零旋钮** —— 物理上最完整的配置。

推理（`joint_contact.py:449` 逐字）：
    # --- 帽的 BC：只放开 y（上/下压），x = z = 0 ------------------------------
    bcs.append(fbc.BCZeroDisplacement(node_set="talus_cap_set", x_dof=1, z_dof=1))
⇒ 距骨帽的 **x/z 本来就全固定**，唯一自由的是 y。
⇒ 只要**距下接触在 t=0 就啮合**（`joint_contact_gap_mm <= 0`），帽的 y 也被定住
⇒ 跟骨被**夹在上下两个刚性面之间**（距下帽 + 地面板）
⇒ **不需要 `hold_spring_k` / `joint_cap_hold_k` 任何旋钮**。

本跑：joint_bc=contact + plantar_bc=contact，双 gap ∈ {0.5/0.5(原), 0/0, -0.02/-0.02}
      hold 两个都 None；penalty=5000；pad_E=0.185（跟垫）；pad 只固定底面。
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
WORK = ROOT / "temp" / "opensim_fe" / "v1_dual"
F0, REF = 26383.133775894516, 186.7563689967737
GRID = [(0.5, 0.5), (0.0, 0.0), (-0.02, -0.02), (-0.02, 0.0), (0.0, -0.02)]


def run(mesh, gj, gp, tag):
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    feb = (wd / f"{tag}.feb").resolve()
    for ext in (".xplt", ".hdf5", ".log"):
        p = feb.with_suffix(ext)
        if p.exists():
            p.unlink()
    T.build_thums_feb(mesh, feb, load_n=F0, use_rigid=False, time_steps=1,
                      plantar_bc="contact", contact_gap_mm=gp,
                      contact_hold_spring_k=None,
                      contact_pad_e_mpa=0.185, contact_pad_thickness_mm=15.4,
                      contact_penalty=5000.0, contact_auto_penalty=0,
                      contact_pad_fix_bottom_only=True,
                      joint_bc="contact", joint_contact_gap_mm=gj,
                      joint_cap_hold_k=None)
    d = feb.read_bytes().decode("iso-8859-1", errors="replace")
    print(f"    deck: <Contact>×{d.count('<Contact')}  "
          f"(plantar + talus_cap)" if d.count('<Contact') else "    deck: 无接触!", flush=True)
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
    out = {"gap_joint": gj, "gap_plantar": gp, "rc": rc, "s": el,
           "n_neg": len(neg), "neg_first": neg[0] if neg else None,
           "steps": len(steps), "last_t": steps[-1] if steps else None}
    if rc == 0:
        vm = _read_domain_vm(feb.with_suffix(".xplt"), "calcaneus", feb.with_suffix(".hdf5"))
        out["sigma_max"] = float(vm.max())
        out["rel"] = float(vm.max()) / REF
    return out


if __name__ == "__main__":
    mesh = T.load_thums_mesh(MESH)
    print(f"h=5 m  **双接触（距下帽 + 跖面板）**  hold=None×2  penalty=5000  pad_E=0.185")
    print(f"参照 fixed = {REF:.4f} MPa\n")
    rows = []
    for gj, gp in GRID:
        tag = f"gj{gj:+g}_gp{gp:+g}"
        print(f"### joint_gap={gj:+g}  plantar_gap={gp:+g}", flush=True)
        r = run(mesh, gj, gp, tag)
        rows.append(r)
        s = (f"✅ 收敛  σ_max={r['sigma_max']:.4f} (fixed×{r['rel']:.4f})"
             if "sigma_max" in r else
             f"❌ 负J×{r['n_neg']} 首={r['neg_first']} 步={r['steps']} t={r['last_t']}")
        print(f"    {s}  ({r['s']:.0f}s)", flush=True)
    out = ROOT / "results" / "opensim_fe" / "v1_dual_contact_scan.json"
    out.write_text(json.dumps({"h_m": 5.0, "rows": rows}, ensure_ascii=False,
                              indent=2), encoding="utf-8")
    print("\n->", out)
