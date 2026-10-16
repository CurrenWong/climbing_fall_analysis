"""★★ 在**真 deck** 上验证 MRE 指出的机制：反转 `pad_top` 的面片绕向。

MRE 结论：`c2_bot` 绕向反转 ⇒ 接触从「惰性/飞走」变成「正常工作」（穿透 0.0038 mm）。
真 deck 现状：`pad_top` 法向 (0,+1,0)、`plantar` 法向 Y 中位 −0.796 ⇒ **两面相向**
—— 正是 MRE 里**失败**的那个组合。

检验：反转 `pad_top` 的面片绕向后，
  A) `contact` + hold_k=100 是否**不再**等于 `spring k=100`（⇒ 接触开始起作用）
  B) `contact` + hold_k=None 是否能收敛（⇒ 不需要旋钮）
"""
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT / "src"), str(ROOT / "scripts" / "opensim_fe")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402
import thums_feb as T  # noqa: E402
from s2_thums_sweep import _read_domain_vm  # noqa: E402
from climbing.coupling.febio_run import run_febio  # noqa: E402

MESH = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_anatframe.npz"
WORK = ROOT / "temp" / "opensim_fe" / "v1_flip"
F0 = 26383.133775894516


def flip_padtop(txt: str) -> str:
    m = re.search(r'(<Surface name="pad_top">)(.*?)(</Surface>)', txt, re.S)
    if not m:
        return txt
    body = re.sub(r'<quad4 id="(\d+)">(.*?)</quad4>',
                  lambda mm: f'<quad4 id="{mm.group(1)}">'
                  + ",".join(reversed(mm.group(2).split(','))) + '</quad4>',
                  m.group(2))
    return txt[:m.start(2)] + body + txt[m.end(2):]


def solve(mesh, tag, flip, **kw):
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    feb = (wd / f"{tag}.feb").resolve()
    for ext in (".xplt", ".hdf5", ".log"):
        p = feb.with_suffix(ext)
        if p.exists():
            p.unlink()
    T.build_thums_feb(mesh, feb, load_n=F0, use_rigid=False, time_steps=1, **kw)
    if flip:
        feb.write_bytes(flip_padtop(
            feb.read_bytes().decode("iso-8859-1")).encode("iso-8859-1"))
    try:
        run_febio(feb, workdir=wd, timeout=1200)
    except Exception:                              # noqa: BLE001
        return None
    return _read_domain_vm(feb.with_suffix(".xplt"), "calcaneus",
                           feb.with_suffix(".hdf5"))


mesh = T.load_thums_mesh(MESH)
print("在真 deck 上验证绕向机制（h=5 m, F=26383.1 N）\n")
C = solve(mesh, "spring_100", False, plantar_bc="spring", spring_k=100.0)
A = solve(mesh, "flip_hold100", True, plantar_bc="contact", contact_gap_mm=0.5,
          contact_hold_spring_k=100.0)
B = solve(mesh, "flip_nohold", True, plantar_bc="contact", contact_gap_mm=0.5,
          contact_hold_spring_k=None)
N = solve(mesh, "noflip_nohold", False, plantar_bc="contact", contact_gap_mm=0.5,
          contact_hold_spring_k=None)


def sh(nm, v):
    print(f"  {nm:34s} " + ("❌ 不收敛" if v is None else
                            f"max={v.max():.6f} 中位={np.median(v):.6f}"))


sh("C  spring k=100（对照）", C)
sh("A  flip_pad_top + hold100", A)
sh("B  flip_pad_top + **无 hold**", B)
sh("N  原绕向 + 无 hold", N)
print()
if A is not None and C is not None:
    same = np.allclose(A, C, rtol=0, atol=0)
    print(f"  A（反转绕向 + hold100）== C（spring）？ **{same}**")
    if not same:
        print(f"     最大差 {np.abs(A-C).max():.6f} 相对 {np.abs(A-C).max()/max(C.max(),1e-30):.3%}")
        print("     ⇒ **接触开始起作用了！绕向就是机制。**")
print(f"  B（反转绕向 + 无 hold）收敛？ **{B is not None}**"
      + ("" if B is None else "  ← 若 True ⇒ 终于不需要旋钮"))
