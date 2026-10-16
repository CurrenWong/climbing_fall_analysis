"""★★ 检验我自己的「接触惰性」结论是否成立。

上一轮：`contact`+`hold_k=100`（**gap=+0.5 初始分离**）与 `spring k=100` **逐位相同**
=> 我判「接触惰性」。**但更简约的解释是：gap=+0.5 时骨被弹簧撑住、根本没碰到板。**

决定性检验：把 gap 改成**负**（强制啮合），再与 spring 对照。
  A: contact, gap=-0.05, hold_k=100   ← 强制啮合
  B: contact, gap=+0.5,  hold_k=100   ← 上轮（= spring 的那个）
  C: spring,  k=100                   ← 对照
若 A != C ⇒ **接触是真实的**，上轮「惰性」结论**错误，必须撤回**。
若 A == C ⇒ 接触确实惰性。
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
WORK = ROOT / "temp" / "opensim_fe" / "v1_inert"
F0 = 26383.133775894516


def solve(mesh, tag, **kw):
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    feb = (wd / f"{tag}.feb").resolve()
    for ext in (".xplt", ".hdf5", ".log"):
        p = feb.with_suffix(ext)
        if p.exists():
            p.unlink()
    T.build_thums_feb(mesh, feb, load_n=F0, use_rigid=False, time_steps=1, **kw)
    try:
        run_febio(feb, workdir=wd, timeout=1200)
    except Exception:                              # noqa: BLE001
        return None
    return _read_domain_vm(feb.with_suffix(".xplt"), "calcaneus",
                           feb.with_suffix(".hdf5"))


mesh = T.load_thums_mesh(MESH)
print("检验「接触惰性」结论（h=5 m, F=26383.1 N）\n")
A = solve(mesh, "A_contact_gap-0.05_hold100", plantar_bc="contact",
          contact_gap_mm=-0.05, contact_hold_spring_k=100.0)
B = solve(mesh, "B_contact_gap+0.5_hold100", plantar_bc="contact",
          contact_gap_mm=0.5, contact_hold_spring_k=100.0)
C = solve(mesh, "C_spring_100", plantar_bc="spring", spring_k=100.0)


def show(nm, v):
    print(f"  {nm:34s} " + ("None（不收敛）" if v is None else
                            f"max={v.max():.6f} 中位={np.median(v):.6f}"))


show("A contact gap=-0.05 hold_k=100", A)
show("B contact gap=+0.5  hold_k=100", B)
show("C spring  k=100", C)
print()
if A is not None and C is not None:
    same = np.allclose(A, C, rtol=0, atol=0)
    print(f"  **A(强制啮合) == C(spring)？ {same}**")
    if not same:
        d = np.abs(A - C)
        print(f"     最大差 {d.max():.6f}  相对 {d.max()/max(C.max(),1e-30):.4%}  "
              f"相关 {np.corrcoef(A, C)[0,1]:.6f}")
        print("     ⇒ **接触是真实的！上轮「惰性」结论错误，须撤回。**")
    else:
        print("     ⇒ 接触确实惰性。")
if B is not None and C is not None:
    print(f"  B(gap=+0.5) == C(spring)？ {np.allclose(B, C, rtol=0, atol=0)}  "
          f"（上轮的对照条件）")
