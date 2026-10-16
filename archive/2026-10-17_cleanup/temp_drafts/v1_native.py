"""验证：原生 flip_primary/flip_secondary 与我上一轮的**文本反转节点序**是否等价。

上一轮（文本后处理）结果：
  原绕向 + hold100              -> 151.896519（= spring k=100，惰性）
  两个都反转 + hold100          -> **152.491512**（接触生效）
本跑（原生 flip_*）应复现同一对数。
"""
import re
import sys
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
WORK = ROOT / "temp" / "opensim_fe" / "v1_native"
F0 = 26383.133775894516
EXP_ORIG, EXP_BOTH = 151.896519, 152.491512


def solve(mesh, tag, fp, fs, hold=100.0, gap=0.5):
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    feb = (wd / f"{tag}.feb").resolve()
    for ext in (".xplt", ".hdf5", ".log"):
        p = feb.with_suffix(ext)
        if p.exists():
            p.unlink()
    kw = {}
    if fp is not None:
        kw["contact_flip_primary"] = fp
    if fs is not None:
        kw["contact_flip_secondary"] = fs
    T.build_thums_feb(mesh, feb, load_n=F0, use_rigid=False, time_steps=1,
                      plantar_bc="contact", contact_gap_mm=gap,
                      contact_hold_spring_k=hold, **kw)
    d = feb.read_bytes().decode("iso-8859-1", errors="replace")
    got = re.findall(r'<flip_(primary|secondary)>(\d)</flip_(?:primary|secondary)>', d)
    print(f"    [{tag}] deck flip 标签: {got}", flush=True)
    try:
        run_febio(feb, workdir=wd, timeout=1200)
    except Exception:                              # noqa: BLE001
        return None
    return _read_domain_vm(feb.with_suffix(".xplt"), "calcaneus",
                           feb.with_suffix(".hdf5"))


mesh = T.load_thums_mesh(MESH)
print("验证原生 flip_* 与文本反转等价（h=5 m, hold_k=100, gap=+0.5）\n")
cases = [("orig", None, None), ("fp1", 1, 0), ("fs1", 0, 1), ("fp1fs1", 1, 1)]
res = {}
for tag, fp, fs in cases:
    v = solve(mesh, tag, fp, fs)
    res[tag] = v
    print(f"  {tag:8s} " + ("❌ 不收敛" if v is None else
                            f"max={v.max():.6f} 中位={np.median(v):.6f}"), flush=True)
print()
if res["orig"] is not None:
    print(f"  orig  == 151.896519（文本原绕向）？ "
          f"{abs(res['orig'].max()-EXP_ORIG) < 1e-4}")
if res["fp1fs1"] is not None:
    print(f"  fp1fs1 == 152.491512（文本两反转）？ "
          f"{abs(res['fp1fs1'].max()-EXP_BOTH) < 1e-4}")
    print(f"  fp1fs1 与 orig 逐位相同？ "
          f"{np.allclose(res['fp1fs1'], res['orig'], rtol=0, atol=0)}"
          "   （应为 False ⇒ 接触生效）")
