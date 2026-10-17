# -*- coding: utf-8 -*-
"""一致性审计（第 2 步）：**同一份 deck** 过两套路径，逐层定位 7.6% 差异。

第 1 步已证：探针路径与 sweep 路径生成的 deck **逐位相同**（只差 plotfile 文件名）
⇒ 差异必在 **求解** 或 **后处理**。本件把同一 deck 分别用两套调用求解+后处理，
并**逐层打印中间量**（vm.max / vm.mean / gauge_max / 域名 / 单元数），定位到具体一层。

用法::

    .venv\\Scripts\\python.exe scripts\\opensim_fe\\consistency_probe_vs_sweep.py --deep
"""
from __future__ import annotations

import hashlib
import io
import subprocess
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np                                     # noqa: E402
import thums_feb as T                                  # noqa: E402
import s2_thums_sweep as S2                            # noqa: E402
from climbing.coupling.load_transfer import transfer   # noqa: E402

H = 2.0
ACT = 0.0
OUT = ROOT / "temp" / "opensim_fe" / "consistency"


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16] if p.is_file() else '(missing)'


def _solve(deck: Path, exe) -> float:
    for s in (".xplt", ".hdf5", ".log"):
        p = deck.with_suffix(s)
        if p.exists():
            p.unlink()
    t0 = time.time()
    pr = subprocess.run([str(exe), "-i", str(deck)], capture_output=True,
                        text=True, cwd=str(deck.parent), timeout=1800)
    return pr.returncode, time.time() - t0


def main() -> int:
    from climbing.coupling.febio_run import find_febio
    from climbing.coupling import fe_post
    from pyfebio import xplt as xplt_mod

    OUT.mkdir(parents=True, exist_ok=True)
    mesh = T.load_thums_mesh()
    f_sub, t1 = S2._load_at(H)
    spec = transfer(f_sub, achilles_activation=ACT)
    cort_c, cort_v = T.hex_centroids_volumes(mesh["nodes"], mesh["elements_cort"])
    gauge_mm = float(getattr(S2, "REG_LEN_MM", 4.0))
    print('[ctx] mesh=%s  cort=%d  gauge_mm=%g  subtalar=%.6f N' % (
        Path(mesh['source_npz']).name, len(cort_c), gauge_mm, spec.subtalar_n))

    deck = OUT / "C_common.feb"
    if deck.exists():
        deck.unlink()
    T.build_thums_feb(mesh, deck, load_n=spec.subtalar_n,
                      load_dir=spec.subtalar_dir, achilles_n=spec.achilles_n,
                      achilles_dir=spec.achilles_dir, use_rigid=False,
                      time_steps=1, plantar_bc="fixed")
    print('[deck] %s  %d B  sha16=%s' % (deck.name, deck.stat().st_size, _sha(deck)))

    exe = find_febio()
    rc, el = _solve(deck, exe)
    xplt, hdf5 = deck.with_suffix(".xplt"), deck.with_suffix(".hdf5")
    print('[solve] rc=%s  %.1fs  xplt=%d B' % (rc, el, xplt.stat().st_size))

    print()
    print('--- 层 1：两套 to_hdf5 调用 → h5 结构 ---')
    xplt_mod.to_hdf5(str(xplt), str(hdf5))
    import h5py
    with h5py.File(hdf5, "r") as h:
        st = sorted(h["states"].keys(), key=int)
        doms = list(h["states"][st[-1]]["element_data"]["stress"].keys())
        print('   states=%s  stress 域=%s' % (st, doms))

    print()
    print('--- 层 2：两套读场函数 → vm 统计 ---')
    vm_sweep = S2._read_domain_vm(xplt, "calcaneus", hdf5)
    stress, used = fe_post._read_last_stress(h5py.File(hdf5, "r"), "calcaneus")
    vm_direct = fe_post.von_mises_from_voigt(stress)
    print('   S2._read_domain_vm : n=%d  max=%.6f  mean=%.6f  used=%s' % (
        len(vm_sweep), vm_sweep.max(), vm_sweep.mean(), used))
    print('   fe_post 直接读      : n=%d  max=%.6f  mean=%.6f  used=%s' % (
        len(vm_direct), vm_direct.max(), vm_direct.mean(), used))
    print('   两路 vm 逐位相同:', bool(np.array_equal(vm_sweep, vm_direct)))

    print()
    print('--- 层 3：gauge（含 ℓ 正则化） ---')
    g1 = T.gauge_von_mises_mixed(vm_sweep, cort_c, cort_v, gauge_mm)
    g2 = T.gauge_von_mises_mixed(vm_direct, cort_c, cort_v, gauge_mm)
    print('   gauge_max(vm_sweep)  = %.6f' % g1.max())
    print('   gauge_max(vm_direct) = %.6f' % g2.max())
    c1 = int(np.argmax(g1))
    print('   gauge 峰值单元 idx=%d  质心=%s' % (
        c1, [round(float(x), 3) for x in cort_c[c1]]))
    print('   vm 峰值单元 idx=%d  质心=%s' % (
        int(np.argmax(vm_sweep)), [round(float(x), 3) for x in cort_c[int(np.argmax(vm_sweep))]]))

    print()
    print('--- 层 4：peak_von_mises（公共入口，会自己 to_hdf5）---')
    pk = fe_post.peak_von_mises(xplt, hdf5_path=deck.with_suffix(".pk.hdf5"),
                                domain="calcaneus")
    print('   peak_von_mises = %.6f   （vm.max = %.6f）' % (pk, vm_sweep.max()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
