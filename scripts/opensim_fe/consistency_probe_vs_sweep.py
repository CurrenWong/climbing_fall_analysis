# -*- coding: utf-8 -*-
"""一致性审计：探针路径 vs sweep 路径，同一 (h=2, act=0, fixed) 的 deck 是否逐位相同？

**动因**（`V1_PLANTAR_BC_SENSITIVITY.md` §4 的「已知未解」）：同一网格/高度下
探针 `ref_fixed` 给 **362.56 MPa**，sweep 流水线给 **336.92 MPa**，相差 7.6% 未归因。
**先判定差异在「deck」还是在「后处理」** —— 这是决定性的一刀。

用法::

    .venv\\Scripts\\python.exe scripts\\opensim_fe\\consistency_probe_vs_sweep.py
"""
from __future__ import annotations

import hashlib
import io
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(ROOT / "src"))

import thums_feb as T                                  # noqa: E402
import s2_thums_sweep as S2                            # noqa: E402
from climbing.coupling.load_transfer import transfer   # noqa: E402

H = 2.0
ACT = 0.0
OUT = ROOT / "temp" / "opensim_fe" / "consistency"


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16] if p.is_file() else '(missing)'


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    mesh = T.load_thums_mesh()
    print('[mesh] source=%s  nodes=%d cort=%d' % (
        mesh.get('source_npz'), len(mesh['nodes']), len(mesh['elements_cort'])))
    f_sub, t1 = S2._load_at(H)
    spec = transfer(f_sub, achilles_activation=ACT)
    print('[load] h=%g act=%g  subtalar=%.6f N  achilles=%.6f N  dir=%s' % (
        H, ACT, spec.subtalar_n, spec.achilles_n,
        tuple(round(float(x), 6) for x in spec.subtalar_dir)))

    # ---- A：sweep 路径（逐字复刻 `s2_thums_sweep._solve_case` 的 build 调用）----
    fa = OUT / "A_sweep.feb"
    if fa.exists():
        fa.unlink()
    T.build_thums_feb(
        mesh, fa,
        load_n=spec.subtalar_n, load_dir=spec.subtalar_dir,
        achilles_n=spec.achilles_n, achilles_dir=spec.achilles_dir,
        use_rigid=False, time_steps=1,
    )

    # ---- B：探针路径（`plantar_contact_fix_probe._build_one` 的 build 调用）----
    fb = OUT / "B_probe.feb"
    if fb.exists():
        fb.unlink()
    T.build_thums_feb(
        mesh, fb,
        load_n=spec.subtalar_n, load_dir=spec.subtalar_dir,
        achilles_n=spec.achilles_n, achilles_dir=spec.achilles_dir,
        use_rigid=False, time_steps=1,
        plantar_bc="fixed",
    )

    sa, sb = _sha(fa), _sha(fb)
    print()
    print('[A_sweep] %s  %d B  sha16=%s' % (fa, fa.stat().st_size, sa))
    print('[B_probe] %s  %d B  sha16=%s' % (fb, fb.stat().st_size, sb))
    print('逐位相同:', sa == sb)
    if sa != sb:
        ta = fa.read_bytes().decode('iso-8859-1', errors='replace').split('\r\n')
        tb = fb.read_bytes().decode('iso-8859-1', errors='replace').split('\r\n')
        print('\n--- 首个不同行 ---')
        for i in range(max(len(ta), len(tb))):
            x = ta[i] if i < len(ta) else '<EOF>'
            y = tb[i] if i < len(tb) else '<EOF>'
            if x != y:
                print('  line %d' % (i + 1))
                print('   A| %s' % x[:150])
                print('   B| %s' % y[:150])
                break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
