"""A4 · V6-A 机制判据：滑动卸载是否随内翻角单调增强。

不依赖 ANKLE_AVULSION_LOAD_N（215 N）、不依赖 [T22] 时标。
判据：同一高度下，`ankle_slide_frac`(30°) > `ankle_slide_frac`(10°)，且 30° 处 mode=sliding。
"""
import math
import sys

sys.path.insert(0, r'D:\Project\climbing_fall_analysis\src')
import climbing.pad as P
from climbing.pad2d import Posture2D, simulate_boulder_fall_2d

ANKLE_D = dict(mu_slide=0.4, ankle_d_lateral_m=0.030,
               ankle_d_medial_m=0.045, ligament_r_mm=22.0)
BETAS = (10.0, 20.0, 30.0)


def run(h, beta, on_pad, musc):
    old = dict(P.TWO_SEG_LEG)
    P.TWO_SEG_LEG = {} if musc is None else {'m_ecc_max_nm': musc, 'ecc_rate_ref_ms': 1.0}
    try:
        p2 = Posture2D(ankle_beta0_rad=math.radians(beta), **ANKLE_D)
        return simulate_boulder_fall_2d(height_m=h, mass_kg=80.0,
                                        posture='controlled-drop',
                                        on_pad=on_pad, posture2d=p2)
    finally:
        P.TWO_SEG_LEG.clear(); P.TWO_SEG_LEG.update(old)


print('V6-A：滑动卸载份额 slide_frac（越高＝卸载越强）')
print(f'{"地面":<4}{"m_ecc":>7}{"h":>5}  ' + '  '.join(f'{b:>4.0f}°' for b in BETAS) + '   判据')
for on_pad in (False, True):
    for musc in (None, 210.0):
        for h in (5.0, 10.0, 20.0):
            vals, modes = [], []
            for b in BETAS:
                r = run(h, b, on_pad, musc)
                vals.append(r.ankle_slide_frac)
                modes.append(r.ankle_mode)
            ok = vals[2] > vals[0] and modes[2] == 'sliding'
            tag = '垫' if on_pad else '刚'
            print(f'{tag:<4}{str(musc or 400):>7}{h:>5.0f}  '
                  + '  '.join(f'{v:5.3f}' for v in vals)
                  + f'   {"PASS" if ok else "FAIL"}  modes={modes}')
