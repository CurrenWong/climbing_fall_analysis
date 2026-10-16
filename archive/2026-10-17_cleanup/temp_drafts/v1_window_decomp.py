"""第 2 步：定性软接触下 `F_flex(h)` 的**收敛非单调** —— 在**冲击窗内**把 f_flex 拆成弹簧/阻尼项。
（上一次分解用了全轨迹，与 f_flex_peak_n 不同口径 ⇒ 不可比。）只读 + 临时门。"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

G, LEG0, C, KSTOP = 9.80665, 0.85, 400.0, 3.0e5
post = P.POSTURES["anchored-self"]
XM = post.x_flex_max
K = post.k_flex
print(f"anchored-self: k_flex={K/1e3:.1f} kN/m, x_flex_max={XM} m, "
      f"k_stop={KSTOP:.1e}, c_flex={C}")
print("（k_stop 段起点力 = k·xm = %.1f kN）\n" % (K * XM / 1e3))

o = P.RIGID_CONTACT["k_n_per_m"]
try:
    for kn in (None, 2.1e6):
        P.RIGID_CONTACT["k_n_per_m"] = kn
        tag = "None" if kn is None else f"{kn:.1e}"
        print(f"=== k_n = {tag} ===")
        print(f"  {'h(m)':>5} {'窗':>16} {'f_flex峰':>9} {'弹簧峰':>9} {'阻尼峰':>9} "
              f"{'t*峰(s)':>9} {'xf@峰':>8} {'>xm?':>6} {'窗口maxXf':>10}")
        for h in (1.0, 2.0, 3.0, 5.0):
            r = P.simulate_boulder_fall(height_m=h, mass_kg=80.0, on_pad=False,
                                        posture="anchored-self")
            t = np.asarray(r.t_s, float)
            zf, zt = np.asarray(r.z_foot_m, float), np.asarray(r.z_torso_m, float)
            t0, t1 = r.meta["impact_window_s"]
            i0, i1 = int(np.searchsorted(t, t0)), int(np.searchsorted(t, t1))
            i1 = max(i1, i0 + 2)
            xf = np.maximum(0.0, LEG0 - (zf - zt))
            rate = np.maximum(0.0, np.gradient(xf, t))
            f_spr = np.where(xf <= XM, K * xf, K * XM + KSTOP * (xf - XM))
            f_dmp = C * rate
            f_tot = np.where(xf > 0, f_spr + f_dmp, 0.0)
            sl = slice(i0, i1)
            i_star = i0 + int(np.argmax(f_tot[sl]))
            print(f"  {h:>5.1f} {t0:>7.4f}-{t1:<8.4f} {f_tot[sl].max()/1e3:>9.3f} "
                  f"{f_spr[sl].max()/1e3:>9.3f} {f_dmp[sl].max()/1e3:>9.3f} "
                  f"{t[i_star]:>9.4f} {xf[i_star]:>8.4f} "
                  f"{('是' if xf[i_star] > XM else '否'):>6} {xf[sl].max():>10.4f}")
            print(f"        pad.py 自报 f_flex_peak_n = {r.meta['f_flex_peak_n']/1e3:.3f} kN"
                  f"  （与上面 f_flex峰 应一致）")
        print()
finally:
    P.RIGID_CONTACT["k_n_per_m"] = o
    print(f"(已恢复 {o!r})")
