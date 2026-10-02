import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from climbing import G
from climbing.rope import simulate_rope_fall, default_rope
from climbing.belay import get_device
from climbing.pad import simulate_boulder_fall, POSTURES
from climbing.metrics import hic_from_result, rope_risk, pad_risk

print("=== 绳索: UIAA 条件 80kg, FF=1.77, L0=2m, ATC-braked ===")
r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.77,
                       device="ATC-braked", t_hold_n=get_device("ATC-braked").t_hold_n)
r.hic = hic_from_result(r)
print(rope_risk(r))
# 绳绷紧时的速度应等于**自由落差**那一段的自由落体速度，
# 不是总落差。FF=1.77 时只有 (1.77-1)*2 = 1.54 m 是自由落体。
v_taut = np.sqrt(2 * G * r.free_fall_m)
print(f"  绷紧时速度: 期望 {v_taut:.2f} m/s, 仿真 {np.abs(r.v_ms).max():.2f} m/s")

print()
print("=== 绳索: 保护器对比 (80kg, FF=1.77, L0=2m) ===")
for name in ["ATC-braked", "ATC-nobrake", "Reverso", "Grigri", "Figure8-passive"]:
    d = get_device(name)
    rr = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.77,
                            device=name, t_hold_n=d.t_hold_n)
    rr.hic = hic_from_result(rr)
    print(f"  {name:18s} peak={rr.peak_tension_kn:6.2f} kN  g={rr.peak_accel_g:6.1f}  "
          f"drop={rr.max_drop_m:5.2f}m  payout={rr.max_payout_m:5.2f}m  slide={rr.total_rope_slide_m:5.2f}m")

print()
print("=== 抱石: 3m / 80kg / 20cm pad ===")
for name in POSTURES:
    b = simulate_boulder_fall(height_m=3.0, mass_kg=80, posture=name)
    b.hic = hic_from_result(b)
    d = pad_risk(b)
    print(f"  {name:20s} F={d['peak_kN']:7.2f} kN  torso={d['torso_g']:6.1f}g  "
          f"leg={d['leg_g']:6.1f}g  comp={d['comp_cm']:5.1f}cm  bottom={d['bottomed']}  "
          f"absorbed={d['absorbed']*100:5.1f}%  resid={d['residual_kJ']:.2f}kJ  HIC={b.hic:7.1f}")
