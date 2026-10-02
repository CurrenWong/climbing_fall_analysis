"""绳索模型验证：FF 单调性 / UIAA 12 kN 判据 / 保护器对比。"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from climbing.rope import simulate_rope_fall, default_rope
from climbing.belay import get_device, DEVICES, belayer_side_force
from climbing.metrics import hic_from_result

print("=" * 88)
print("1) UIAA 101 条件 (80 kg, FF=1.77, L0=2 m) 标定结果")
print("=" * 88)
r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.77,
                       device="ATC-braked", t_hold_n=1e9)
print(f"  峰值力      {r.peak_tension_kn:6.2f} kN   (UIAA 实测 8-9 kN)")
print(f"  峰值加速度  {r.peak_accel_g:6.1f} g")
print(f"  最大伸长    {r.max_elongation_pct:6.1f} %   (UIAA 实测 30-40%)")
print(f"  冲击位移    {r.max_drop_m:6.2f} m")
print(f"  UIAA 12kN  {'通过' if r.is_uiaa_pass else '不通过'}")
print(f"  绷紧时绳速  {np.sqrt(2*9.80665*r.free_fall_m):.2f} m/s")

print()
print("=" * 88)
print("2) FF 单调性（锁死保护器，80 kg, L0=2 m）")
print("=" * 88)
print(f"{'FF':>6s} {'自由落差m':>9s} {'峰值kN':>8s} {'峰值g':>7s} {'位移m':>7s} {'伸长%':>7s} {'UIAA':>6s}")
prev = -1.0
mono = True
for ff in [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.77, 2.0]:
    rr = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=ff,
                            t_hold_n=1e9, device="locked")
    ok = "OK" if rr.peak_tension_kn >= prev - 1e-6 else "非单调!"
    if rr.peak_tension_kn < prev - 1e-6:
        mono = False
    prev = rr.peak_tension_kn
    print(f"{ff:6.2f} {rr.free_fall_m:9.2f} {rr.peak_tension_kn:8.2f} {rr.peak_accel_g:7.1f} "
          f"{rr.max_drop_m:7.2f} {rr.max_elongation_pct:7.1f} {ok:>6s}")
print(f"  -> 峰值对 FF 单调递增: {'是' if mono else '否'}")

print()
print("=" * 88)
print("3) 保护器对比（80 kg, FF=1.77, L0=2 m）—— 标定后应能区分开")
print("=" * 88)
print(f"{'保护器':<18s} {'保持力kN':>9s} {'峰值kN':>8s} {'峰值g':>7s} {'放绳m':>7s} "
      f"{'滑动态s':>8s} {'绳总滑m':>8s}")
for name, d in DEVICES.items():
    rr = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.77,
                            device=name, t_hold_n=d.t_hold_n)
    print(f"{d.name:<18s} {d.t_hold_n/1e3:9.1f} {rr.peak_tension_kn:8.2f} {rr.peak_accel_g:7.1f} "
          f"{rr.max_payout_m:7.2f} {rr.device_sat_time_s:8.4f} {rr.total_rope_slide_m:8.2f}")

print()
print("=" * 88)
print("4) 保护者侧张力（欧拉绞盘）—— 多挂点为何最后一挂受力最大")
print("=" * 88)
print(f"{'系统':<22s} {'包角rad':>8s} {'坠落者kN':>9s} {'保护者kN':>9s} {'放大':>6s}")
for label, th in [("顶绳/单段", 0.30), ("运动式两挂", 1.10), ("传统多段", 2.20)]:
    f = belayer_side_force(8500.0, mu=0.20, theta_rad=th)
    print(f"{label:<22s} {th:8.2f} {8.5:9.2f} {f/1e3:9.2f} {f/8500:6.2f}x")

print()
print("=" * 88)
print("5) 冲击质量敏感性（FF=1.77, L0=2 m, 锁死）")
print("=" * 88)
print(f"{'质量kg':>7s} {'峰值kN':>8s} {'峰值g':>7s} {'伸长%':>7s} {'HIC':>7s}")
for mk in [40, 60, 80, 100, 120]:
    rr = simulate_rope_fall(mass_kg=mk, rope_length_m=2.0, fall_factor=1.77, t_hold_n=1e9)
    rr.hic = hic_from_result(rr)
    print(f"{mk:7d} {rr.peak_tension_kn:8.2f} {rr.peak_accel_g:7.1f} "
          f"{rr.max_elongation_pct:7.1f} {rr.hic:7.0f}")
