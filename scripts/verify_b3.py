"""B-3 链模型验证：是否真的修掉了 B-2 的两个结构性问题？

判据（写死，不允许事后调参凑）：
  A. 屈膝缓冲的**总腿屈曲行程**应落在 0.20-0.50 m（B-2 只有 0.087 m）
  B. 躯干峰值加速度应显著低于 B-2 的 37 g（3 m）
  C. 损伤判读在低落差下应能出现"低"档（B-2 全部顶"极高"）
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import time

from climbing.bodies import simulate_chain_fall, CHAIN_POSTURES, default_chain
from climbing.pad import simulate_boulder_fall, POSTURES
from climbing.metrics import hic_from_result, hic, classify_hic
from climbing.injury import assess_boulder_fall

H, M = 3.0, 80.0

print("=" * 92)
print("B-3 链模型 vs B-2 双质点模型   3.0 m / 80 kg / 20 cm 软垫")
print("=" * 92)
print(f"{'姿势':<22s} {'B2力kN':>7s} {'B3力kN':>7s} {'B2压力':>7s} {'B3压力':>7s} "
      f"{'B2接触g':>8s} {'B3接触g':>8s} {'B2屈曲cm':>9s} {'B3屈曲cm':>9s}")
print("-" * 92)
rows = []
for name in CHAIN_POSTURES:
    t0 = time.time()
    b3 = simulate_chain_fall(height_m=H, mass_kg=M, posture=name)
    b2 = simulate_boulder_fall(height_m=H, mass_kg=M, posture=name)
    p2 = float(np.max(b2.pad_force_n / np.maximum(b2.contact_area_m2, 1e-9))) / 1e3
    rows.append((name, b2, b3, p2))
    print(f"{CHAIN_POSTURES[name].name:<22s} {b2.peak_force_kn:7.1f} {b3.peak_force_kn:7.1f} "
          f"{p2:7.0f} {b3.peak_pressure_kpa:7.0f} "
          f"{b2.peak_primary_g:8.1f} {b3.peak_accel_foot_g:8.1f} "
          f"{100*np.max(np.maximum(0,0.85-(b2.z_foot_m-b2.z_torso_m))):9.1f} "
          f"{100*b3.total_leg_flex_m:9.1f}   ({time.time()-t0:.1f}s)")

print()
print("=" * 92)
print("判据 A: 屈膝缓冲的总腿屈曲行程应在 0.20-0.50 m")
print("=" * 92)
b3 = simulate_chain_fall(height_m=H, mass_kg=M, posture="controlled-drop")
b2 = simulate_boulder_fall(height_m=H, mass_kg=M, posture="controlled-drop")
x2 = 0.85 - (b2.z_foot_m - b2.z_torso_m)
print(f"  B-2 双质点: {100*np.max(x2):.1f} cm   <- 躯干无第二条通路，行程被截断")
print(f"  B-3 链模型: {100*b3.total_leg_flex_m:.1f} cm  "
      f"(踝{100*b3.total_leg_flex_m*0+np.max(b3.ankle_x):.0f} "
      f"膝{100*np.max(b3.knee_x):.0f} 髋{100*np.max(b3.hip_x):.0f})")
ok_a = 0.20 <= b3.total_leg_flex_m <= 0.50
print(f"  -> {'通过' if ok_a else '未通过'}（判据 0.20-0.50 m）")

print()
print("=" * 92)
print("判据 B: 躯干峰值加速度应低于 B-2")
print("=" * 92)
print(f"  B-2 躯干 {b2.peak_accel_torso_g:.1f} g  ->  B-3 躯干 {b3.peak_accel_trunk_g:.1f} g"
      f"  (降低 {100*(1-b3.peak_accel_trunk_g/b2.peak_accel_torso_g):.0f}%)")
print(f"  B-2 头   {b2.peak_primary_g:.1f} g  ->  B-3 头   {b3.peak_accel_head_g:.1f} g")
ok_b = b3.peak_accel_trunk_g < 0.75 * b2.peak_accel_torso_g
print(f"  -> {'通过' if ok_b else '未通过'}（判据：至少降低 25%）")

print()
print("=" * 92)
print("判据 C: 低落差下应能出现'低'档（B-2 全部顶'极高'）")
print("=" * 92)
print(f"{'落差m':>6s} {'B3躯干g':>8s} {'B3脊柱kN':>9s} {'B3压力kPa':>10s} {'B3 HIC':>7s} {'B2判读':>7s}")
for h in [0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
    r3 = simulate_chain_fall(height_m=h, mass_kg=M, posture="controlled-drop")
    a_head = r3.primary_accel_g()
    h3 = hic(r3.t_s, a_head)
    r2 = simulate_boulder_fall(height_m=h, mass_kg=M, posture="controlled-drop")
    r2.hic = hic_from_result(r2)
    v2 = assess_boulder_fall(r2, "controlled-drop").verdict
    band = classify_hic(h3)
    print(f"{h:6.1f} {r3.peak_accel_trunk_g:8.1f} {r3.max_spine_load_kn:9.1f} "
          f"{r3.peak_pressure_kpa:10.0f} {h3:7.0f} {v2:>7s}   <- B3 HIC 分档: {band}")

print()
print("=" * 92)
print("关节峰值载荷 (kN) —— 膝/髋是否有合理的生理上限")
print("=" * 92)
print(f"{'姿势':<22s} {'脊柱kN':>8s} {'膝kN':>7s} {'髋kN':>7s} {'首触节点':>9s}")
for name, b2, b3, p2 in rows:
    print(f"{CHAIN_POSTURES[name].name:<22s} {b3.max_spine_load_kn:8.1f} "
          f"{b3.max_knee_load_kn:7.1f} {b3.max_hip_load_kn:7.1f} {b3.first_contact_node:9d}")
print("  (活体膝关节轴向失效量级约 8-16 kN；超过说明肌肉/行程参数需重标)")
