"""Phase B-2 分析：找出"落在合规软垫上仍然受伤"的姿势与原因。"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from climbing.pad import simulate_boulder_fall, POSTURES, CrashPad, hard_surface
from climbing.metrics import hic_from_result
from climbing.injury import assess_boulder_fall, PRESSURE_BANDS

H, M = 3.0, 80.0

def rep(name, height=H, mass=M, pad=None, on_pad=True):
    b = simulate_boulder_fall(height_m=height, mass_kg=mass, posture=name,
                              pad=pad, on_pad=on_pad)
    b.hic = hic_from_result(b)
    return b, assess_boulder_fall(b, name)

print("#" * 92)
print(f"# Phase B-2：合规软垫（20 cm 标准 PU）上仍然受伤的姿势   落差 {H} m / {M:.0f} kg")
print("#" * 92)
print()
print(f"{'姿势':<22s} {'峰值力kN':>9s} {'面积cm2':>8s} {'压力kPa':>8s} {'接触g':>7s} "
      f"{'HIC':>6s} {'脊柱kN':>7s} {'判读':>5s}")
print("-" * 92)
data = []
for name, p in POSTURES.items():
    b, r = rep(name)
    p_at = r.parts[0]
    spine = [x for x in r.parts if x.name == "脊柱"][0]
    data.append((name, p, b, r))
    print(f"{p.name:<22s} {b.peak_force_kn:9.1f} {p_at.detail['area_cm2_at_ppeak']:8.0f} "
          f"{r.peak_pressure_kpa:8.0f} {b.peak_primary_g:7.1f} {r.hic:6.0f} "
          f"{spine.metric:7.1f} {r.verdict:>5s}")

print()
print("=" * 92)
print("仍然受伤的姿势（判读 != 低），按严重度排序")
print("=" * 92)
rank = {"低": 0, "中": 1, "高": 2, "极高": 3}
hurt = [d for d in data if rank[d[3].verdict] >= 1]
for name, p, b, r in sorted(hurt, key=lambda d: -rank[d[3].verdict]):
    print(f"\n【{p.name}】 判读={r.verdict}   压力={r.peak_pressure_kpa:.0f} kPa  "
          f"力={r.peak_force_kn:.1f} kN  HIC={r.hic:.0f}  接触g={b.peak_primary_g:.1f}")
    for part in r.parts:
        if rank.get(part.band, 0) >= 1:
            print(f"   - {part.name}: {part.metric_name}={part.metric:.1f}{part.unit} [{part.band}]")
            print(f"     原因: {part.mechanism}")

print()
print("=" * 92)
print("判读随高度变化（20 cm 标准垫）—— 看分档有没有区分度")
print("=" * 92)
heights = [0.3, 0.5, 0.8, 1.2, 1.8, 2.5, 3.0, 4.0, 5.0]
SYM = {"低": ".", "中": "中", "高": "高", "极高": "极"}
print(f"{'姿势':<22s}" + "".join(f"{h:>6.1f}m" for h in heights))
print("-" * 92)
for name, p in POSTURES.items():
    cells = []
    for h in heights:
        b = simulate_boulder_fall(height_m=h, mass_kg=M, posture=name)
        b.hic = hic_from_result(b)
        cells.append(SYM[assess_boulder_fall(b, name).verdict])
    print(f"{p.name:<22s}" + "".join(f"{c:>6s}" for c in cells))

print()
print("=" * 92)
print("判读随质量变化（3.0 m，20 cm 标准垫）")
print("=" * 92)
masses = [50, 60, 70, 80, 90, 100, 120]
print(f"{'姿势':<22s}" + "".join(f"{m:>6d}kg" for m in masses))
print("-" * 92)
for name, p in POSTURES.items():
    cells = []
    for mk in masses:
        b = simulate_boulder_fall(height_m=H, mass_kg=mk, posture=name)
        b.hic = hic_from_result(b)
        cells.append(SYM[assess_boulder_fall(b, name).verdict])
    print(f"{p.name:<22s}" + "".join(f"{c:>6s}" for c in cells))

print()
print("=" * 92)
print("每个姿势多高开始受伤？（20 cm 标准垫，判读 >= 中 即视为受伤）")
print("=" * 92)
for name, p in POSTURES.items():
    hit = None
    for h in [x / 20 for x in range(4, 121)]:
        b = simulate_boulder_fall(height_m=h, mass_kg=M, posture=name)
        b.hic = hic_from_result(b)
        r = assess_boulder_fall(b, name)
        if rank[r.verdict] >= 1:
            hit = (h, r.verdict, r)
            break
    if hit:
        h, v, r = hit
        print(f"  {p.name:<22s} {h:4.2f} m  (判读={v}, 压力={r.peak_pressure_kpa:.0f} kPa, "
              f"HIC={r.hic:.0f})")
    else:
        print(f"  {p.name:<22s}  >6.0 m 未见受伤")

print()
print("=" * 92)
print("软垫有效性：同样姿势，峰值接触压力（软垫 vs 硬地面）")
print("=" * 92)
print(f"{'姿势':<22s} {'垫kPa':>7s} {'硬地kPa':>9s} {'降幅':>7s} {'垫判读':>7s} {'硬地判读':>9s}")
for name, p in POSTURES.items():
    _, ron = rep(name, on_pad=True)
    _, roff = rep(name, on_pad=False)
    red = 100 * (1 - ron.peak_pressure_kpa / roff.peak_pressure_kpa)
    print(f"{p.name:<22s} {ron.peak_pressure_kpa:7.0f} {roff.peak_pressure_kpa:9.0f} "
          f"{red:6.0f}% {ron.verdict:>7s} {roff.verdict:>9s}")

print()
print("=" * 92)
print("软垫厚度的影响（3.0 m，姿势=前脚掌点地 —— 压力最敏感的姿势）")
print("=" * 92)
print(f"{'垫厚cm':>7s} {'F kN':>7s} {'压力kPa':>9s} {'压缩cm':>7s} {'HIC':>6s} {'判读':>5s}")
for th in [0.05, 0.10, 0.15, 0.20, 0.30, 0.40]:
    b, r = rep("toe-point", pad=CrashPad(thickness_m=th))
    print(f"{100*th:7.0f} {b.peak_force_kn:7.1f} {r.peak_pressure_kpa:9.0f} "
          f"{100*b.max_compression_m:7.1f} {r.hic:6.0f} {r.verdict:>5s}")
