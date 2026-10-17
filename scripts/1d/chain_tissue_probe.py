"""G5 · v2 柔顺层探针（判据跑前冻结，逐条如实判）。

对照：chain2d.TissueLayer (E=2 MPa) vs chain2d_tissue (E=185 kPa in-vivo)。
口径：{2,3,5,10,20} m 刚地面，80 kg，feet-first-stiff。

判据（G5 预注册 §4）
  1 [主] 20m 足节点峰值变形 ∈ Bennett&Ker in-vivo **8.5–11.3 mm**  ← 唯一独立文献靶
  2 [暴露] σ_胫骨峰值 vs 皮质骨抗压强度 130–180 MPa；> 180 须标注超骨强度
  3 自洽  E_现/E_v2 = 10.8×（由 §1 直接算出，写死）
  4 [仅记录] 吸能占比（能量闸门已声明坠落能量下不可作判据）
  5 [仅记录] V8b 峰时（已声明不可能改善）
  6 [方向] 20m ξ_r 终值 > 0（永久压实生效 ⇒ 消除整身回弹）
"""
import sys
import numpy as np
sys.path.insert(0, "src")
from climbing.chain2d import PlanarChainConfig, simulate_planar_chain_fall
from climbing.chain2d_crush import CrushChainConfig, simulate_crush_chain_fall
from climbing.chain2d_tissue import (
    TissueChainConfig, simulate_tissue_chain_fall,
    tissue_dissipated_energy_j, tissue_peak_deformation_mm,
    tibia_peak_stress_mpa, tibia_local_peak,
)

HEIGHTS = (2.0, 3.0, 5.0, 10.0, 20.0)
BK_LO, BK_HI = 8.5, 11.3      # Bennett & Ker in-vivo 变形量 (mm)
BONE_LO, BONE_HI = 130.0, 180.0   # 皮质骨抗压强度 (MPa)
MASS, G = 80.0, 9.81


def head_rise(r, horizon=0.25):
    t0, _ = r.meta["impact_window_s"]
    m = (r.t_s >= t0) & (r.t_s <= t0 + horizon)
    zh = r.node_z[m, -1]
    return float(max(0.0, zh[0] - zh.min())) * 1e3


print("=" * 78)
print("G5 · v2 柔顺层探针  (E1=185 kPa in-vivo, t=15.4 mm, xi_ref=0.70, Ed=20 MPa)")
print("=" * 78)
print(f"{'h(m)':>5} | {'模型':>9} | {'变形mm':>8} {'σ_胫骨MPa':>10} {'峰力kN':>8} "
      f"{'吸能J':>8} {'吸能%':>7} {'V8b_ms':>7} {'headrise':>9} {'ξr终':>7}")
print("-" * 78)

rows = {}
for h in HEIGHTS:
    tag = f"{h:.0f}m"
    r_g1 = simulate_planar_chain_fall(h, MASS, "feet-first-stiff", False,
                                      PlanarChainConfig(mass_kg=MASS))
    r_g4 = simulate_crush_chain_fall(h, MASS, "feet-first-stiff", False,
                                     CrushChainConfig(mass_kg=MASS))
    r_v2 = simulate_tissue_chain_fall(h, MASS, "feet-first-stiff", False,
                                      TissueChainConfig(mass_kg=MASS))
    rows[h] = (r_g1, r_g4, r_v2)
    for lbl, r, diss in (("G1 2MPa", r_g1, 0.0),
                         ("G4 crush", r_g4, None),
                         ("v2 in-vivo", r_v2, tissue_dissipated_energy_j(r_v2))):
        defo = tissue_peak_deformation_mm(r)
        sg, _ = tibia_peak_stress_mpa(r)
        tp, _ = tibia_local_peak(r)
        d = tissue_dissipated_energy_j(r) if diss is None else diss
        print(f"{tag:>5} | {lbl:>9} | {defo:8.2f} {sg:10.1f} "
              f"{r.peak_force_n/1e3:8.2f} {d:8.1f} {100*d/r.impact_energy_j:6.2f}% "
              f"{tp*1e3:7.2f} {head_rise(r):9.1f} "
              f"{r.meta['xi_r'][-1,0] if r.meta.get('xi_r') is not None else 0.0:7.4f}")
    print("-" * 78)

# ---- 判据逐条判定（只用 20 m） ----
r_v2 = rows[20.0][2]
defo = tissue_peak_deformation_mm(r_v2)
sg, tms = tibia_peak_stress_mpa(r_v2)
er = r_v2.meta["xi_r"][-1, 0]
print("\n===== 判据判定 @20 m =====")
c1 = BK_LO <= defo <= BK_HI
print(f"1 [主] 峰值变形 {defo:.2f} mm ∈ [{BK_LO}, {BK_HI}] (Bennett&Ker in-vivo) : "
      f"{'PASS' if c1 else 'FAIL'}")
print(f"2 [暴露] σ_胫骨 {sg:.1f} MPa vs 骨强度 130–180 MPa : "
      f"{'OK' if sg <= BONE_HI else '超骨强度 %.1f× (标注)' % (sg/BONE_HI)}"
      f"  [峰时 {tms:.2f} ms]")
print(f"3 [自洽] E_现/E_v2 = 2e6/0.185e6 = {2e6/0.185e6:.2f}× （写死 10.8×）")
print(f"6 [方向] ξ_r 终值 = {er:.5f} > 0（永久压实生效）：{'PASS' if er > 0 else 'FAIL'}")
