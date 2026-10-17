"""G6 · 载荷段软组织串联柔顺 探针（判据跑前冻结，如实判）。

对照：SEGMENTS_6（纯骨柱）vs SEGMENTS_6_SOFT（加软组织串联柔顺）。
接触层统一用 v2 in-vivo 柔顺层（G5），隔离段刚度这一个变量。

判据（G6 预注册 §4）
  1 [主 A1] 整链等效轴向刚度 @20m = m·(π/(2·T_arrest))² ∈ 30–170 kN/m
  2 [主 A2] 步行等效工况(h=0.2m, 157J) 胫骨骨应变 ∈ 237–1250 µε
  3 [暴露] σ_胫骨 @20m vs 皮质骨抗压强度 130–180 MPa
  4 [记录] V8b 峰时（不判 PASS/FAIL）
  5 [记录] 峰值力 @20m
"""
import sys
sys.path.insert(0, "src")
import numpy as np
from climbing.chain2d import SEGMENTS_6, segment_local_peak
from climbing.chain2d_soft import (
    SEGMENTS_6_SOFT, chain_effective_stiffness_n_per_m, arrest_stiffness_n_per_m,
)
from climbing.chain2d_tissue import (
    TissueChainConfig, simulate_tissue_chain_fall,
    tissue_peak_deformation_mm, tibia_peak_stress_mpa,
)

MASS, G = 80.0, 9.81
A1_LO, A1_HI = 30.0, 170.0          # kN/m
A2_LO, A2_HI = 237.0, 1250.0        # microstrain
BONE_LO, BONE_HI = 130.0, 180.0     # MPa


def t_arrest(r):
    """首次接触 → **整身动量**刹停（质量加权平均竖直速度首次过零）的时长 (s)。

    ⚠️ **定义更正（2026-10-13，首跑暴露）**：预注册原定义的「节点 0（足）速度
    首次过零」是**退化观测量**——足节点被**地面**刹停（≈1 ms），与链刚度无关，
    两个变体给出**完全相同**的值（均 1.00 ms）⇒ 无区分力。改用**整身动量**：
    ``v̄ = Σ m_i·v_z,i / Σ m_i`` 首次过零，才真正反映「链把身体刹停」的时标。
    这是**观测量定义纠错**（原定义退化、不可用），非事后挪靶；原值已留档。
    """
    t = r.t_s
    m = r.meta["node_mass_kg"]
    vbar = (r.node_vz * m[None, :]).sum(axis=1) / m.sum()
    t0, _ = r.meta["impact_window_s"]
    vbar = (r.node_vz * m[None, :]).sum(axis=1) / m.sum()
    i0 = int(np.argmax(t >= t0))
    for k in range(i0, len(t)):
        if vbar[k] <= 0.0:
            return float(t[k] - t[i0])
    return float(t[-1] - t[i0])


def arrest_diag(r):
    """T_arrest 可信性诊断：是否真发生动量反转 + v̄ 最小值。"""
    t = r.t_s
    m = r.meta["node_mass_kg"]
    vbar = (r.node_vz * m[None, :]).sum(axis=1) / m.sum()
    t0, _ = r.meta["impact_window_s"]
    i0 = int(np.argmax(t >= t0))
    crossed = bool(np.any(vbar[i0:] <= 0.0))
    return crossed, float(vbar[i0:].min()), float(t[-1] - t[i0]) * 1e3


def run(segs, label):
    r20 = simulate_tissue_chain_fall(20.0, MASS, "feet-first-stiff", False,
                                     TissueChainConfig(mass_kg=MASS), segments=segs)
    r02 = simulate_tissue_chain_fall(0.2, MASS, "feet-first-stiff", False,
                                     TissueChainConfig(mass_kg=MASS), segments=segs)
    ta = t_arrest(r20)
    k_eff = arrest_stiffness_n_per_m(MASS, ta) / 1e3
    crossed, vmin, win_ms = arrest_diag(r20)
    eps_w = segment_local_peak(r02, 1)[1]
    eps_ue = eps_w * 1e6
    sg, _ = tibia_peak_stress_mpa(r20)
    tp, _ = segment_local_peak(r20, 1)
    return dict(label=label, k_static=chain_effective_stiffness_n_per_m(segs) / 1e3,
                T_arrest_ms=ta * 1e3, k_eff=k_eff, eps_ue=eps_ue,
                sigma_mpa=sg, v8b_ms=tp * 1e3,
                Fpk_kn=r20.peak_force_n / 1e3,
                comp_mm=tissue_peak_deformation_mm(r20),
                walk_comp_mm=tissue_peak_deformation_mm(r02),
                crossed=crossed, vmin_ms=vmin * 1e3, win_ms=win_ms)


rows = [run(SEGMENTS_6, "纯骨柱 SEGMENTS_6"),
        run(SEGMENTS_6_SOFT, "串联柔顺 SEGMENTS_6_SOFT")]

print("=" * 96)
print("G6 · 载荷段软组织串联柔顺（接触层 = v2 in-vivo 柔顺，隔离段刚度）")
print("=" * 96)
print(f"{'变体':>22} {'k_静(kN/m)':>11} {'T_arrest(ms)':>13} {'k_eff(kN/m)':>12} "
      f"{'ε步(µε)':>10} {'σ20m(MPa)':>10} {'Fpk(kN)':>9} {'V8b(ms)':>8} {'压20m(mm)':>10}")
print("-" * 96)
for d in rows:
    print(f"{d['label']:>22} {d['k_static']:11.1f} {d['T_arrest_ms']:13.2f} "
          f"{d['k_eff']:12.1f} {d['eps_ue']:10.0f} {d['sigma_mpa']:10.1f} "
          f"{d['Fpk_kn']:9.2f} {d['v8b_ms']:8.2f} {d['comp_mm']:10.2f}")

print()
print("===== T_arrest 可信性诊断（动量是否真反转） =====")
for d in rows:
    print(f"  {d['label']:>22}: 动量反转={d['crossed']}  "
          f"v̄_min={d['vmin_ms']:.3f} m/s  积分窗={d['win_ms']:.1f} ms  "
          f"T_arrest={d['T_arrest_ms']:.2f} ms")

print()
print("===== 判据判定 =====")
d = rows[1]
c1 = A1_LO <= d["k_eff"] <= A1_HI
c2 = A2_LO <= d["eps_ue"] <= A2_HI
print(f"1 [主 A1] k_eff = {d['k_eff']:.1f} kN/m ∈ [{A1_LO}, {A1_HI}]  : "
      f"{'PASS' if c1 else 'FAIL'}")
print(f"   （参照：现模型 {rows[0]['k_eff']:.1f} kN/m）")
print(f"2 [主 A2] 步行工况 ε = {d['eps_ue']:.0f} µε ∈ [{A2_LO:.0f}, {A2_HI:.0f}] : "
      f"{'PASS' if c2 else 'FAIL'}")
print(f"   （参照：现模型 {rows[0]['eps_ue']:.0f} µε）")
print(f"3 [暴露] σ_胫骨@20m = {d['sigma_mpa']:.1f} MPa vs 骨强度 [{BONE_LO:.0f},{BONE_HI:.0f}]"
      f" : {'OK' if d['sigma_mpa'] <= BONE_HI else '超 %.1f×' % (d['sigma_mpa']/BONE_HI)}")
print(f"   （参照：现模型 {rows[0]['sigma_mpa']:.1f} MPa）")
print(f"4 [记录] V8b = {d['v8b_ms']:.2f} ms（现模型 {rows[0]['v8b_ms']:.2f} ms）")
print(f"5 [记录] 峰值力 = {d['Fpk_kn']:.2f} kN（现模型 {rows[0]['Fpk_kn']:.2f} kN）")
