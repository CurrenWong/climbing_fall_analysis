"""诊断：把耗散按**时间序列**累计，看能量是否闭合。

假设：`tissue_dissipated_energy_j` 返回标量 ⇒ 加进 e(t) 是常数 ⇒ 对 drift 零贡献。
真正应闭合的是  KE + PE_g + PE_j + W_into_pad。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from climbing.chain2d_tissue import (
    POSTURES, TissueChainConfig, TissueParams, simulate_tissue_chain_fall,
    _build_tissue_ctx,
)

_P = TissueParams()
T, E1 = _P.thickness_m, _P.E1_pa
CFG = TissueChainConfig(
    mass_kg=80.0, zeta_joint=0.0, zeta_contact=0.0,
    c_axial_n_s_per_m=0.0, mu_foot=0.0, rtol=1e-8, atol=1e-11,
)


def go(h, cfg=CFG):
    r = simulate_tissue_chain_fall(h, 80.0, "feet-first-stiff", False, config=cfg)
    ctx = _build_tissue_ctx(POSTURES["feet-first-stiff"], False, cfg)
    m = r.meta["node_mass_kg"]
    ke = 0.5 * np.sum(m[None, :] * (r.node_vx**2 + r.node_vz**2), axis=1)
    pe_g = -np.sum(m[None, :] * cfg.g * r.node_z, axis=1)
    dx = r.node_x[:, 1:] - r.node_x[:, :-1]
    dz = r.node_z[:, 1:] - r.node_z[:, :-1]
    ell = np.maximum(np.hypot(dx, dz), 1e-9)
    comp, ext = ctx.rest[None, :] - ell, ell - ctx.rest[None, :]
    pe_j = 0.5 * np.sum(ctx.k_ax[None, :] * (ell - ctx.rest[None, :]) ** 2, axis=1)
    pe_j += np.sum(0.5 * ctx.k_comp_stop[None, :]
                   * np.maximum(comp - ctx.comp_limit[None, :], 0.0) ** 2, axis=1)
    pe_j += np.sum(0.5 * ctx.k_ext_stop[None, :]
                   * np.maximum(ext - ctx.ext_limit[None, :], 0.0) ** 2, axis=1)

    # ---- W_into_pad 的时间序列（有符号，梯形） ----
    xi = np.maximum(r.node_z, 0.0) / T
    d_xi_t = np.diff(xi, axis=0) * T              # m  (n_t-1, n_node)
    F = r.contact_force_n                          # (n_t, n_node)
    w_inc = 0.5 * (F[:-1] + F[1:]) * d_xi_t        # J  per step per node
    w_inc = w_inc.sum(axis=1)
    w_into = np.concatenate([[0.0], np.cumsum(w_inc)])

    u = np.maximum(xi - r.meta["xi_r"], 0.0)
    pe_c = 0.5 * E1 * np.sum(u**2 * r.contact_area_m2 * T, axis=1)

    def dev(e):
        return (e - e[0]) / abs(e[0]) * 100.0

    print(f"\n=== h={h} m   n_t={len(ke)} ===")
    for name, e in (
        ("KE+PEg+PEj+PEc      (现测试量)", ke + pe_g + pe_j + pe_c),
        ("KE+PEg+PEj+W_into   (真闭合量)", ke + pe_g + pe_j + w_into),
    ):
        d = dev(e)
        print(f"  {name}: dev∈[{d.min():+.4f}%, {d.max():+.4f}%]  "
              f"max|dev|={np.abs(d).max():.4f}%")
    print(f"  W_into[end]={w_into[-1]:9.4f} J   "
          f"PEc[end]={pe_c[-1]:9.4f} J   "
          f"diss(标量,新)={float(__import__('climbing.chain2d_tissue',fromlist=['x']).tissue_dissipated_energy_j(r)):9.4f} J")
    return r


if __name__ == "__main__":
    for h, cfg in ((0.5, CFG), (1.5, CFG), (2.0, CFG)):
        go(h, cfg)
