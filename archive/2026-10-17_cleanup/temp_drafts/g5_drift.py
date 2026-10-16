"""量致密化段能量账目漂移（h 扫描）—— L4(b) 的改前/改后对照。

保守配置（无关节/接触阻尼、无摩擦、无轴向阻尼），机械能闭合漂移。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from climbing.chain2d_tissue import (
    POSTURES, TissueChainConfig, TissueParams, simulate_tissue_chain_fall,
    tissue_dissipated_energy_j, _build_tissue_ctx,
)

_P = TissueParams()
T = _P.thickness_m
E1 = _P.E1_pa

CFG = TissueChainConfig(
    mass_kg=80.0, zeta_joint=0.0, zeta_contact=0.0,
    c_axial_n_s_per_m=0.0, mu_foot=0.0, rtol=1e-8, atol=1e-11,
)


def drift(h, tol="tight"):
    cfg = CFG
    if tol != "tight":
        cfg = TissueChainConfig(
            mass_kg=80.0, zeta_joint=0.0, zeta_contact=0.0,
            c_axial_n_s_per_m=0.0, mu_foot=0.0, rtol=1e-4, atol=1e-6,
        )
    r = simulate_tissue_chain_fall(h, 80.0, "feet-first-stiff", False, config=cfg)
    ctx = _build_tissue_ctx(POSTURES["feet-first-stiff"], False, cfg)
    m = r.meta["node_mass_kg"]
    ke = 0.5 * np.sum(m[None, :] * (r.node_vx**2 + r.node_vz**2), axis=1)
    pe_g = -np.sum(m[None, :] * cfg.g * r.node_z, axis=1)
    dx = r.node_x[:, 1:] - r.node_x[:, :-1]
    dz = r.node_z[:, 1:] - r.node_z[:, :-1]
    ell = np.maximum(np.hypot(dx, dz), 1e-9)
    comp = ctx.rest[None, :] - ell
    ext = ell - ctx.rest[None, :]
    pe_j = 0.5 * np.sum(ctx.k_ax[None, :] * (ell - ctx.rest[None, :]) ** 2, axis=1)
    pe_j += np.sum(0.5 * ctx.k_comp_stop[None, :]
                   * np.maximum(comp - ctx.comp_limit[None, :], 0.0) ** 2, axis=1)
    pe_j += np.sum(0.5 * ctx.k_ext_stop[None, :]
                   * np.maximum(ext - ctx.ext_limit[None, :], 0.0) ** 2, axis=1)
    xi = np.maximum(r.node_z, 0.0) / T
    u = np.maximum(xi - r.meta["xi_r"], 0.0)
    pe_c = 0.5 * E1 * np.sum(u**2 * r.contact_area_m2 * T, axis=1)
    diss = tissue_dissipated_energy_j(r)
    e = ke + pe_g + pe_j + pe_c + diss
    d = max(abs(e.max() - e[0]), abs(e.min() - e[0])) / abs(e[0])
    return d, float(diss), float(np.max(np.maximum(r.node_z, 0.0)) * 1e3)


if __name__ == "__main__":
    print("  h(m)  tol       drift%    diss_J   peak_mm")
    for h in (0.5, 1.0, 1.5, 2.0, 3.0):
        for tol in ("tight", "loose"):
            try:
                d, dj, pk = drift(h, tol)
            except RuntimeError:
                print(f"  {h:4.1f}  {tol:5s}    --- 不收敛 ---")
                continue
            print(f"  {h:4.1f}  {tol:5s}  {d*100:7.3f}  {dj:8.2f}   {pk:6.2f}")
            break
