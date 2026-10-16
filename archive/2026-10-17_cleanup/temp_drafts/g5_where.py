"""定位能量漂移发生在哪个相位（L4(b) 诊断）。

报告 e(t) 的 min/max 出现在何时、以及各分量在同一时刻的占比。
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
T, E1 = _P.thickness_m, _P.E1_pa
CFG = TissueChainConfig(
    mass_kg=80.0, zeta_joint=0.0, zeta_contact=0.0,
    c_axial_n_s_per_m=0.0, mu_foot=0.0, rtol=1e-8, atol=1e-11,
)


def parts(h, cfg=CFG):
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
    xi = np.maximum(r.node_z, 0.0) / T
    u = np.maximum(xi - r.meta["xi_r"], 0.0)
    pe_c = 0.5 * E1 * np.sum(u**2 * r.contact_area_m2 * T, axis=1)
    diss = tissue_dissipated_energy_j(r)
    a = result_ = r
    return r, dict(ke=ke, pe_g=pe_g, pe_j=pe_j, pe_c=pe_c, diss=diss,
                   e=ke + pe_g + pe_j + pe_c + diss)


if __name__ == "__main__":
    for h, cfg in ((0.5, CFG), (1.5, CFG), (2.0, CFG)):
        r, p = parts(h, cfg)
        e = p["e"]
        i0 = 0
        i_min, i_max = int(np.argmin(e)), int(np.argmax(e))
        print(f"\n=== h={h} m   n_t={len(e)}  t_end={r.t_s[-1]*1e3:.1f} ms ===")
        print(f"  e[0]      = {e[0]:10.4f}   (t={r.t_s[0]*1e3:.3f} ms)")
        print(f"  min e     = {e.min():10.4f} @ idx {i_min} (t={r.t_s[i_min]*1e3:.3f} ms)"
              f"  dev {(e.min()-e[0])/abs(e[0])*100:+.3f}%")
        print(f"  max e     = {e.max():10.4f} @ idx {i_max} (t={r.t_s[i_max]*1e3:.3f} ms)"
              f"  dev {(e.max()-e[0])/abs(e[0])*100:+.3f}%")
        print(f"  e[end]    = {e[-1]:10.4f}  dev {(e[-1]-e[0])/abs(e[0])*100:+.3f}%")
        # 相位：以最大压缩时刻为界
        ic = int(np.argmax(np.max(np.maximum(r.node_z, 0.0), axis=1)))
        print(f"  最大压缩 idx {ic} (t={r.t_s[ic]*1e3:.3f} ms)")
        print(f"  早于压缩相位（idx<{ic}）: dev ∈ "
              f"[{(e[:ic].min()-e[0])/abs(e[0])*100:+.3f}%, "
              f"{(e[:ic].max()-e[0])/abs(e[0])*100:+.3f}%]")
        print(f"  晚于压缩相位（idx>{ic}）: dev ∈ "
              f"[{(e[ic:].min()-e[0])/abs(e[0])*100:+.3f}%, "
              f"{(e[ic:].max()-e[0])/abs(e[0])*100:+.3f}%]")
        print(f"  分量 @e_min : ke={p['ke'][i_min]:9.4f} pe_g={p['pe_g'][i_min]:9.4f} "
              f"pe_j={p['pe_j'][i_min]:8.4f} pe_c={p['pe_c'][i_min]:8.4f} "
              f"diss={p['diss'][i_min]:7.4f}")
        print(f"  分量 @e_max : ke={p['ke'][i_max]:9.4f} pe_g={p['pe_g'][i_max]:9.4f} "
              f"pe_j={p['pe_j'][i_max]:8.4f} pe_c={p['pe_c'][i_max]:8.4f} "
              f"diss={p['diss'][i_max]:7.4f}")
