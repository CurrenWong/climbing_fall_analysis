import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from climbing.chain2d_tissue import (
    POSTURES, TissueChainConfig, TissueParams, simulate_tissue_chain_fall,
    tissue_dissipated_energy_t_j, tissue_pad_work_t_j, tissue_dissipated_energy_j,
    _build_tissue_ctx,
)
P = TissueParams()
cfg = TissueChainConfig(mass_kg=80.0, zeta_joint=0.0, zeta_contact=0.0,
                        c_axial_n_s_per_m=0.0, mu_foot=0.0, rtol=1e-8, atol=1e-11)
r = simulate_tissue_chain_fall(0.5, 80.0, "feet-first-stiff", False, config=cfg)
w = tissue_pad_work_t_j(r)
d = tissue_dissipated_energy_t_j(r)
dd = np.diff(d)
print("min(diff(d)) =", dd.min(), " n_neg =", int((dd < 0).sum()), "/", len(dd))
print("w[0],w[-1] =", w[0], w[-1])
print("d[0],d[-1] =", d[0], d[-1])
print("scalar tissue_dissipated_energy_j =", tissue_dissipated_energy_j(r))
i = int(np.argmin(dd))
print("worst idx", i, "t=", r.t_s[i]*1e3, "dd=", dd[i],
      " 附近 w:", w[i-1:i+3], " d:", d[i-1:i+3])
