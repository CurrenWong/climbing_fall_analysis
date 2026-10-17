"""P1 v2：裕度表（修域后缀映射 + 外推断裂载荷 + 奇点剔除）。"""
import io
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
PY = R / ".venv/Scripts/python.exe"
W = R / "temp/thums/febio_lig"
sys.path.insert(0, str(R / "scripts/opensim_fe"))

# ── 1) 用 acceptance 跑 +1.6（它会打印合力 + NORMAL）──
r = subprocess.run([str(PY), "-u", str(R / "scripts/ankle_fe/acceptance.py"),
                    "DRIVER=dispx", "DISP_MM=+1.6", "DET_MIN=0.0",
                    "TRUNC_MM=60", "AXIS_SRC=full", "ALIGN_AXIS=1",
                    "FIX_TALUS=transverse", "CART_NEW=1", "TIE=1",
                    "CART_DIR=temp/thums/cart_m1", "POSE_IE=0"],
                   cwd=str(R), capture_output=True, text=True,
                   encoding="utf-8", errors="replace", timeout=3600)
out = (r.stdout or "") + (r.stderr or "")
(io.open(R / "temp/pyfebio_demo/_p1v2.log", "w", encoding="utf-8")
 .write(out))
mF = re.search(r"两顶面合力.*?\|F\|=([\d.]+)", out, re.S)
F = float(mF.group(1)) if mF else float("nan")
normal = ("N O R M A L" in out) or ("NORMAL TERMINATION" in out)
print(f"[跑模型] F = {F:.4f} N   NORMAL={normal}")

# ── 2) xplt -> h5 ──
xplt = W / "ankle_lig_WK.xplt"
h5 = W / "ankle_lig_WK_stress.h5"
import pyfebio.xplt as xp_mod
_orig = xp_mod._parse_domain_section


def _patched(buf):
    doms = _orig(buf)
    for i, d in enumerate(doms):
        d["id"] = np.array([i + 1], dtype=np.int32)
        nm = d.get("name")
        if (not nm) or (isinstance(nm, bytes) and not nm.strip(b"\x00")) \
                or (isinstance(nm, str) and nm.strip() == ""):
            d["name"] = f"__unnamed_{i+1}"
    return doms


xp_mod._parse_domain_section = _patched
try:
    xp_mod.to_hdf5(str(xplt), str(h5))
finally:
    xp_mod._parse_domain_section = _orig


def von_mises(s):
    s = np.asarray(s, dtype=np.float64)
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(0.5*((xx-yy)**2 + (yy-zz)**2 + (zz-xx)**2)
                   + 3.0*(xy**2 + yz**2 + xz**2))


import h5py
with h5py.File(h5, "r") as f:
    sids = sorted(f["states"].keys(), key=int)
    sg = f["states"][sids[-1]]["element_data"]["stress"]
    data = {d: von_mises(np.asarray(sg[d][:])) for d in sg.keys()}

# ── 3) 强度映射（剥掉 __<elem> 后缀）──
STR = {
    "tibia_cort": ("tibia_ends", 30, 70), "tibia_spon": ("tibia_ends", 30, 70),
    "tibia_marrow": ("tibia_ends", 30, 70),
    "fibula_cort": ("fibula_ends", 30, 70), "fibula_spon": ("fibula_ends", 30, 70),
    "talus_cort": ("talus*", 30, 70), "talus_spon": ("talus*", 30, 70),
    "calcaneus_cort": ("calcaneus", 100, 150),
    "calcaneus_spon": ("calcaneus", 100, 150),
    "navicular_spon": ("navicular*", 30, 70),
    "tibiotalar_tibia": ("软骨", None, None),
    "tibiotalar_talus": ("软骨", None, None),
}

print("\n" + "=" * 108)
print(f"🦴 P1 骨折裕度表 @ 压缩 {F:.1f} N（DISP_MM=+1.6，NORMAL={normal}）")
print("   裕度 = σ_vm / σ_strength（<1 未越阈）；外推断裂载荷 = F × σ_strength / σ_vm_max")
print("=" * 108)
print(f"{'域':<34}{'单元':>6}{'σ_max':>9}{'p99':>8}{'强度拉/压':>11}"
      f"{'裕度拉':>8}{'裕度压':>8}{'外推F_断(拉)':>13}")
print("-" * 108)
rows = []
for d in sorted(data.keys()):
    vm = data[d]
    base = re.sub(r"__(hex8|tet4|penta6|tri3|quad4)$", "", d)
    lab, st_, sc_ = STR.get(base, ("—", None, None))
    n = len(vm)
    mx, p99 = float(vm.max()), float(np.percentile(vm, 99))
    if st_ is None:
        print(f"{d:<34}{n:>6}{mx:>9.2f}{p99:>8.2f}{'—':>11}{'—':>8}{'—':>8}{'—':>13}")
        continue
    mr, mc = mx/st_, mx/sc_
    ffrac = F * st_ / mx
    rows.append((d, n, mx, p99, st_, sc_, mr, mc, ffrac))
    print(f"{d:<34}{n:>6}{mx:>9.2f}{p99:>8.2f}{f'{st_}/{sc_}':>11}"
          f"{mr:>8.2f}{mc:>8.2f}{ffrac:>13.0f}")
print("-" * 108)
print("❗ 已知伪影：`tibiotalar_talus__penta6` 的 σ_max=104.85 是**近退化软骨楔形单元**"
      "（detJ 2.4e-3）的奇点 ⇒ 看 p99 而非 max。")
print("\n【只看承重骨，按裕度排序】")
for d, n, mx, p99, st_, sc_, mr, mc, ff in sorted(rows, key=lambda x: -x[6]):
    print(f"  {d:<30} σ_max={mx:6.2f}  p99={p99:5.2f}  "
          f"裕度(拉)={mr:.2f}  外推断裂载荷 ≈ {ff:6.0f} N")
print("\n⚠️ caveat（bone.py 自带）：材料强度 ≠ 整体骨失效载荷；此表**系统性低估**骨折风险")
print("   （材料强度×截面路线对跟骨高估约 k 倍）。独立路线：跟骨用『足跟接触力』判据。")
