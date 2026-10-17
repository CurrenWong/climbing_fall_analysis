"""P1 v3（旗舰交付物）：不同姿势的骨折裕度表。

固定 DISP_MM=+1.0（此前 9 姿势全 NORMAL），扫 POSE_IE，逐姿势读骨 σ_vm。
"""
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
import pyfebio.xplt as xp_mod
import h5py

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


def von_mises(s):
    s = np.asarray(s, dtype=np.float64)
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(0.5*((xx-yy)**2 + (yy-zz)**2 + (zz-xx)**2)
                   + 3.0*(xy**2 + yz**2 + xz**2))


STR = {"tibia_cort": 30, "tibia_spon": 30, "fibula_cort": 30, "fibula_spon": 30,
       "talus_cort": 30, "talus_spon": 30, "calcaneus_cort": 100,
       "calcaneus_spon": 100, "navicular_spon": 30}
WATCH = ["fibula_cort", "tibia_cort", "talus_cort", "calcaneus_cort",
         "fibula_spon", "tibia_spon"]

res = {}
for ie in (-3.0, -1.0, 0.0, 1.0, 3.0, 5.0):
    env = dict(os.environ)
    env.update({"DRIVER": "dispx", "DISP_MM": "+1.0", "DET_MIN": "0.0",
                "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
                "FIX_TALUS": "transverse", "CART_NEW": "1", "TIE": "1",
                "CART_DIR": str(R / "temp/thums/cart_m1"),
                "POSE_IE": f"{ie:+.1f}", "PYTHONIOENCODING": "utf-8"})
    r = subprocess.run([str(PY), "-u", str(R / "scripts/ankle_fe/acceptance.py"),
                        "DRIVER=dispx", "DISP_MM=+1.0", "DET_MIN=0.0",
                        "TRUNC_MM=60", "AXIS_SRC=full", "ALIGN_AXIS=1",
                        "FIX_TALUS=transverse", "CART_NEW=1", "TIE=1",
                        "CART_DIR=temp/thums/cart_m1", f"POSE_IE={ie:+.1f}"],
                       cwd=str(R), env=env, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    o = (r.stdout or "") + (r.stderr or "")
    (io.open(R / f"temp/pyfebio_demo/_pose_{ie:+.0f}.log", "w",
             encoding="utf-8").write(o))
    mF = re.search(r"两顶面合力.*?\|F\|=([\d.]+)", o, re.S)
    F = float(mF.group(1)) if mF else float("nan")
    normal = ("N O R M A L" in o) or ("NORMAL TERMINATION" in o)
    xplt = W / "ankle_lig_WK.xplt"
    h5 = W / f"_pose_{ie:+.0f}.h5"
    xp_mod._parse_domain_section = _patched
    try:
        xp_mod.to_hdf5(str(xplt), str(h5))
    finally:
        xp_mod._parse_domain_section = _orig
    with h5py.File(h5, "r") as f:
        sids = sorted(f["states"].keys(), key=int)
        sg = f["states"][sids[-1]]["element_data"]["stress"]
        d0 = {d: von_mises(np.asarray(sg[d][:])) for d in sg.keys()}
    vm = {}
    for base in WATCH:
        ks = [k for k in d0 if re.sub(r"__(hex8|tet4|penta6)$", "", k) == base]
        if ks:
            vm[base] = float(d0[ks[0]].max())
    res[ie] = (F, normal, vm)
    print(f"  POSE_IE={ie:+5.1f}  F={F:7.2f} N  NORMAL={str(normal):<5}  "
          + "  ".join(f"{k}={v:6.2f}" for k, v in vm.items()))

print("\n" + "=" * 104)
print("🦴 不同姿势骨折裕度表（σ_vm_max / σ_strength；σ=30 MPa 骨 / 100 MPa 跟骨）")
print("=" * 104)
print(f"{'POSE_IE':>9}{'F(N)':>9}" + "".join(f"{k:>14}" for k in WATCH))
print("-" * 104)
for ie, (F, ok, vm) in res.items():
    print(f"{ie:>+9.1f}{F:>9.1f}"
          + "".join(f"{vm.get(k, float('nan')):>14.2f}" for k in WATCH))
print("-" * 104)
print("⇒ 裕度 = σ/30（跟骨 σ/100）")
print(f"{'POSE_IE':>9}{'F(N)':>9}" + "".join(f"{k:>14}" for k in WATCH))
print("-" * 104)
for ie, (F, ok, vm) in res.items():
    cells = []
    for k in WATCH:
        v = vm.get(k)
        cells.append("—" if v is None
                     else f"{v/STR[k]:.2f}")
    print(f"{ie:>+9.1f}{F:>9.1f}" + "".join(f"{c:>14}" for c in cells))
print("-" * 104)
mx = max(res.items(), key=lambda kv: (kv[1][2].get("fibula_cort") or 0))
print(f"⚠️ 最大裕度姿势 = POSE_IE={mx[0]:+.1f}（腓骨 {mx[1][2].get('fibula_cort', 0):.2f} MPa "
      f"⇒ 裕度 {mx[1][2].get('fibula_cort', 0)/30:.2f}）")
print("⚠️ caveat：材料强度≠整体骨失效（bone.py 自带）⇒ 本表**低估**骨折风险；")
print("   且软骨域有退化单元奇点（σ_max 104.85 vs p99 6.14）⇒ 软骨应力暂不可用。")
