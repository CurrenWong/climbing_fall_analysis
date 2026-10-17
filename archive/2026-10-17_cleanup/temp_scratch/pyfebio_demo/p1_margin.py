"""P1 v1：踝模型骨折裕度表。

链路：跑 +1.6mm（592.6 N）→ xplt 转 h5（复用 bone_feb 补丁）
     → 逐域 von Mises → 除以 bone.py 的骨强度 ⇒ 裕度。

⚠️ caveat（bone.py 自带）：材料强度 ≠ 整体骨失效载荷；
   材料强度×截面路线对跟骨"高估约 k 倍"（即低估骨折风险）。
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

# ── 1) 跑模型到 +1.6mm ──
env = dict(os.environ)
env.update({"DRIVER": "dispx", "DISP_MM": "+1.6", "DET_MIN": "0.0",
            "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
            "FIX_TALUS": "transverse", "CART_NEW": "1", "TIE": "1",
            "CART_DIR": str(R / "temp/thums/cart_m1"), "POSE_IE": "0",
            "PYTHONIOENCODING": "utf-8"})
r = subprocess.run([str(PY), "-u", str(R / "scripts/ankle_fe/ankle_lig_feb.py")],
                   cwd=str(R), env=env, capture_output=True, text=True,
                   encoding="utf-8", errors="replace", timeout=3600)
out = (r.stdout or "") + (r.stderr or "")
m = re.findall(r"两顶面合力 = \[[^\]]*\] N\s+\|F\|=([\d.]+)", out)
F = float(m[0]) if m else float("nan")
print(f"[跑模型] rc={r.returncode}  两顶面合力 = {F:.4f} N")
if "NORMAL TERMINATION" not in out and "N O R M A L" not in out:
    print("  ⚠ 未确认 NORMAL TERMINATION，尾部：")
    print("  " + out.strip().splitlines()[-1][:150])

# ── 2) xplt -> h5（带补丁）──
xplt = W / "ankle_lig_WK.xplt"
h5 = W / "ankle_lig_WK_stress.h5"
if not xplt.exists():
    raise SystemExit(f"✗ 缺 {xplt}")

import pyfebio.xplt as xp_mod

_orig = xp_mod._parse_domain_section


def _patched(buf):
    doms = _orig(buf)
    for i, d in enumerate(doms):
        d["id"] = np.array([i + 1], dtype=np.int32)
        nm = d.get("name")
        empty = (not nm) or (isinstance(nm, bytes) and not nm.strip(b"\x00")) \
            or (isinstance(nm, str) and nm.strip() == "")
        if empty:
            d["name"] = f"__unnamed_{i+1}"
    return doms


xp_mod._parse_domain_section = _patched
try:
    xp_mod.to_hdf5(str(xplt), str(h5))
finally:
    xp_mod._parse_domain_section = _orig
print(f"[转换] {h5.name} {h5.stat().st_size/1e6:.1f} MB")

# ── 3) 逐域 von Mises ──
import h5py


def von_mises(s):
    s = np.asarray(s, dtype=np.float64)
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(0.5 * ((xx-yy)**2 + (yy-zz)**2 + (zz-xx)**2)
                   + 3.0 * (xy**2 + yz**2 + xz**2))


with h5py.File(h5, "r") as f:
    sids = sorted(f["states"].keys(), key=int)
    st = f["states"][sids[-1]]
    sg = st.get("element_data", {}).get("stress")
    print(f"[h5] {len(sids)} 个 state；末 state={sids[-1]}；"
          f"stress 域 = {list(sg.keys())}")
    doms = list(sg.keys())
    stats = {}
    for d in doms:
        vm = von_mises(np.asarray(sg[d][:]))
        stats[d] = (len(vm), vm.max(), np.percentile(vm, 99),
                    np.percentile(vm, 95), vm.mean())

# ── 4) 强度映射 + 裕度 ──
STR = {"tibia_cort": ("tibia_ends", 30, 70), "tibia_spon": ("tibia_ends", 30, 70),
       "tibia_marrow": ("tibia_ends", 30, 70),
       "fibula_cort": ("fibula_ends", 30, 70),
       "fibula_spon": ("fibula_ends", 30, 70),
       "talus_cort": ("talus(代理=tibia_ends)", 30, 70),
       "talus_spon": ("talus(代理=tibia_ends)", 30, 70),
       "calcaneus_cort": ("calcaneus", 100, 150),
       "calcaneus_spon": ("calcaneus", 100, 150)}

print("\n" + "=" * 96)
print(f"骨折裕度表 @ 压缩 {F:.1f} N（DISP_MM=+1.6）  判据 = σ_vm_max / σ_strength")
print("=" * 96)
print(f"{'域':<18}{'单元':>7}{'σ_vm_max':>10}{'p99':>9}{'p95':>9}{'均值':>9}"
      f"{'强度(拉/压)':>14}{'裕度(拉/压)':>14}")
print("-" * 96)
for d in doms:
    n, mx, p99, p95, mu = stats[d]
    key = d.lower()
    lab, st_, sc_ = STR.get(key, ("?（未映射）", None, None))
    if st_ is not None:
        mr = f"{mx/st_:.2f} / {mx/sc_:.2f}"
        sl = f"{st_}/{sc_}"
    else:
        mr, sl = "—", "—"
    print(f"{d:<18}{n:>7}{mx:>10.2f}{p99:>9.2f}{p95:>9.2f}{mu:>9.2f}"
          f"{sl:>14}{mr:>14}")
print("-" * 96)
print("解读：裕度 <1 = 未越阈；>1 = 该域最大应力已超材料强度（**不必然是整体骨折**）。")
print("⚠️ caveat：bone.py 自带警告 —— 材料强度 ≠ 整体骨失效载荷；")
print("   材料强度路线对跟骨高估约 k 倍 ⇒ 此表系统性**低估**骨折风险。")
print("   独立佐证路线：跟骨用『足跟接触力』判据（Voo HIPC / Barnes 2019）。")
