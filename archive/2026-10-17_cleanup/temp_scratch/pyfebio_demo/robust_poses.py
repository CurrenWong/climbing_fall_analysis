"""第1项：相对结论稳健性 —— 旧/新软骨网格 × 同组姿势。

每格：跑 acceptance(DISP=+1.0, POSE_IE) → 读骨 σ_vm → 记录 F 与关键骨应力。
最后比较「外−内跨度」与姿势排序是否在两种网格下一致。
"""
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
import pyfebio.xplt as xp_mod       # noqa: E402
import h5py                          # noqa: E402

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


def vm(s):
    s = np.asarray(s, float)
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(0.5*((xx-yy)**2+(yy-zz)**2+(zz-xx)**2)
                   + 3.0*(xy**2+yz**2+xz**2))


POSES = (-3.0, -1.0, 0.0, 1.0, 3.0, 5.0)
MESHES = {"old": "temp/thums/cart_m1", "fix": "temp/thums/cart_m1_fix"}
WATCH = ("fibula_cort", "tibia_cort", "talus_cort", "calcaneus_cort")

RES = {}
for mtag, cdir in MESHES.items():
    for ie in POSES:
        env = dict(os.environ)
        env.update({"DRIVER": "dispx", "DISP_MM": "+1.0", "DET_MIN": "0.0",
                    "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
                    "FIX_TALUS": "transverse", "CART_NEW": "1", "TIE": "1",
                    "CART_DIR": cdir, "POSE_IE": f"{ie:+.1f}",
                    "PYTHONIOENCODING": "utf-8"})
        r = subprocess.run(
            [str(PY), "-u", str(R / "scripts/ankle_fe/acceptance.py"),
             "DRIVER=dispx", "DISP_MM=+1.0", "DET_MIN=0.0", "TRUNC_MM=60",
             "AXIS_SRC=full", "ALIGN_AXIS=1", "FIX_TALUS=transverse",
             "CART_NEW=1", "TIE=1", "CART_DIR=" + cdir, f"POSE_IE={ie:+.1f}"],
            cwd=str(R), env=env, capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=3600)
        o = (r.stdout or "") + (r.stderr or "")
        mF = re.search(r"两顶面合力.*?\|F\|=([\d.]+)", o, re.S)
        F = float(mF.group(1)) if mF else float("nan")
        ok = ("N O R M A L" in o) or ("NORMAL TERMINATION" in o)
        # 转 h5 读应力
        h5 = W / f"_rob_{mtag}_{ie:+.0f}.h5"
        xplt = W / "ankle_lig_WK.xplt"
        st = {}
        try:
            xp_mod._parse_domain_section = _patched
            xp_mod.to_hdf5(str(xplt), str(h5))
        finally:
            xp_mod._parse_domain_section = _orig
        with h5py.File(h5, "r") as f:
            sids = sorted(f["states"].keys(), key=int)
            sg = f["states"][sids[-1]]["element_data"]["stress"]
            for base in WATCH:
                ks = [k for k in sg.keys()
                      if re.sub(r"__(hex8|tet4|penta6)$", "", k) == base]
                if ks:
                    st[base] = float(vm(np.asarray(sg[ks[0]][:])).max())
        RES[(mtag, ie)] = (F, ok, st)
        print(f"  {mtag:<4} POSE_IE={ie:+5.1f}  F={F:7.2f} N  NORMAL={str(ok):<5}"
              + "".join(f"  {k}={st.get(k, float('nan')):6.2f}" for k in WATCH))

print("\n" + "=" * 108)
print("对照：同一姿势下 旧 vs 新 网格")
print("=" * 108)
print(f"{'POSE_IE':>9}{'F_old':>9}{'F_fix':>9}{'F比':>7}"
      + "".join(f"{k+'_old':>13}{k+'_fix':>13}{'比':>6}" for k in WATCH))
for ie in POSES:
    Fo, _, so = RES[("old", ie)]
    Fn, _, sn = RES[("fix", ie)]
    row = f"{ie:>+9.1f}{Fo:>9.2f}{Fn:>9.2f}{Fn/Fo:>7.2f}"
    for k in WATCH:
        a, b = so.get(k, np.nan), sn.get(k, np.nan)
        row += f"{a:>13.2f}{b:>13.2f}{b/a:>6.2f}"
    print(row)

print("\n" + "=" * 108)
print("⭐ 相对结论：外翻侧(−) 与 内翻侧(+) 的跨度  σ(+5°) − σ(−3°)")
print("=" * 108)
for k in WATCH:
    o = RES[("old", 5.0)][2].get(k, np.nan) - RES[("old", -3.0)][2].get(k, np.nan)
    n = RES[("fix", 5.0)][2].get(k, np.nan) - RES[("fix", -3.0)][2].get(k, np.nan)
    print(f"  {k:<18} 旧 {o:>8.3f}   新 {n:>8.3f}   "
          f"{'同号 ✓' if o*n > 0 else '**异号 🚨**'}   新/旧 = {n/o if o else float('nan'):>6.2f}")

print("\n" + "=" * 108)
print("⭐ 姿势风险排序（按各骨 σ 从大到小）")
print("=" * 108)
for mt in ("old", "fix"):
    for k in WATCH:
        rank = sorted(POSES, key=lambda ie: -RES[(mt, ie)][2].get(k, -1))
        print(f"  {mt:<4} {k:<18} " + " > ".join(f"{x:+.0f}" for x in rank))
