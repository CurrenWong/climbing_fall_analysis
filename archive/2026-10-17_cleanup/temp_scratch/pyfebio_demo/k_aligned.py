"""关节实际压缩量 δ + 真实关节刚度 k=F/δ，对照物理估计 E·A/t。"""
import inspect
import os
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X

from closest_vec import closest_vec

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio"
PY = ROOT / ".venv/Scripts/python.exe"
AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)
_src = textwrap.dedent(inspect.getsource(X.parse_mesh))
for _o, _n in (('domain["name"]', 'domain.get("name", str(domain["id"][0]))'),
               ('surface["name"]', 'surface.get("name", str(surface["id"][0]))'),
               ("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
               ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
               ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
               ("{nodeset['name']}", "{nodeset.get('name', str(nodeset['id'][0]))}"),
               ("{elementset['name']}",
                "{elementset.get('name', str(elementset['id'][0]))}")):
    _src = _src.replace(_o, _n)
_ns = X.__dict__
exec(compile(_src, "<p>", "exec"), _ns)
X.parse_mesh = _ns["parse_mesh"]
_ps = textwrap.dedent(inspect.getsource(X.parse_state))
for _a, _b in (
    ('set_name = mesh_dict[set_lut[key]][set_id + 1]',
     'set_name = mesh_dict[set_lut[key]].get(set_id + 1) or '
     'mesh_dict[set_lut[key]].get(int(set_id) + 1) or f"set{set_id}"'),
    ('set_name = mesh_dict[set_lut[key]][set_id]',
     'set_name = mesh_dict[set_lut[key]].get(set_id) or '
     'mesh_dict[set_lut[key]].get(int(set_id)) or f"set{set_id}"'),
):
    _ps = _ps.replace(_a, _b)
_ns2 = X.__dict__
exec(compile(_ps, "<p2>", "exec"), _ns2)
X.parse_state = _ns2["parse_state"]

K = 8
g = np.linspace(0, 1, K + 1)
tbc = []
for i in range(K):
    for j in range(K - i):
        a, b, c = (i / K, j / K), ((i + 1) / K, j / K), (i / K, (j + 1) / K)
        d = ((i + 1) / K, (j + 1) / K)
        if i + j + 1 <= K:
            tbc.append((a, b, c))
        if i + j + 2 <= K:
            tbc.append((b, d, c))
tbc = np.array(tbc)


def subdiv(T):
    def s(k):
        return (tbc[:, k, 0][:, None, None] * T[None, :, 0, :]
                + tbc[:, k, 1][:, None, None] * T[None, :, 1, :]
                + (1 - tbc[:, k].sum(1))[:, None, None] * T[None, :, 2, :])
    return np.stack([s(0), s(1), s(2)], 2).reshape(-1, 3, 3)


base = dict(os.environ)
# AUGLAG 分支下 build_feb.py:696 自动写 auto_penalty=0 => 天然符合项目铁律
base.update({"TIE": "0", "LAUGON": os.environ.get("LAUGON", "AUGLAG"),
             "TIME_STEPS": "20", "DRIVER": "disp", "DISP_MM": "0.1472"})
PENS = [p for p in os.environ.get("PENS", "2,8,32").split(",") if p]
rows = []
for pen in PENS:
    env = dict(base)
    env["PENALTY"] = pen
    subprocess.run([str(PY), "-u", "scripts/ankle_fe/build_feb.py"], cwd=str(ROOT),
                   env=env, capture_output=True, text=True,
                   timeout=int(os.environ.get("FEB_TIMEOUT", "1800")))
    fp, xp = W / "ankle_tibiotalar.feb", W / "ankle_tibiotalar.xplt"
    # ── AUTO0=1（2026-10-09 E1）：关掉 auto_penalty 重跑同一 deck ─────────────
    # FEBio 官方语义：auto 开启时 <penalty> 只是乘性 scale factor（基值=primary
    # surface 刚度+网格尺度自动估计）⇒ 我们扫的不是物理 k。项目铁律见
    # ankle_lig_feb.py:127-128（罚因子扫描必须关 auto_penalty）。
    if os.environ.get("AUTO0") == "1":
        _t = fp.read_text(encoding="latin-1")
        assert _t.count("<auto_penalty>1</auto_penalty>") == 1, "deck 无 auto_penalty=1"
        fp2 = W / "ankle_tibiotalar_auto0.feb"
        fp2.write_text(_t.replace("<auto_penalty>1</auto_penalty>",
                                  "<auto_penalty>0</auto_penalty>"),
                       encoding="latin-1")
        r2 = subprocess.run([r"D:\Program\FEBioStudio\bin\febio4.exe", "-i", fp2.name],
                            cwd=str(W), capture_output=True, text=True, timeout=1800)
        xp = W / "ankle_tibiotalar_auto0.xplt"
        print(f"    [AUTO0] PENALTY={pen} FEBio rc={r2.returncode} "
              f"xplt={xp.exists()}", flush=True)
    txt = fp.read_text(encoding="latin-1")
    XYZ = {}
    for m in re.finditer(
            r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
        XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                         float(m.group(4))])
    sf = {}
    for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
        sf[m.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                          re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))]
    X.to_hdf5(str(xp), str(W / "_k.h5"))
    with h5py.File(W / "_k.h5", "r") as h:
        ns = sorted((int(q) for q in h["states"]), key=int)
        u0 = np.asarray(h[f"states/{ns[0]}/node_data/displacement/all"][()], float)
        u1 = np.asarray(h[f"states/{ns[-1]}/node_data/displacement/all"][()], float)
        R = np.asarray(h[f"states/{ns[-1]}/node_data/reaction forces/all"][()], float)
    NX = max(XYZ) + 1
    U0 = np.zeros((NX, 3))
    U0[:len(u0)] = u0
    U1 = np.zeros((NX, 3))
    U1[:len(u1)] = u1
    Tti0 = np.array([np.array([XYZ[i] + U0[i - 1] for i in f])
                     for f in sf["surf_tibiotalar_tibia"]])
    Tta = np.array([np.array([XYZ[i] for i in f]) for f in sf["surf_tibiotalar_talus"]])
    sub = subdiv(Tta)
    ar = np.linalg.norm(np.cross(sub[:, 1] - sub[:, 0], sub[:, 2] - sub[:, 0]),
                        axis=1) / 2
    cen = sub.mean(1)
    Tti1 = np.array([np.array([XYZ[i] + U1[i - 1] for i in f])
                     for f in sf["surf_tibiotalar_tibia"]])
    g0 = closest_vec(cen, Tti0)[0]
    g1 = closest_vec(cen, Tti1)[0]
    F = float((R @ AX).sum())
    sel = g0 < 0.5
    A0 = ar[sel].sum()
    dl = g0[sel] - g1[sel]
    A1 = ar[g1 < 0.5].sum()
    # ── 口径修复（R11，2026-10-09）────────────────────────────────────────
    # 旧口径 k = F_total / δmax：分子是**全接触面**反力积分、分母是**单点最大**
    # 闭合 ⇒ 分子分母不对齐，系统性偏高。
    # 对齐口径（均匀压缩）：k = (F/A) / (δ̄/t) = F·t/(A·δ̄)，
    # 其中 A/δ̄ 只统计**真正受压**的点（dl>0），与 k_phys=E·A/t 同源。
    eng = dl > 0
    A_eng = float(ar[sel][eng].sum())
    dbar = float(dl[eng].mean()) if eng.any() else 0.0
    rows.append((pen, F, float(dl.max()) if sel.any() else 0.0,
                 float(np.median(dl)) if sel.any() else 0.0, A0, A1,
                 A_eng, dbar, int(eng.sum())))
    print(f"  PENALTY={pen:>6s}  F={F:8.2f} N  δmax={rows[-1][2]:.5f} "
          f"δ中位={rows[-1][3]:.5f} mm  闭合面积 {A0:.0f}→{A1:.0f} mm²  "
          f"受压面积={A_eng:.0f} mm²({eng.sum()}点) δ̄={dbar:.5f} mm", flush=True)

print("\n" + "=" * 92)
print("关节刚度：旧口径(F/δmax) vs 对齐口径(F·t/(A·δ̄))，t=3.0mm, E=12MPa")
print("=" * 92)
print(f"  {'PENALTY':>8} {'F(N)':>9} {'δmax':>9} {'δ̄':>9} {'A_eng':>8} "
      f"{'k_old':>8} {'k_aligned':>11} {'E·A/t(本模型)':>14} {'A/t×E vs 1059':>15}")
for pen, F, dm, dmed, A0, A1, A_eng, dbar, npts in rows:
    k_old = F / dm if dm > 0 else float("nan")
    k_new = F * 3.0 / (A_eng * dbar) if (A_eng > 0 and dbar > 0) else float("nan")
    kphys = 12.0 * A_eng / 3.0
    print(f"  {pen:>8} {F:>9.2f} {dm:>9.5f} {dbar:>9.5f} {A_eng:>8.0f} "
          f"{k_old:>8.1f} {k_new:>11.1f} {kphys:>14.1f} "
          f"{k_new / 1059 * 100:>14.1f}%")
print(f"\n  Anderson §2.4 反推估计：k ≈ E·A/t = 12×300/3.4 ≈ 1059 N/mm"
      f"（A=300mm²/t=3.4mm 均为 **Anderson 的几何**，非本模型）")
