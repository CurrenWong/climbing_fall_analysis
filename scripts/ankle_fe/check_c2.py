"""C2 重定义：姿势输入保真性 + 解剖语义 + 韧带方向（纯几何，NO_RUN=1，**不求解**）。

背景 —— 原 C2「末态出现非零内外翻角」为什么废：
    原实现测的是**距骨 t=0→末态**的转动，而 `POSE_IE` 是把足部组**烘焙进参考构型**
    ⇒ 结构上测不到（实测仅 0.0014~0.0033°）。**架构错位，不是模型错。**

新 C2 = 验"姿势作为**输入**被忠实实现、且语义方向正确"：

  C2a-1 幅值：|实测相对转角 − 请求角| ≤ 0.5°      （距骨 vs 胫骨，逐骨 Kabsch）
  C2a-2 轴  ：实测轴与 AP（足前后轴）夹角 ≤ 5°
  C2a-3 语义：**三个独立锚定判据必须一致**
              ① 足底法向：外侧分量减小 ⇒ 内翻（内翻的定义 = 足底朝向内侧）
              ② 足缘高度：外侧缘下沉 ⇒ 内翻（内翻 = 站在外侧缘上）
              ③ 韧带长度：外侧韧带(ATFL/CFL/PTFL)被拉长 ⇒ 内翻
                 （解剖学：内翻 = 外侧韧带扭伤机制，ATFL/CFL 必受拉）
              ★ ①② 是纯几何、③ 是几何+解剖部件号（韧带左右侧来自 THUMS 部件号）

  ⚠️ 若 ①② 一致但 ③ 相反 ⇒ **韧带插入点几何错了**（不是读数错）
     已知：ATFL 建模长 36.99 mm vs 解剖 15~20 mm（≈2× 过长）、CFL 34.33 vs ~28.5
  ⚠️ 跨姿势必须用**中性网格的 ID 列表**索引（见 `[ERR-20261006-001]`）

用法：
    python check_c2.py                    # IE = 0,±1,±3,±5,+10
    python check_c2.py IE=0,5,-5          # 自定义
    python check_c2.py AXIS=DF
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
C2D = WORK / "c2"
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
C2D.mkdir(parents=True, exist_ok=True)

BASE = {
    "DET_MIN": "0.0", "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
    "FIX_TALUS": "transverse", "DRIVER": "dispx", "DISP_MM": "-1.0",
    "LOAD_SIGN": "+1", "CART_NEW": "1", "TIE": "1",
    "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "LOG_DATA": "1",
    "NO_RUN": "1",
}
# 韧带分侧（右侧 = 连腓骨；依据 THUMS 部件号，见 scripts/ankle_fe/ligaments.py）
LAT = ("ATFL", "CFL", "PTFL")
MED = ("MD1", "MD2", "MD3", "MD4")


def parse_feb(p: Path):
    t = p.read_text(encoding="latin-1", errors="replace")
    xyz = {}
    for m in re.finditer(r'<node id="(\d+)">([-\d.eE+,]+)</node>', t):
        xyz[int(m.group(1))] = np.array(
            [float(x) for x in m.group(2).split(",")], float)
    sets = {}
    for m in re.finditer(r'<NodeSet name="([^"]+)">([\d,]*)</NodeSet>', t):
        sets[m.group(1)] = [int(x) for x in m.group(2).split(",") if x]
    return xyz, sets


def rigid_rot(P0, P1):
    c0, c1 = P0.mean(0), P1.mean(0)
    H = (P0 - c0).T @ (P1 - c1)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    res = np.linalg.norm((P1 - c1) - (P0 - c0) @ R.T, axis=1)
    tot = np.linalg.norm(P0 - c0, axis=1)
    return R, float(res.mean() / max(tot.mean(), 1e-12))


def rot_angle_axis(R):
    c = np.clip((np.trace(R) - 1) / 2, -1, 1)
    ang = float(np.degrees(np.arccos(c)))
    v = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    n = np.linalg.norm(v)
    return ang, (v / n if n > 1e-12 else np.zeros(3))


def median_plane_normal(P):
    Q = P - P.mean(0)
    _, _, Vt = np.linalg.svd(Q, full_matrices=False)
    return Vt[2] / np.linalg.norm(Vt[2])


def build(varname: str, val: float, tag_out: str):
    env = os.environ.copy()
    env.update(BASE)
    env[varname] = str(val)
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    src = WORK / f"ankle_lig_{env.get('MATSET', 'WK')}.feb"
    dst = C2D / f"{tag_out}.feb"
    if not src.exists():
        print(f"  ✗ 建失败 {tag_out}\n{out[-1200:]}")
        return None, out
    shutil.copy2(src, dst)
    shutil.copy2(WORK / f"lig_read_{env.get('MATSET', 'WK')}.npz",
                 C2D / f"{tag_out}_lig.npz")
    return dst, out


def load_axis_from(out: str, xyz: dict, sets: dict):
    m = re.search(r"加载轴（胫骨->距骨）= \[([^\]]+)\]", out)
    if m:
        a = np.array([float(x) for x in m.group(1).split()], float)
        return a / np.linalg.norm(a)
    a = (np.array([xyz[i] for i in sets["talus_low_ns"]]).mean(0)
         - np.array([xyz[i] for i in sets["tibia_top_ns"]]).mean(0))
    return a / np.linalg.norm(a)


# ─────────────────────────── 主流程 ───────────────────────────
args = {a.split("=", 1)[0]: a.split("=", 1)[1]
        for a in sys.argv[1:] if "=" in a}
VAR = args.get("AXIS", "IE").upper()
vals = [float(x) for x in args.get(VAR, "0,1,-1,3,-3,5,-5,10").split(",")]

print(f"=== C2 重定义检查：{VAR} 扫描 {vals}（只建不解 NO_RUN=1）===")
ref_feb, ref_out = build(f"POSE_{VAR}", 0.0, f"{VAR}_0")
if ref_feb is None:
    sys.exit(1)
xyz0, sets0 = parse_feb(ref_feb)
axis = load_axis_from(ref_out, xyz0, sets0)
print(f"参考构型（{VAR}=0）：{ref_feb.name}  {len(xyz0)} 节点")
print(f"加载轴 axis = {np.round(axis, 4)}")

fib = np.array([xyz0[i] for i in sets0["fibula_top_ns"]]).mean(0)
tib = np.array([xyz0[i] for i in sets0["tibia_top_ns"]]).mean(0)
lat = fib - tib
lat = lat - (lat @ axis) * axis
lat = lat / np.linalg.norm(lat)
print(f"外侧方向 lat（腓骨−胫骨，⊥加载轴）= {np.round(lat, 4)}"
      f"   |lat·axis|={abs(float(lat @ axis)):.4f}")

# AP（足前后轴）：`POSE_IE=0` 时脚本的姿势段**不执行** ⇒ 拿不到它算的 `_ap`。
#   ⇒ 用一个 ε 姿势（0.001°，位移量 ~1e-5 mm 可忽略）建一次，只为采 stdout 里的 AP/ML。
#   不要自己近似（先前用 NodeSet 质心估 AP，与脚本口径差 12° ⇒ C2a-2 假 FAIL）。
_eps_feb, _eps_out = build(f"POSE_{VAR}", 0.001, f"{VAR}_eps")
ap = None
_m = re.search(r"前后轴 AP = \[([^\]]+)\]", _eps_out)
if _m:
    ap = np.array([float(x) for x in _m.group(1).split()], float)
    ap = ap / np.linalg.norm(ap)
print(f"前后轴 AP（脚本口径，ε 构建采集）= {np.round(ap, 4)}"
      f"   |AP·axis|={abs(float(ap @ axis)):.4f}" if ap is not None else
      "⚠ 未采到 AP（C2a-2 将跳过）")

# ★ `[ERR-20261006-001]`：这些 NodeSet 由 `extreme_faces()` 选出，**姿势一变集合就变**
#   （实测 calcaneus_bottom_ns 中性 194 → 某姿势 227 节点）⇒ 跨姿势必须固定用 r=0 的 ID。
NS_CAL = list(sets0["calcaneus_bottom_ns"])
NS_TAL = list(sets0["talus_low_ns"])
NS_TIB = list(sets0["tibia_top_ns"])

# 韧带杆对：拓扑与姿势无关，只从中性构建取一次
_z = np.load(C2D / f"{VAR}_0_lig.npz", allow_pickle=True)
LIG_NAMES = [str(x) for x in _z["lig_names"]]
LIG_PAIRS = _z["lig_pairs"]
LIG_NP = _z["lig_npair"]
print(f"韧带 {len(LIG_NAMES)} 条：" + "，".join(LIG_NAMES))

cal0 = np.array([xyz0[i] for i in NS_CAL])
n0 = median_plane_normal(cal0)
if n0 @ axis > 0:
    n0 = -n0
n0_lat = float(n0 @ lat)
c0_lat = cal0 @ lat
half = np.median(c0_lat)
sel_L, sel_M = c0_lat > half, c0_lat <= half
print(f"足底中性法向 n0·lat = {n0_lat:+.4f}")


def lig_dL(xyz):
    """每条韧带的平均长度变化（相对中性，mm）。"""
    out = {}
    for i, nm in enumerate(LIG_NAMES):
        n = int(LIG_NP[i])
        if n <= 0:
            out[nm] = float("nan")
            continue
        d = []
        for j in range(n):
            a, b = int(LIG_PAIRS[i, j, 0]), int(LIG_PAIRS[i, j, 1])
            if a < 0 or b < 0:
                continue
            d.append(np.linalg.norm(xyz[b] - xyz[a]) - np.linalg.norm(xyz0[b] - xyz0[a]))
        out[nm] = float(np.mean(d)) if d else float("nan")
    return out


print(f"\n{'═'*104}")
hdr = (f"{VAR:>4} {'法向Δlat':>9} {'外缘Δh':>8} {'内缘Δh':>8} "
       f"{'距骨转':>8} {'胫骨转':>8} {'相对转':>8} {'|轴·AP|':>8}  "
       f"{'外侧韧带ΔL':>26} {'内侧韧带ΔL':>18}")
print(hdr)
print("─" * 104)

rows = []
dl_rows = []
for v in vals:
    if abs(v) < 1e-12:
        rows.append((v, 0, 0, 0, 0, 0, 0, 1.0, "中性"))
        dl_rows.append((v, {nm: 0.0 for nm in LIG_NAMES}))
        print(f"{v:>4.0f} {'—':>9} {'—':>8} {'—':>8} {'—':>8} {'—':>8} "
              f"{'—':>8} {'—':>8}  {'（基准）':>26} {'':>18}")
        continue
    tag = f"{VAR}_{int(round(v*100))}"
    p, out = build(f"POSE_{VAR}", v, tag)
    if p is None:
        continue
    xyz, sets = parse_feb(p)

    def P(xy, ids):
        return np.array([xy[i] for i in ids])

    Rtal, _ = rigid_rot(P(xyz0, NS_TAL), P(xyz, NS_TAL))
    Rtib, _ = rigid_rot(P(xyz0, NS_TIB), P(xyz, NS_TIB))
    Rrel = Rtib.T @ Rtal
    a_rel, ax_rel = rot_angle_axis(Rrel)
    a_tal, _ = rot_angle_axis(Rtal)
    a_tib, _ = rot_angle_axis(Rtib)
    axap = abs(float(ax_rel @ ap)) if ap is not None else float("nan")

    cal = P(xyz, NS_CAL)
    n = median_plane_normal(cal)
    if n @ axis > 0:
        n = -n
    d_nlat = float(n @ lat) - n0_lat
    dh = (cal - cal0) @ axis
    dh_L, dh_M = float(dh[sel_L].mean()), float(dh[sel_M].mean())

    dL = lig_dL(xyz)
    dl_rows.append((v, dL))
    lat_avg = float(np.nanmean([dL[n] for n in LIG_NAMES
                                if any(k in n for k in LAT)]))
    med_avg = float(np.nanmean([dL[n] for n in LIG_NAMES
                                if any(k in n for k in MED)]))

    # 语义：①②几何 + ③韧带
    g_inv = (d_nlat < 0) and (dh_L < dh_M)
    g_ev = (d_nlat > 0) and (dh_L > dh_M)
    l_inv = lat_avg > 0 and lat_avg > med_avg      # 外侧被拉长 ⇒ 内翻
    l_ev = med_avg > 0 and med_avg > lat_avg       # 内侧被拉长 ⇒ 外翻
    sem_g = "内翻" if g_inv else ("外翻" if g_ev else "⚠几何不一致")
    sem_l = "内翻" if l_inv else ("外翻" if l_ev else "⚠韧带不一致")
    sem = sem_g if sem_g == sem_l else f"⚠冲突({sem_g}/韧带{sem_l})"
    rows.append((v, d_nlat, dh_L, dh_M, a_tal, a_tib, a_rel, axap, sem))

    lt = " ".join(f"{n.split('_')[0]}:{dL[n]:+.2f}" for n in LIG_NAMES
                  if any(k in n for k in LAT))
    md = " ".join(f"{n.split('_')[0]}:{dL[n]:+.2f}" for n in LIG_NAMES
                  if any(k in n for k in MED))
    print(f"{v:>4.0f} {d_nlat:>9.4f} {dh_L:>8.3f} {dh_M:>8.3f} {a_tal:>8.4f} "
          f"{a_tib:>8.4f} {a_rel:>8.4f} {axap:>8.4f}  {lt:>26} {md:>18}  {sem}")

print(f"\n{'='*104}\nC2a 判据汇总\n{'='*104}")
ok1 = ok2 = ok3 = True
for v, d_nlat, dh_L, dh_M, a_tal, a_tib, a_rel, axap, sem in rows:
    if abs(v) < 1e-12:
        continue
    e1 = abs(a_rel - abs(v))
    c1 = e1 <= 0.5
    # ★ 旋转是**绕 AP** 做的 ⇒ 实测轴应与 AP **平行**（|轴·AP| → 1），不是垂直。
    #   （先前把判据写反过一次：拿"垂直"当期望 ⇒ 0.9787 被误判 FAIL。）
    c2 = (axap >= float(np.cos(np.radians(5)))) if ap is not None else True
    c3 = sem in ("内翻", "外翻")
    ok1 &= c1
    ok2 &= c2
    ok3 &= c3
    print(f"  {VAR}={v:>+5.1f}°  C2a-1 幅值误差 {e1:>6.3f}° {'✓' if c1 else '✗'}"
          f"   C2a-2 |轴·AP|={axap:.4f} {'✓' if c2 else '✗'}"
          f"   C2a-3 {sem} {'✓' if c3 else '✗'}")
print(f"\n  C2a-1 幅值保真 : {'PASS ✓' if ok1 else 'FAIL ✗'}")
print(f"  C2a-2 轴向保真 : {'PASS ✓' if ok2 else 'FAIL ✗'}")
print(f"  C2a-3 语义自洽 : {'PASS ✓' if ok3 else 'FAIL ✗（见上方 ⚠）'}")

print(f"\n{'='*104}\n⭐ 符号语义核对（本项目最易错处）\n{'='*104}")
print("  代码注释：POSE_IE(内翻+ / 外翻−)")
neg = sorted({r[8] for r in rows if r[0] < 0})
pos = sorted({r[8] for r in rows if r[0] > 0})
print(f"  实测 负值 −IE ⇒ {neg}")
print(f"  实测 正值 +IE ⇒ {pos}")
if pos == ["内翻"] and neg == ["外翻"]:
    print("  ⇒ **注释正确**：+IE = 内翻 / −IE = 外翻")
elif pos == ["外翻"] and neg == ["内翻"]:
    print("  ⇒ **注释反了**：+IE = 外翻 / −IE = 内翻")
else:
    print("  ⇒ ⚠ 不单调或三判据冲突，需人工看表")

print(f"\n{'='*104}\n📏 韧带建模长度 vs 解剖（`ligaments.npz`）\n{'='*104}")
_z2 = np.load(ROOT / "temp/thums/ligaments.npz", allow_pickle=True)
_nm2, _L2 = [str(x) for x in _z2["names"]], _z2["length"]
_ANAT = {"ATFL": "15~20", "CFL": "25~32", "PTFL": "20~25"}
for nm, L in zip(_nm2, _L2):
    k = nm.split("_")[0]
    tag = f"  ← 解剖 {_ANAT[k]} mm" if k in _ANAT else ""
    print(f"  {nm:<30} {L:>7.2f} mm{tag}")
print("\n（注：C2b「载荷下姿势保留」需两次真求解，见 ANKLE_FIX_TODO.md 待办）")
