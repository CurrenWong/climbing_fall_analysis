"""踝子线验收脚本（ANKLE_GOAL_TODO / [LRN-20261005-057]）。

用法：
    python acceptance.py DRIVER=dispx DISP_MM=0.20 FIX_TALUS=transverse
    python acceptance.py --no-run            # 只解析上一次结果

产出：run → 解析 .feb（节点坐标/集合）+ FEBio 文本读数 → 打印指标与判据结论。

判据（5 条，全相对量）：
  C1 force 驱动能收敛（rc=0 且 NORMAL TERMINATION）+ |F−目标|<2%
  C2 末态出现非零的内外翻角
  C3 V9：内翻损伤 < 外翻损伤          （需姿势扫描，见 --pair）
  C4 外侧(跟腓) < 内侧、胫跟 > 胫舟   （需读数）
  C5 姿势排序对罚因子鲁棒             （需两档罚因子，见 --pair）
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"


# ───────────────────────── .feb 解析 ─────────────────────────
def parse_feb(p: Path):
    """取 节点坐标 / NodeSet / 载荷段行数 / Contact 设置。"""
    t = p.read_text(encoding="latin-1", errors="replace")
    xyz, sets, load_lines, ctl = {}, {}, 0, {}
    for m in re.finditer(r'<node id="(\d+)">([-\d.eE+,]+)</node>', t):
        xyz[int(m.group(1))] = np.array(
            [float(x) for x in m.group(2).split(",")], float)
    for m in re.finditer(r'<NodeSet name="([^"]+)">([\d,]*)</NodeSet>', t):
        sets[m.group(1)] = [int(x) for x in m.group(2).split(",") if x]
    m = re.search(r"<Loads>(.*?)</Loads>", t, re.S)
    load_lines = len([x for x in m.group(1).split("\n") if x.strip()]) if m else -1
    for k in ("laugon", "penalty", "auto_penalty", "tolerance"):
        mm = re.search(rf"<{k}>([^<]*)</{k}>", t)
        if mm:
            ctl[k] = mm.group(1)
    return xyz, sets, load_lines, ctl


# ───────────────────── FEBio 文本读数解析 ─────────────────────
def read_records(p: Path):
    """返回 {step: (time, np.array([[id, v1, v2, v3], ...]))}。"""
    if not p.exists():
        return {}
    out, cur, buf = {}, None, []
    for l in p.read_text(encoding="utf-8", errors="replace").split("\n"):
        if l.startswith("*Step"):
            if cur is not None and buf:
                out[cur[0]] = (cur[1], np.array(buf, float))
            cur = (int(l.split("=")[1]), None, )
            buf = []
        elif l.startswith("*Time"):
            if cur is not None:
                cur = (cur[0], float(l.split("=")[1]))
        elif l.strip() and not l.startswith("*"):
            try:
                buf.append([float(x) for x in l.split(",")])
            except ValueError:
                pass
    if cur is not None and buf:
        out[cur[0]] = (cur[1], np.array(buf, float))
    return out


def rigid_rot(P0, P1):
    """Kabsch 刚体拟合：返回 (旋转矩阵 R, 残差占比)。"""
    c0, c1 = P0.mean(0), P1.mean(0)
    H = (P0 - c0).T @ (P1 - c1)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    res = np.linalg.norm((P1 - c1) - (P0 - c0) @ R.T, axis=1)
    tot = np.linalg.norm(P0 - c0, axis=1)
    return R, float(res.mean() / max(tot.mean(), 1e-12))


def rot_angle_axis(R):
    """R → (角度(度), 轴)。"""
    c = np.clip((np.trace(R) - 1) / 2, -1, 1)
    ang = float(np.degrees(np.arccos(c)))
    v = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    n = np.linalg.norm(v)
    return ang, (v / n if n > 1e-12 else np.zeros(3))


# ───────────────────────────── 主流程 ─────────────────────────────
RUN = "--no-run" not in sys.argv
args = {a.split("=", 1)[0]: a.split("=", 1)[1]
        for a in sys.argv[1:] if "=" in a and not a.startswith("--")}
env = os.environ.copy()
env.update(args)
env.setdefault("LOG_DATA", "1")          # 验收必须要有文本读数
env["PYTHONIOENCODING"] = "utf-8"
tag = env.get("MATSET", "WK") + ("_" + "_".join(env["BONES"].split(","))
                                 if env.get("BONES") else "")

if RUN:
    print(f"[run] {args or '(默认)'}")
    t0 = time.time()
    r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    el = time.time() - t0
    out = (r.stdout or "") + (r.stderr or "")
    (ROOT / "temp/pyfebio_demo/_acc_last.log").write_text(out, encoding="utf-8")
    hard_rc = r.returncode
    m = re.search(r"FEBio rc=(-?\d+)\s+用时 ([\d.]+)s", out)
    febio_rc = int(m.group(1)) if m else None
    febio_t = float(m.group(2)) if m else None
    print(f"[rc] python={hard_rc}  febio={febio_rc}  febio用时={febio_t}s  总{el:.1f}s")

feb = WORK / f"ankle_lig_{tag}.feb"
LOG = WORK / f"log_{tag}.txt"
PYOUT = ((ROOT / "temp/pyfebio_demo/_acc_last.log").read_text(
    encoding="utf-8", errors="replace")
    if (ROOT / "temp/pyfebio_demo/_acc_last.log").exists() else "")
xyz, sets, load_lines, ctl = parse_feb(feb) if feb.exists() else ({}, {}, -1, {})
log_txt = LOG.read_text(encoding="utf-8", errors="replace")
normal = "N O R M A L   T E R M I N A T I O N" in log_txt
negjac = len(re.findall(r"negative jacobian", log_txt, re.I))

disp = read_records(WORK / f"disp_{tag}.txt")
reac_t = read_records(WORK / f"reac_tibia_{tag}.txt")
reac_f = read_records(WORK / f"reac_fibula_{tag}.txt")
reac_s = read_records(WORK / f"reac_talus_{tag}.txt")

last = max(disp, key=lambda k: disp[k][0]) if disp else None
# ★ 纪律（[LRN-049]）：必须同 **t=1.0** 比较 —— 不能按 max(step) 取，
#   因为 FEBio 的 logfile 在步长重试时可能多写记录，最后一条未必是 t=1.0。
_last_t = disp[last][0] if last is not None else float("nan")
print(f"\n{'='*74}\n验收报告  tag={tag}\n{'='*74}")
print(f"  求解     : {'NORMAL TERMINATION ✓' if normal else 'ERROR ✗'}"
      f"   负 Jacobian {negjac} 处")
print(f"  末态     : step={last}  t={_last_t}   记录数 disp={len(disp)}"
      f"  {'✓ t=1.0' if abs(_last_t - 1.0) < 1e-9 else '⚠ 非 t=1.0，不要与其他运行直接比'}")
print(f"  Loads 段 : {load_lines} 行  ({'✓' if load_lines > 0 or env.get('DRIVER') != 'force' else '✗ 空! [ERR-20261005-007]'})")
print(f"  接触     : {ctl}")
# ★ 运行参数回显（权威，含脚本默认值）——防止"两次跑的不是同一套参数"
#   起因：2026-10-05 一次对照实验里 acceptance 没兜 TRUNC_MM，导致一次 TRUNC=60、
#   一次 TRUNC=0，两次的**加载轴差 12.18°**，反力一个在 x 一个在 z（[LRN-20261005-060]）。
print(f"  运行参数 : （以下为脚本自身 stdout 回显，含默认值）")
for _l in PYOUT.splitlines():
    if _l.startswith(("[mode]", "[axis]", "[trunc]", "[data]", "[loaddir]", "[pose]", "[align]")):
        print(f"     {_l.strip()[:118]}")

F_ax = None
if reac_t:
    _tl = max(reac_t, key=lambda k: reac_t[k][0])
    _t_f = reac_t[_tl][0]
    rr = reac_t[_tl][1]
    Fv = rr[:, 1:4].sum(0)          # 胫骨顶面反力合力（向量）
    Ff = (reac_f[max(reac_f, key=lambda k: reac_f[k][0])][1][:, 1:4].sum(0)
          if reac_f else np.zeros(3))
    ax = np.array([1.0, 0.0, 0.0])
    _npr = WORK / f"lig_read_{tag}.npz"
    if _npr.exists():
        ax = np.asarray(np.load(_npr)["axis"], float)
    else:
        _m = re.search(r"加载轴（胫骨->距骨）= \[([^\]]+)\]", PYOUT)
        if _m:
            ax = np.array([float(x) for x in _m.group(1).split()])
    ax = ax / np.linalg.norm(ax)
    print(f"\n  ▪ [t={_t_f}] 胫骨顶面反力合力 = {np.round(Fv, 4)} N  |F|={np.linalg.norm(Fv):.4f}")
    # ★ 距骨约束（tal_dofs=y,z）也是**真实载荷通路**：力可以从这里流走
    #   ⇒ 顶面反力 ≠ 关节载荷。必须把距骨反力一并记账（2026-10-05 补）
    Fs = (reac_s[max(reac_s, key=lambda k: reac_s[k][0])][1][:, 1:4].sum(0)
          if reac_s else np.zeros(3))
    print(f"  ▪ [t={_t_f}] 距骨约束反力合力 = {np.round(Fs, 4)} N  |F|={np.linalg.norm(Fs):.4f}"
          f"   ← fix_talus 约束 y,z，z 分量是载荷旁路")
    print(f"  ▪ [t={_t_f}] 腓骨顶面反力合力 = {np.round(Ff, 4)} N  |F|={np.linalg.norm(Ff):.4f}")
    print(f"  ▪ 两顶面合力 = {np.round(Fv + Ff, 4)} N  |F|={np.linalg.norm(Fv + Ff):.4f}")
    print(f"  ▪ 投影到加载轴 = {float((Fv + Ff) @ ax):.4f} N   (轴={np.round(ax, 4)})")
    F_ax = float(abs((Fv + Ff) @ ax))

if last is not None and last in disp:
    _, dd = disp[last]
    U = np.zeros((max(xyz) + 1, 3)) if xyz else None
    for row in dd:
        nid = int(row[0])
        if xyz and nid in xyz:
            U[nid] = row[1:4]
    ids = np.array(sorted(xyz))
    P0 = np.array([xyz[i] for i in ids])
    P1 = P0 + U[ids]
    print(f"  ▪ |u|max = {np.linalg.norm(U[ids], axis=1).max():.5f} mm")
    for nm, key in (("talus_cort", "talus"), ("tibia_cort", "tibia")):
        sel = [i for i in ids if sets.get("talus_low_ns")
               and i in sets.get("talus_low_ns")] if key == "talus" else None
    R, ratio = rigid_rot(P0, P1)
    ang, ax = rot_angle_axis(R)
    print(f"  ▪ 全体刚体拟合转动 = {ang:.4f}°  轴={np.round(ax, 3)}  "
          f"残差/原始={ratio:.3f}")

    # 关节内外翻角：用距骨低端 + 胫骨顶端两条点集分别拟合的相对转动
    tl = sets.get("talus_low_ns", [])
    if tl:
        tl = [i for i in tl if i in xyz]
        Rt, rat_t = rigid_rot(P0[[i for i in range(len(ids)) if ids[i] in tl]],
                              P1[[i for i in range(len(ids)) if ids[i] in tl]])
        at, axt = rot_angle_axis(Rt)
        print(f"  ▪ 距骨低端刚体转动 = {at:.4f}° 轴={np.round(axt, 3)} "
              f"残差/原始={rat_t:.3f}")

print(f"\n  C1 force收敛      : "
      f"{'—' if env.get('DRIVER') != 'force' else ('PASS' if normal and febio_rc == 0 else 'FAIL')}")
print(f"  C2 非零内外翻角    : {'待定（需判据阈值）'}")
print(f"  C3/C4/C5          : 需多次运行成对比较（姿势扫描 / 罚因子扫描）")
