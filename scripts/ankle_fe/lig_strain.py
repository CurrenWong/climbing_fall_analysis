"""权威韧带应变读数（从几何直接算，不依赖任何旧读数）。

为什么需要它 —— 旧 §4「逐韧带应变」有三个问题：
 ① **读数来源不明**：`ankle_lig_feb.py` 根本不计算应变，`LOG_DATA` 也只写位移/反力
    ⇒ 那批数字没有可追溯的产生路径（推测是已删除的临时脚本）。
 ② **与几何矛盾**：纯几何量出「内翻(+IE) 时 ATFL 长度 +1.75 mm ⇒ 应变 +4.7%」，
    而旧读数说 ATFL「几乎不动，0.17~0.20%」⇒ 差 ~24×。
 ③ **标签被循环论证污染**：旧结论把 `POSE_IE=−3` 当内翻，而这个判断源自韧带数据本身，
    再用它宣布韧带数据"物理正确"。`check_c2.py` 的两个**解剖锚定**判据证明
    **+IE = 内翻** ⇒ 标签其实没反，是**读数**反了。

本脚本给三种应变（都要，因为它们回答不同问题）：
  (a) **姿势预应变** = (L_posed − L_neutral)/L_neutral
      ⇒ 姿势本身把韧带拉了多少。**这是内翻扭伤的机制所在**，不能省。
  (b) **载荷增量应变** = (L_loaded − L_posed)/L_posed
      ⇒ 等价于 FEBio 在 t=1 报的元素应变（其参考构型 = posed 网格）。
  (c) **总应变** = (L_loaded − L_neutral)/L_neutral
      ⇒ **应对它设损伤阈值**，因为韧带的生理参考态是中性姿势。
      ⚠️ 只报 (b) 会把姿势预应变整段丢掉 —— 那正是要判的东西。

用法：
    python lig_strain.py NEUTRAL=c2/IE_0.feb TAG=WK
    python lig_strain.py NEUTRAL=c2/IE_0.feb TAG=WK POSED=c2/IE_500.feb
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
C2D = WORK / "c2"

LAT_KEYS = ("ATFL", "CFL", "PTFL")
MED_KEYS = ("MD1", "MD2", "MD3", "MD4")


def parse_feb(p: Path):
    t = p.read_text(encoding="latin-1", errors="replace")
    xyz = {}
    for m in re.finditer(r'<node id="(\d+)">([-\d.eE+,]+)</node>', t):
        xyz[int(m.group(1))] = np.array(
            [float(x) for x in m.group(2).split(",")], float)
    return xyz


def read_disp_last(p: Path):
    """读 FEBio node_data txt，返回 {t: (time, disp数组)}，取 t→1 最大者。"""
    if not p.exists():
        return None, None
    cur, best = None, None
    rows, tval = [], None
    for l in p.read_text(encoding="utf-8", errors="replace").split("\n"):
        if l.startswith("*Step"):
            if rows and tval is not None:
                arr = np.array(rows, float)
                if best is None or abs(tval - 1.0) < abs(best[0] - 1.0):
                    best = (tval, arr)
            rows, tval = [], None
        elif l.startswith("*Time"):
            tval = float(l.split("=")[1])
        elif l.strip() and not l.startswith("*"):
            try:
                rows.append([float(x) for x in l.split(",")])
            except ValueError:
                pass
    if rows and tval is not None:
        arr = np.array(rows, float)
        if best is None or abs(tval - 1.0) < abs(best[0] - 1.0):
            best = (tval, arr)
    if best is None:
        return None, None
    t0, arr = best
    return t0, {int(r[0]): r[1:4] for r in arr}


def lig_lengths(xyz: dict, names, pairs, npair):
    """每条韧带每根杆的长度（mm）。"""
    out = {}
    for i, nm in enumerate(names):
        ls = []
        for j in range(int(npair[i])):
            a, b = int(pairs[i, j, 0]), int(pairs[i, j, 1])
            if a < 0 or b < 0 or a not in xyz or b not in xyz:
                continue
            ls.append(float(np.linalg.norm(xyz[b] - xyz[a])))
        out[nm] = np.array(ls)
    return out


def main():
    a = {x.split("=", 1)[0]: x.split("=", 1)[1] for x in sys.argv[1:] if "=" in x}
    neutral_f = C2D / "IE_0.feb"
    if "NEUTRAL" in a:
        neutral_f = WORK.parent / a["NEUTRAL"] if not a["NEUTRAL"].startswith(
            ("c2/", "c2\\")) else WORK / a["NEUTRAL"]
    tag = a.get("TAG", "WK")
    posed_f = Path(a["POSED"]) if "POSED" in a else None

    z = np.load(Path(a["LIGNPZ"]) if "LIGNPZ" in a
                else WORK / f"lig_read_{tag}.npz", allow_pickle=True)
    names = [str(x) for x in z["lig_names"]]
    pairs, npair = z["lig_pairs"], z["lig_npair"]

    xyz_n = parse_feb(neutral_f)
    xyz_p = parse_feb(posed_f) if posed_f else xyz_n

    # 变形坐标 = posed 参考 + 位移
    disp_f = Path(a["DISP"]) if "DISP" in a else WORK / f"disp_{tag}.txt"
    t0, disp = read_disp_last(disp_f)
    print(f"中性网格   : {neutral_f.name}  {len(xyz_n)} 节点")
    print(f"姿势参考   : {'（=中性）' if posed_f is None else posed_f.name}"
          f"  {len(xyz_p)} 节点")
    print(f"载荷位移   : {disp_f.name}  取到 t={t0}"
          f"  {len(disp) if disp else 0} 条记录")
    if disp:
        xyz_l = {i: (xyz_p[i] + disp[i] if i in disp else xyz_p[i])
                 for i in xyz_p}
        # 检查覆盖率
        n_miss = sum(1 for i in xyz_p if i not in disp)
        print(f"位移覆盖   : {len(xyz_p)-n_miss}/{len(xyz_p)}"
              f"（缺 {n_miss} 个节点按零位移处理）")

    L_n = lig_lengths(xyz_n, names, pairs, npair)
    L_p = lig_lengths(xyz_p, names, pairs, npair)
    L_l = lig_lengths(xyz_l, names, pairs, npair) if disp else None

    print(f"\n{'韧带':<30} {'n杆':>4} {'L中性':>7} {'L姿势':>7} "
          f"{'L载荷':>7} {'预应变%':>9} {'增量%':>8} {'总应变%':>9}")
    print("─" * 92)
    res = {}
    for nm in names:
        n0, p0 = L_n[nm], L_p[nm]
        if len(n0) == 0:
            continue
        # ★ 用**逐杆**应变再取统计（长度不同不能先平均长度）
        pre = (p0 - n0) / n0 * 100.0
        if L_l is not None and len(L_l[nm]) == len(p0):
            inc = (L_l[nm] - p0) / p0 * 100.0
            tot = (L_l[nm] - n0) / n0 * 100.0
            inc_s, tot_s = f"{inc.mean():>8.3f}", f"{tot.mean():>9.3f}"
        else:
            inc_s, tot_s = f"{'—':>8}", f"{'—':>9}"
        print(f"{nm:<30} {len(n0):>4} {n0.mean():>7.2f} {p0.mean():>7.2f} "
              f"{(L_l[nm].mean() if L_l is not None else float('nan')):>7.2f} "
              f"{pre.mean():>9.3f} {inc_s} {tot_s}")
        res[nm] = dict(pre=pre, inc=(L_l[nm] - p0) / p0 * 100.0
                       if L_l is not None else None)

    # 分侧汇总
    print(f"\n{'─'*92}")
    for lbl, keys in (("外侧（连腓骨）", LAT_KEYS), ("内侧（三角韧带）", MED_KEYS)):
        sel = [n for n in names if any(k in n for k in keys)]
        if not sel:
            continue
        pre = np.nanmean([res[n]["pre"].mean() for n in sel if n in res])
        tot = np.nanmean([res[n]["inc"].mean() for n in sel
                          if n in res and res[n]["inc"] is not None])
        print(f"  {lbl:<16} 预应变 {pre:>+7.3f}%   载荷增量 {tot:>+7.3f}%")
    print("\n★ 判据提示：内翻(+IE) 应满足『外侧预应变 > 0 且 > 内侧』；外翻相反。")
    print("★ 阈值应用 (c) 总应变（相对中性姿势），不是 (b) 载荷增量。")


if __name__ == "__main__":
    main()
