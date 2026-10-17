"""P2③ 成因：坏单元处的厚度 vs 全局，判断是"薄"还是"抬回折叠"。"""
import re
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"
txt = FEB.read_text(encoding="latin-1", errors="replace")
GP2 = 1 / np.sqrt(3)


def dNw(r, s, t):
    return np.array([-(1-t), -(1-t), -(1-r-s), (1-t), 0, -r, 0, (1-t), -s,
                     -t, -t, (1-r-s), t, 0, r, 0, t, s]).reshape(6, 3)


PTS = [(1/6, 1/6, .5), (2/3, 1/6, .5), (1/6, 2/3, .5),
       (1/3, 1/3, 0.), (1/3, 1/3, 1.), (1/3, 1/3, .5),
       (0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
       (0., 0., 1.), (1., 0., 1.), (0., 1., 1.)]

XYZ, ELEMS = {}, {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,'
                     r'\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])
for m in re.finditer(r'<Elements\s+type="(\w+)"[^>]*?name="([^"]*)"[^>]*>'
                     r'(.*?)</Elements>', txt, re.S):
    kind, nm, body = m.group(1), m.group(2), m.group(3)
    conn = [[int(x) for x in c.split(",")] for c in
            re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*</elem>', body)]
    ELEMS[nm] = (kind, conn)


def metrics(P):
    w = min(float(np.linalg.det(dNw(*p).T @ P)) for p in PTS)
    t = float(np.linalg.norm(P[[3, 4, 5]].mean(0) - P[[0, 1, 2]].mean(0)))
    # 底/顶三角的"扁平度"：底三角最长边 vs 高
    def tri_shape(Q):
        L = [np.linalg.norm(Q[1]-Q[0]), np.linalg.norm(Q[2]-Q[1]),
             np.linalg.norm(Q[0]-Q[2])]
        s = sum(L) / 2
        A = float(np.linalg.norm(np.cross(Q[1]-Q[0], Q[2]-Q[0]))) / 2
        h = 2 * A / max(L) if max(L) > 0 else 0
        mn = 4 * A / (sum(L)**2) * 4 * np.sqrt(3) if sum(L) > 0 else 0
        return A, h, min(L) / max(L) if max(L) > 0 else 0
    A, h, r = tri_shape(P[[0, 1, 2]])
    # 底/顶三角平面的夹角（扭角）
    n1 = np.cross(P[1]-P[0], P[2]-P[0])
    n2 = np.cross(P[4]-P[3], P[5]-P[3])
    ang = np.degrees(np.arccos(np.clip(
        abs(n1 @ n2) / (np.linalg.norm(n1) * np.linalg.norm(n2) + 1e-12), -1, 1)))
    return w, t, A, min(L_h := [np.linalg.norm(P[1]-P[0]),
                                np.linalg.norm(P[2]-P[1]),
                                np.linalg.norm(P[0]-P[2])]) / max(L_h), ang


print("=" * 104)
print("坏单元 vs 全体：厚度 / 底三角面积 / 最瘦边比 / 上下三角扭角")
print("=" * 104)
print(f"{'域/单元':<34}{'worst-detJ':>13}{'厚度':>9}{'底三角面积':>11}"
      f"{'边比min':>9}{'扭角(°)':>9}")
allm = []
for nm in ("tibiotalar_talus__penta6", "tibiotalar_tibia__penta6"):
    kind, conn = ELEMS[nm]
    rows = []
    for idx, c in enumerate(conn):
        P = np.array([XYZ[n] for n in c], float)
        rows.append((idx, *metrics(P)))
    allm += [(nm, r) for r in rows]
    bad = [r for r in rows if r[1] <= 1e-2]
    good = [r for r in rows if r[1] > 1e-2]
    print(f"\n--- {nm}（坏 {len(bad)} / 共 {len(rows)}）---")
    for idx, w, t, A, ratio, ang in sorted(bad, key=lambda x: x[1]):
        print(f"  #{idx:<4}                          {w:>13.4e}{t:>9.4f}"
              f"{A:>11.4f}{ratio:>9.3f}{ang:>9.2f}")
    if good:
        gt = [r[2] for r in good]
        gA = [r[3] for r in good]
        gr = [r[4] for r in good]
        ga = [r[5] for r in good]
        print(f"  基准（好的 {len(good)} 个）中位：厚度 {np.median(gt):.4f}  "
              f"面积 {np.median(gA):.4f}  边比 {np.median(gr):.3f}  "
              f"扭角 {np.median(ga):.2f}")

print("\n" + "=" * 104)
print("👉 判定")
print("=" * 104)
bads = [(nm, r) for nm, r in allm if r[1] <= 1e-2]
goods = [(nm, r) for nm, r in allm if r[1] > 1e-2]
print(f"  坏单元厚度 中位 {np.median([r[2] for _, r in bads]):.4f} mm  vs  "
      f"好的 {np.median([r[2] for _, r in goods]):.4f} mm")
print(f"  坏单元面积 中位 {np.median([r[3] for _, r in bads]):.4f} mm² vs  "
      f"好的 {np.median([r[3] for _, r in goods]):.4f} mm²")
print(f"  坏单元边比 中位 {np.median([r[4] for _, r in bads]):.3f}   vs  "
      f"好的 {np.median([r[4] for _, r in goods]):.3f}")
print(f"  坏单元扭角 中位 {np.median([r[5] for _, r in bads]):.2f}°   vs  "
      f"好的 {np.median([r[5] for _, r in goods]):.2f}°")
print("\n  若坏单元面积/边比显著更差 ⇒ 是「瘦长/畸形底三角」，来自参数化网格")
print("  若扭角显著更大 ⇒ 是「上下三角扭转」，来自抬回时方向不一致")
