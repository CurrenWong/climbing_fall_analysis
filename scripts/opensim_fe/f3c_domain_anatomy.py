# -*- coding: utf-8 -*-
"""F-3′c · 域-解剖对应核对（正式版，零解算）。

结论（2026-10-14，双解剖锚点 + 变换复现，证据不可翻转）：
- top_load（+Z 端 30%）= 髁端（远端，膝关节侧）；bottom_fix（-Z 端 30%）= 股骨头端（近端）
- partition_ends_shaft 的 ends = 头端25%（头+颈+大转子）∪ 髁端25%（两髁）合并带
  ⇒ A2-股骨需换解剖域（近端/骨干/远端三分）并重做判定（并入 FE BC 大件 ②）

证据链：
1. .feb 节点 id 已重编号 ⇒ 按坐标变换复现 align_pca（bone_feb._align_pca_rotate verbatim）
2. 旋转后股骨点集 ↔ .feb 节点集最近邻 RMS=0.0000 mm ⇒ 同一几何
3. 髋骨/胫骨锚点绕股骨质心同 R 旋转 ⇒ .feb 系内判定
4. THUMS 原始系独立交叉复核（髋骨/胫骨双锚点一致，差距 4×+）

运行：.venv\\Scripts\\python.exe scripts/opensim_fe/f3c_domain_anatomy.py
"""
from pathlib import Path
import re
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
FEB = ROOT / "temp" / "opensim_fe" / "fe_bc_perheight" / "femur_r" / "h50" / "femur_r.feb"
K = ROOT / "model" / "AM50_V71_Occupant" / "THUMS_model" / "THUMS_AM50_V71_Occupant_202411.k"

FEMUR_PIDS = [81000000, 81000001, 81000100]
HIP_PIDS = [83500200, 83500201, 83500202, 83500203, 83500204, 83500205,
            83500206, 83500207, 83500208, 83500209]
TIBIA_PIDS = [81000600, 81000601, 81000700]


def fields8(line):
    """LS-DYNA 定宽 8 字符/字段，数字可无空格连写。"""
    return [line[i:i + 8].strip() for i in range(0, len(line), 8)]


def align_pca(nodes):
    """verbatim: bone_feb._align_pca_rotate。返回 (旋转后, R, 质心)。"""
    nodes = np.array(nodes, dtype=np.float64, copy=True)
    c = nodes.mean(0)
    X = nodes - c
    w, v = np.linalg.eigh(X.T @ X)
    ax = v[:, -1]
    if ax[2] < 0:
        ax = -ax
    z = np.array([0.0, 0.0, 1.0])
    vv = np.cross(ax, z)
    s = np.linalg.norm(vv)
    if s < 1e-9:
        R = np.eye(3)
    else:
        vv = vv / s
        ang = np.arctan2(s, float(np.dot(ax, z)))
        K = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
        R = np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * (K @ K)
    return (X @ R.T) + c, R, c


def parse_feb(path: Path):
    text = path.read_text(encoding="ISO-8859-1")
    m = re.search(r"<Nodes name=\"([^\"]+)\">(.*?)</Nodes>", text, re.S)
    nodes = []
    for mm in re.finditer(r'<node id="(\d+)">([^<]+)</node>', m.group(2)):
        nodes.append([float(x) for x in mm.group(2).split(",")])
    ns = {}
    for mm in re.finditer(r'<NodeSet name="([^"]+)">([^<]+)</NodeSet>', text):
        ns[mm.group(1)] = [int(x) for x in mm.group(2).split(",")]
    return np.array(nodes), ns


def load_thums_parts():
    with K.open("r", encoding="ISO-8859-1", errors="replace") as f:
        lines = f.readlines()
    node_map = {}
    in_node = False
    for s in lines:
        t = s.strip()
        if t.startswith("*"):
            in_node = t.startswith("*NODE")
            continue
        if in_node and t and not t.startswith("$"):
            f = t.split()
            if len(f) >= 4:
                node_map[int(f[0])] = (float(f[1]), float(f[2]), float(f[3]))
    pid_nodes = {}
    in_es = False
    for s in lines:
        t = s.strip()
        if t.startswith("*ELEMENT_SOLID"):
            in_es = True
            continue
        if t.startswith("*") and not t.startswith("*ELEMENT_SOLID"):
            in_es = False
            continue
        if in_es and t and not t.startswith("$"):
            f = fields8(s)
            if len(f) >= 6:
                try:
                    pid = int(f[1])
                except ValueError:
                    continue
                pid_nodes.setdefault(pid, []).extend(int(x) for x in f[2:] if x)
    return node_map, pid_nodes


def part_pts(pids, node_map, pid_nodes):
    ids = set()
    for p in pids:
        ids.update(pid_nodes.get(p, []))
    return np.array([node_map[i] for i in ids if i in node_map])


def main() -> int:
    node_map, pid_nodes = load_thums_parts()
    femur_orig = part_pts(FEMUR_PIDS, node_map, pid_nodes)
    hip_orig = part_pts(HIP_PIDS, node_map, pid_nodes)
    tib_orig = part_pts(TIBIA_PIDS, node_map, pid_nodes)
    feb_pts, feb_ns = parse_feb(FEB)
    print(f"[f3c] THUMS 股骨 {len(femur_orig)} / 髋骨 {len(hip_orig)} / 胫骨 {len(tib_orig)}")
    print(f"[f3c] .feb 节点 {len(feb_pts)}  bbox span=[{np.ptp(feb_pts, 0).round(2)}]")

    rot_pts, R, c = align_pca(femur_orig)
    d, _ = cKDTree(feb_pts).query(rot_pts, k=1)
    rms = float(np.sqrt((d ** 2).mean()))
    print(f"[f3c] 变换复现一致性：旋转后股骨 ↔ .feb 最近邻 RMS={rms:.4f} mm")

    def rot_about(pts):
        return (pts - c) @ R.T + c

    hip_c = rot_about(hip_orig).mean(0)
    tz = rot_about(tib_orig)[:, 2]
    tib_top_c = rot_about(tib_orig)[tz >= tz.max() - 0.10 * (tz.max() - tz.min())].mean(0)
    print(f"[f3c] .feb 系锚点 髋骨质心=({hip_c[0]:.1f},{hip_c[1]:.1f},{hip_c[2]:.1f})  "
          f"胫骨上端=({tib_top_c[0]:.1f},{tib_top_c[1]:.1f},{tib_top_c[2]:.1f})")

    verdicts = []
    for name in ("top_load", "bottom_fix", "lat_support_0"):
        ids = feb_ns.get(name, [])
        P = feb_pts[[int(i) - 1 for i in ids]]
        ctr = P.mean(0)
        d_hip = float(np.linalg.norm(ctr - hip_c))
        d_tib = float(np.linalg.norm(ctr - tib_top_c))
        v = "股骨头端(近端)" if d_hip < d_tib else "髁端(远端)"
        verdicts.append((name, v))
        print(f"[f3c] {name:14s} n={len(ids):5d}  .feb质心=({ctr[0]:7.1f},{ctr[1]:7.1f},{ctr[2]:7.1f})  "
              f"距髋骨 {d_hip:6.1f} mm  距胫骨上端 {d_tib:6.1f} mm  ⇒ {v}")

    # 原始系独立交叉复核
    z = femur_orig[:, 2]
    zlo, zhi = z.min(), z.max()
    band = 0.30 * (zhi - zlo)
    tc = femur_orig[z >= zhi - band].mean(0)
    bc = femur_orig[z <= zlo + band].mean(0)
    hip_c0 = hip_orig.mean(0)
    tib_pts0 = tib_orig
    tz0 = tib_pts0[:, 2]
    tib_top0 = tib_pts0[tz0 >= tz0.max() - 0.10 * (tz0.max() - tz0.min())].mean(0)
    print(f"\n[f3c] 原始系交叉复核（THUMS 全局坐标，独立于 R）：")
    print(f"[f3c]   股骨 -Z 端30% 质心距髋骨 {np.linalg.norm(bc - hip_c0):.2f} mm（+Z 端 {np.linalg.norm(tc - hip_c0):.2f} mm）")
    print(f"[f3c]   胫骨上端 距股骨 +Z 端 {np.linalg.norm(tib_top0 - tc):.2f} mm（-Z 端 {np.linalg.norm(tib_top0 - bc):.2f} mm）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
