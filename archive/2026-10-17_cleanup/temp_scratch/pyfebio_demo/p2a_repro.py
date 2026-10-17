"""P2'a 判定：在**同一批面**上逐式复现，找出 47.8 vs 555.5 的差从哪来。"""
import re
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"
txt = FEB.read_text(encoding="latin-1", errors="replace")

XYZ = {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,'
                     r'\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])

AXIS = np.array([0.0, 0.0, 1.0])


def faces_of(name):
    m = re.search(r'<Surface name="' + name + r'"[^>]*>(.*?)</Surface>', txt, re.S)
    if not m:
        return []
    out = []
    for mm in re.finditer(r'<(quad4|tri3) id="\d+">\s*([\d,]+)\s*</\1>', m.group(1)):
        out.append((mm.group(1), [int(x) for x in mm.group(2).split(",")]))
    return out


def area_proj_model(faces):
    """完全照抄 ankle_lig_feb.py:837 的 area_proj。"""
    a = 0.0
    for kind, n in faces:
        P = np.array([XYZ[c] for c in n], float)
        if len(P) == 3:
            a += 0.5 * abs(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), AXIS))
        else:
            v = np.cross(P[1] - P[0], P[3] - P[2])
            a += 0.5 * abs(np.dot(v, AXIS))
    return a


def area_proj_correct(faces):
    """正解：quad 拆 2 三角，各自 0.5|n·axis|。"""
    a = 0.0
    for kind, n in faces:
        P = np.array([XYZ[c] for c in n], float)
        tris = [(P[0], P[1], P[2]), (P[0], P[2], P[3])] if len(P) == 4 else \
               [(P[0], P[1], P[2])]
        for t in tris:
            cr = np.cross(t[1] - t[0], t[2] - t[0])
            a += 0.5 * abs(float(cr @ AXIS))
    return a


def real_area(faces):
    a = 0.0
    for kind, n in faces:
        P = np.array([XYZ[c] for c in n], float)
        tris = [(P[0], P[1], P[2]), (P[0], P[2], P[3])] if len(P) == 4 else \
               [(P[0], P[1], P[2])]
        for t in tris:
            a += 0.5 * np.linalg.norm(np.cross(t[1] - t[0], t[2] - t[0]))
    return a


def face_normal(P):
    if len(P) == 3:
        n = np.cross(P[1] - P[0], P[2] - P[0])
    else:
        n = np.cross(P[1] - P[0], P[3] - P[2])
    return n


for nm in ("tibia_top", "fibula_top"):
    fs = faces_of(nm)
    if not fs:
        print(f"✗ {nm} 未找到")
        continue
    q = sum(1 for k, _ in fs if k == "quad4")
    t3 = sum(1 for k, _ in fs if k == "tri3")
    print("=" * 88)
    print(f"{nm}: {len(fs)} 面（quad4 {q} / tri3 {t3}）")
    print("=" * 88)
    print(f"  模型 area_proj 公式      = {area_proj_model(fs):>10.2f} mm²")
    print(f"  正解（quad 拆2三角）投影  = {area_proj_correct(fs):>10.2f} mm²")
    print(f"  真实面积（未投影）        = {real_area(fs):>10.2f} mm²")
    # 法向·轴 分布
    cos = []
    for k, n in fs:
        P = np.array([XYZ[c] for c in n], float)
        nn = face_normal(P)
        L = np.linalg.norm(nn)
        if L > 0:
            cos.append(abs(float(nn @ AXIS)) / L)
    cos = np.array(cos)
    print(f"  |n·axis| 分布: min {cos.min():.3f} 中位 {np.median(cos):.3f} "
          f"max {cos.max():.3f}")
    print(f"    其中 >0.6 的面数 = {(cos > 0.6).sum()} / {len(cos)}"
          f"（extreme_faces 的 cos_min 阈值）")
    m6 = [(k, n) for (k, n), c in zip(fs, cos) if c > 0.6]
    print(f"    仅取 |n·axis|>0.6 时：模型式 {area_proj_model(m6):.2f} / "
          f"正解 {area_proj_correct(m6):.2f} mm²")

print("\n" + "=" * 88)
print("模型打印 = 胫骨 47.8 / 腓骨 12.7 mm²")
print("=" * 88)
