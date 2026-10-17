"""对比三个 OpenSim 模型：Maharaj2021 / Rajagopal FullBody / IITD 踝足。

重点回答：Maharaj2021 换掉 IITD 是否更好？(对 FE 几何 / 对多体 两个维度)
"""
import os
import sys
from pathlib import Path

VENV = Path(r"D:\Project\climbing_fall_analysis\.venv")
SP = VENV / "Lib" / "site-packages"
for d in (SP, SP / "opensim"):
    if d.exists():
        try:
            os.add_dll_directory(str(d))
        except Exception:
            pass

import opensim as osim

MODELS = [
    ("Maharaj2021_BothLegs",
     r"D:\Project\climbing_fall_analysis\model\FootwithinLowerlimbModel\Maharaj2021_BothLegs.osim"),
    ("Maharaj2021_V2",
     r"D:\Project\climbing_fall_analysis\model\FootwithinLowerlimbModel\Maharaj2021_BothLegs_V2.osim"),
    ("Maharaj_SampleData_SCALED",
     r"D:\Project\climbing_fall_analysis\model\FootwithinLowerlimbModel\SampleData_SCALED.osim"),
    ("Rajagopal FullBodyModel-4.0",
     r"D:\Project\climbing_fall_analysis\model\opensim\FullBodyModel-4.0\Rajagopal2015.osim"),
    ("IITD 踝足",
     r"D:\Project\climbing_fall_analysis\model\Model\IITD_Ankle_foot_model_NLSRclass.OSIM"),
]

print(f"{'model':28s} {'bodies':>6s} {'joints':>6s} {'coords':>6s} {'muscles':>7s} "
      f"{'mass(kg)':>8s} {'DOF':>4s}")
print("-" * 78)

loaded = {}
for label, p in MODELS:
    path = Path(p)
    if not path.exists():
        print(f"{label:28s}  MISSING")
        continue
    try:
        m = osim.Model(p)
        ms = m.getMuscles()
        cs = m.getCoordinateSet()
        bs = m.getBodySet()
        loaded[label] = m
        # OpenSim 4.6 py 绑定 getTotalMass() 签名异常，直接累加 body 质量
        mass = sum(bs.get(i).getMass() for i in range(bs.getSize()))
        print(f"{label:28s} {bs.getSize():6d} {m.getJointSet().getSize():6d} "
              f"{cs.getSize():6d} {ms.getSize():7d} {mass:8.2f} "
              f"{cs.getSize():4d}")
    except Exception as exc:
        msg = str(exc).split("\n")[0][:90]
        print(f"{label:28s}  LOAD ERROR {type(exc).__name__}: {msg}")

print()
print("=== 各模型 body 名称（判断解剖范围）===")
for label in ("Maharaj2021_BothLegs", "IITD 踝足"):
    m = loaded.get(label)
    if m is None:
        continue
    names = [m.getBodySet().get(i).getName() for i in range(m.getBodySet().getSize())]
    print(f"{label} ({len(names)}): {names}")

print()
print("=== 关节名（右腿）===")
for label in ("Maharaj2021_BothLegs", "IITD 踝足"):
    m = loaded.get(label)
    if m is None:
        continue
    js = m.getJointSet()
    names = [js.get(i).getName() for i in range(js.getSize())]
    print(f"{label}: {names}")

print()
print("=== geometry 文件解析测试 (Maharaj) ===")
mdir = Path(r"D:\Project\climbing_fall_analysis\model\FootwithinLowerlimbModel")
refs = ["sacrum.vtp", "pelvis.vtp", "l_pelvis.vtp", "femur_r.vtp", "tibia_r.vtp",
        "fibula_r.vtp", "pat.vtp", "FootModel_talus_r.stl", "FootModel_calcn_r.stl",
        "FootModel_midfoot_r.stl", "FootModel_forefoot_r.stl", "FootModel_digits_r.stl",
        "femur_l.vtp", "tibia_l.vtp", "l_fibula.vtp", "l_pat.vtp",
        "FootModel_talus_l.stl", "FootModel_calcn_l.stl", "FootModel_midfoot_l.stl",
        "FootModel_forefoot_l.stl", "FootModel_digits_l.stl"]
raj_geo = Path(r"D:\Project\climbing_fall_analysis\model\opensim\FullBodyModel-4.0\Geometry")
found_local, found_raj, missing = [], [], []
for r in refs:
    if (mdir / r).exists():
        found_local.append(r)
    elif (raj_geo / r).exists():
        found_raj.append(r)
    else:
        missing.append(r)

print(f"本目录直接存在   ({len(found_local):2d}/{len(refs)}): {found_local}")
print(f"仅 Rajagopal 目录有 ({len(found_raj):2d}): {found_raj}")
print(f"两处都缺         ({len(missing):2d}): {missing}")
