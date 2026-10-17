"""导出踝部韧带/肌肉的完整参数表 (含 PathPoint 附着点坐标) 到 JSON。"""

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GEOM = ROOT / "temp" / "ankle_geom"
osim = (GEOM / "IITD_Ankle_foot_model_NLSRclass.OSIM").read_text(encoding="utf-8", errors="replace")

i = osim.find("<ForceSet>")
j = osim.find("</ForceSet>")
fs = osim[i : j + 11]

OUT = {}


def get_num(blob, key, cast=float, default=None):
    m = re.search(rf"<{key}>([^<]*)</{key}>", blob)
    if not m:
        return default
    v = m.group(1).strip()
    try:
        return cast(v)
    except (ValueError, TypeError):
        return v


def get_pathpoints(blob):
    pts = []
    for pm in re.finditer(
        r"<PathPoint name=\"([^\"]*)\">\s*<location>([^<]*)</location>\s*<body>([^<]*)</body>",
        blob,
    ):
        loc = [float(x) for x in pm.group(2).split()]
        pts.append({"name": pm.group(1), "body": pm.group(3).strip(), "location_m": loc})
    return pts


def get_wraps(blob):
    return re.findall(r"<wrap_object>([^<]*)</wrap_object>", blob)


# ---- NonLinearLigament ----
nll = re.findall(r'<NonLinearLigament name="([^"]+)">(.*?)</NonLinearLigament>', fs, re.S)
lig = re.findall(r'<Ligament name="([^"]+)">(.*?)</Ligament>', fs, re.S)
mus = re.findall(r'<Millard2012EquilibriumMuscle name="([^"]+)">(.*?)</Millard2012EquilibriumMuscle>', fs, re.S)

# 论文对标的四条踝韧带 (命名对应关系)
WUKAI4 = {
    "pTTL_r": "下胫腓联合前韧带 (ATFL of tibiofibular syndesmosis)",
    "aTFL_r": "胫距前韧带 -> 论文'下胫腓联合前韧带'的近亲/论文称前下胫腓",
    "CFL_r": "跟腓韧带",
    "TNL_r": "胫舟韧带 (tibionavicular)",
    "TCL_r": "胫距后韧带",
    "TPS_r": "距腓前/后突韧带 (talofibular)",
    "FPS_r": "腓距后韧带",
    "aTTL_r": "胫距前韧带 (anterior talofibular) <- 论文最可能对应这条",
    "pTTL_r2": "下胫腓前韧带(别名)",
    "aTTL": "胫距前韧带",
}

out_lig = []
for name, blob in nll + lig:
    if not re.search(r"_(r|l)$", name):
        continue
    side = name.rsplit("_", 1)[-1]
    if side != "r":
        continue
    base = name[:-2]
    if not re.search(r"TTL|TFL|CFL|TNL|TCL|TPS|FPS|TF|PTF|ATFL|PTFL", base, re.I):
        continue
    kind = "NonLinearLigament" if (name, blob) in nll else "Ligament"
    out_lig.append(
        {
            "name": name,
            "type": kind,
            "resting_length_m": get_num(blob, "resting_length"),
            "pcsa_force_N": get_num(blob, "pcsa_force"),
            "eta_a": get_num(blob, "eta_a"),
            "alpha_a": get_num(blob, "alphaa_a"),
            "beta_a": get_num(blob, "beta_a"),
            "k_a": get_num(blob, "ka_a"),
            "area_a": get_num(blob, "areaa_a"),
            "max_strain": get_num(blob, "max_strain"),
            "isDisabled": get_num(blob, "isDisabled", cast=str),
            "attachment_points": get_pathpoints(blob),
            "wrap_objects": get_wraps(blob),
            "wukai_mapping": WUKAI4.get(name, ""),
        }
    )

out_mus = []
for name, blob in mus:
    if not name.endswith("_r"):
        continue
    if not re.search(
        r"ankle|tibial|soleus|gastroc|perone|fibul|plantar|dorsi|flex|extensor|evert",
        name,
        re.I,
    ):
        continue
    out_mus.append(
        {
            "name": name,
            "max_isometric_force_N": get_num(blob, "max_isometric_force"),
            "optimal_fiber_length_m": get_num(blob, "optimal_fiber_length"),
            "tendon_slack_length_m": get_num(blob, "tendon_slack_length"),
            "pennation_deg": get_num(blob, "pennation_angle_at_optimal"),
            "fiber_damping": get_num(blob, "fiber_damping"),
            "attachment_points": get_pathpoints(blob),
        }
    )

# 踝/距下关节运动学
kin = {}
for jn in ["Ankle_r", "Subtalar_r", "Chopart_r"]:
    m = re.search(rf'<CustomJoint name="{jn}">(.*?)</CustomJoint>', osim, re.S)
    if not m:
        continue
    blob = m.group(1)
    axes = re.findall(
        r'<TransformAxis name="(\w+)">\s*<coordinates>([^<]*)</coordinates>\s*<axis>([^<]*)</axis>',
        blob,
    )
    coord = re.search(r'<Coordinate name="([^"]+)"', blob)
    rng = re.search(r"<range>([^<]*)</range>", blob)
    kin[jn] = {
        "coordinate": coord.group(1) if coord else None,
        "range_rad": [float(x) for x in rng.group(1).split()] if rng else None,
        "axes": [{"name": a, "coords": c.strip(), "axis": [float(x) for x in v.split()]} for a, c, v in axes],
    }

OUT = {
    "source": "IITD_Ankle_foot_model_NLSRclass.OSIM (inside paper/ankleModel-latest.zip)",
    "units": {"length": "m", "force": "N"},
    "ankle_ligaments_right": out_lig,
    "ankle_muscles_right": out_mus,
    "joint_kinematics": kin,
}

dst = ROOT / "temp" / "ankle_geom" / "osim_params.json"
dst.write_text(json.dumps(OUT, indent=2, ensure_ascii=False), encoding="utf-8")

print(f"写出 {dst} ({dst.stat().st_size:,} bytes)")
print(f"  踝韧带 {len(out_lig)} 条, 踝肌肉 {len(out_mus)} 条, 关节 {len(kin)} 个")
print("\n=== 韧带参数表 ===")
print(f"{'name':<10}{'L0(m)':>9}{'F0(N)':>9}{'eta':>7}{'alpha':>8}{'beta':>7}{'k_a':>9}{'maxε':>7}{'Pts':>5}")
for d in out_lig:
    print(
        f"  {d['name']:<8}{d['resting_length_m'] or 0:>9.4f}{d['pcsa_force_N'] or 0:>9.1f}"
        f"{d['eta_a'] or 0:>7.0f}{d['alpha_a'] or 0:>8.3f}{d['beta_a'] or 0:>7.1f}"
        f"{d['k_a'] or 0:>9.0f}{d['max_strain'] or 0:>7.2f}{len(d['attachment_points']):>5}"
    )
print("\n=== 肌肉 ===")
print(f"{'name':<24}{'Fmax(N)':>10}{'Lf(m)':>10}{'Lt(m)':>10}{'pen°':>7}")
for d in out_mus:
    print(
        f"  {d['name']:<22}{d['max_isometric_force_N'] or 0:>10.1f}"
        f"{d['optimal_fiber_length_m'] or 0:>10.4f}{d['tendon_slack_length_m'] or 0:>10.4f}"
        f"{d['pennation_deg'] or 0:>7.2f}"
    )
print("\n=== 关节运动学 ===")
for jn, d in kin.items():
    print(f"  {jn:<14} coord={d['coordinate']:<18} range={d['range_rad']}")
    for a in d["axes"][:3]:
        print(f"      {a['name']:<12} [{a['coords']}] axis={a['axis']}")
