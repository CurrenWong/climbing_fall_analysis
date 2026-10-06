"""C2-C5: 解析 THUMS 骨材料，建立「骨部件 -> 材料 -> 参数」对照表。

★ 实测：THUMS 的骨相关材料用了 5 种不同的卡，字段布局各不相同 ——
  按卡类型分派，并用 **ν ∈ (0, 0.5)** 作硬闸门（ν=4.0 这类不可能值必须报错，
  这个闸门正是抓出「把 G0 当成 ν」那个 bug 的东西）。

  MAT_ELASTIC (001)                     : MID RO **E** **PR** DA DB K
  MAT_PIECEWISE_LINEAR_PLASTICITY (024) : MID RO **E** **PR** SIGY ETAN FAIL TDEL
                                          + card3/4 = 塑性应变/应力曲线
  MAT_PLASTICITY_WITH_DAMAGE (081)      : MID RO **E** **PR** SIGY ETAN ...
  MAT_DAMAGE_2 (105)                    : MID RO **E** **PR** SIGY ETAN EPPF ...
  MAT_VISCOELASTIC (006)                : MID RO **BULK** **G0** GI BETA
                                          -> E = 9KG/(3K+G), ν = (3K-2G)/(2(3K+G))

★ 单位（脚本内硬校验，非假设）：mm-s-tonne-N-MPa
   ρ = 2.0E-9 tonne/mm³ = 2000 kg/m³（皮质骨）；E = 18000 MPa = 18 GPa
   ⇒ **E 已是 MPa，与 FEBio mm-N-MPa 一致，零换算**
"""
from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
MAIN = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"
MATS = {
    "bone_fracture": "AM50_V71_Occupant/THUMS_model/mat_bone_AM50_V71_Occ_fracture.k",
    "bone_no_fracture": "AM50_V71_Occupant/THUMS_model/mat_bone_AM50_V71_Occ_no_fracture.k",
}
OUT = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
RE_KW = re.compile(r"^\*([A-Z_0-9]+)")

# 卡类型 -> 字段映射。("E", iE, iNu) = 弹性模量/泊松比在 card1 的字段序号
ELASTIC_LIKE = {
    "MAT_ELASTIC": (2, 3),
    "MAT_PIECEWISE_LINEAR_PLASTICITY": (2, 3),
    "MAT_PLASTICITY_WITH_DAMAGE": (2, 3),
    "MAT_DAMAGE_2": (2, 3),
    "MAT_ISOTROPIC_ELASTIC_PLASTIC": (2, 3),
}


def f10(s: str, i: int):
    seg = s[i * 10 : i * 10 + 10].strip()
    if not seg:
        return None
    try:
        return float(seg)
    except ValueError:
        return None


def parse_cards(txt: str):
    out, kw, data = [], None, []
    for raw in txt.split("\n"):
        s = raw.rstrip("\r")
        st = s.strip()
        if not st or st[0] == "$":
            continue
        if st[0] == "*":
            if kw:
                out.append((kw, data))
            m = RE_KW.match(st.upper())
            kw, data = (m.group(1) if m else ""), []
            continue
        if kw:
            data.append(s)
    if kw:
        out.append((kw, data))
    return out


def build_entry(kw: str, data: list[str]):
    """按卡类型解析一条材料。无法处理返回 None。"""
    if not data:
        return None
    mid = f10(data[0], 0)
    if mid is None:
        return None
    e: dict = {"card": kw, "MID": int(mid), "RO": f10(data[0], 1)}

    if kw == "MAT_VISCOELASTIC":
        bulk, g0, gi, beta = (f10(data[0], 2), f10(data[0], 3),
                              f10(data[0], 4), f10(data[0], 5))
        e.update({"BULK": bulk, "G0": g0, "GI": gi, "BETA": beta})
        if bulk and g0 and (3 * bulk + g0) != 0:
            e["E"] = 9 * bulk * g0 / (3 * bulk + g0)
            e["PR"] = (3 * bulk - 2 * g0) / (2 * (3 * bulk + g0))
            e["_derived"] = "E,nu from BULK/G0"
    elif kw in ELASTIC_LIKE:
        iE, iNu = ELASTIC_LIKE[kw]
        e.update({"E": f10(data[0], iE), "PR": f10(data[0], iNu),
                  "SIGY": f10(data[0], 4), "ETAN": f10(data[0], 5),
                  "FAIL": f10(data[0], 6)})
        if len(data) >= 4:
            pts = [(f10(data[2], i), f10(data[3], i)) for i in range(8)]
            pts = [(a, b) for a, b in pts if a is not None and b is not None and b]
            e["plastic_curve"] = pts
            e["max_stress_MPa"] = max((b for _, b in pts), default=None)
    else:
        e["raw_card1"] = data[0]
        e["_unhandled_card"] = True
    return e


# ── 1. 材料库：骨材料文件 + 主 .k 内联 ────────────────────────
matlib: dict[str, dict[int, dict]] = {}
for tag, inner in MATS.items():
    with zipfile.ZipFile(ZIP) as z:
        txt = z.read(inner).decode("latin-1")
    lib = {}
    for kw, data in parse_cards(txt):
        ent = build_entry(kw, data)
        if ent:
            lib[ent["MID"]] = ent
    matlib[tag] = lib
    print(f"[{tag}] {len(lib)} 条材料")

# ── 2. 主 .k：*PART(pids/secid/mid/title) + 内联材料 ──────────
print("\n扫描主 .k …")
part_info: dict[int, dict] = {}
main_cards: list[tuple[str, list[str]]] = []
pending = None
cur_kw = None
cur_data: list[str] = []
buf = ""


def flush_main():
    global cur_data
    if cur_kw and cur_data:
        main_cards.append((cur_kw, list(cur_data)))
    cur_data = []


with zipfile.ZipFile(ZIP) as z, z.open(MAIN) as f:
    while True:
        chunk = f.read(8 * 1024 * 1024)
        if not chunk:
            break
        buf += chunk.decode("latin-1")
        lines = buf.split("\n")
        buf = lines.pop()
        for raw in lines:
            s = raw.rstrip("\r")
            st = s.strip()
            if not st or st[0] == "$":
                continue
            if st[0] == "*":
                flush_main()
                m = RE_KW.match(st.upper())
                cur_kw = m.group(1) if m else None
                if cur_kw == "PART":
                    pending = None
                continue
            if cur_kw == "PART":
                if pending is None:
                    pending = st
                else:
                    pid, secid, mid = f10(st, 0), f10(st, 1), f10(st, 2)
                    if pid is not None:
                        part_info[int(pid)] = {
                            "title": pending,
                            "secid": int(secid) if secid is not None else None,
                            "mid": int(mid) if mid is not None else None,
                        }
                    pending = None
            elif cur_kw and cur_kw.startswith("MAT"):
                cur_data.append(s)
    flush_main()

main_lib: dict[int, dict] = {}
for kw, data in main_cards:
    ent = build_entry(kw, data)
    if ent:
        main_lib[ent["MID"]] = ent
matlib["main_k"] = main_lib
print(f"  *PART {len(part_info):,} 个；主 .k 材料卡 {len(main_cards):,} 条 -> "
      f"解析 {len(main_lib)} 个")
from collections import Counter
print(f"  主 .k 材料卡类型分布: {dict(Counter(k for k, _ in main_cards))}")


def lookup(mid: int):
    for tag in ("bone_fracture", "bone_no_fracture", "main_k"):
        if mid in matlib.get(tag, {}):
            d = dict(matlib[tag][mid])
            d["_source"] = tag
            return d
    return None


# ── 3. 目标骨部件 ─────────────────────────────────────────────
BONES = {
    81000900: "R_FIBULA_CORT", 81000800: "right_fibula_end_spon",
    81000801: "right_fibula_center_spon",
    81000700: "R_TIBIA_CORT", 81000600: "right_tibia_end_spon",
    81000601: "right_tibia_center_spon",
    81001100: "R_TALUS_CORT", 81001000: "R_TALUS_SPON",
    81001300: "R_CALCANEUS_CORT", 81001200: "R_CALCANEUS_SPON",
    81000100: "R_FEMUR_CORT", 81000000: "right_femur_end_spon",
    81000001: "right_femur_center_spon",
    81000300: "R_PATELLA_CORT", 81000200: "R_PATELLA_SPON",
    81001401: "R_NAVICULAR_CORT", 81001501: "R_CUBOID_CORT",
    81001901: "R_1ST_METATARSAL_CORT", 81001900: "R_1ST_METATARSAL_SPON",
}

print("\n=== 骨部件 -> 材料 -> 参数 ===")
print(f"{'pid':>10s} {'mid':>10s} {'source':>16s} {'card':>32s} "
      f"{'ρ(kg/m³)':>9s} {'E(MPa)':>9s} {'ν':>6s} {'SIGY':>7s} {'σmax':>6s}")
table: dict = {}
bad: list[str] = []
for pid in sorted(BONES):
    info = part_info.get(pid)
    if not info:
        continue
    ent = lookup(info["mid"])
    if not ent:
        print(f"{pid:>10d} {info['mid']:>10d} {'—':>16s} {'未找到':>32s}  {info['title']}")
        continue
    ro = ent.get("RO")
    ro_kgm3 = ro * 1e12 if ro else None
    E, nu = ent.get("E"), ent.get("PR")
    # ★ 硬闸门：泊松比必须在 (0, 0.5)
    if nu is not None and not (0.0 < nu < 0.5):
        bad.append(f"pid {pid} ({info['title']}): ν={nu} 不在 (0,0.5) —— 卡解析有误")
    if E is not None and E <= 0:
        bad.append(f"pid {pid}: E={E} <= 0")
    s_ro = f"{ro_kgm3:.0f}" if ro_kgm3 else "—"
    s_e = f"{E:.1f}" if E else "—"
    s_nu = f"{nu:.3f}" if nu is not None else "—"
    print(f"{pid:>10d} {info['mid']:>10d} {ent['_source']:>16s} {ent['card']:>32s} "
          f"{s_ro:>9s} {s_e:>9s} {s_nu:>6s} "
          f"{str(ent.get('SIGY')):>7s} {str(ent.get('max_stress_MPa')):>6s}"
          + (f"  {info['title']}" if ent.get("_derived") else ""))
    table[pid] = {"name": BONES[pid], "title": info["title"], "mid": info["mid"],
                  "source": ent["_source"], "card": ent["card"],
                  "RO_kg_m3": ro_kgm3, "E_MPa": E, "nu": nu,
                  "SIGY_MPa": ent.get("SIGY"), "ETAN": ent.get("ETAN"),
                  "FAIL": ent.get("FAIL"), "max_stress_MPa": ent.get("max_stress_MPa"),
                  "plastic_curve": ent.get("plastic_curve"),
                  "BULK": ent.get("BULK"), "G0": ent.get("G0"),
                  "GI": ent.get("GI"), "BETA": ent.get("BETA"),
                  "derived": ent.get("_derived")}

# ── 4. 单位与物理硬校验 ───────────────────────────────────────
print("\n=== 硬校验（单位 + 物理合理性）===")
checks: list[tuple[str, object, bool]] = []
for pid, lo, hi, label in (
    (81000700, 1700, 2100, "皮质骨密度"), (81000900, 1700, 2100, "皮质骨密度"),
    (81001100, 1700, 2100, "皮质骨密度"),
):
    t = table.get(pid)
    if t and t["RO_kg_m3"]:
        checks.append((f"{label} 1700–2100 (pid {pid})", t["RO_kg_m3"],
                       1700 <= t["RO_kg_m3"] <= 2100))
for pid in (81000600, 81000800, 81001000, 81001200):
    t = table.get(pid)
    if t and t["RO_kg_m3"]:
        ok = 300 <= t["RO_kg_m3"] <= 1200
        checks.append((f"松质骨密度 300–1200 (pid {pid})", t["RO_kg_m3"], ok))
for pid in (81000700, 81000900, 81001100, 81001300):
    t = table.get(pid)
    if t and t["E_MPa"]:
        checks.append((f"皮质 E 15–22 GPa (pid {pid})", t["E_MPa"],
                       15000 <= t["E_MPa"] <= 22000))
for pid in (81000600, 81000800, 81001000, 81001200):
    t = table.get(pid)
    if t and t["E_MPa"]:
        checks.append((f"松质 E 50–2000 MPa (pid {pid})", t["E_MPa"],
                       50 <= t["E_MPa"] <= 2000))
for pid in (81000601, 81000801, 81000001):
    t = table.get(pid)
    if t and t["RO_kg_m3"]:
        checks.append((f"骨髓密度 ≈1000 (pid {pid})", t["RO_kg_m3"],
                       900 <= t["RO_kg_m3"] <= 1100))
ok_all = all(c[2] for c in checks) and not bad
for desc, val, ok in checks:
    print(f"  {'✅' if ok else '❌'} {desc:36s} = {val}")
if bad:
    print("\n  ⚠️ 物理不可能值：")
    for b in bad:
        print(f"    ❌ {b}")
print(f"\n  总判定：{'✅ 全部通过 —— mm-s-tonne-N-MPa，E 已是 MPa，FEBio 零换算' if ok_all else '❌ 有不通过项'}")

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "bone_materials.json").write_text(json.dumps({
    "units": "mm-s-tonne-N-MPa",
    "all_checks_passed": bool(ok_all),
    "hard_gate": "0 < nu < 0.5",
    "source": "THUMS AM50 V7.1 (mat_bone_*_fracture.k + 主 .k 内联)",
    "card_field_maps": {
        "MAT_ELASTIC": "MID RO E PR DA DB K",
        "MAT_PIECEWISE_LINEAR_PLASTICITY": "MID RO E PR SIGY ETAN FAIL TDEL + 曲线",
        "MAT_PLASTICITY_WITH_DAMAGE": "MID RO E PR SIGY ETAN",
        "MAT_DAMAGE_2": "MID RO E PR SIGY ETAN EPPF",
        "MAT_VISCOELASTIC": "MID RO BULK G0 GI BETA -> E=9KG/(3K+G), nu=(3K-2G)/(2(3K+G))",
    },
    "parts": table}, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"\n写出 {OUT/'bone_materials.json'}")
