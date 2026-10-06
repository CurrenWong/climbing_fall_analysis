"""从 THUMS V7.1 主 .k（zip 内，219 MB）流式提取下肢骨 / 韧带的 FE 网格。

★ 实测格式（务必记住 —— 同一个文件里三种字段宽度并存）★
  *NODE                    : I8 + 3*F16     有空格，空白分隔即可
  *ELEMENT_SOLID           : **8 字符定宽，无空格**  eid,pid,n1..n8   行=80 字符
  *ELEMENT_SHELL           : **8 字符定宽，无空格**  eid,pid,n1..n4   行=48 字符
  *ELEMENT_BEAM            : **8 字符定宽**          eid,pid,n1,...
  *ELEMENT_SHELL_THICKNESS : 厚度值（不是单元），10 字符浮点
  *PART 数值卡             : **10 字符定宽**         pid,secid,mid... 行=90 字符

踩过的坑：按空白 split 会让 1,696,785 行（占 99.998%）只解析出 1 个 token，
静默全废 —— 必须定宽切片，并用「引用节点必须存在于 *NODE」做交叉校验。

输出：temp/thums/thums_lowerlimb.npz + summary.json
"""
from __future__ import annotations

import json
import re
import sys
import zipfile
from pathlib import Path

import numpy as np

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"
OUT = Path(r"D:\Project\climbing_fall_analysis\temp\thums")

# pid -> (名称, 期望类别)   —— 来自 diag_k_fmt4.py 的实测发现
TARGETS: dict[int, tuple[str, str]] = {
    # ── 右腿 腓骨（本项目缺的就是它）──────────────
    81000900: ("R_FIBULA_CORT", "solid"),
    81000800: ("right_fibula_end_spon", "solid"),
    81000801: ("right_fibula_center_spon", "solid"),
    # ── 右腿 胫骨 ────────────────────────────
    81000700: ("R_TIBIA_CORT", "solid"),
    81000600: ("right_tibia_end_spon", "solid"),
    81000601: ("right_tibia_center_spon", "solid"),
    # ── 右腿 距骨 / 跟骨 ──────────────────────
    81001100: ("R_TALUS_CORT", "solid"),
    81001000: ("R_TALUS_SPON", "solid"),
    81001300: ("R_CALCANEUS_CORT", "solid"),
    81001200: ("R_CALCANEUS_SPON", "solid"),
    # ── 右腿 股骨 / 髌骨 ──────────────────────
    81000100: ("R_FEMUR_CORT", "solid"),
    81000000: ("right_femur_end_spon", "solid"),
    81000001: ("right_femur_center_spon", "solid"),
    81000300: ("R_PATELLA_CORT", "solid"),
    81000200: ("R_PATELLA_SPON", "solid"),
    # ── 踝部韧带（壳单元）──────────────────────
    81100801: ("L_LIG_CALCANEOFIBULARE", "shell"),      # 跟腓
    81100701: ("L_LIG_TALOFIBULARE_POSTERIUS", "shell"),
    81100901: ("L_LIG_TALOFIBULARE_ANTERIUS", "shell"),
    81101101: ("L_LIG_TIBIOFIBULARE", "shell"),          # 下胫腓
    81101001: ("L_LIG_TALOCALCANEUM_POSTERIUS", "shell"),
    81101401: ("L_LIG_LONG_PLANTAR", "shell"),
    81101201: ("L_LIG_MEDIALE1", "shell"),
    81101202: ("L_LIG_MEDIALE2", "shell"),
    81101203: ("L_LIG_MEDIALE3", "shell"),
    81101204: ("L_LIG_MEDIALE4", "shell"),
    81100100: ("R_LIG_COLLATERAL_FIBULAR_TISSUES", "solid"),
    81100200: ("R_LIG_COLLATERAL_TIBIAL_TISSUES", "solid"),
    # ── 远端跗骨（★ 2026-10-04 补：胫舟韧带需要舟骨锚定）────────────
    #    实测：跗骨 SPON=solid(hex/tet) / CORT=shell（与距跟骨的 CORT=solid 不同）
    81001400: ("R_NAVICULAR_SPON", "solid"),
    81001401: ("R_NAVICULAR_CORT", "shell"),
    81001500: ("R_CUBOID_SPON", "solid"),
    81001501: ("R_CUBOID_CORT", "shell"),
    81001600: ("R_MEDIAL_CUNEIFORM_SPON", "solid"),
    81001601: ("R_MEDIAL_CUNEIFORM_CORT", "shell"),
    81001700: ("R_INTERMED_CUNEIFORM_SPON", "solid"),
    81001701: ("R_INTERMED_CUNEIFORM_CORT", "shell"),
    81001800: ("R_LATERAL_CUNEIFORM_SPON", "solid"),
    81001801: ("R_LATERAL_CUNEIFORM_CORT", "shell"),
    # ── 足部韧带（★ 补：胫舟 / 弹簧韧带可能在其中）──────────────────
    81101301: ("L_LIG_FOOT_1", "shell"),
    81101302: ("L_LIG_FOOT_2", "shell"),
    81101303: ("L_LIG_FOOT_3", "shell"),
    81101304: ("L_LIG_FOOT_4", "shell"),
    81101305: ("L_LIG_FOOT_5", "shell"),
    81101306: ("L_LIG_FOOT_6", "shell"),
    81101307: ("L_LIG_FOOT_7", "shell"),
    81101308: ("L_LIG_FOOT_8", "shell"),
    81101408: ("L_LIG_FOOT_10", "shell"),
    81101409: ("L_LIG_FOOT_11", "shell"),
    81101410: ("L_LIG_CAPITIS_FIBULAE", "shell"),
}

RE_KW = re.compile(r"^\*([A-Z_0-9]+)")
ELEM_W8 = {"ELEMENT_SOLID", "ELEMENT_SHELL", "ELEMENT_BEAM"}


def slice_fields(s: str, width: int) -> list[int] | None:
    """按定宽切分并转 int；任何一段非整数则返回 None（不静默吞掉）。"""
    vals: list[int] = []
    for i in range(0, len(s) - len(s) % width, width):
        f = s[i : i + width].strip()
        if not f:
            continue
        try:
            vals.append(int(f))
        except ValueError:
            return None
    return vals


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    node_ids: list[int] = []
    node_xyz: list[tuple[float, float, float]] = []
    part_titles: dict[int, str] = {}
    # (pid, kw) -> list[list[node ids]]
    elems: dict[tuple[int, str], list[list[int]]] = {}
    fails: dict[str, int] = {}

    cur: str | None = None
    pending_title: str | None = None
    buf = ""
    nbytes = 0

    with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
        while True:
            chunk = f.read(8 * 1024 * 1024)
            if not chunk:
                break
            nbytes += len(chunk)
            buf += chunk.decode("latin-1")
            lines = buf.split("\n")
            buf = lines.pop()

            for raw in lines:
                s = raw.rstrip("\r")
                st = s.strip()
                if not st or st[0] == "$":
                    continue
                if st[0] == "*":
                    m = RE_KW.match(st.upper())
                    cur = m.group(1) if m else None
                    if cur == "PART":
                        pending_title = None
                    continue

                if cur == "NODE":
                    t = st.split()
                    if len(t) >= 4:
                        try:
                            node_ids.append(int(t[0]))
                            node_xyz.append((float(t[1]), float(t[2]), float(t[3])))
                            continue
                        except ValueError:
                            pass
                    fails["NODE"] = fails.get("NODE", 0) + 1

                elif cur == "PART":
                    if pending_title is None:
                        pending_title = st
                    else:
                        v = slice_fields(st, 10)          # ★ 10 字符
                        if v:
                            part_titles[v[0]] = pending_title
                        else:
                            fails["PART"] = fails.get("PART", 0) + 1
                        pending_title = None

                elif cur in ELEM_W8:
                    v = slice_fields(s, 8)                # ★ 8 字符
                    if v is None or len(v) < 3:
                        fails[cur] = fails.get(cur, 0) + 1
                        continue
                    pid = v[1]
                    if pid in TARGETS and TARGETS[pid][1] == cur.split("_")[1].lower():
                        elems.setdefault((pid, cur), []).append(v[2:])

            if nbytes % (64 * 1024 * 1024) < 8 * 1024 * 1024:
                print(f"  ... {nbytes/1024/1024:6.1f} MB  nodes={len(node_ids):,}",
                      flush=True)

    print(f"\n扫描完成 {nbytes/1024/1024:.1f} MB")
    print(f"节点 {len(node_ids):,}  |  *PART 名 {len(part_titles):,}  |  解析失败 {fails}")

    nid = np.asarray(node_ids, dtype=np.int64)
    xyz = np.asarray(node_xyz, dtype=np.float64)
    idx = {int(v): i for i, v in enumerate(nid)}

    print(f"\n{'pid':>10s} {'kw':16s} {'type':6s} {'elems':>7s} {'nodes':>8s} "
          f"{'miss':>5s}  name")
    summary: dict[str, dict] = {}
    keep: set[int] = set()
    for (pid, kw), els in sorted(elems.items()):
        name, want = TARGETS[pid]
        refd = {n for e in els for n in e}
        missing = [n for n in refd if n not in idx]
        keep |= (refd - set(missing))
        print(f"{pid:>10d} {kw:16s} {want:6s} {len(els):>7,} {len(refd):>8,} "
              f"{len(missing):>5,}  {name}"
              + (f"  ⚠️ 缺 {sorted(missing)[:4]}" if missing else ""))
        summary[f"{pid}_{want}"] = {
            "pid": pid, "name": name, "kind": want, "block": kw,
            "n_elements": len(els), "n_nodes": len(refd),
            "n_missing_nodes": len(missing),
            "nodes_per_element": sorted({len(e) for e in els}),
            "part_title": part_titles.get(pid, ""),
        }

    # 目标缺失检查
    got = {pid for pid, _ in elems}
    for pid, (name, kind) in TARGETS.items():
        if pid not in got:
            print(f"  ⚠️ 未提取到: {pid} {name} ({kind})")

    kn = np.asarray(sorted(keep), dtype=np.int64)
    kmap = np.asarray([idx[int(v)] for v in kn], dtype=np.int64)
    arrays = {
        "node_ids": kn,
        "node_xyz": xyz[kmap],
    }
    for (pid, kw), els in elems.items():
        name, want = TARGETS[pid]
        w = max(len(e) for e in els)
        arr = np.full((len(els), w), -1, dtype=np.int64)
        for i, e in enumerate(els):
            arr[i, : len(e)] = e
        arrays[f"{pid}_{want}"] = arr

    np.savez_compressed(OUT / "thums_lowerlimb.npz", **arrays)
    (OUT / "summary.json").write_text(
        json.dumps({"n_nodes_total": len(node_ids), "n_nodes_kept": int(len(kn)),
                    "nbytes": nbytes, "parse_fails": fails,
                    "n_parts_named": len(part_titles), "parts": summary},
                   ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n写出 {OUT/'thums_lowerlimb.npz'}（保留 {len(kn):,} 节点）+ summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
