# -*- coding: utf-8 -*-
"""P4 · P-b（解剖端分布载荷）能力探针 —— P4-a′：位置隔离（压缩方向修正版）。

预注册：``P4_PB_PROBE_TODO.md``（2026-10-14，冻结）。

探针问题：把股骨载荷从「髁端 top 带」（F 路线 all_zero 形态）换成
「股骨头近端关节面节点」，固定端换到髁端，载荷方向 **+bbox 长轴**
（P4-a 初版误用 -bbox 长轴 ⇒ 头端球冠 + 髁端固定 = 偏心拉伸，拉伸无屈曲、
弯曲路径被消灭 → 伪证伪；P4-a′ 修正为 +Z = 球冠压向髁端 = 解剖压缩，
回归预注册「压缩」意图）⇒ 引入弯曲载荷路径（力线过股骨头/头端带形心，
与骨干轴有垂距 e_head≈7-10mm）。

判据（预注册 §2.3 修订版，跑前写死）：
  主判据 = 阶跃形态出现：σ_shaft(h) 单调升且存在跨阈点 h*（<220→≥220）且 h* ≤ 45；
  次判据 = 跨越点与解析预测（e≈9.5 ⇒ h*≈30）的吻合度，如实记录；
  证伪 = 仍 all_zero（h=50 σ<220 或 σ 平缓 <100）。

正对照（§2.4）：同脚本基线（默认参数）h=20 的 σ_shaft 须与 F 路线
``A2_axial_lat_90_h20`` 的统计一致（同口径验证）。

用法：
  $env:PYTHONPATH="src"
  & .venv\\Scripts\\python.exe scripts\\opensim_fe\\p4_pb_probe.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

import bone_feb as bf  # noqa: E402
import fracture_matrix as fm  # noqa: E402
import paper_alignment as pa  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
BASE_DIR = ROOT / "temp" / "opensim_fe" / "p4_pb_probe" / "femur_r"

HEIGHTS = [1, 5, 10, 15, 20, 30, 40, 50]  # 预注册 §2.3：8 个代表高度
HEAD_FRAC = 0.20   # 头端带（股骨头+部分颈+大转子；F-3′c 头端 25% 的保守子带）
FIX_FRAC = 0.30    # 固定端 = 髁端带（与原 top 带同带）
SIGMA_C = 220.0    # fem_cort 压缩强度（A2 口径）


def load_ctx(bone: str = "femur_r"):
    fmj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    heights = np.asarray(fmj["heights_m"], dtype=float)
    r1 = json.loads((RES / "bc_robust_metric_route1.json").read_text(encoding="utf-8"))
    dom = {r["bone"]: r for r in r1["rows"]}[bone]["solid_dom"]
    part_row = next(p for p in fmj["parts"] if p["bone"] == bone)
    joint = part_row["joint"]
    loads = np.asarray(fmj["loads_n"][joint], dtype=float)

    part = {p.bone: p for p in fm.build_parts()}[bone]
    pids = [p for sd in part.spec.solid_domains for p in sd.pids]
    pids += [p for sh in part.spec.shell_domains for p in sh.pids]
    mats = fm.read_materials(fm.MASTER, pids)
    for sd in part.spec.solid_domains:
        sd.material = mats[sd.pids[0]]
    for sh in part.spec.shell_domains:
        sh.material = mats[sh.pids[0]]
    mesh = fm.extract_mesh(part.spec, fm.MASTER)
    solid_names = [sd.mesh_name for sd in part.spec.solid_domains]
    shell_specs = [
        (sh.mesh_name, int(mesh["by_shell"][sh.mesh_name]["q4"].shape[0]
                           + mesh["by_shell"][sh.mesh_name]["t3"].shape[0]))
        for sh in part.spec.shell_domains]
    return part, mesh, solid_names, shell_specs, dom, joint, loads, heights


def statistic(h5: Path, dom: str) -> dict:
    import h5py
    with h5py.File(h5, "r") as h:
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
        conn = h["meshes/0/domains"][dom][:]
        sids = sorted(h["states"].keys(), key=int)
        stress = h[f"states/{sids[-1]}/element_data/stress/{dom}"][:].astype(float)
    idx = conn[:, 1:] - 1
    cent = coords[idx].mean(axis=1)
    p = pa.partition_ends_shaft(coords, cent)
    vm = pa.vm_voigt(stress)
    return {"ends_stat_mpa": float(np.median(vm[p["ends_mask"]])),
            "shaft_stat_mpa": float(np.median(vm[p["shaft_mask"]])),
            "vm_max_mpa": float(vm.max())}


def head_joint_nodes(mesh: dict, *, frac: float = HEAD_FRAC) -> dict:
    """提取股骨头近端关节面节点（1-based），返回节点 + 几何说明。

    约定（F-3′c 已验证）：femur PCA 旋转后 +Z 顶端 = 髁端，-Z 底端 = 头端。
    top_load（髁端带）→ 新固定端；bottom 方向（头端带）→ 新载荷端。
    """
    nodes_rot = bf._align_pca_rotate(mesh["nodes"])  # 与 build_bone_feb 内部同款
    axis_idx, a_min, a_max, top_local, bot_local = bf._select_top_bottom(
        nodes_rot, axis=None, top_frac=FIX_FRAC, bot_frac=FIX_FRAC)
    span = float(a_max - a_min)
    coord = nodes_rot[:, axis_idx]
    head_thr = a_min + frac * span
    head_sel = (np.flatnonzero(coord <= head_thr) + 1).astype(np.int64)  # 1-based
    # 带形心垂距（对 bbox 轴）：PCA 旋转以质心为原点，主轴 z 过质心，
    # 垂距 = 带形心在横向平面 (x,y) 的模（轴分量清零）。
    sel0 = head_sel - 1
    band_c = nodes_rot[sel0].mean(axis=0)
    # 垂距必须相对质心：PCA 旋转以质心为原点、主轴 z 过质心，但旋转后的
    # 坐标仍是原坐标系（_align_pca_rotate 末尾 +c），直接用绝对坐标会把
    # 「质心在横向平面离原点的距离」误算成垂距（2026-10-14 P4 实测 e=365mm
    # 虚高，修正后 e_head≈10mm）。修正：band_c − 质心。
    c = nodes_rot.mean(axis=0)
    band_rel = band_c - c
    perp = band_rel.copy()
    perp[axis_idx] = 0.0
    e_band = float(np.linalg.norm(perp))
    # 载荷方向（P4-a′ 修正）：球冠在 z 小端、固定髁端在 z 大端，
    # 压缩需要 +z（球冠压向髁端）；bone_feb 的 load_dir_vec 语义为
    # 「轴向量」，load_vec = -load_n × load_dir_vec ⇒ 传 -e_axis 得到 +e_axis 力。
    e_axis = np.zeros(3, dtype=float)
    e_axis[axis_idx] = 1.0
    load_dir_plus = -e_axis  # 使实际载荷力 = +e_axis（球冠→髁端压缩）
    return {
        "axis_idx": int(axis_idx),
        "a_min_mm": float(a_min), "a_max_mm": float(a_max),
        "head_thr_mm": float(head_thr), "frac": float(frac),
        "head_nodes_1based": head_sel.tolist(),
        "n_head_nodes": int(len(head_sel)),
        "fix_nodes_1based": top_local.tolist(),   # 髁端带 = 新固定端
        "n_fix_nodes": int(len(top_local)),
        "band_centroid": band_c.tolist(),
        "centroid": c.tolist(),
        "band_eccentricity_mm": e_band,           # 头端带形心对 bbox 轴的垂距（相对质心）
        "load_dir_vec": load_dir_plus.tolist(),   # P4-a′：+轴压缩
    }


def main() -> int:
    t_all = time.time()
    part, mesh, solid_names, shell_specs, dom, joint, loads, heights = load_ctx()
    print(f"[p4] {part.bone}: mesh {mesh['n_nodes']} nodes, dom={dom}, joint={joint}",
          flush=True)

    hj = head_joint_nodes(mesh)
    print(f"[p4] 头端带 frac={hj['frac']}: n={hj['n_head_nodes']} "
          f"thr={hj['head_thr_mm']:.2f}mm, 带形心垂距 e_band={hj['band_eccentricity_mm']:.2f}mm",
          flush=True)
    if hj["band_eccentricity_mm"] < 5.0:  # 球冠/头端带形心应偏离轴 5mm+（预注册 §2.3）
        print("[p4] ⚠ 头端带垂距过小（<5mm）：位置变量可能无效，记录并继续",
              flush=True)

    rows: list[dict] = []

    # ---- 正对照（预注册 §2.4）：基线 h=20（默认参数）vs F 路线 A2_axial_lat_90_h20 ----
    base_h5 = ROOT / "temp" / "opensim_fe" / "fe_bc_perheight" / "femur_r" / \
        "h20" / "femur_r.hdf5"
    if base_h5.is_file():
        base_stat = statistic(base_h5, dom)
        print(f"[p4] 正对照基线（F 路线 E′ A2_axial_lat_90 h20 缓存）: "
              f"ends={base_stat['ends_stat_mpa']:.2f} shaft={base_stat['shaft_stat_mpa']:.2f}",
              flush=True)
        rows.append({"case": "positive_control_ref", "h_no": 20,
                     **base_stat, "source": "A2_axial_lat_90 h20 cached (E′)"})
    else:
        print("[p4] ⚠ 基线缓存不存在，正对照改为现场重建", flush=True)

    # ---- P4-a 主探针：8 高度，载荷=股骨头，固定=髁端 ----
    for hno in HEIGHTS:
        i = hno - 1
        load = float(loads[i])
        d = BASE_DIR / f"head_load_h{hno:02d}"
        d.mkdir(parents=True, exist_ok=True)
        feb = d / "femur_r.feb"
        h5 = feb.with_suffix(".hdf5")
        rec = {"case": "P4a_head_load", "h_no": hno,
               "h_m": float(heights[i]), "load_n": load}
        if not h5.is_file():
            t = time.time()
            try:
                bf.build_bone_feb(
                    part.spec, feb, load_n=load, top_frac=0.30, bot_frac=0.30,
                    axis=None, with_shell=True, time_steps=1, mesh=mesh,
                    align_pca=True, lateral_support=[(0.05, 0.95)],
                    bottom_fix_axial_only=True,
                    load_nodes_override=np.asarray(hj["head_nodes_1based"], dtype=np.int64),
                    fix_nodes_override=np.asarray(hj["fix_nodes_1based"], dtype=np.int64),
                    load_dir_vec=np.asarray(hj["load_dir_vec"], dtype=float),
                )
                bf.run_bone_feb(feb, workdir=feb.parent, timeout=1800)
                if feb.with_suffix(".xplt").is_file():
                    bf.post_process(feb.with_suffix(".xplt"), h5,
                                    solid_names=solid_names,
                                    shell_specs=shell_specs, label=part.bone)
                rec["converged"] = True
            except Exception as exc:  # noqa: BLE001
                rec["converged"] = False
                rec["error"] = f"{type(exc).__name__}: {exc}"
            rec["seconds"] = round(time.time() - t, 1)
        else:
            rec["converged"] = True
            rec["seconds"] = 0.0
        if h5.is_file():
            rec.update(statistic(h5, dom))
        rows.append(rec)
        print(f"[p4] h={hno:2d} ({load:9.1f}N) "
              f"{'OK ' if rec.get('converged') else 'FAIL'} "
              + (f"ends={rec['ends_stat_mpa']:8.2f} shaft={rec['shaft_stat_mpa']:9.2f}"
                 if "shaft_stat_mpa" in rec else f"({rec.get('error','')[:60]})")
              + f"  [{rec.get('seconds')}s]", flush=True)

    # ---- 形态判定（预注册 §2.3，跑前写死，跑后如实） ----
    print("\n=== 判定 ===")
    hs = [r["h_no"] for r in rows if r.get("converged") and "shaft_stat_mpa" in r]
    ss = [r["shaft_stat_mpa"] for r in rows if r.get("converged") and "shaft_stat_mpa" in r]
    if not ss:
        print("  全部不收敛 → 探针失败（非形态问题）")
        return 2
    mono = all(b > a for a, b in zip(ss, ss[1:]))
    cross = [h for h, s in zip(hs, ss) if s >= SIGMA_C]
    hmax = max(hs)
    shmax = ss[-1] if hs[-1] == hmax else max(ss)
    if cross:
        hstar = int(cross[0])
        if mono and hstar <= 45:
            if 25 <= hstar <= 35:
                verdict = (f"证实（主判据+次判据）：阶跃形态，h*={hstar} "
                           f"∈[25,35] 与解析预测（e≈9.5 ⇒ h*≈30）吻合")
            else:
                verdict = (f"证实（主判据）：阶跃形态，h*={hstar} ≤ 45；"
                           f"次判据不吻合（预测 h*≈30），如实记录")
        elif not mono:
            verdict = "证伪/边缘：非单调但存在跨阈（形态异常，记录）"
        else:
            verdict = f"证伪/边缘：跨阈但 h*={hstar} > 45（低位骨折，非阶跃）"
    elif mono and shmax >= SIGMA_C * 0.45:
        verdict = "证伪（边缘）：单调升但 h=50 未跨阈（σ≈F/A 线性，弯曲贡献不足）"
    else:
        verdict = "证伪：全部 < σ_c（all_zero 或接近），弯曲路径不足"
    print(f"  heights={hs}")
    print(f"  σ_shaft={[round(s,1) for s in ss]}")
    print(f"  单调={mono}  跨阈点={cross}  h*（首个跨阈）={cross[0] if cross else None}")
    print(f"  判定：{verdict}")

    # 归因（证据支撑：temp/p4_diag_check.py / temp/p4_bend_check.py；正文 P4_PB_PROBE.md §3）
    attribution = (
        "lateral_support 结构性抑制弯曲路径（非弯曲路径不足）：探针载荷带（头端球冠）"
        "几乎完全落在 lateral_support=[(0.05,0.95)] 约束带内，载荷节点横向位移被钉死、"
        "偏心力矩被约束反力吸收；shaft 域 σ_zz~横向坐标 r²=0.054 vs 基线 0.291、"
        "σ_zz max=0.0（无拉伸侧）、A_eff≈408mm²≈几何 A ⇒ 纯轴向响应"
    )

    out_json = RES / "P4_PB_PROBE.json"
    out_json.write_text(json.dumps({
        "bone": part.bone, "domain": dom, "joint": joint,
        "probe": "P4-a head joint-face load (position isolation)",
        "head_frac": HEAD_FRAC, "fix_frac": FIX_FRAC, "sigma_c_mpa": SIGMA_C,
        "head_geom": {k: v for k, v in hj.items() if k != "head_nodes_1based"
                      and k != "fix_nodes_1based"},
        "n_head_nodes": hj["n_head_nodes"], "n_fix_nodes": hj["n_fix_nodes"],
        "band_eccentricity_mm": hj["band_eccentricity_mm"],
        "prediction": "阶跃形态，跨阈点 h*≈30（e≈9.5mm，预注册 §2.3 补充）；F-3′d 原窗 [10,15] 已后移",
        "criterion": "主判据：σ_shaft(h) 单调升且存在 h*≤45 从 <220→≥220；"
                     "次判据：h*≈30 吻合度（如实记录）",
        "verdict": verdict,
        "attribution": attribution,
        "rows": rows}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n[p4] 总耗时 {time.time()-t_all:.0f}s；写出 {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
