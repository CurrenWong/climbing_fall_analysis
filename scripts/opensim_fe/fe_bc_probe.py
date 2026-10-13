# -*- coding: utf-8 -*-
"""FE_BC 大件 · W2 能力探针：端部 BC 能不能把股骨从 `all_one` 里救出来？

背景（10-13 诊断，冻结件 `results/opensim_fe/A2_MISS_DIAGNOSIS.md`）：
  股骨 `ends_stat` = 14.94 MPa（阈 80，伪低）、`shaft_stat` = 426.77 MPa（阈 220，比值 1.94）。
  S4 判据对长骨是 **OR**：`(ends_h>=sig_e) | (shaft_h>=sig_s)` ⇒ 只要 `shaft_h>=220` 对
  所有 h>=1 成立就是 `all_one`（无过渡 ⇒ 不可能是「阶跃」）。
  ⇒ **要翻案，`shaft_stat` 必须掉 >= 1.94x**（或 ends/shaft 的相对关系改变）。

代码病根（`bone_feb.py:1027-1051`）：
  * P-a  `BCZeroDisplacement(node_set="bottom_fix", x_dof=1, y_dof=1, z_dof=1)`
        覆盖 `bot_frac=0.30` 的**全 3 DOF 刚性夹持带**——而 `partition_ends_shaft`
        的 `END_FRACTION=0.25` ⇒ **整个远端 ends 落在夹持带内** ⇒ ends 必然趋零。
  * P-b  载荷 = 沿骨长轴/PCA 轴的**均匀 NodalForce** 铺在 `top_frac=0.30` 带上。

本脚本**只做能力探针**（不建全链、不改判据、不写任何冻结件）：
  对**股骨单骨**、**与 BC2 逐字相同的弹性条件**，只动一个旋钮，看
  `ends_stat` / `shaft_stat` 能移动多少 ⇒ 回答「F 路线值不值得开」这个可判定性问题。

纪律
----
* 冻结输入**只读**；输出全部在 `temp/opensim_fe/fe_bc_probe/`。
* 统计量与 `s4_alignment_bc2.py:159-171` **逐字相同**（median of von Mises over
  `pa.partition_ends_shaft` 的 ends/shaft mask，main solid domain）。
* `V0_control` 用**与 BC2 完全相同的参数**重建重解 ⇒ 它必须复现冻结的
  14.94 / 426.77；**复现不了就说明本探针测的不是同一个量，后续变体结论作废**。
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

import bone_feb as bf  # noqa: E402
import fracture_matrix as fm  # noqa: E402
import paper_alignment as pa  # noqa: E402

BONE = "femur_r"
REF_IDX = 4                                   # s4: sigma_h_at_ref = vm_h[4]
BC2_TEMP = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2"
PROBE_TEMP = ROOT / "temp" / "opensim_fe" / "fe_bc_probe"
OUT_JSON = ROOT / "results" / "opensim_fe" / "FE_BC_PROBE.json"
OUT_MD = ROOT / "results" / "opensim_fe" / "FE_BC_PROBE.md"

#: (name, top_frac, bot_frac, note)
VARIANTS: list[tuple[str, float, float, str]] = [
    ("V0_control", 0.30, 0.30, "与 BC2 逐字相同（正对照：必须复现 14.94 / 426.77）"),
    ("V1_bot015",  0.30, 0.15, "只缩夹持带：ends 里还剩多少是伪低？"),
    ("V2_bot005",  0.30, 0.05, "极限缩夹持带"),
    ("V3_wide",    0.50, 0.05, "夹持带最小 + 加载带最大（两端都远离解剖端理想化）"),
]


class Log:
    """同时往 stdout 和探针目录写（`solve_bone_fe` 要一个有 .write 的 log）。"""

    def __init__(self, path: Path) -> None:
        self.path = path
        path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = path.open("a", encoding="utf-8")

    def __call__(self, msg: str) -> None:
        print(msg, flush=True)
        self.fh.write(str(msg) + "\n")
        self.fh.flush()

    def write(self, msg: str) -> None:
        self.__call__(msg)


def statistic(h5_path: Path, dom: str) -> dict:
    """`s4_alignment_bc2.py:159-171` 的逐字复刻（median 口径）。"""
    import h5py

    with h5py.File(h5_path, "r") as h:
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
        conn = h["meshes/0/domains"][dom][:]
        sids = sorted(h["states"].keys(), key=int)
        stress = h[f"states/{sids[-1]}/element_data/stress/{dom}"][:].astype(float)

    idx = conn[:, 1:] - 1
    cent = coords[idx].mean(axis=1)
    part = pa.partition_ends_shaft(coords, cent)
    vm = pa.vm_voigt(stress)
    return {
        "domain": dom,
        "n_elem": int(vm.shape[0]),
        "n_ends": int(part["n_ends"]),
        "load_axis": part["load_axis"],
        "ends_stat_mpa": float(np.median(vm[part["ends_mask"]])),
        "shaft_stat_mpa": float(np.median(vm[part["shaft_mask"]])),
        "vm_max_mpa": float(vm.max()),
        "vm_p99_mpa": float(np.percentile(vm, 99)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="只跑名字含该子串的变体")
    ap.add_argument("--time-steps", type=int, default=1)
    args = ap.parse_args()

    PROBE_TEMP.mkdir(parents=True, exist_ok=True)
    log = Log(PROBE_TEMP / "run.log")

    # 输入（冻结）只读；输出根改到探针目录
    fm.OUT_TEMP = PROBE_TEMP
    src = BC2_TEMP / "load_chain.npz"
    if not (PROBE_TEMP / "load_chain.npz").is_file() and src.is_file():
        shutil.copy2(src, PROBE_TEMP / "load_chain.npz")

    frozen = json.loads(
        (ROOT / "results" / "opensim_fe" / "fracture_matrix_bc2.json")
        .read_text(encoding="utf-8"))
    ref = next(p for p in frozen["parts"] if p["bone"] == BONE)
    rl = float(ref["fe_solve_load_n"])
    joint = ref["joint"]
    r1 = json.loads(
        (ROOT / "results" / "opensim_fe" / "bc_robust_metric_route1.json")
        .read_text(encoding="utf-8"))
    dom = {r["bone"]: r for r in r1["rows"]}[BONE]["solid_dom"]

    log(f"[probe] {BONE}: 冻结 fe_solve_load_n={rl:g} N (joint={joint}), "
        f"main domain={dom}")

    # 冻结基线（正对照的靶）
    fz = statistic(BC2_TEMP / BONE / f"{BONE}.hdf5", dom)
    log(f"[probe] 冻结 BC2 统计量: ends={fz['ends_stat_mpa']:.2f} "
        f"shaft={fz['shaft_stat_mpa']:.2f} (n_ends={fz['n_ends']})")

    parts = {p.bone: p for p in fm.build_parts()}
    part = parts[BONE]

    rows: list[dict] = []
    for name, top_frac, bot_frac, note in VARIANTS:
        if args.only and args.only not in name:
            continue
        d = PROBE_TEMP / name
        d.mkdir(parents=True, exist_ok=True)
        h5 = d / f"{BONE}.hdf5"
        rec: dict = {"variant": name, "top_frac": top_frac, "bot_frac": bot_frac,
                     "note": note}

        def _build(part_, mesh, mats, feb, load, time_steps,
                   _t=top_frac, _b=bot_frac) -> None:
            for sd in part_.spec.solid_domains:
                sd.material = mats[sd.pids[0]]
            for sh in part_.spec.shell_domains:
                sh.material = mats[sh.pids[0]]
            bf.build_bone_feb(part_.spec, feb, load_n=load, top_frac=_t,
                              bot_frac=_b, axis=None, with_shell=True,
                              time_steps=time_steps, mesh=mesh, align_pca=True)
            bf.run_bone_feb(feb, workdir=feb.parent, timeout=1800)

        if h5.is_file():
            log(f"[probe] {name}: 已存在，跳过解算")
        else:
            fm._build_run = _build
            t0 = time.time()
            try:
                fm.solve_bone_fe(part, rl, log, reuse=False,
                                 time_steps=args.time_steps)
            except Exception as exc:  # noqa: BLE001
                rec["error"] = f"{type(exc).__name__}: {exc}"
                rec["traceback"] = traceback.format_exc()[-1200:]
                log(f"[probe] {name}: FAIL {rec['error']}")
                rows.append(rec)
                continue
            # solve_bone_fe 写到 OUT_TEMP/<bone>/；搬到变体目录
            produced = PROBE_TEMP / BONE / f"{BONE}.hdf5"
            if produced.is_file() and produced.resolve() != h5.resolve():
                shutil.move(str(produced), str(h5))
            rec["seconds"] = round(time.time() - t0, 1)
            log(f"[probe] {name}: 解算完成 {rec['seconds']}s")

        if not h5.is_file():
            rec["error"] = "hdf5 未产出"
            rows.append(rec)
            continue

        st = statistic(h5, dom)
        rec.update(st)
        fe_json = PROBE_TEMP / BONE / f"{BONE}_fe.json"
        if fe_json.is_file():
            try:
                rec["fe_record"] = json.loads(fe_json.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
                pass
        rec["ends_ratio_vs_80"] = rec["ends_stat_mpa"] / 80.0
        rec["shaft_ratio_vs_220"] = rec["shaft_stat_mpa"] / 220.0
        rec["d_ends_vs_frozen"] = rec["ends_stat_mpa"] - fz["ends_stat_mpa"]
        rec["d_shaft_vs_frozen"] = rec["shaft_stat_mpa"] - fz["shaft_stat_mpa"]
        log(f"[probe] {name}: ends={rec['ends_stat_mpa']:8.2f} "
            f"(x{rec['ends_ratio_vs_80']:.2f} vs 80)  "
            f"shaft={rec['shaft_stat_mpa']:8.2f} "
            f"(x{rec['shaft_ratio_vs_220']:.2f} vs 220)")
        rows.append(rec)

    payload = {
        "bone": BONE, "joint": joint, "domain": dom, "ref_load_n": rl,
        "frozen_baseline": fz, "bc2_contract": {
            "top_frac": 0.30, "bot_frac": 0.30, "align_pca": True,
            "axis": None, "with_shell": True, "material": "linear elastic",
            "source": "scripts/opensim_fe/rerun_matrix_bc2.py::_build_run_bc2",
        },
        "criteria": {"ends_sigma_c_mpa": 80.0, "shaft_sigma_c_mpa": 220.0,
                     "long_bone_rule": "(ends_h>=80) | (shaft_h>=220)"},
        "variants": rows,
        "reads_frozen_only": True,
        "writes": [str(OUT_JSON), str(OUT_MD)],
    }
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                        encoding="utf-8")
    log(f"[probe] 写出 {OUT_JSON}")
    print(json.dumps(rows, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
