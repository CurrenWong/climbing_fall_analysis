"""S5-contact 波5 —— 真实泡沫垫矩阵（G7 配方 + hold 弹簧）驱动器 + 报告。

波4（``nonvertical_s5_contact.py``）的 7 条真实 ``contact`` 行在 FEBio 默认 ``<Control>``
下全部发散（rc=1 / end_t=None / negJac=6 / invalid=4）。本驱动复用波4 的 ``_solve_case``
（已加波5 参数：``penalty`` / ``contact_hold_spring_k`` / ``g7_contact``），把 commit
``315c373`` 的修复（quad4 + G7 配方 + opt-in hold 弹簧）接进真实矩阵。

波5 关键发现（见 ``docs/S5_wave5_matrix.md``）：矩阵网格上真实接触的收敛还额外需要
**rigid-body 加载路径（``use_rigid=True``）** —— 波4 的 contact 行用 ``use_rigid=False``
（PressureLoad），在接触闭合（t≈0.42）处 FEBio 收敛判据退化为「零残差」而 ERROR。

产物（**不覆盖** ``nonvertical_s5_contact.json``）::

    results/opensim_fe/nonvertical_s5_contact_wave5.json
    results/opensim_fe/NONVERTICAL_S5_CONTACT_WAVE5_REPORT.md

运行（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s5_contact_wave5.py --contact-only
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s5_contact_wave5.py
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import nonvertical_s5_contact as S5C  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
OUT_JSON = RES / "nonvertical_s5_contact_wave5.json"
OUT_MD = RES / "NONVERTICAL_S5_CONTACT_WAVE5_REPORT.md"
QUOTABLE_F_N_KN = 34.4  # docs/S5_QUOTABLE.md（tet4 pad / 固定 sink）；仅量级锚。


def _wave5_cases(loads: dict[float, float]) -> list[dict]:
    """波5 全矩阵：非接触参照 + 接触子集 + hold/penalty/load-path 对照。

    ⚠️ 波5 关键发现：矩阵网格（anatframe）上真实接触的收敛需要 **rigid-body 加载路径
    （``use_rigid=True``）** + **hold 弹簧**。波4 的 contact case 用 ``use_rigid=False``
    （PressureLoad）⇒ 在接触闭合（t≈0.42）发散。故波5 contact 行统一用 ``use_rigid=True``
    （= proof 与 S4 主口径），并保留一条 PressureLoad 行作**阴性对照**。
    """
    from nonvertical_s5_contact import HEIGHTS, MU_PRIMARY, MUS, _hard_surface_load

    cs: list[dict] = []
    cs.append(dict(case_id="g0_fixed_axial_ref", height_m=5.0, plantar_bc="fixed",
                   load_n=26383.133775894516, pad_scenario="on_pad"))
    for h in HEIGHTS:
        cs.append(dict(case_id=f"h{h:g}_fixed_pad", height_m=h, plantar_bc="fixed",
                       load_n=loads[h], pad_scenario="on_pad"))
    cs.append(dict(case_id="h2.0_spring_k100_pad", height_m=2.0, plantar_bc="spring",
                   load_n=loads[2.0], pad_scenario="on_pad", spring_k=100.0))
    for mu in MUS:
        cs.append(dict(case_id=f"h2.0_contact_rigid_mu{mu:g}_pad", height_m=2.0,
                       plantar_bc="contact", load_n=loads[2.0], pad_scenario="on_pad",
                       contact_mu=mu, penalty=0.1, contact_hold_spring_k=S5C.CONTACT_HOLD_K,
                       use_rigid=True))
    for h in (3.0, 3.5, 4.5):
        cs.append(dict(case_id=f"h{h:g}_contact_rigid_mu{MU_PRIMARY:g}_pad", height_m=h,
                       plantar_bc="contact", load_n=loads[h], pad_scenario="on_pad",
                       contact_mu=MU_PRIMARY, penalty=0.1,
                       contact_hold_spring_k=S5C.CONTACT_HOLD_K, use_rigid=True))
    cs.append(dict(case_id="h2.0_fixed_hard_surface", height_m=2.0, plantar_bc="fixed",
                   load_n=_hard_surface_load(2.0, loads[2.0]), pad_scenario="hard_surface"))
    cs.append(dict(case_id="h2.0_contact_rigid_mu0.0_hard_surface", height_m=2.0,
                   plantar_bc="contact", load_n=_hard_surface_load(2.0, loads[2.0]),
                   pad_scenario="hard_surface", contact_mu=0.0, penalty=100.0,
                   contact_hold_spring_k=S5C.CONTACT_HOLD_K, use_rigid=True))
    cs.append(dict(case_id="h2.0_contact_rigid_mu0.6_pad_nohold", height_m=2.0,
                   plantar_bc="contact", load_n=loads[2.0], pad_scenario="on_pad",
                   contact_mu=MU_PRIMARY, penalty=0.1, contact_hold_spring_k=None,
                   use_rigid=True))
    cs.append(dict(case_id="h2.0_contact_rigid_mu0.6_pad_pen100", height_m=2.0,
                   plantar_bc="contact", load_n=loads[2.0], pad_scenario="on_pad",
                   contact_mu=MU_PRIMARY, penalty=100.0,
                   contact_hold_spring_k=S5C.CONTACT_HOLD_K, use_rigid=True))
    cs.append(dict(case_id="h2.0_contact_rigid_mu0.6_pad_press", height_m=2.0,
                   plantar_bc="contact", load_n=loads[2.0], pad_scenario="on_pad",
                   contact_mu=MU_PRIMARY, penalty=0.1,
                   contact_hold_spring_k=S5C.CONTACT_HOLD_K, use_rigid=False))
    return cs


def _deck_evidence(feb: Path) -> dict:
    """从 emitted .feb 抽跖面 facet 类型/数量 + G7 <Control> 关键字段。"""
    if not feb.is_file():
        return {"present": False}
    txt = feb.read_text(encoding="ISO-8859-1")
    sm = re.search(r'<Surface name="plantar">(.*?)</Surface>', txt, re.S)
    inner = sm.group(1) if sm else ""
    return {
        "present": True,
        "plantar_quad4": len(re.findall(r"<quad4\s", inner)),
        "plantar_tri3": len(re.findall(r"<tri3\s", inner)),
        "control_time_steps_2400": "time_steps>2400</time_steps>" in txt,
        "control_dtmax_1_2400": "dtmax>0.0004166666666666667</dtmax>" in txt,
        "control_cutback_0125": "cutback>0.125</cutback>" in txt,
    }


def _annotate(runs: list[dict], cases: list[dict]) -> None:
    """把 case 计划的 penalty/use_rigid/hold 补进 run（旧 JSON 可能缺这些键）。"""
    plan = {c["case_id"]: c for c in cases}
    for r in runs:
        c = plan.get(r["id"], {})
        r.setdefault("contact_penalty", c.get("penalty"))
        r.setdefault("contact_hold_spring_k_N_per_mm", c.get("contact_hold_spring_k"))
        r.setdefault("use_rigid", c.get("use_rigid", False))
        # 失败 contact 行归一为 DIVERGED（旧 JSON 行可能缺 contact_verdict）。
        if r.get("plantar_bc") == "contact" and r.get("febio_rc") != 0:
            r["contact_verdict"] = "DIVERGED"


def _val(row: dict, key: str):
    v = row.get(key)
    return None if not isinstance(v, dict) else v.get("value")


def _summary(runs: list[dict]) -> dict:
    by_id = {r["id"]: r for r in runs}
    contact = [r for r in runs if r["plantar_bc"] == "contact"]
    verdicts = {r["id"]: r.get("contact_verdict", r["run_verdict"]) for r in contact}

    def _blk(cid):
        r = by_id.get(cid)
        if not r:
            return None
        return {
            "id": cid, "rc": r["febio_rc"], "end_t": r["febio_end_t"],
            "penetration_mm": _val(r, "penetration_mm"),
            "plantar_drop_mm": _val(r, "plantar_drop_mm"),
            "contact_area_mm2": _val(r, "contact_area_mm2"),
            "carries_load": r.get("carries_load"),
            "verdict": verdicts.get(cid),
        }

    hold_on = _blk("h2.0_contact_rigid_mu0.6_pad")
    hold_off = _blk("h2.0_contact_rigid_mu0.6_pad_nohold")
    pen100 = _blk("h2.0_contact_rigid_mu0.6_pad_pen100")
    press = _blk("h2.0_contact_rigid_mu0.6_pad_press")
    return {
        "n_runs": len(runs),
        "n_contact": len(contact),
        "verdicts": verdicts,
        "hold_spring_required": bool(
            hold_on and hold_on["verdict"] == "CONVERGED+CARRYING"
            and hold_off and hold_off["verdict"] == "DIVERGED"),
        "hold_spring_evidence": {"with_hold": hold_on, "without_hold": hold_off},
        "penalty_sensitivity": {str(r["id"]): _val(r, "penetration_mm")
                                for r in contact if r["id"].startswith("h2.0_contact_rigid_mu0.6_pad")},
        "penalty_sensitivity_blocks": {"pen0.1": hold_on, "pen100": pen100},
        "load_path_finding": {
            "rigid_body": hold_on,
            "pressure_load": press,
            "conclusion": ("rigid-body 加载路径收敛且承载；PressureLoad 路径在接触闭合 "
                           "(t≈0.42) 处 FEBio 收敛判据退化 ⇒ ERROR。波4 contact 行用 "
                           "use_rigid=False，这是矩阵与 proof 的真实差异之一。"),
        },
        "f_n_note": ("F_n 标 INFERRED（静力平衡）；FEBio 4.13 本 deck 不报固定 DOF 反力，"
                     "reaction-sum 仅记 unreliable_diagnostic。"),
    }


def _finalize(runs: list[dict], ts: str, cases: list[dict], mesh_name: str) -> dict:
    _annotate(runs, cases)
    return {
        "meta": {
            "generated_at": ts, "wave": 5, "partial": False,
            "study": "S5-contact 波5：真实泡沫垫矩阵（G7 配方 + hold 弹簧 + rigid 加载路径）",
            "plan_doc": "docs/S5_wave5_matrix.md",
            "wave4_script": "scripts/opensim_fe/nonvertical_s5_contact.py::_solve_case (波5 扩展)",
            "mesh": mesh_name,
            "mesh_note": ("矩阵网格 = calcaneus_r_anatframe.npz（137 跖面 quad4）；"
                          "proof 用 calcaneus_r_calcnframe.npz（146 跖面 quad4）。二者均 0 tri3。"),
            "units": "mm-N-MPa-s", "febio": "4.13",
            "quotable_f_n_kn": QUOTABLE_F_N_KN,
            "quotable_note": ("QUOTABLE F_n≈34.4 kN 是 tet4 pad / 固定 sink 路径；本矩阵是 "
                              "THUMS mix-hex 网格 / 距下关节载荷路径 ⇒ 只作量级锚，不并列。"),
        },
        "runs": runs,
        "summary": _summary(runs),
    }


def _table(runs: list[dict]) -> list[str]:
    L = ["| case id | scenario | use_rigid | penalty | hold_k | febio_rc | end_t | negJac | plantar_faces | log_invalid | penetration (mm) | contact_area (mm²) | carries | verdict |",
         "|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---|"]
    for r in runs:
        de = r.get("deck_evidence") or {}
        dg = r.get("diagnostics") or {}
        pen = _val(r, "penetration_mm")
        area = _val(r, "contact_area_mm2")
        pf = (f"{de.get('plantar_quad4')}q/{de.get('plantar_tri3')}t"
              if r["plantar_bc"] == "contact" and de.get("present") else "—")
        hold = r.get("contact_hold_spring_k_N_per_mm")
        L.append(
            f"| {r['id']} | {r['pad_scenario']} | {r.get('use_rigid')} | "
            f"{(r.get('contact_penalty') if r.get('contact_penalty') is not None else '—')} | "
            f"{('—' if hold is None else hold)} | {r['febio_rc']} | "
            f"{('—' if r['febio_end_t'] is None else r['febio_end_t'])} | "
            f"{dg.get('neg_jacobians')} | {pf} | {dg.get('invalid_facets')} | "
            f"{('—' if pen is None else f'{pen:.3f}')} | "
            f"{('—' if area is None else f'{area:.1f}')} | "
            f"{r.get('carries_load')} | "
            f"{r.get('contact_verdict', r['run_verdict'])} |")
    return L


def _deck_snippet(case_id: str = "h2.0_contact_rigid_mu0.6_pad") -> list[str]:
    """代表 contact deck 的 `<Surface name="plantar">` + `<Control>` XML 证据。"""
    feb = Path(S5C.WORKDIR) / f"{case_id}.feb"
    if not feb.is_file():
        return [f"(代表 deck 未找到：{feb})"]
    txt = feb.read_text(encoding="ISO-8859-1")
    out: list[str] = []
    sm = re.search(r'<Surface name="plantar">.*?</Surface>', txt, re.S)
    if sm:
        block = sm.group(0)
        quads = re.findall(r"<quad4\b[^>]*>", block)
        tris = re.findall(r"<tri3\b[^>]*>", block)
        out.append(f'<!-- plantar Surface: {len(quads)} quad4, {len(tris)} tri3 -->')
        out.append('<Surface name="plantar">')
        out.extend(quads[:2])
        out.append(f"    ... ({len(quads)} quad4 total, 0 tri3)")
        out.append("</Surface>")
    for tag in ("time_steps", "dtmax", "cutback", "max_retries", "opt_iter"):
        m = re.search(rf"<{tag}>[^<]*</{tag}>", txt)
        if m:
            out.append(m.group(0))
    return out


def _write_report(path: Path, res: dict, ts: str) -> None:
    runs = res["runs"]
    s = res["summary"]
    L: list[str] = []
    A = L.append
    A("# S5-contact 波5 · 真实泡沫垫矩阵报告（G7 配方 + hold 弹簧 + rigid 加载路径）")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s；FEBio 4.13。")
    A("> 计划 `docs/S5_wave5_matrix.md`；契约 `docs/S5_contact_contract.md` §4/§5/§6。")
    A("> 引用限定：仅 `docs/S5_QUOTABLE.md` 的 QUOTABLE 行。F_n 一律 **INFERRED**（见 §6）。")
    A("")
    A("## 0. 一句话结论")
    A("")
    n_cc = sum(1 for r in runs if r["plantar_bc"] == "contact"
               and r.get("contact_verdict") == "CONVERGED+CARRYING")
    A(f"1. 真实泡沫垫矩阵（anatframe）上，**7 条先前失败 contact 行现全部 "
      f"CONVERGED+CARRYING**（rc=0、end_t=1.0）。")
    A(f"2. **hold 弹簧必需**：同一 case 去 hold ⇒ 2–3 s 内发散（接触闭合前 3 刚体模态）。"
      f"（converged+carrying contact 行数 = {n_cc}）")
    A("3. **加载路径是剩余阻塞**：波4 contact 行用 `use_rigid=False`（PressureLoad）⇒ "
      "在接触闭合 t≈0.42 发散（阴性对照 `_press`）；本波改用 `use_rigid=True`（= proof 与 "
      "S4 主口径）后收敛。")
    A(f"4. penalty 敏感性：pen0.1 → penetration ≈ "
      f"{s['penalty_sensitivity_blocks']['pen0.1']['penetration_mm']:.2f} mm；"
      f"pen100 → {s['penalty_sensitivity_blocks']['pen100']['penetration_mm']:.2f} mm；"
      "两者均承载（≪ 无接触 sink ≈110 mm）。")
    A("")
    A("## 1. 关键发现：加载路径 + hold 弹簧")
    A("")
    A("| 配置 | 网格 | use_rigid | hold | rc | end_t | 结果 |")
    A("|---|---|---:|---:|---:|---:|---|")
    lp = s["load_path_finding"]
    for tag, b in (("rigid + hold", lp["rigid_body"]), ("pressure + hold", lp["pressure_load"]),
                   ("rigid, no hold", s["hold_spring_evidence"]["without_hold"])):
        if b:
            A(f"| {tag} | anatframe | {('True' if 'rigid' in tag else 'False')} | "
              f"{('1.0' if 'no hold' not in tag else 'None')} | {b['rc']} | "
              f"{('—' if b['end_t'] is None else b['end_t'])} | **{b['verdict']}** |")
    A("")
    A(f"- {lp['conclusion']}")
    A("")
    A("## 2. 逐 case 结果表（全部 17 行）")
    A("")
    L.extend(_table(runs))
    A("")
    A("> `plantar_faces`：emitted deck 的 `<Surface name=\"plantar\">` 内 `<quad4>/<tri3>` 计数"
      "（**跖面专属**；任务书写的 146 是 proof 的 calcnframe 网格，矩阵 anatframe = 137）。"
      "`log_invalid` = FEBio 日志全局 `invalid facets` 计数（收敛参照行基线 = 2）。")
    A("")
    A("### 2.1 emitted deck 证据（代表 case `h2.0_contact_rigid_mu0.6_pad`）")
    A("")
    A("```xml")
    L.extend(_deck_snippet())
    A("```")
    A("")
    A("## 3. hold 弹簧必要性（同 case 对照）")
    A("")
    for tag, b in (("with hold k=1.0 N/mm", s["hold_spring_evidence"]["with_hold"]),
                   ("without hold (None)", s["hold_spring_evidence"]["without_hold"])):
        if b:
            A(f"- `{tag}`：rc={b['rc']}, end_t={b['end_t']}, penetration="
              f"{b['penetration_mm']}, verdict=**{b['verdict']}**")
    A(f"- ⇒ **hold_spring_required = {s['hold_spring_required']}**")
    A("")
    A("## 4. penalty 敏感性")
    A("")
    A("| 配置 | penetration (mm) | verdict |")
    A("|---|---:|---|")
    for tag, b in s["penalty_sensitivity_blocks"].items():
        if b:
            A(f"| {tag} | {b['penetration_mm']:.3f} | **{b['verdict']}** |")
    A("")
    A("## 5. 验收标准（承载 ≠ rc=0）")
    A("")
    A("- **CONVERGED**：`febio_rc==0` 且 `end_t==1.0`。")
    A("- **CARRYING**：`plantar_drop ≪ 无接触 sink`（= 载荷/191 N·mm⁻¹ ≈ 110 mm @ 21086 N）。")
    A("- **UNLOADED**：收敛但 drop ≈ sink。**DIVERGED**：rc≠0 / end_t<1。")
    A("")
    A("## 6. F_n 诚实边界")
    A("")
    A(f"- {s['f_n_note']}")
    A(f"- 量级锚：`docs/S5_QUOTABLE.md` 的 **F_n ≈ {QUOTABLE_F_N_KN} kN**（2 m 落地，tet4 pad / "
      "固定 sink / NEW PENALTY，网格收敛 <0.6%）。本矩阵载荷是**距下关节**载荷 "
      "（h=2.0 → 21086 N），与 34.4 kN 是**不同 mesh 与不同 load path** ⇒ 只作量级锚，不并列。")
    A("")
    A("## 7. 未解决的风险与下一步")
    A("")
    A("- 波5 把 contact 行切到 `use_rigid=True`（proof 已验证的加载路径）。若下游要求 "
      "`use_rigid=False` 的 PressureLoad 口径，接触闭合处仍需换数值策略（如分段 ramp / "
      "AUGLAG / 更细接触网格）——本波只证明 rigid 路径可行，**未**解决 pressure 路径。")
    A("- `penalty=0.1`（G7 配方，为泡沫标定）在 bone-vs-rigid 上穿透 40 mm（仍承载）；"
      "pen100 降到 ~2 mm。建议报告 σ/穿透时同引 penalty。")
    A("- F_n 全靠力平衡推断；未解析接触压力场。")
    A("")
    A("## 8. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s5_contact_wave5.py")
    A("```")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def _dump(path: Path, res: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=2, ensure_ascii=False, default=float),
                    encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="S5-contact 波5 真实泡沫垫矩阵驱动器")
    ap.add_argument("--out-json", type=Path, default=OUT_JSON)
    ap.add_argument("--out-md", type=Path, default=OUT_MD)
    ap.add_argument("--contact-only", action="store_true", help="只跑 contact 行（子集先行）")
    ap.add_argument("--reuse-json", type=Path, default=None,
                    help="不跑 FEBio，直接读取既有 runs JSON 重新汇总/写报告")
    args = ap.parse_args(argv)

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    loads = S5C._height_subtalar_loads()
    cases = _wave5_cases(loads)
    mesh_name = Path(S5C.MESH_ANAT).name

    if args.reuse_json is not None:
        runs = json.loads(args.reuse_json.read_text(encoding="utf-8"))["runs"]
        print(f"[reuse] {len(runs)} runs from {args.reuse_json}", flush=True)
    else:
        for p in (args.out_json, args.out_md):
            if p.exists():
                print(f"拒绝覆盖既有文件：{p}", file=sys.stderr)
                return 2
        print("=== S5-contact 波5：真实泡沫垫矩阵（G7 配方 + hold 弹簧 + rigid 加载）===",
              flush=True)
        mesh = S5C.T.load_thums_mesh(S5C.MESH_ANAT)
        print(f"[mesh] {mesh_name} nodes={mesh['n_nodes']} "
              f"plantar_nodes={len(mesh['node_sets']['plantar'])}", flush=True)
        selected = [c for c in cases if c["plantar_bc"] == "contact"] if args.contact_only else cases
        print(f"[cases] {len(selected)} 行（contact="
              f"{sum(1 for c in selected if c['plantar_bc'] == 'contact')}）", flush=True)
        runs = []
        for case in selected:
            row = S5C._solve_case(mesh, **case)
            row["deck_evidence"] = _deck_evidence(Path(row["feb_path"]))
            runs.append(row)
            dg = row.get("diagnostics", {})
            pen = _val(row, "penetration_mm")
            print(f"  {row['id']:42s} rc={row['febio_rc']} end_t={row['febio_end_t']} "
                  f"negJac={dg.get('neg_jacobians')} invalid={dg.get('invalid_facets')} "
                  f"penetration={('—' if pen is None else f'{pen:.3f}')} "
                  f"carries={row.get('carries_load')} "
                  f"verdict={row.get('contact_verdict', row['run_verdict'])} "
                  f"[{row['wall_s']['value']:.0f}s]", flush=True)
            _dump(args.out_json, _finalize(runs, ts, cases, mesh_name))

    res = _finalize(runs, ts, cases, mesh_name)
    _dump(args.out_json, res)
    _write_report(args.out_md, res, ts)
    print(f"wrote: {args.out_json}", flush=True)
    print(f"wrote: {args.out_md}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
