"""接触刚度**均匀压缩正对照**（R11 后续 / R5 规格落地）。

问题背景：`k` 这个量在本项目从未被定义清楚（旧口径 `F/δmax`=880 vs 对齐口径
`F·t/(A·δ̄)`=97.1 N/mm，差 9.1×），而"目标 1059"用的是第三方几何 ⇒ 任何 k 数字
都没有判据意义。本脚本提供**解析解已知的正对照**，把"测量机器 + 接触表述"
放到有真值的场景里验。

两个算例（均为 mm-N-MPa，与项目一致）：

A. **纯材料 patch test（无接触）** —— 验证测量机器
   10×10×3 mm 软骨块（E=12 MPa, ν=0.42, hex8 2×2×2），底面固定、
   顶面只约束 z（x/y 自由 ⇒ 侧向可外伸 ⇒ 单轴应力态）
   ⇒ **解析 k = E·A/t = 12×100/3 = 400.00 N/mm**（精确）

B. **接触串联正对照（主用例）** —— 同软骨块 + 上面加刚性压板
   （10×10×2 mm, E=200 GPa），界面 `sliding-elastic` **frictionless**（μ=0）
   ⇒ 仍是单轴应力态；压板顶面只约束 z（可侧移）
   解析：δ_total = δ_contact + δ_slab，δ_slab = F·t/(E·A)（**已知**）
        ⇒ **k_contact = F / (δ_total − δ_slab)**（把接触柔度单独隔离出来）
   扫不同接触表述（laugon / penalty / auto_penalty）⇒ 直接测出各表述的
   **有效接触刚度**，与"罚参数就是 k"的假设直接对撞。

用法:
  .venv\\Scripts\\python.exe -u scripts/ankle_fe/contact_patch_test.py --case A
  .venv\\Scripts\\python.exe -u scripts/ankle_fe/contact_patch_test.py --case B
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import textwrap
import time
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
OUT = ROOT / "temp/ankle_fe/patch"
RES = ROOT / "results/opensim_fe/contact_patch_test.json"

# ---- 几何/材料（与踝线一致：E=12 MPa 软骨、ν=0.42、总厚 3 mm）----------------
LX = LY = 10.0
T_SLAB = 3.0
T_PLAT = 2.0
E_CART, NU_CART = 12.0, 0.42
E_PLAT = 200000.0
DISP = 0.05          # mm，压板顶面 z 位移
K_ANALYTIC = E_CART * (LX * LY) / T_SLAB          # = 400.00 N/mm

# 接触表述扫描：(tag, laugon, penalty, auto_penalty)
CONTACT_VARIANTS = [
    ("pen_auto1_p1", "PENALTY", 1.0, 1),      # 踝线生产语义
    ("pen_auto0_p1", "PENALTY", 1.0, 0),      # 项目铁律语义
    ("pen_auto0_p8", "PENALTY", 8.0, 0),
    ("pen_auto0_p32", "PENALTY", 32.0, 0),
    ("aug_auto0_p2", "AUGLAG", 2.0, 0),
    ("aug_auto0_p8", "AUGLAG", 8.0, 0),
    ("aug_auto0_p16", "AUGLAG", 16.0, 0),
]


def hex_block(nx: int, ny: int, nz: int, lx: float, ly: float, lz: float,
              z0: float, nid0: int, eid0: int):
    """结构化 hex8 块。返回 (nodes[(id,x,y,z)], elems[(id,n1..n8)],
    top_quads, bot_quads, top_ns, bot_ns)。"""
    xs = np.linspace(0, lx, nx + 1)
    ys = np.linspace(0, ly, ny + 1)
    zs = np.linspace(z0, z0 + lz, nz + 1)
    idx = {}
    nodes = []
    n = nid0
    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                idx[(i, j, k)] = n
                nodes.append((n, float(x), float(y), float(z)))
                n += 1
    elems = []
    e = eid0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                elems.append((e, idx[(i, j, k)], idx[(i + 1, j, k)],
                              idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)],
                              idx[(i, j, k + 1)], idx[(i + 1, j, k + 1)],
                              idx[(i + 1, j + 1, k + 1)], idx[(i, j + 1, k + 1)]))
                e += 1

    def face(k, rev: bool):
        quads = []
        q = 1
        for j in range(ny):
            for i in range(nx):
                c = [idx[(i, j, k)], idx[(i + 1, j, k)],
                     idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)]]
                if rev:
                    c = c[::-1]
                quads.append((q, c))
                q += 1
        return quads

    top_q, bot_q = face(nz, False), face(0, True)
    top_ns = [idx[(i, j, nz)] for j in range(ny + 1) for i in range(nx + 1)]
    bot_ns = [idx[(i, j, 0)] for j in range(ny + 1) for i in range(nx + 1)]
    return nodes, elems, top_q, bot_q, top_ns, bot_ns, idx


GLOBALS = ("  <Globals><Constants><T>0</T><R>0</R><Fc>0</Fc></Constants>"
            "</Globals>")
EMPTY_LOADS = "  <Loads>\n  </Loads>"


def _lc(points=((0.0, 0.0), (1.0, 1.0))) -> str:
    pts = "\n".join(f"        <pt>{a},{b}</pt>" for a, b in points)
    return (f'    <load_controller id="1" type="loadcurve">\n'
            f"      <interpolate>LINEAR</interpolate>\n"
            f"      <extend>CONSTANT</extend>\n      <points>\n{pts}\n"
            f"      </points>\n    </load_controller>")


def emit_case_a(path: Path, nu: float = NU_CART,
                fix_dofs: str = "xyz") -> dict:
    """纯材料 patch test。
    解析（**底面只固定 z、顶面只约束 z ⇒ 侧向自由 ⇒ 单轴应力态**）：
        F = E·A·(δ/t)，k = E·A/t = 400 N/mm（精确，与 ν 无关）
    若底面把 x/y 也夹住，则底面 Poisson 外伸被抑制 ⇒ k 偏高 ⇒ 用作**机制判别**。
    """
    nodes, elems, top_q, bot_q, top_ns, bot_ns, _ = hex_block(
        2, 2, 2, LX, LY, T_SLAB, 0.0, 1, 1)
    L = [f"<?xml version='1.0' encoding='ISO-8859-1'?>", "<febio_spec version='4.0'>",
         "  <Module type='solid'/>", GLOBALS, "  <Material>",
         "    <material name='cart' id='1' type='isotropic elastic'>",
         f"      <density>1e-06</density><E>{E_CART}</E><v>{nu}</v>",
         "    </material>", "  </Material>", "  <Mesh>",
         "    <Nodes name='nodes'>"]
    L += [f"      <node id='{i}'>{x},{y},{z}</node>" for i, x, y, z in nodes]
    L += ["    </Nodes>", "    <Elements name='slab' type='hex8'>"]
    L += [f"      <elem id='{e}'>{','.join(str(v) for v in el)}</elem>"
          for e, *el in elems]
    L += ["    </Elements>",
          f"    <NodeSet name='top'>{','.join(map(str, top_ns))}</NodeSet>",
          f"    <NodeSet name='bot'>{','.join(map(str, bot_ns))}</NodeSet>",
          "  </Mesh>", "  <MeshDomains>",
          "    <SolidDomain name='slab' mat='cart'/>", "  </MeshDomains>",
          "  <LoadData>", _lc(), "  </LoadData>", EMPTY_LOADS, "  <Boundary>",
          "    <bc name='fix' type='zero displacement' node_set='bot'>",
          f"      <x_dof>{1 if 'x' in fix_dofs else 0}</x_dof>"
          f"<y_dof>{1 if 'y' in fix_dofs else 0}</y_dof><z_dof>1</z_dof>",
          "    </bc>",
          "    <bc name='push' type='prescribed displacement' node_set='top'>",
          f"      <dof>z</dof><value lc='1'>{-DISP}</value><relative>0</relative>",
          "    </bc>", "  </Boundary>",
          "  <Output>",
          "    <plotfile type='febio'><var type='displacement'/>"
          "<var type='reaction forces'/></plotfile>", "  </Output>",
          _control(), "</febio_spec>"]
    path.write_text("\n".join(L), encoding="latin-1")
    return {"nodes": len(nodes), "elems": len(elems)}


def emit_case_b(path: Path, laugon: str, penalty: float, auto: int,
                nu: float = NU_CART) -> dict:
    """软骨块 + 刚性压板 + frictionless 接触（串联 ⇒ 可隔离接触柔度）。"""
    sn, se, s_top_q, s_bot_q, s_top_ns, s_bot_ns, sidx = hex_block(
        2, 2, 2, LX, LY, T_SLAB, 0.0, 1, 1)
    nid = max(n for n, *_ in sn) + 1
    eid = max(e for e, *_ in se) + 1
    pn, pe, p_bot_q, p_top_q, p_bot_ns, p_top_ns, pidx = hex_block(
        2, 2, 1, LX, LY, T_PLAT, T_SLAB, nid, eid)

    L = [f"<?xml version='1.0' encoding='ISO-8859-1'?>", "<febio_spec version='4.0'>",
         "  <Module type='solid'/>", GLOBALS, "  <Material>",
         "    <material name='cart' id='1' type='isotropic elastic'>",
         f"      <density>1e-06</density><E>{E_CART}</E><v>{nu}</v>",
         "    </material>",
         "    <material name='platen' id='2' type='isotropic elastic'>",
         f"      <density>1e-06</density><E>{E_PLAT}</E><v>0.3</v>",
         "    </material>", "  </Material>", "  <Mesh>",
         "    <Nodes name='nodes'>"]
    L += [f"      <node id='{i}'>{x},{y},{z}</node>" for i, x, y, z in sn + pn]
    L += ["    </Nodes>"]
    L += ["    <Elements name='slab' type='hex8'>"]
    L += [f"      <elem id='{e}'>{','.join(str(v) for v in el)}</elem>"
          for e, *el in se]
    L += ["    </Elements>", "    <Elements name='platen' type='hex8'>"]
    L += [f"      <elem id='{e}'>{','.join(str(v) for v in el)}</elem>"
          for e, *el in pe]
    L += ["    </Elements>",
          f"    <NodeSet name='s_top'>{','.join(map(str, s_top_ns))}</NodeSet>",
          f"    <NodeSet name='s_bot'>{','.join(map(str, s_bot_ns))}</NodeSet>",
          f"    <NodeSet name='p_top'>{','.join(map(str, p_top_ns))}</NodeSet>",
          "    <Surface name='slab_top'>"]
    L += [f"      <quad4 id='{q}'>{','.join(map(str, c))}</quad4>"
          for q, c in s_top_q]
    L += ["    </Surface>", "    <Surface name='platen_bot'>"]
    L += [f"      <quad4 id='{q}'>{','.join(map(str, c))}</quad4>"
          for q, c in p_bot_q]
    L += ["    </Surface>",
          "    <SurfacePair name='cart_platen_pair'>",
          "      <primary>platen_bot</primary>",
          "      <secondary>slab_top</secondary>", "    </SurfacePair>",
          "  </Mesh>", "  <MeshDomains>",
          "    <SolidDomain name='slab' mat='cart'/>",
          "    <SolidDomain name='platen' mat='platen'/>", "  </MeshDomains>",
          "  <LoadData>", _lc(), "  </LoadData>", EMPTY_LOADS, "  <Boundary>",
          "    <bc name='fix' type='zero displacement' node_set='s_bot'>",
          "      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>", "    </bc>",
          "    <bc name='push' type='prescribed displacement' node_set='p_top'>",
          f"      <dof>z</dof><value lc='1'>{-DISP}</value><relative>0</relative>",
          "    </bc>", "  </Boundary>", "  <Contact>",
          "    <contact name='cp' surface_pair='cart_platen_pair' "
          "type='sliding-elastic'>",
          f"      <laugon>{laugon}</laugon>", f"      <penalty>{penalty:g}</penalty>",
          f"      <auto_penalty>{auto}</auto_penalty>", "      <two_pass>0</two_pass>",
          "      <node_reloc>0</node_reloc>", "      <symmetric_stiffness>0</symmetric_stiffness>",
          "      <tolerance>0.005</tolerance>", "      <search_tol>0.01</search_tol>",
          "      <search_radius>20.0</search_radius>", "      <fric_coeff>0.0</fric_coeff>",
          "      <minaug>0</minaug>", "      <maxaug>20</maxaug>", "    </contact>",
          "  </Contact>", "  <Output>",
          "    <plotfile type='febio'><var type='displacement'/>"
          "<var type='reaction forces'/></plotfile>", "  </Output>",
          _control(), "</febio_spec>"]
    path.write_text("\n".join(L), encoding="latin-1")
    return {"nodes": len(sn) + len(pn), "elems": len(se) + len(pe),
            "ns_slab_top": s_top_ns, "ns_platen_top": p_top_ns}


def _control(time_steps: int = 20) -> str:
    # ⚠️ 逐字照抄**已知可跑通**的踝线 deck（temp/thums/febio/ankle_tibiotalar.feb）
    # 的 Control 段——自拟的 solver 标签会让 FEBio 4.13 直接访问违例崩溃
    # (0xC0000005，无可读报错)，故此处不做"优化"，只改步数。
    return textwrap.dedent(f"""\
      <Control>
        <analysis>STATIC</analysis><time_steps>{time_steps}</time_steps><step_size>{1.0/time_steps:.6f}</step_size>
        <plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range>
        <plot_level>PLOT_MAJOR_ITRS</plot_level><output_level>OUTPUT_MAJOR_ITRS</output_level>
        <plot_stride>1</plot_stride><output_stride>1</output_stride>
        <adaptor_re_solve>1</adaptor_re_solve>
        <time_stepper type="default"><max_retries>10</max_retries><opt_iter>15</opt_iter><cutback>0.25</cutback></time_stepper>
        <solver type="solid"><linear_solver type="pardiso"/><max_refs>25</max_refs><reform_each_time_step>1</reform_each_time_step></solver>
      </Control>""")


# ---- pyfebio 读取（沿用项目已验证的 patch 手法）--------------------------
_PYFEBIO_PATCHED = False


def _patch_pyfebio():
    """pyfebio 的 parse_mesh/parse_state 对缺 name 的集合会 KeyError（项目已踩过）。
    patch 必须**幂等**：一旦替换过，`inspect.getsource` 就再也拿不到源码 ⇒ 只做一次。
    """
    global _PYFEBIO_PATCHED
    import inspect

    import pyfebio.xplt as X
    if _PYFEBIO_PATCHED:
        return X
    _PYFEBIO_PATCHED = True
    src = textwrap.dedent(inspect.getsource(X.parse_mesh))
    for o, n in (('domain["name"]', 'domain.get("name", str(domain["id"][0]))'),
                 ('surface["name"]', 'surface.get("name", str(surface["id"][0]))'),
                 ("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
                 ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
                 ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
                 ("{nodeset['name']}", "{nodeset.get('name', str(nodeset['id'][0]))}"),
                 ("{elementset['name']}",
                  "{elementset.get('name', str(elementset['id'][0]))}")):
        src = src.replace(o, n)
    ns = X.__dict__
    exec(compile(src, "<p>", "exec"), ns)
    X.parse_mesh = ns["parse_mesh"]
    src2 = textwrap.dedent(inspect.getsource(X.parse_state))
    for a, b in (
            ('set_name = mesh_dict[set_lut[key]][set_id + 1]',
             'set_name = mesh_dict[set_lut[key]].get(set_id + 1) or '
             'mesh_dict[set_lut[key]].get(int(set_id) + 1) or f"set{set_id}"'),
            ('set_name = mesh_dict[set_lut[key]][set_id]',
             'set_name = mesh_dict[set_lut[key]].get(set_id) or '
             'mesh_dict[set_lut[key]].get(int(set_id)) or f"set{set_id}"')):
        src2 = src2.replace(a, b)
    ns2 = X.__dict__
    exec(compile(src2, "<p2>", "exec"), ns2)
    X.parse_state = ns2["parse_state"]
    return X


def read_forces(case_dir: Path, ns_slab_top: list[int] | None = None,
                ns_platen_top: list[int] | None = None) -> tuple[float, float, float]:
    """返回 (t=1 反力 z 合力, 压板顶面平均 uz, 软骨顶面平均 uz)。

    ⚠️ 接触柔度用**实测**分离：δ_contact = uz_platen_top − uz_slab_top
    （不用 F·t/(E·A) 假定期ointe，底面横向夹持会让解析值偏离实测 ~35%，见 Case A）。
    """
    import h5py

    X = _patch_pyfebio()
    xp = case_dir / "model.xplt"
    h5 = case_dir / "model.h5"
    if not xp.exists():
        return float("nan"), float("nan"), float("nan")
    X.to_hdf5(str(xp), str(h5))
    with h5py.File(h5, "r") as f:
        ns = sorted((int(q) for q in f["states"]), key=int)
        last = ns[-1]
        nd = f[f"states/{last}/node_data"]

        def grab(qty: str) -> np.ndarray:
            g = nd[qty]
            key = "nodes" if "nodes" in g else "all"
            return np.asarray(g[key][()], float)

        R, U = grab("reaction forces"), grab("displacement")
    rz = float(R[:, 2].sum())                     # 底面约束反力 = 施加载荷

    def mean_uz(ids: list[int] | None) -> float:
        if not ids:
            return float("nan")
        return float(np.mean([U[i - 1, 2] for i in ids]))   # node id → 0-based

    return rz, mean_uz(ns_platen_top), mean_uz(ns_slab_top)


def run(tag: str, emit) -> dict:
    d = OUT / tag
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    deck = d / "model.feb"
    meta = emit(deck) or {}
    t0 = time.time()
    r = subprocess.run([str(FEBIO), "-i", "model.feb"], cwd=str(d),
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    (d / "febio.log").write_text(out, encoding="utf-8")
    log = out.replace("  ", " ")
    bad = re.search(r'tag "([^"]+)"[^\n]*unrecognized tag', out)
    F = uz_p = uz_s = float("nan")
    if r.returncode == 0 and "N O R M A L" in log and not bad:
        try:
            F, uz_p, uz_s = read_forces(d, meta.get("ns_slab_top"),
                                        meta.get("ns_platen_top"))
        except Exception as exc:                      # noqa: BLE001
            (d / "read_error.txt").write_text(repr(exc), encoding="utf-8")
    return {"tag": tag, "rc": r.returncode,
            "normal": "N O R M A L" in log, "bad_tag": bad.group(1) if bad else None,
            "log_kb": round(len(out) / 1024, 1),
            "F_z": F, "uz_platen": uz_p, "uz_slab_top": uz_s,
            "wall_s": round(time.time() - t0, 1)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="AB", choices=["A", "B", "AB"])
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    if args.case in ("A", "AB"):
        # 三个**解析已知**的机制判别变体（底面夹持方式 / ν）：
        #   A_xyz_nu42 : 底面 xyz 全夹 + ν=0.42 → 底面 Poisson 外伸被抑制 ⇒ k 偏高
        #   A_xyz_nu0  : 同上但 ν=0（无横向耦合）⇒ 仍应精确 = E·A/t = 20 N
        #   A_z_nu42   : 底面只夹 z（侧向全自由）⇒ **精确** = E·A/t = 20 N
        rows.append(run("caseA_xyz_nu0.42",
                        lambda p: emit_case_a(p, 0.42, "xyz")))
        rows.append(run("caseA_xyz_nu0", lambda p: emit_case_a(p, 0.0, "xyz")))
        rows.append(run("caseA_zonly_nu0.42",
                        lambda p: emit_case_a(p, 0.42, "z")))
    if args.case in ("B", "AB"):
        for tag, lg, pen, au in CONTACT_VARIANTS:
            rows.append(run(f"caseB_{tag}",
                            lambda p, lg=lg, pen=pen, au=au: emit_case_b(p, lg, pen, au)))

    for r in rows:
        r["k_meas"] = (abs(r["F_z"]) / DISP) if r["F_z"] == r["F_z"] else float("nan")
        if r["tag"].startswith("caseB") and r["k_meas"] == r["k_meas"] \
                and r["uz_slab_top"] == r["uz_slab_top"]:
            # 接触穿透 = 压板顶面位移 − 软骨顶面位移（**实测之差**，不用解析减法）
            r["d_contact"] = abs(r["uz_platen"]) - abs(r["uz_slab_top"])
            r["k_contact"] = (abs(r["F_z"]) / r["d_contact"]
                              if r["d_contact"] > 1e-9 else float("nan"))
            r["k_slab_fe"] = (abs(r["F_z"]) / abs(r["uz_slab_top"])
                              if r["uz_slab_top"] else float("nan"))

    print("\n" + "=" * 92)
    print(f"解析：k_slab = E·A/t = {E_CART}×{LX*LY:.0f}/{T_SLAB} = {K_ANALYTIC:.2f} N/mm"
          f"   （Case A 应精确命中）")
    print("=" * 92)
    print(f"{'tag':>18} {'rc':>3} {'NORM':>5} {'F_z(N)':>11} {'uz_plat':>9} "
          f"{'uz_slab':>9} {'k_meas':>9} {'vs 400':>7} {'k_contact':>11} {'logKB':>6}")
    for r in rows:
        kc = r.get("k_contact", float("nan"))
        print(f"{r['tag']:>18} {r['rc']:>3} {str(r['normal']):>5} "
              f"{r['F_z']:>11.4f} {r['uz_platen']:>9.5f} {r['uz_slab_top']:>9.5f} "
              f"{r['k_meas']:>9.2f} {r['k_meas']/K_ANALYTIC*100:>6.1f}% "
              f"{kc:>11.2f} {r['log_kb']:>6.0f}")
    RES.parent.mkdir(parents=True, exist_ok=True)
    RES.write_text(json.dumps({"disp_mm": DISP, "k_analytic": K_ANALYTIC,
                               "rows": rows,
                               "generated": time.strftime("%Y-%m-%d %H:%M:%S")},
                              indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\njson -> {RES}")


if __name__ == "__main__":
    main()