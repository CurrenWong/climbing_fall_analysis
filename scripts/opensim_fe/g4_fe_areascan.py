"""G4 ② · A：载荷面积扫描（恒定 1000 N，只变加载面积）

──────────────────────────────────────────────────────────────────────────
预注册预测（**先写后跑**；不得在看过结果后改判据或调面积阶梯）
──────────────────────────────────────────────────────────────────────────
目标量：``σ_vm_max``（主域 ``fib_cort``，MPa）。**固定** F=1000 N、BC（bottom 30%
固定）、方向（沿骨长轴 PCA、压缩）、网格、材料；**只变加载面积 A**。

物理参照：皮质截面 A≈102 mm² ⇒ 名义 F/A ≈ **9.8 MPa**（§9.8）。判据目标带 =
**10–30 MPa**（§9.11 的探针判据）。

两个互斥假设（本扫描即判别实验）：
  * **H_pressure**：应力由「载荷面积」主导（接触/压载压力 σ ≈ C/A）
    ⇒ σ(A) 单调陡降，在 A ≳ 100 mm² 处可落入 10–30 MPa；
  * **H_bend**：应力由「细长骨弯曲 + BC 过渡」主导
    ⇒ σ(A) ≈ C/A + **σ₀**（面积无关的地板）⇒ 即使 A 到 5000 mm² 也停在 σ₀。

判据方式：对 traction 阶梯最小二乘拟合 ``σ = C/A + σ₀``，报告 **C 与 σ₀**
（以及仅看陡降/饱和的形状），再回答「10–30 MPa 是否可达」。

⚠️ 若结果指向 H_bend，**不得**「缩面积到好看为止」——那正是 G4 ② 红线
（以调参/换判据掩盖 BC 缺陷）。正确结论是：判据 / 载荷口径 / BC 三者之一需重审。

──────────────────────────────────────────────────────────────────────────
方法
──────────────────────────────────────────────────────────────────────────
1. 几何：``facet_geom(top_frac=0.30, align_pca=True)`` 取骨端外表面（band 面片）,
   apex = 沿轴最高节点（对应刚体压板接触的那个圆顶）；
2. 面片按「到 apex 的**测地距离**（面片图 Dijkstra、边权=边长）均值」自近及远排序；
   k 个面片 = 一个绕 apex 的同心补丁 ⇒ 面积 A = Σ 面片面积（真实曲面面积）；
3. 载荷：``<surface_load type="traction">``，矢量 = ``-t·axis``（**轴向压缩**）、
   ``t = 1000/A`` ⇒ **轴向合力恰为 1000 N**（traction 是按参考面积给的矢量，
   与面片朝向无关；不像 pressure 会因斜端面产生横向分量 = §9.8 的 V1 发散嫌疑）；
4. 另有 1 个 **点载荷锚点**（单节点 1000 N 节点力，= P1-a/V0 的极端机制），
   用于把「机制（点 vs 均布）」与「面积」分离；
5. 每级都加 **ramp**（``<scale lc="1">``，0→1）与 reaction 输出；
6. 产物记账：``R_axis``（**必须 ≈ ±1000 N** —— 静力平衡；这是「载荷被静默丢弃」
   的第一判据）、``n_loaded``、``max|u|``（空空间/发散哨兵）、峰值单元位置
   （在载荷补丁内？在 bottom_fix？在骨干？）+ 补丁内局部峰值 ``vm_patch_max``。

索引约定（**已用正对照验证**，2026-10-12）：
  * Python 侧 ``facet_geom()['cap'/'band'][i]['poly']`` 是 **1-based 全局 id**；
  * h5 ``meshes/0/domains/<dom>`` 连接表是 **0-based 位置**（= id−1）；
  * h5 ``meshes/0/nodes[id-1]`` = 该 id 的坐标；``nodesets/*`` 也是 0-based 位置。
  脚本启动时做 **集合级正对照**（``conn+1`` 的 id 集合必须等于 Python 侧 hex 单元
  的 id 集合），不一致立即报错。

复现（repo root）::

    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_areascan.py

只 ADD 本模块；复用 ``g4_fe_p1b`` / ``bone_feb`` / ``g4_fe_stab``，不改既有模块。
"""
from __future__ import annotations

import argparse
import heapq
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import bone_feb as bf  # noqa: E402
import g4_fe_p1b as p1b  # noqa: E402
import g4_fe_stab as gs  # noqa: E402

OUT_TEMP = ROOT / "temp" / "opensim_fe" / "g4_fe_areascan"
OUT_RESULTS = ROOT / "results" / "opensim_fe"

LOAD_N = 1000.0
TOP_FRAC = 0.30
BOT_FRAC = 0.30
MAIN_DOM = "fib_cort"

#: 默认阶梯：``"point"`` = 单节点 1000 N（点载荷锚点）；整数 = 补丁面片数 k。
DEFAULT_LEVELS = ("point", 1, 3, 5, 8, 13, 17, 28, 52, 95, 173, 291)


# ===========================================================================
# 几何：测地距离 / 支面积 / 面片阶梯
# ===========================================================================
def geodesic_from(nodes: np.ndarray, polys, src_id: int) -> dict[int, float]:
    """面片图上的测地距离（Dijkstra，边权 = 边长）。polys 为 1-based id 元组。"""
    adj: dict[int, list[tuple[int, float]]] = {}
    for poly in polys:
        ids = [int(x) for x in poly]
        m = len(ids)
        for k in range(m):
            a, b = ids[k], ids[(k + 1) % m]
            d = float(np.linalg.norm(nodes[a - 1] - nodes[b - 1]))
            adj.setdefault(a, []).append((b, d))
            adj.setdefault(b, []).append((a, d))
    dist = {src_id: 0.0}
    pq = [(0.0, src_id)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, np.inf) + 1e-12:
            continue
        for v, w in adj.get(u, ()):
            nd = d + w
            if nd < dist.get(v, np.inf):
                dist[v] = nd
                heapq.heappush(pq, (nd, v))
    return dist


def tributary_areas(facets: list[dict]) -> dict[int, float]:
    """节点 → 支面积（Σ 相邻面片面积 / 该面片顶点数）。Σ = 面片总面积。"""
    trib: dict[int, float] = {}
    for f in facets:
        ids = [int(n) for n in f["poly"]]
        for n in ids:
            trib[n] = trib.get(n, 0.0) + float(f["area"]) / len(ids)
    return trib


def facet_ladder(nodes: np.ndarray, facets: list[dict], apex: int):
    """按「面片节点测地距离均值」排序 ⇒ 绕 apex 的同心补丁序列 + 累积面积。"""
    dist = geodesic_from(nodes, [f["poly"] for f in facets], apex)
    order = sorted(
        facets,
        key=lambda f: float(np.mean([dist.get(int(n), np.inf) for n in f["poly"]])),
    )
    cum = np.cumsum([float(f["area"]) for f in order])
    cum_proj = np.cumsum([float(f["area"]) * float(f["axdiv"]) for f in order])
    return order, cum, cum_proj, dist


def k_for_target(cum: np.ndarray, target: float) -> int:
    k = int(np.searchsorted(cum, target) + 1)
    return min(max(k, 1), len(cum))


# ===========================================================================
# .feb 注入
# ===========================================================================
def inject_fragment(text: str, block: str) -> str:
    """把一段 XML 注入 <Mesh> 尾部（与 p1b.inject_surface 用同一唯一标记）。"""
    if text.count("  </Mesh>") != 1:
        raise RuntimeError("inject_fragment: 找不到唯一的 '  </Mesh>'")
    return text.replace("  </Mesh>", block + "  </Mesh>", 1)


def inject_boundary(text: str, block: str) -> str:
    """把一段 BC 注入 step 的 <Boundary> 尾部（C：横向支撑）。"""
    marker = "      </Boundary>"
    if text.count(marker) != 1:
        raise RuntimeError("inject_boundary: 找不到唯一的 '</Boundary>'")
    return text.replace(marker, block + marker, 1)


def lateral_bc_xml(node_set: str) -> str:
    """横向（x/y）零位移支撑 —— 骨间膜 / 胫骨侧向约束的**刚性上界**理想化。"""
    return (f'        <bc type="zero displacement" node_set="{node_set}">\n'
            '          <x_dof>1</x_dof>\n'
            '          <y_dof>1</y_dof>\n'
            '        </bc>\n')


def brace_ids_for(nodes: np.ndarray, geom: dict, lateral, patch_ids,
                  bot_frac: float) -> list[int]:
    """按 z_frac 区间（或 "patch"）取横向支撑节点；始终排除载荷补丁节点。"""
    ax_i = int(geom["axis_idx"])
    z_lo, z_hi = float(geom["a_min"]), float(geom["a_max"])
    zf = (nodes[:, ax_i] - z_lo) / (z_hi - z_lo)
    ps = {int(x) for x in patch_ids}
    if lateral == "patch":
        return sorted(ps)
    lo, hi = float(lateral[0]), float(lateral[1])
    sel = [i + 1 for i in np.flatnonzero((zf >= lo) & (zf <= hi))]
    return sorted(int(i) for i in sel if int(i) not in ps)


def nodeset_xml(name: str, ids) -> str:
    return ('    <NodeSet name="%s">%s</NodeSet>\n'
            % (name, ",".join(str(int(i)) for i in ids)))


def nodal_load_xml(node_set: str, vec) -> str:
    v = f"{float(vec[0]):.10g},{float(vec[1]):.10g},{float(vec[2]):.10g}"
    return ('      <Loads>\n'
            f'        <nodal_load type="nodal_force" node_set="{node_set}">\n'
            '          <relative>0</relative>\n'
            f'          <value>{v}</value>\n'
            '        </nodal_load>\n'
            '      </Loads>')


#: 弹簧"接地面"的横向偏移（mm）。**必须远大于位移**：地面节点在 p+(d,0,0)，
#: 骨节点移动 (ux,uy,uz) 后弹簧长度 ≈ d − ux + (uy²+uz²)/(2d) ⇒
#: ① 第一项 = 想要的 x 向刚度；② 第二项 = **横向位移平方引起的张力**
#: （力 = k·u²/(2d)）。2026-10-12 实测：d=1e-3 mm 时 k=1e6 的弹簧变成刚性系绳，
#: 把骨轴向夹死（vm 0.1 MPa、|u|≈0，正对照失败）。取 d=1000 mm ⇒
#: k·u²/(2d) 在 k=1e6、u=0.3 mm 时仅 ≈0.05 N，可忽略。
SPRING_OFF_MM = 1000.0


def spring_blocks(brace_ids, nodes: np.ndarray, k: float, first_node_id: int):
    """分布式侧向**弹簧**支撑（软支撑）：每节点 2 根线性弹簧 → 地面节点。

    FEBio 弹簧是 2 节点单元（``<delem>n1,n2</delem>``，沿两节点连线），所以
    「接地面」用「远处的地面节点」实现：地面节点位置 = 骨节点 + (d,0,0) /
    (0,d,0)，d = ``SPRING_OFF_MM`` ⇒ 分别提供 x / y 向刚度 k（见上面的偏移分析）。
    地面节点全部零位移约束。返回 (nodes_xml, mesh_xml, discrete_xml, ground_ids)。
    """
    off = float(SPRING_OFF_MM)
    ax_i = None
    node_lines = ['    <Nodes name="lat_ground">']
    fix_ids: list[int] = []
    delems: list[str] = []
    nid = int(first_node_id)
    for b in brace_ids:
        p = np.asarray(nodes[int(b) - 1], dtype=float)
        for d in (0, 1):
            q = p.copy(); q[d] += off
            node_lines.append(f'      <node id="{nid}">'
                              f'{q[0]:.10e},{q[1]:.10e},{q[2]:.10e}</node>')
            fix_ids.append(nid)
            delems.append(f'      <delem>{int(b)},{nid}</delem>')
            nid += 1
    node_lines.append("    </Nodes>")
    mesh_xml = (
        "    <NodeSet name=\"lat_ground_fix\">"
        + ",".join(str(i) for i in fix_ids) + "</NodeSet>\n"
        + "    <DiscreteSet name=\"lat_springs\">\n"
        + "\n".join(delems) + "\n    </DiscreteSet>\n")
    discrete_xml = (
        "  <Discrete>\n"
        "    <discrete_material id=\"1\" name=\"lat_spring\" "
        "type=\"linear spring\">\n"
        f"      <E>{float(k):.10g}</E>\n"
        "    </discrete_material>\n"
        "    <discrete dmat=\"1\" discrete_set=\"lat_springs\"/>\n"
        "  </Discrete>\n")
    return "\n".join(node_lines) + "\n", mesh_xml, discrete_xml, fix_ids


def ground_bc_xml(node_set: str = "lat_ground_fix") -> str:
    return (f'        <bc type="zero displacement" node_set="{node_set}">\n'
            '          <x_dof>1</x_dof>\n'
            '          <y_dof>1</y_dof>\n'
            '          <z_dof>1</z_dof>\n'
            '        </bc>\n')


def build_level_feb(setup, out_path: Path, *, level, geom, order, cum, cum_proj,
                    apex: int, trib: dict, steps: int = 5,
                    load_n: float = LOAD_N, bot_frac: float = BOT_FRAC,
                    nodes: np.ndarray | None = None, lateral=None,
                    spring_k: float | None = None) -> dict:
    """为一级（"point" 或 k 个面片）构建 .feb，返回该级的载荷记账。

    ``bot_frac`` 只改**固定带长度**（⇒ 自由长度 L = span·(1−bot_frac)），
    载荷补丁仍由 ``facet_geom(top_frac=0.30)`` 决定 ⇒ 可把「自由长度」与
    「加载区」两个因子分离开（B：BC 灵敏度扫描）。

    ``lateral``（C：横向支撑）：``None`` = 无；``"patch"`` = 只约束载荷补丁的 x/y
    （= 导向压头）；``(lo, hi)`` = 约束 z_frac ∈ [lo, hi] 的节点 x/y（= 骨间膜/胫骨
    侧向约束的**刚性上界**理想化；载荷补丁节点始终排除）。
    """
    setup.fill_materials()
    meta = bf.build_bone_feb(
        setup.spec, out_path, load_n=1.0, top_frac=TOP_FRAC, bot_frac=bot_frac,
        axis=None, with_shell=True, time_steps=steps, mesh=setup.mesh,
        align_pca=True)
    text = out_path.read_text(encoding="ISO-8859-1")
    axis = np.asarray(geom["axis"], dtype=float)

    info: dict = {"level": level, "axis_idx": int(geom["axis_idx"]),
                  "load_n": float(load_n), "bot_frac": float(bot_frac),
                  "n_bot_fix": int(meta["n_bot_nodes"]),
                  "lateral": (list(lateral) if isinstance(lateral, tuple)
                              else lateral)}
    if level == "point":
        ids = [int(apex)]
        a_patch = float(trib[apex])
        text = inject_fragment(text, nodeset_xml("a_patch", ids))
        text = p1b.replace_loads(
            text, load_xml=nodal_load_xml("a_patch", -load_n * axis))
        text = p1b.add_load_ramp(text, kind="nodal")
        info.update({"mode": "nodal_point", "k": 0, "a_patch_mm2": a_patch,
                     "a_proj_mm2": float("nan"), "t_mpa": None,
                     "n_loaded": 1})
    else:
        k = int(level)
        patch = order[:k]
        ids = sorted({int(n) for f in patch for n in f["poly"]})
        a_patch = float(cum[k - 1])
        a_proj = float(cum_proj[k - 1])
        t = load_n / a_patch
        text = p1b.inject_surface(text, name="a_patch",
                                  polys=[f["poly"] for f in patch])
        text = inject_fragment(text, nodeset_xml("a_patch", ids))
        text = p1b.replace_loads(
            text, load_xml=p1b.traction_load_xml(surface="a_patch", vec=-t * axis))
        text = p1b.add_load_ramp(text, kind="pressure")
        info.update({"mode": "traction_patch", "k": k, "a_patch_mm2": a_patch,
                     "a_proj_mm2": a_proj, "t_mpa": t, "n_loaded": len(ids)})

    if lateral is not None:
        if nodes is None:
            raise RuntimeError("lateral 需要传入 nodes=")
        br = brace_ids_for(nodes, geom, lateral, ids, bot_frac)
        if not br:
            raise RuntimeError(f"lateral={lateral} 选到 0 个支撑节点")
        z_lo2, z_hi2 = float(geom["a_min"]), float(geom["a_max"])
        zfz = [(float(nodes[n - 1, geom["axis_idx"]]) - z_lo2) / (z_hi2 - z_lo2)
               for n in (br[0], br[-1])]
        if spring_k is None:
            text = inject_fragment(text, nodeset_xml("lat_brace", br))
            text = inject_boundary(text, lateral_bc_xml("lat_brace"))
            info["support"] = "hard"
        else:
            g_n, g_m, g_d, g_ids = spring_blocks(
                br, nodes, spring_k, int(setup.mesh["n_nodes"]) + 1)
            text = text.replace("    </Nodes>", "    </Nodes>\n" + g_n, 1)
            text = inject_fragment(text, g_m)
            if text.count("  <Step>") != 1:
                raise RuntimeError("spring: 找不到唯一的 '  <Step>'")
            text = text.replace("  <Step>", g_d + "  <Step>", 1)
            text = inject_boundary(text, ground_bc_xml())
            info["support"] = "spring"
            info["spring_k_n_per_mm"] = float(spring_k)
            info["n_ground_nodes"] = len(g_ids)
        info["n_brace"] = len(br)
        info["brace_z_frac_ends"] = [float(zfz[0]), float(zfz[1])]
    else:
        info["n_brace"] = 0
        info["support"] = "none"

    text = p1b.add_reaction_var(text)
    out_path.write_text(text, encoding="ISO-8859-1")
    info["patch_node_ids"] = ids
    return info


# ===========================================================================
# 后处理：反力 / 空空间哨兵 / 峰值定位
# ===========================================================================
def peak_locate(hdf5: Path, *, patch_ids, setup, geom, hex_rows,
                bot_frac: float = BOT_FRAC) -> dict:
    """峰值单元位置 + 补丁内局部峰值；含索引基数**正对照**。

    ⚠️ h5 ``meshes/0/domains/<dom>`` 的真实布局是 **(n_elem, 1 + n_nodes)**：
    第 0 列 = **单元 id（1-based）**，第 1..n 列 = **节点位置（0-based）**。
    （2026-10-12 实证：fib_cort 行0 = ``[1, 44,72,88,85,192,193,194,195]``
    ⇒ 单元 1 的节点位置 = 44..195 ⇒ id = 45,73,89,86,193,194,195,196
    = Python 侧 ``hex_rows[0]`` 逐字相同。把它当成「8 节点 + 材料 id」会错位。）
    """
    import h5py

    out: dict = {}
    with h5py.File(hdf5, "r") as h5:
        m = h5["meshes/0"]
        c = np.asarray(m["domains"][MAIN_DOM])
        n_nodes_elem = c.shape[1] - 1
        if n_nodes_elem != 8:
            raise RuntimeError(f"{MAIN_DOM} 期望 hex8（9 列），得到 {c.shape}")
        conn1 = c[:, 1:] + 1                                  # 1-based 全局 id

        # --- 正对照（集合级 + 逐行） ---
        a = set(int(x) for x in np.unique(conn1))
        b = set(int(x) for e in hex_rows for x in e)
        out["index_basis_ok"] = bool(a == b)
        out["index_rowwise_equal"] = int(sum(
            np.array_equal(conn1[i], np.asarray(hex_rows[i], dtype=int))
            for i in range(min(len(hex_rows), conn1.shape[0]))))
        out["index_n_elem"] = int(conn1.shape[0])
        if not out["index_basis_ok"]:
            raise RuntimeError(
                f"索引基数正对照失败：h5 侧 {len(a)} 个 id vs Python 侧 {len(b)} 个"
                f"（交集 {len(a & b)}）")

        nid = np.asarray(m["nodes"]["id"])
        xyz = np.column_stack([np.asarray(m["nodes"][c_]) for c_ in "xyz"])
        out["node_id_identity"] = bool(
            np.array_equal(nid, np.arange(1, nid.size + 1)))
        out["apex_coord_ok"] = bool(abs(
            float(xyz[geom["apex_node"] - 1, geom["axis_idx"]])
            - float(geom["a_max"])) < 1e-3)

        sid = sorted(h5["states"].keys(), key=int)[-1]
        s = np.asarray(h5[f"states/{sid}/element_data/stress/{MAIN_DOM}"][:],
                       dtype=np.float64)
        if s.shape[0] != conn1.shape[0]:
            raise RuntimeError(f"stress {s.shape} vs conn {conn1.shape}")
        vm = bf._von_mises_from_voigt(s)
        out["vm_mean_mpa"] = float(vm.mean())

        kpk = int(np.argmax(vm))
        pk_nodes = conn1[kpk]
        cen = xyz[[int(n) - 1 for n in pk_nodes]].mean(0)
        ps = np.array(sorted({int(x) for x in patch_ids}))
        z = xyz[:, geom["axis_idx"]]
        z_lo, z_hi = float(geom["a_min"]), float(geom["a_max"])
        bot_span = bot_frac * (z_hi - z_lo)
        out.update({
            "vm_max_elem": kpk,
            "peak_centroid": [float(x) for x in cen],
            "peak_z_frac": float((cen[geom["axis_idx"]] - z_lo) / (z_hi - z_lo)),
            "peak_nodes_in_patch_frac": float(
                np.mean([int(n) in set(ps.tolist()) for n in pk_nodes])),
            "peak_nodes_in_botfix_frac": float(
                np.mean(z[[int(n) - 1 for n in pk_nodes]] <= z_lo + bot_span)),
        })

        # 补丁内局部峰值：补丁是**表面**补丁，一个 hex 只在部分节点落入补丁里，
        # 故口径取「至少 1 个节点在补丁内」的单元（= 承载表层）。
        in_patch = np.any(np.isin(conn1, ps), axis=1)
        out["n_elem_touch_patch"] = int(in_patch.sum())
        out["vm_patch_max_mpa"] = (float(vm[in_patch].max())
                                   if in_patch.any() else None)
        # bottom_fix 区内峰值（判「峰值是否落在固定端」）
        zc = z[[[int(n) - 1 for n in row] for row in conn1]].mean(1)
        in_bot = zc <= (z_lo + bot_span)
        out["vm_botfix_max_mpa"] = float(vm[in_bot].max()) if in_bot.any() else None

        # 沿轴峰值剖面（相对骨长 0=远端, 1=apex 端）
        ecen = np.stack([xyz[[int(n) - 1 for n in row]].mean(0) for row in conn1])
        ez = (ecen[:, geom["axis_idx"]] - z_lo) / (z_hi - z_lo)
        bins = np.linspace(0.0, 1.0, 11)
        prof = []
        for i in range(10):
            sel = (ez >= bins[i]) & ((ez < bins[i + 1]) if i < 9 else (ez <= 1.0))
            prof.append(float(vm[sel].max()) if sel.any() else None)
        out["vm_profile_by_z_frac"] = prof

        # 空空间/发散哨兵 + 位移
        nd = h5[f"states/{sid}/node_data"]
        u = gs._node_var(nd, "displacement")
        out["max_abs_u_mm"] = float(np.abs(u).max())
        out["u_axial_patch_mean"] = float(u[ps - 1, geom["axis_idx"]].mean())
        out["u_axial_apex"] = float(u[geom["apex_node"] - 1, geom["axis_idx"]])
    return out


def run_level(setup, *, level, geom, order, cum, cum_proj, apex, trib,
              hex_rows, steps: int, timeout: float,
              load_n: float = LOAD_N, bot_frac: float = BOT_FRAC,
              nodes: np.ndarray | None = None, lateral=None,
              spring_k: float | None = None) -> dict:
    stem = ("scan_l%g_" % load_n) + (("b%.2f_" % bot_frac).replace("0.", "")
                                     ) + ("lat%s_" % (
                                         "patch" if lateral == "patch"
                                         else ("%.2f-%.2f" % lateral
                                               if isinstance(lateral, tuple)
                                               else "none"))
                                     ) + (("k%g_" % spring_k) if spring_k
                                          else "") + ("point" if level == "point"
                                                      else f"k{int(level):03d}")
    feb = OUT_TEMP / f"{stem}.feb"
    info = build_level_feb(setup, feb, level=level, geom=geom, order=order,
                           cum=cum, cum_proj=cum_proj, apex=apex, trib=trib,
                           steps=steps, load_n=load_n, bot_frac=bot_frac,
                           nodes=nodes, lateral=lateral, spring_k=spring_k)
    r = p1b.run_and_post(setup, feb, timeout=timeout)
    row: dict = {"level": level, **{k: info[k] for k in
                                    ("mode", "k", "load_n", "bot_frac",
                                     "lateral", "support", "n_brace",
                                     "a_patch_mm2", "a_proj_mm2", "t_mpa",
                                     "n_loaded")},
                 "spring_k": spring_k,
                 "ok": bool(r.get("ok")), "rc": r.get("rc"),
                 "wall_s": r.get("wall_clock_s")}
    if not r.get("ok"):
        row["err"] = (r.get("error_reason") or r.get("post_error")
                      or str(r.get("error"))[:300])
        return row
    row.update({"vm_max_mpa": r.get("peak_vm_mpa"),
                "vm_p95_mpa": r.get("p95_vm_mpa"),
                "vm_mean_mpa": r.get("vm_mean_mpa"),
                "R_axis_n": r.get("R_axis_n"), "R_net_n": r.get("R_net_n")})
    try:
        row.update(peak_locate(feb.with_suffix(".hdf5"),
                               patch_ids=info["patch_node_ids"], setup=setup,
                               geom=geom, hex_rows=hex_rows,
                               bot_frac=bot_frac))
    except Exception as exc:  # noqa: BLE001
        row["locate_error"] = f"{type(exc).__name__}: {exc}"
    return row


def fit_two_term(rows: list[dict], *, load_n: float = LOAD_N,
                 bot_frac: float = BOT_FRAC) -> dict:
    """σ = C/A + σ₀ 最小二乘（只用同一载荷、同一 BC 下的 traction 阶梯）。"""
    pts = [(r["a_patch_mm2"], r["vm_max_mpa"]) for r in rows
           if r.get("ok") and r.get("mode") == "traction_patch"
           and r.get("vm_max_mpa") and r.get("load_n") == load_n
           and r.get("bot_frac") == bot_frac]
    if len(pts) < 3:
        return {"ok": False, "n": len(pts)}
    A = np.array([p[0] for p in pts], dtype=float)
    S = np.array([p[1] for p in pts], dtype=float)
    X = np.column_stack([1.0 / A, np.ones_like(A)])
    coef, res, *_ = np.linalg.lstsq(X, S, rcond=None)
    pred = X @ coef
    ss_tot = float(((S - S.mean()) ** 2).sum())
    r2 = 1.0 - float(((S - pred) ** 2).sum()) / ss_tot if ss_tot > 0 else None
    return {"ok": True, "n": len(pts), "C_mpa_mm2": float(coef[0]),
            "sigma0_mpa": float(coef[1]), "r2": r2,
            "sigma0_over_30": float(coef[1]) / 30.0}


# ===========================================================================
# main
# ===========================================================================
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--levels", default=",".join(str(x) for x in DEFAULT_LEVELS),
                    help="逗号分隔；point = 单节点点载荷，整数 = 面片数 k")
    ap.add_argument("--steps", type=int, default=5)
    ap.add_argument("--timeout", type=float, default=900.0)
    ap.add_argument("--load", type=float, default=LOAD_N,
                    help="总载荷 N（线性度正对照用；默认 1000）")
    ap.add_argument("--loads", default="",
                    help="逗号分隔的多载荷（一个进程内扫；用于定位线性失效点）")
    ap.add_argument("--tag", default="")
    ap.add_argument("--bot-fracs", default="",
                    help="逗号分隔的固定带比例（B：自由长度扫描）；空 = 默认 0.30")
    ap.add_argument("--laterals", default="",
                    help="逗号分隔的横向支撑（C）：none/patch/mid/long/full 或 lo:hi")
    ap.add_argument("--spring-ks", default="",
                    help="逗号分隔的侧向弹簧刚度 k (N/mm，每节点)；空 = 用硬 BC")
    args = ap.parse_args()

    levels = []
    for x in args.levels.split(","):
        x = x.strip()
        if not x:
            continue
        levels.append("point" if x.lower().startswith("p") else int(x))

    t0 = time.time()
    setup = p1b.Setup()
    nodes = setup.nodes_for(True)
    geom = setup.facet_geom(top_frac=TOP_FRAC, align_pca=True)
    ax_i = int(geom["axis_idx"])
    band = geom["band"]
    trib = tributary_areas(band)
    cap_nodes = sorted({int(n) for f in geom["cap"] for n in f["poly"]})
    pz = nodes[[c - 1 for c in cap_nodes]]
    apex = int(cap_nodes[int(np.argmax(pz[:, ax_i]))])

    order, cum, cum_proj, dist = facet_ladder(nodes, band, apex)
    hex_rows = [e for e in setup.elems if len(e) == 8]
    geom["apex_node"] = int(apex)
    print(f"[setup] {time.time() - t0:.0f}s  axis={'XYZ'[ax_i]}"
          f"  apex={apex}  n_band_facets={len(band)}  a_band={cum[-1]:.1f} mm^2"
          f"  trib_apex={trib[apex]:.2f}  hex_rows={len(hex_rows)}", flush=True)
    print(f"[nominal] 皮质截面 F/A = {LOAD_N / 102.0:.2f} MPa（§9.8: A≈102 mm²）")

    rows: list[dict] = []
    loads = ([float(x) for x in args.loads.split(",") if x.strip()]
             if args.loads.strip() else [args.load])
    bot_fracs = ([float(x) for x in args.bot_fracs.split(",") if x.strip()]
                 if args.bot_fracs.strip() else [BOT_FRAC])
    laterals = ([parse_lateral(x) for x in args.laterals.split(",") if x.strip()]
                if args.laterals.strip() else [None])
    spring_ks = ([float(x) for x in args.spring_ks.split(",") if x.strip()]
                 if args.spring_ks.strip() else [None])
    span = float(geom["a_max"] - geom["a_min"])
    # Euler 定标：以 bot_frac=0.30 实测拐点 ≈250 N 为锚，P_cr ∝ 1/L_free²
    if len(bot_fracs) > 1:
        print("[Euler 预测] P_cr ≈ 250 N × (0.70/(1−bf))²  （锚=bf0.30 的实测拐点；"
              "L_free=span·(1−bf)）")
        for bf in bot_fracs:
            L = span * (1.0 - bf)
            print(f"    bf={bf:.2f}  L_free={L:6.1f} mm  P_cr≈{250*(0.70/(1-bf))**2:7.1f} N")
    for bf in bot_fracs:
        for lat in laterals:
            for sk in spring_ks:
                for load_n in loads:
                    for lv in levels:
                        t1 = time.time()
                        row = run_level(setup, level=lv, geom=geom, order=order,
                                        cum=cum, cum_proj=cum_proj, apex=apex,
                                        trib=trib, hex_rows=hex_rows,
                                        steps=args.steps, timeout=args.timeout,
                                        load_n=load_n, bot_frac=bf,
                                        nodes=nodes, lateral=lat,
                                        spring_k=sk)
                        rows.append(row)
                        tag = "point " if lv == "point" else f"k={lv:<4d}"
                        print(f"  [bf={bf:.2f} lat={lateral_label(lat):<9s}"
                              f" k={sk if sk else 'hard':<6} F={load_n:<6g} {tag}]"
                              f" A={row['a_patch_mm2']:8.1f}  n_br={row['n_brace']:4d}"
                              f"  ok={row['ok']}  R={_f(row.get('R_axis_n'))}"
                              f"  vm={_f(row.get('vm_max_mpa'), 1)}"
                              f"  |u|max={_f(row.get('max_abs_u_mm'), 2)}"
                              f"  u_ax={_f(row.get('u_axial_apex'), 2)}"
                              f"  z={_f(row.get('peak_z_frac'), 2)}"
                              f"  ({time.time() - t1:.0f}s)", flush=True)
                        if not row["ok"]:
                            print(f"         err: {row.get('err')}", flush=True)

    fit = fit_two_term(rows, load_n=loads[-1], bot_frac=bot_fracs[0])
    out = {"load_n": loads[-1], "loads": loads, "bot_fracs": bot_fracs,
           "laterals": [lateral_label(x) for x in laterals],
           "top_frac": TOP_FRAC, "bot_frac": bot_fracs[0],
           "apex_node": apex, "axis_idx": ax_i,
           "a_band_mm2": float(cum[-1]),
           "nominal_f_over_a_mpa": loads[-1] / 102.0, "fit": fit, "rows": rows}
    tag = f"_{args.tag}" if args.tag else ""
    p = OUT_RESULTS / f"G4_FE_AREASCAN{tag}.json"
    p.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n[fit] σ = C/A + σ₀ : C={fit.get('C_mpa_mm2'):.1f} MPa·mm²"
          f"  σ₀={fit.get('sigma0_mpa'):.1f} MPa  R²={fit.get('r2')}"
          f"  (σ₀/30 = {fit.get('sigma0_over_30'):.2f})" if fit.get("ok")
          else f"[fit] 级数不足: {fit}")
    print(f"[out] {p}")
    print(f"[done] {time.time() - t0:.0f}s")


def _f(x, nd: int = 3) -> str:
    return "?" if x is None else f"{x:.{nd}f}"


def _per_n(v, f) -> float | None:
    return None if (v is None or not f) else float(v) / float(f)


#: C：横向支撑的命名区间（z_frac；0=远端, 1=apex 端；bottom_fix 占 [0,0.30]）
LATERAL_PRESETS = {"none": None, "patch": "patch",
                   "mid": (0.45, 0.65), "long": (0.35, 0.80),
                   "full": (0.32, 0.85)}


def parse_lateral(s: str):
    s = s.strip().lower()
    if s in LATERAL_PRESETS:
        return LATERAL_PRESETS[s]
    if ":" in s:
        lo, hi = s.split(":")
        return (float(lo), float(hi))
    raise ValueError(f"无法解析 lateral={s!r}（用 none/patch/mid/long/full 或 lo:hi）")


def lateral_label(lat) -> str:
    if lat is None:
        return "none"
    if lat == "patch":
        return "patch"
    return "%.2f-%.2f" % lat


if __name__ == "__main__":
    main()
