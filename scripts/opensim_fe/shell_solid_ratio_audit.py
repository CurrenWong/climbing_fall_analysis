# -*- coding: utf-8 -*-
"""闸门 7 · 皮质**壳厚** vs 壳下**实体单元厚度**审计（零 FE，带同款架构对照）。

依据：FEBio 手册 §8.3.1.4 + L1981 —— 「It is the user's responsibility to ensure that the
**shell thickness is less than that of the solid element**; otherwise a negative Jacobian
error will be **automatically** generated.」

问：T6 的皮质壳 1.5 mm 是否真的 ≥ 壳下实体单元厚度？多少单元越界？
对照：股骨（皮质是 solid，按项目记录无壳）——需实测确认。

用法：.venv\\Scripts\\python.exe scripts\\opensim_fe\\shell_solid_ratio_audit.py

立此件的**口径教训**：旧记 T6「特征边长 1.72 mm」系**中位单元体积的立方根**（等效立方体），
对四面体**系统性偏小**；按手册要求的量（**壳下实体单元高**）实测**中位 2.965 mm** ⇒ 比 **0.506**。
⇒ **报「特征长度」必须给定义；与判据阈值比时必须用判据自己定义的那个量**（LRN-20261015-114）。
对照部位 **L3**（同款「松质 tet + 皮质壳」架构）比值 **0.487** ⇒ **T6 不特殊**。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import p4_pb_probe as p4                                            # noqa: E402


def tet_vol(p: np.ndarray) -> np.ndarray:
    """p: (n,4,3) -> (n,) 体积（无符号）。"""
    a, b, c, d = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0], None
    del d
    return np.abs(np.einsum("ij,ij->i", np.cross(b, c), a)) / 6.0


def tri_area(q: np.ndarray) -> np.ndarray:
    """q: (n,3,3) -> (n,) 面积。"""
    return 0.5 * np.linalg.norm(np.cross(q[:, 1] - q[:, 0], q[:, 2] - q[:, 0]), axis=1)


def stats(x: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return {}
    return {"n": int(x.size), "min": float(x.min()), "p05": float(np.percentile(x, 5)),
            "p25": float(np.percentile(x, 25)), "median": float(np.median(x)),
            "mean": float(x.mean()), "p75": float(np.percentile(x, 75)),
            "p95": float(np.percentile(x, 95)), "max": float(x.max())}


def probe(bone: str) -> dict:
    part, mesh, solid_names, shell_specs, dom, joint, loads, h_arr = p4.load_ctx(bone)
    N = mesh["nodes"]
    out = {"bone": bone, "n_nodes": int(mesh["n_nodes"]),
           "n_id_solid": int(mesh["n_id_solid"]), "n_id_shell": int(mesh["n_id_shell"]),
           "solid_domains": {}, "shell_domains": {}}

    # ---- 实体单元：全局边长 ----
    all_edges = []
    solid_tets = []
    for mn, b in mesh["by_solid"].items():
        t = np.asarray(b["tets"], int)
        hx = np.asarray(b["hexes"], int)
        out["solid_domains"][mn] = {"n_tet": int(t.shape[0]), "n_hex": int(hx.shape[0])}
        if t.size:
            solid_tets.append((mn, t))
            P = N[t - 1]                                   # (n,4,3)
            for i in range(4):
                for j in range(i + 1, 4):
                    all_edges.append(np.linalg.norm(P[:, i] - P[:, j], axis=1))
    if all_edges:
        E = np.concatenate(all_edges)
        out["solid_edge_len"] = stats(E)

    # ---- 壳单元：厚度 + 节点集 ----
    shell_elems: list[frozenset] = []
    for mn, b in mesh["by_shell"].items():
        q4 = np.asarray(b["q4"], int)
        t3 = np.asarray(b["t3"], int)
        out["shell_domains"][mn] = {
            "n_q4": int(q4.shape[0]), "n_t3": int(t3.shape[0]),
            "thickness_mm": float(b["thickness_mm"]), "mat": b["mat_name"]}
        for row in q4:
            shell_elems.append(frozenset(int(x) for x in row))
        for row in t3:
            shell_elems.append(frozenset(int(x) for x in row))

    if not shell_elems or not solid_tets:
        out["note"] = "无壳 或 无实体 —— 本部位不适用壳/实体厚度冲突检验"
        return out

    # ---- 壳下实体单元厚度：tet 面 ⊂ 壳单元节点集 ⇒ 高度 = 3V/A ----
    key: dict[frozenset, list] = {}
    for mn, t in solid_tets:
        P = N[t - 1]
        V = tet_vol(P)
        for ti in range(t.shape[0]):
            for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
                key.setdefault(frozenset(int(x) for x in t[ti, list(f)]), []).append(
                    (mn, ti, tuple(int(x) for x in t[ti, list(f)]), V[ti]))

    heights: list[float] = []
    n_shell_elems_matched = 0
    n_face_hits = 0
    for S in shell_elems:
        hit = False
        for k, ent in key.items():
            if not k <= S:                      # tet 面必须是壳单元的节点子集
                continue
            for _mn, _ti, tri, V in ent:
                A = tri_area(N[np.asarray(tri, int) - 1][None, :, :])[0]
                if A > 1e-12:
                    heights.append(3.0 * V / A)
                    n_face_hits += 1
                    hit = True
        if hit:
            n_shell_elems_matched += 1
    out["shell_face_query"] = {"n_shell_elems": len(shell_elems),
                               "n_shell_elems_matched": n_shell_elems_matched,
                               "n_tet_face_hits": n_face_hits}
    H = np.asarray(heights, float)
    out["solid_height_under_shell"] = stats(H)

    # ---- 判据：thickness < height ? ----
    th = max(v["thickness_mm"] for v in out["shell_domains"].values())
    out["verdict"] = {
        "shell_thickness_mm": th,
        "n_violating": int((H <= th).sum()),
        "frac_violating": float((H <= th).mean()) if H.size else None,
        "ratio_th_over_median_h": float(th / np.median(H)) if H.size else None,
        "worst_ratio_th_over_h": float(th / H.min()) if H.size else None,
    }
    return out


def main() -> int:
    res = {}
    for bone in ("T6", "femur_r", "L3"):
        try:
            res[bone] = probe(bone)
        except Exception as exc:                                    # noqa: BLE001
            res[bone] = {"bone": bone, "error": "%s: %s" % (type(exc).__name__, exc)}
    dst = ROOT / "results" / "opensim_fe" / "SHELL_SOLID_RATIO_AUDIT.json"
    dst.write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str),
                   encoding="utf-8")

    for bone, r in res.items():
        print("=" * 88)
        print(bone)
        print("=" * 88)
        if "error" in r:
            print("  ERROR:", r["error"])
            continue
        print("  nodes=%s solid=%s shell=%s" % (r["n_nodes"], r["n_id_solid"], r["n_id_shell"]))
        print("  solid domains:", {k: v for k, v in r["solid_domains"].items()})
        print("  shell domains:", {k: v for k, v in r["shell_domains"].items()})
        if r.get("solid_edge_len"):
            e = r["solid_edge_len"]
            print("  实体单元边长: min=%.3f p05=%.3f median=%.3f mean=%.3f p95=%.3f max=%.3f"
                  % (e["min"], e["p05"], e["median"], e["mean"], e["p95"], e["max"]))
        if r.get("solid_height_under_shell"):
            h = r["solid_height_under_shell"]
            print("  壳下实体厚度: n=%d min=%.3f p05=%.3f median=%.3f mean=%.3f "
                  "p95=%.3f max=%.3f"
                  % (h["n"], h["min"], h["p05"], h["median"], h["mean"], h["p95"], h["max"]))
        if r.get("verdict"):
            v = r["verdict"]
            print("  ⇒ 壳厚 %.2f mm | 违例(厚>=高) %d / %s = %.1f%% | 厚/中位高 %.3f | "
                  "厚/最小高 %.3f"
                  % (v["shell_thickness_mm"], v["n_violating"],
                     r["solid_height_under_shell"]["n"],
                     100 * v["frac_violating"], v["ratio_th_over_median_h"],
                     v["worst_ratio_th_over_h"]))
        if r.get("shell_face_query"):
            q = r["shell_face_query"]
            print("  壳面→实体匹配: 壳单元 %d，匹配上 %d，tet 面命中 %d"
                  % (q["n_shell_elems"], q["n_shell_elems_matched"], q["n_tet_face_hits"]))
        if r.get("note"):
            print("  " + r["note"])
    print("\n产物:", dst.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
