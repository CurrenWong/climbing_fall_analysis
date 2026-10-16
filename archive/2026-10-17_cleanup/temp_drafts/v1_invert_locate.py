"""定位倒置单元：用**最后收敛状态**的位移算每个单元的 det F，看是骨还是板、在哪。

用最优组合跑（gap=-0.05, pad_E=0.185, auto_penalty=1, DYNAMIC+阻尼）。
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT / "src"), str(ROOT / "scripts" / "opensim_fe")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402
import thums_feb as T  # noqa: E402
from climbing.coupling import fe_post  # noqa: E402

MESH = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_anatframe.npz"
CAND = [ROOT / "temp" / "opensim_fe" / "v1_combo" / "g-0.05_e0.185_ap1",
        ROOT / "temp" / "opensim_fe" / "v1_combo" / "g-0.05_e0.185"]

# hex8 在形心的 dN/dxi（±1/8）
S = np.array([[-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1],
              [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1]], float)
DN8 = S / 8.0


def detF_hex8(X, u):
    """X (n,8,3) 参考坐标, u (n,8,3) 位移 ⇒ 形心处的 det F。"""
    J0 = np.einsum('ia,nib->nab', DN8, X)
    J1 = np.einsum('ia,nib->nab', DN8, X + u)
    return np.linalg.det(np.linalg.solve(J0, J1))


def detF_tet4(X, u):
    D = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    J0 = np.einsum('ia,nib->nab', D, X)
    J1 = np.einsum('ia,nib->nab', D, X + u)
    return np.linalg.det(np.linalg.solve(J0, J1))


def main():
    mesh = T.load_thums_mesh(MESH)
    nodes = np.asarray(mesh["nodes"], float)
    ec = np.asarray(mesh["elements_cort"], int)[:, :8]
    es = np.asarray(mesh["elements_spon"], int)[:, :4]
    pn = np.asarray(mesh["node_sets"]["plantar"], int)
    cc = nodes[ec].mean(axis=1)

    for wd in CAND:
        xp = list(wd.glob("*.xplt"))
        if not xp:
            print(f"[skip] {wd.name}: 无 xplt")
            continue
        xp = xp[0]
        try:
            du = fe_post.read_displacement(xp, last=True)
        except Exception as exc:                   # noqa: BLE001
            print(f"[skip] {wd.name}: 读位移失败 {type(exc).__name__}")
            continue
        u = np.asarray(du, float).reshape(-1, 3)
        nb = len(nodes)
        if len(u) < nb:
            print(f"[skip] {wd.name}: 位移节点数 {len(u)} < 网格 {nb}")
            continue
        # 自检：末尾应是**板**（n_bone + n_pad），板底固定 ⇒ 位移应近零
        tail = u[nb:]
        print(f"\n=== {wd.name} ===")
        print(f"  位移节点 {len(u)} = 骨 {nb} + 板 {len(tail)}  "
              f"（板位移模 max={np.linalg.norm(tail, axis=1).max():.4f} mm）")
        u = u[:nb]
        jc = detF_hex8(nodes[ec], u[ec])
        js = detF_tet4(nodes[es], u[es])
        bad_c, bad_s = jc <= 0.0, js <= 0.0
        print(f"\n=== {wd.name} （最后收敛状态）===")
        print(f"  CORT hex8 : {int(bad_c.sum())}/{len(ec)} 倒置  "
              f"det F: min={jc.min():.3e} 中位={np.median(jc):.3f}")
        print(f"  SPON tet4 : {int(bad_s.sum())}/{len(es)} 倒置  "
              f"det F: min={js.min():.3e} 中位={np.median(js):.3f}")
        if bad_c.any():
            d_pl = np.linalg.norm(cc[:, None, :] - nodes[pn][None, :, :],
                                  axis=2).min(axis=1)
            h = float(np.median(np.linalg.norm(nodes[ec[:, 0]] - nodes[ec[:, 1]],
                                               axis=1)))
            print(f"  倒置 CORT 单元 距跖面: 中位 {np.median(d_pl[bad_c]):.2f} mm "
                  f"({np.median(d_pl[bad_c])/h:.2f} h)  "
                  f"范围 [{d_pl[bad_c].min():.2f}, {d_pl[bad_c].max():.2f}]")
            idx = np.where(bad_c)[0][:12]
            print(f"  前 12 个倒置单元 id: {idx.tolist()}")
            print(f"  其 det F: {np.round(jc[idx], 3).tolist()}")
        # 整体位移量级（判断是否大变形）
        mag = np.linalg.norm(u, axis=1)
        print(f"  位移模: max={mag.max():.3f} mm  中位={np.median(mag):.4f} mm")


if __name__ == "__main__":
    main()
