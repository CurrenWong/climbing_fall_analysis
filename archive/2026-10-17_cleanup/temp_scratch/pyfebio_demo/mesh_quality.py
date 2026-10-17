"""自己算 tet 质量（不依赖 gmsh API），避免上轮 API 误用。

gamma = 12 * (3V)^(2/3) / sum(edge^2)   —— 正四面体 = 1
  < 0.1  = sliver，会产生伪应力
  FEBio/gmsh 实践：min gamma > 0.1 ~ 0.2 才适合直接求解
另报：长宽比 max_edge/min_edge、最小边长、负体积
"""
from pathlib import Path
import os, sys
import numpy as np
import gmsh

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"
os.environ["PATH"] = f"{LIB};{SDK/'bin'};{os.environ['PATH']}"
sys.path.insert(0, str(LIB))

MSH = Path(r"D:\Project\climbing_fall_analysis\temp\ankle_mesh")


def tet_gamma(v: np.ndarray) -> np.ndarray:
    """v: (n,4,3) 顶点. 返回每个 tet 的 gamma 质量."""
    a, b, c, d = v[:, 0], v[:, 1], v[:, 2], v[:, 3]
    vol = np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a)) / 6.0
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    s = sum(((v[:, i] - v[:, j]) ** 2).sum(1) for i, j in pairs)
    return 12.0 * (3.0 * vol) ** (2.0 / 3.0) / np.maximum(s, 1e-300)


def report(p: Path):
    gmsh.clear()
    gmsh.open(str(p))
    et, _, en = gmsh.model.mesh.getElements(3)
    el = np.concatenate([np.asarray(n, np.int64).reshape(-1, 4)
                         for e, n in zip(et, en) if e == 4])
    ntags, ncoord, _ = gmsh.model.mesh.getNodes()
    nc = np.asarray(ncoord, float).reshape(-1, 3)
    tag2row = {int(t): i for i, t in enumerate(ntags)}
    v = nc[np.vectorize(tag2row.get)(el)]           # (n,4,3) 单位 m
    g = tet_gamma(v)

    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    L = np.stack([np.linalg.norm(v[:, i] - v[:, j], axis=1) for i, j in pairs], 1)
    Lmax, Lmin = L.max(1), L.min(1)
    vol = np.abs(np.einsum("ij,ij->i",
                           np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]),
                           v[:, 3] - v[:, 0])) / 6.0

    def pct(t):
        return 100.0 * (g < t).sum() / len(g)

    ps = np.percentile(g, [0, 1, 5, 25, 50, 95, 100])
    print(f"\n{p.name}")
    print(f"  tet {len(g):6d} | 体积 {vol.sum()*1e6:9.4f} cm^3 | 负体积 {int((vol<=0).sum())}")
    print(f"  gamma  min={ps[0]:7.4f} p1={ps[1]:7.4f} p5={ps[2]:7.4f} "
          f"p25={ps[3]:7.4f} 中位={ps[4]:7.4f} p95={ps[5]:6.4f} max={ps[6]:7.4f}")
    print(f"  gamma<0.01: {pct(0.01):5.2f}%   <0.05: {pct(0.05):5.2f}%   "
          f"<0.10: {pct(0.10):5.2f}%   <0.20: {pct(0.20):5.2f}%")
    print(f"  边长   min={Lmin.min()*1000:7.3f} mm  中位={np.median(L)*1000:7.3f} mm  "
          f"max={Lmax.max()*1000:8.3f} mm")
    print(f"  长宽比 max/min  中位={np.median(Lmax/Lmin):7.3f}  "
          f"p95={np.percentile(Lmax/Lmin,95):8.3f}  max={(Lmax/Lmin).max():9.2f}")
    print(f"  结论: ", end="")
    if g.min() <= 0:
        print("❌ 有非正体积单元，无法求解")
    elif pct(0.05) > 1.0:
        print(f"❌ {pct(0.05):.1f}% 单元 gamma<0.05 —— 存在大量 sliver，会直接污染应力场")
    elif pct(0.10) > 1.0:
        print(f"⚠️ {pct(0.10):.1f}% 单元 gamma<0.10 —— 需先做质量优化")
    else:
        print("✅ 质量可用（gamma 中位 %.3f）" % np.median(g))


gmsh.initialize(str(LIB))
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)
for f in sorted(MSH.glob("*.msh")):
    report(f)
gmsh.finalize()
