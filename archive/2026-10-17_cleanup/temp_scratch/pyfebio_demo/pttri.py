"""可靠的 point-triangle 最近点（Ericson 标量版，逐点暴力），带自检。"""
import numpy as np


def closest_on_tri(p, a, b, c):
    ab, ac = b - a, c - a
    ap = p - a
    d1, d2 = ap @ ab, ap @ ac
    if d1 <= 0 and d2 <= 0:
        return a
    bp = p - b
    d3, d4 = bp @ ab, bp @ ac
    if d3 >= 0 and d4 <= d3:
        return b
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        v = d1 / (d1 - d3)
        return a + v * ab
    cp = p - c
    d5, d6 = cp @ ab, cp @ ac
    if d6 >= 0 and d5 <= d6:
        return c
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        w = d2 / (d2 - d6)
        return a + w * ac
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return b + w * (c - b)
    den = 1.0 / (va + vb + vc)
    v, w = vb * den, vc * den
    return a + ab * v + ac * w


def dist_pt_tris(P, tris):
    """P:(n,3), tris:(m,3,3) -> (n,) 最短距离"""
    P = np.atleast_2d(np.asarray(P, float))
    out = np.empty(len(P))
    for k, p in enumerate(P):
        best = np.inf
        for t in tris:
            q = closest_on_tri(p, t[0], t[1], t[2])
            dd = float(np.linalg.norm(p - q))
            if dd < best:
                best = dd
                if best == 0.0:
                    break
        out[k] = best
    return out


if __name__ == "__main__":
    A, B, C = np.array([0., 0, 0]), np.array([1., 0, 0]), np.array([0., 1, 0])
    tris = np.array([[A, B, C]])
    tests = [([0.3, 0.3, 2.0], 2.0, "面正上方"),
             ([5.0, 0.0, 0.0], 4.0, "B 外侧"),
             ([-1.0, 0.0, 0.0], 1.0, "A 外侧"),
             ([0.0, 5.0, 0.0], 4.0, "C 外侧"),
             ([0.5, 0.5, 0.0], 0.0, "面内"),
             ([0.5, 0.6, 3.0], 3.0, "边外正上"),
             ([1.0, 1.0, 0.0], np.sqrt(0.5), "斜边外"),
             ([0.25, 0.25, 0.0], 0.0, "面内2")]
    ok = True
    for p, exp, nm in tests:
        got = dist_pt_tris(np.array([p], float), tris)[0]
        good = abs(got - exp) < 1e-9
        ok &= good
        print(f"  [{'OK' if good else 'XX'}] {nm:10s} 期望 {exp:.6f} 实得 {got:.6f} "
              f"Δ={got-exp:+.2e}")
    print(f"\n自检 {'ALL PASS' if ok else 'FAILED'}")
