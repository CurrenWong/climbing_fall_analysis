"""点到三角形精确距离：Ericson《Real-Time Collision Detection》§5.1.5
标准区域算法，按三角形向量化（一个点对所有三角形）。
"""

import numpy as np


def dist_point_to_tris(p, T):
    """p (3,) 单点；T (N,3,3) 三角形。返回 (N,) 距离。"""
    a, b, c = T[:, 0], T[:, 1], T[:, 2]
    ab, ac = b - a, c - a
    ap = p - a
    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)
    bp = p - b
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    cp = p - c
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)

    va_ = d3 * d6 - d5 * d4
    vb_ = d5 * d2 - d1 * d6
    vc_ = d1 * d4 - d3 * d2

    Q = np.empty_like(a)
    done = np.zeros(len(T), bool)

    # 顶点 a
    m = (d1 <= 0) & (d2 <= 0)
    Q[m] = a[m]
    done |= m
    # 顶点 b
    m = (~done) & (d3 >= 0) & (d4 <= d3)
    Q[m] = b[m]
    done |= m
    # 边 ab
    m = (~done) & (vc_ <= 0) & (d1 >= 0) & (d3 <= 0)
    den = d1 - d3
    den = np.where(np.abs(den) < 1e-20, 1e-20, den)
    v = (d1 / den)
    Q[m] = a[m] + v[m, None] * ab[m]
    done |= m
    # 顶点 c
    m = (~done) & (d6 >= 0) & (d5 <= d6)
    Q[m] = c[m]
    done |= m
    # 边 ac
    m = (~done) & (vb_ <= 0) & (d2 >= 0) & (d6 <= 0)
    den = d2 - d6
    den = np.where(np.abs(den) < 1e-20, 1e-20, den)
    w = (d2 / den)
    Q[m] = a[m] + w[m, None] * ac[m]
    done |= m
    # 边 bc
    m = (~done) & (va_ <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    den = (d4 - d3) + (d5 - d6)
    den = np.where(np.abs(den) < 1e-20, 1e-20, den)
    w = ((d4 - d3) / den)
    Q[m] = b[m] + w[m, None] * (c[m] - b[m])
    done |= m
    # 面内
    m = ~done
    if m.any():
        den = va_[m] + vb_[m] + vc_[m]
        den = np.where(np.abs(den) < 1e-20, 1e-20, den)
        v = (vb_[m] / den)[:, None]
        w = (vc_[m] / den)[:, None]
        Q[m] = a[m] + v * ab[m] + w * ac[m]

    return np.linalg.norm(p - Q, axis=1)
