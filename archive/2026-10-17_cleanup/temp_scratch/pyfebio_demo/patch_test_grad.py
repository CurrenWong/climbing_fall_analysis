"""正对照：合成单元 + 已知线性位移场 ⇒ 重构必须还原解析应力。"""
import numpy as np

dN_hex = lambda xi, et, ze: np.array([  # noqa: E731
    [(-1)*(1+(-1)*et)*(1+(-1)*ze), (-1)*(1+(-1)*xi)*(1+(-1)*ze),
     (-1)*(1+(-1)*xi)*(1+(-1)*et)],
    [(1)*(1+(-1)*et)*(1+(-1)*ze), (-1)*(1+(1)*xi)*(1+(-1)*ze),
     (-1)*(1+(1)*xi)*(1+(-1)*et)],
    [(1)*(1+(1)*et)*(1+(-1)*ze), (1)*(1+(1)*xi)*(1+(-1)*ze),
     (-1)*(1+(1)*xi)*(1+(1)*et)],
    [(-1)*(1+(1)*et)*(1+(-1)*ze), (1)*(1+(-1)*xi)*(1+(-1)*ze),
     (-1)*(1+(-1)*xi)*(1+(1)*et)],
    [(-1)*(1+(-1)*et)*(1+(1)*ze), (-1)*(1+(-1)*xi)*(1+(1)*ze),
     (1)*(1+(-1)*xi)*(1+(-1)*et)],
    [(1)*(1+(-1)*et)*(1+(1)*ze), (-1)*(1+(1)*xi)*(1+(1)*ze),
     (1)*(1+(1)*xi)*(1+(-1)*et)],
    [(1)*(1+(1)*et)*(1+(1)*ze), (1)*(1+(1)*xi)*(1+(1)*ze),
     (1)*(1+(1)*xi)*(1+(1)*et)],
    [(-1)*(1+(1)*et)*(1+(1)*ze), (1)*(1+(-1)*xi)*(1+(1)*ze),
     (1)*(1+(-1)*xi)*(1+(1)*et)],
], float) / 8.0

P0 = np.array([[-1., -1., -1.], [1., -1., -1.], [1., 1., -1.], [-1., 1., -1.],
               [-1., -1., 1.], [1., -1., 1.], [1., 1., 1.], [-1., 1., 1.]])

# 解析位移场：u = H X，H 任意常张量 ⇒ du/dX = H（处处）
H = np.array([[0.004, 0.001, -0.0005],
              [0.0002, -0.002, 0.0008],
              [0.001, -0.0003, 0.003]])
P1 = P0 + (H @ P0.T).T
E, nu = 7300.0, 0.3
lam = E * nu / ((1+nu)*(1-2*nu))
mu = E / (2*(1+nu))


def vm_of(sig):
    xx, yy, zz, xy, yz, xz = sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], \
        sig[1, 2], sig[0, 2]
    return np.sqrt(0.5*((xx-yy)**2+(yy-zz)**2+(zz-xx)**2)
                   + 3*(xy**2+yz**2+xz**2))


eps = 0.5*(H + H.T)
sig_true = lam*np.trace(eps)*np.eye(3) + 2*mu*eps
print("=== 解析（真值）===")
print(f"  应变 tr={np.trace(eps):.6e}  σ_vm = {vm_of(sig_true):.4f} MPa")

GAUSS = [tuple(g) for g in np.array(
    [[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],
    float)*(1/np.sqrt(3))]

for label, solve_fn in (
    ("A) 原写法 solve(J,M).T", lambda J, M: np.linalg.solve(J, M).T),
    ("B) 正解 M @ inv(J)", lambda J, M: M @ np.linalg.inv(J)),
    ("C) solve(J.T, M.T).T", lambda J, M: np.linalg.solve(J.T, M.T).T),
):
    print(f"\n=== {label} ===")
    for p in GAUSS:
        dN = dN_hex(*p)
        J = P0.T @ dN
        M = (P1-P0).T @ dN
        g = solve_fn(J, M)
        e = 0.5*(g + g.T)
        s = lam*np.trace(e)*np.eye(3) + 2*mu*e
        print(f"  gauss{tuple(round(x,3) for x in p)}: σ_vm = {vm_of(s):.4f} MPa"
              f"   (dev from true {abs(vm_of(s)-vm_of(sig_true)):.4e})")

print("\n" + "=" * 78)
print("⇒ 只有还原到 %.4f MPa 的那个写法才可信" % vm_of(sig_true))
print("=" * 78)
