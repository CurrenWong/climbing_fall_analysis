"""修正版：① patch test 用**畸变**单元 ② 梯度用正解 ③ tet4 单独采样。"""
import numpy as np

P0c = np.array([[-1., -1., -1.], [1., -1., -1.], [1., 1., -1.], [-1., 1., -1.],
                [-1., -1., 1.], [1., -1., 1.], [1., 1., 1.], [-1., 1., 1.]])


def dN_hex(xi, et, ze):
    sg = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                   [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], float)
    dN = np.zeros((8, 3))
    for i in range(8):
        a, b, c = sg[i]
        dN[i] = [a*(1+b*et)*(1+c*ze)/8, b*(1+a*xi)*(1+c*ze)/8,
                 c*(1+a*xi)*(1+b*et)/8]
    return dN


GH = [tuple(g) for g in np.array(
    [[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],
    float)*(1/np.sqrt(3))]

# 畸变：把 8 个角点各自挪一下（仍是合法 hex8）
rng = np.random.RandomState(0)
P0d = P0c + rng.uniform(-0.25, 0.25, P0c.shape)

H = np.array([[0.004, 0.001, -0.0005],
              [0.0002, -0.002, 0.0008],
              [0.001, -0.0003, 0.003]])   # 非对称
E, nu = 7300.0, 0.3
lam, mu = E*nu/((1+nu)*(1-2*nu)), E/(2*(1+nu))
eps_true = 0.5*(H+H.T)
sig_true = lam*np.trace(eps_true)*np.eye(3) + 2*mu*eps_true


def vm(s):
    xx, yy, zz, xy, yz, xz = s[0,0], s[1,1], s[2,2], s[0,1], s[1,2], s[0,2]
    return np.sqrt(0.5*((xx-yy)**2+(yy-zz)**2+(zz-xx)**2)+3*(xy**2+yz**2+xz**2))


print("=" * 84)
print(f"解析真值 σ_vm = {vm(sig_true):.4f} MPa")
print("=" * 84)
for label, P0 in (("规则单元", P0c), ("畸变单元", P0d)):
    P1 = P0 + (H @ P0.T).T
    print(f"\n--- {label} ---")
    for form, fn in (
        ("A solve(J,M).T（我原用的）", lambda J, M: np.linalg.solve(J, M).T),
        ("B M @ inv(J)（正解）", lambda J, M: M @ np.linalg.inv(J)),
    ):
        errs = []
        for p in GH:
            dN = dN_hex(*p)
            J = P0.T @ dN
            M = (P1-P0).T @ dN
            e = 0.5*(fn(J, M) + fn(J, M).T)
            s = lam*np.trace(e)*np.eye(3) + 2*mu*e
            errs.append(abs(vm(s) - vm(sig_true)))
        print(f"  {form:<30} max|Δσ_vm| = {max(errs):.4e}"
              f"  {'✅ 精确' if max(errs) < 1e-6 else '❌ 不精确'}")
print("\n⇒ 用畸变单元才能区分；正解必须过畸变对照。")
