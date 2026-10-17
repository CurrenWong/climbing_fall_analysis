"""detJ_penta6 改为真实形函数扫描（覆盖三角形内部 + 轴向）。

形状函数：N0=(1-r-s)(1-t) N1=r(1-t) N2=s(1-t)
          N3=(1-r-s)t     N4=r t     N5=s t
J = Σ P_i ⊗ [dN/dr, dN/ds, dN/dt]
在 (r,s,t) 网格上取 min det —— 单点检查会漏掉三角内的退化。
"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# 找到现有 detJ_penta6 整段并替换
start = s.index("def detJ_penta6(P):")
end = s.index("FN = {", start)
new = '''def detJ_penta6(P):
    """楔形 Jacobian 最小值：在 (r,s,t) 采样网格上取 min det(J)。

    形状函数 N0=(1-r-s)(1-t), N1=r(1-t), N2=s(1-t),
             N3=(1-r-s)t,     N4=r t,     N5=s t
    J = Σ_i P_i ⊗ [dN/dr, dN/ds, dN/dt]

    ★ 只在底面形心算会漏掉：三角**内部**或沿 t 的退化会让 FEBio 在积分点
      报 "Negative jacobian during mesh initialization"，而单点检查算出来是正的。
    """
    RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3)
    TS = (0.0, 0.25, 0.5, 0.75, 1.0)
    worst = None
    for r in RS:
        for sq in RS:
            if r + sq > 1.0 + 1e-9:
                continue
            for t in TS:
                d = np.array([
                    [-(1 - t), -(1 - t), -(1 - r - sq)],
                    [1 - t, 0.0, -r],
                    [0.0, 1 - t, -sq],
                    [-t, -t, (1 - r - sq)],
                    [t, 0.0, r],
                    [0.0, t, sq],
                ])
                J = P.T @ d                      # (3,3)
                det = float(np.linalg.det(J))
                if worst is None or det < worst:
                    worst = det
    return worst


'''
s = s[:start] + new + s[end:]
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched，syntax OK")
