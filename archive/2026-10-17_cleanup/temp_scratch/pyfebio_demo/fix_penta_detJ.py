"""修正 detJ_penta6 公式。

楔形 (r,s,t)，r+s<=1 的形状函数：
  N0=(1-r-s)(1-t) N1=r(1-t) N2=s(1-t)
  N3=(1-r-s)t     N4=r t     N5=s t
在底面形心 (1/3,1/3,0)：
  dN/dr = (-1,+1,0, 0,0,0)
  dN/ds = (-1,0,+1, 0,0,0)
  dN/dt = (-1/3,-1/3,-1/3, +1/3,+1/3,+1/3)
⇒ J 的三列 = [P1-P0, P2-P0, (P3+P4+P5-P0-P1-P2)/3]
✗ 旧写成 P3-P0 —— 偏移方向在三角内变化时结果不同，正是 FEBio 报
  "Negative jacobian during mesh initialization" 的原因。
"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = '''def detJ_penta6(P):
    """楔形：底三角 0,1,2 + 顶三角 3,4,5；用 P0 出发的三条边。
    符号 = 底三角绕向 x (P3-P0) 方向；为负时按奇置换翻转。"""
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0],
                                         P[3] - P[0]])))'''
assert s.count(old) == 1, s.count(old)
s = s.replace(old, '''def detJ_penta6(P):
    """楔形在**底面形心**处的 Jacobian 行列式。

    形状函数 N0=(1-r-s)(1-t), N1=r(1-t), N2=s(1-t),
             N3=(1-r-s)t,     N4=r t,     N5=s t
    ⇒ J 的三列 = [P1-P0, P2-P0, (P3+P4+P5-P0-P1-P2)/3]
    ★ dt 列必须是**上下三角形心差**；写成 P3-P0 会在偏移方向沿三角变化时
      算出正值而 FEBio 报 Negative jacobian（实测踩过）。
    """
    dt = (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3.0
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0], dt])))''')

# 朝向修正分支也要求 >DET_MIN，保持口径一致
old2 = """        d2 = fn(XYZ[[i - 1 for i in alt]])
        if d2 > 0:"""
assert s.count(old2) == 1
s = s.replace(old2, """        d2 = fn(XYZ[[i - 1 for i in alt]])
        if d2 > DET_MIN:""")
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched，syntax OK")
