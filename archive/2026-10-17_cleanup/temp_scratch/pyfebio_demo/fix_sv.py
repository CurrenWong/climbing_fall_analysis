"""修 detJ_penta6：把「朝向判定」与「质量判定」分开。

★ 根因：采样版 detJ（Σ xᵢ⊗∇Nᵢ）**对节点置换不变** ⇒ 闸门里
  「交换上下三角翻符号」对该量**永远无效** ⇒ 负朝向的楔形只能被剔除。
  实测：m=4 的 122/150 个楔形 det 为负（四边形底的绕向与偏移方向相反），
  它们本可被修正，却全被剔除 ⇒ 剔除率 tet4 2% vs penta6 14%。

修法：
  · **朝向判定**用「带符号体积」`det([P1-P0, P2-P0, dt])` —— 它**依赖节点顺序**，
    交换上下三角会正确变号
  · **质量判定**仍用采样版（置换不变，正好只反映形状畸变）
  · 两者都过才保留：朝向已修正 且 采样最小 det > 0（FEBio 的硬要求）
"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# ① 加一个「带符号体积」函数
anchor = "def detJ_penta15(P):"
assert s.count(anchor) == 1
s = s.replace(anchor, '''def sv_penta6(P):
    """楔形**带符号体积**（依赖节点顺序，用于朝向判定/修正）。

    J = [P1-P0, P2-P0, dt]，dt = (P3+P4+P5-P0-P1-P2)/3
    ★ 交换上下三角（3 次对换，奇置换）会正确变号 —— 采样版做不到这点。
    """
    dt = (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3.0
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0], dt])))


''' + anchor)

# ② 闸门：penta6/penta15 用「朝向 + 采样质量」双判据
old = """        P = XYZ[[i - 1 for i in ids]]
        d = fn(P)
        if d > DET_MIN:"""
assert s.count(old) == 1
s = s.replace(old, """        P = XYZ[[i - 1 for i in ids]]
        d = fn(P)
        if et in ("penta6", "penta15"):
            # ★ 采样版对置换不变 ⇒ 必须另用带符号体积判朝向
            if sv_penta6(P[:6]) < 0:
                d = -1e9          # 朝向反 ⇒ 走修正分支
        if d > DET_MIN:""")

# ③ 修正分支：penta6 交换上下三角后，用带符号体积复核符号
old2 = """        if et == "tet4":
            alt = [ids[1], ids[0]] + ids[2:]                 # 1 次对换 = 奇
        elif et == "penta15":
            # ★ 15 节点的奇置换不好构造：
            #   交换上下三角 = 角 3 对换 + 边 3 对换 = 6 次 = **偶**，符号不变；
            #   再加对换又会把中点节点挪到别的棱上（几何就变了）。
            #   penta6 的构造本身已保证楔形朝向正确 ⇒ 继承角点朝向，直接剔除。
            alt = None"""
assert s.count(old2) == 1
s = s.replace(old2, """        if et == "tet4":
            alt = [ids[1], ids[0]] + ids[2:]                 # 1 次对换 = 奇
        elif et == "penta6":
            alt = ids[3:] + ids[:3]                          # 换上下三角 = 奇
        elif et == "penta15":
            # 15 节点的奇置换会把中点节点挪到别的棱上（几何就变了）
            # ⇒ 不修正，直接剔除
            alt = None""")

# ④ 修正分支的判定：penta6 用带符号体积 + 采样质量
old3 = """        d2 = (fn(XYZ[[i - 1 for i in alt]])
              if alt is not None else -1.0)"""
assert s.count(old3) == 1
s = s.replace(old3, """        if alt is None:
            d2 = -1.0
        else:
            P2 = XYZ[[i - 1 for i in alt]]
            d2 = fn(P2)
            if et == "penta6":
                d2 = (min(d2, DET_MIN + 1e-6)
                      if sv_penta6(P2[:6]) > 0 else -1e9)""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched 4 处，syntax OK")
