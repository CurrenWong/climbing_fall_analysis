"""补上缺失的 sv_penta6（fix_sv.py 的 patch ① 实际没落地）。"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

if "def sv_penta6(" in s:
    print("已存在，跳过")
else:
    anchor = "def detJ_penta15(P):"
    assert s.count(anchor) == 1, s.count(anchor)
    s = s.replace(anchor, '''def sv_penta6(P):
    """楔形**带符号体积**（依赖节点顺序，用于朝向判定/修正）。

    J = [P1-P0, P2-P0, dt]，dt = (P3+P4+P5-P0-P1-P2)/3
    ★ 交换上下三角（3 次对换 = 奇置换）会正确变号。
      采样版 detJ_penta6 是按形状函数求和、**不能**可靠地靠调换节点来变号，
      所以朝向判定必须另用本函数 —— 这正是此前负朝向楔形「修不动、只能剔」
      的原因（剔除率 tet4 2% vs penta6 14%）。
    """
    dt = (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3.0
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0], dt])))


''' + anchor)
    p.write_text(s, encoding="utf-8")
    ast.parse(s)
    print("已插入 sv_penta6，syntax OK")

# 复核
s = p.read_text(encoding="utf-8")
print("sv_penta6 def 出现次数:", s.count("def sv_penta6("))
print("sv_penta6 调用次数:", s.count("sv_penta6("))
