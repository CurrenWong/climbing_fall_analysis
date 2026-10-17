"""重写闸门块（干净版）：朝向与质量两判据分开。

  朝向  sv_penta6(P) > 0     —— 依赖节点顺序，交换上下三角可修正
  质量  fn(P) > 0            —— 采样版，反映形状畸变（FEBio 硬要求）

对 penta15：只用采样值 + 不修正（奇置换不好构造）。
"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

i0 = s.index("    for ids in blocks[key]:")
i1 = s.index("    blocks[key] = keep")
new = '''    for ids in blocks[key]:
        P = XYZ[[i - 1 for i in ids]]

        def _ok(PP):
            """朝向 + 质量双判据。"""
            if et in ("penta6", "penta15"):
                if sv_penta6(PP[:6]) <= 0:
                    return None
            return FN[et](PP)

        d = _ok(P)
        if d is not None and d > DET_MIN:
            keep.append(ids); dmins.append(d); continue
        # ★ 朝向修正必须用**奇置换**：
        #   tet4 = 交换两个节点（1 次对换）
        #   penta6 = 交换上下三角（3 次对换）
        #   ✗ 整体倒序是偶置换，行列式符号不变
        if et == "tet4":
            alt = [ids[1], ids[0]] + ids[2:]
        elif et == "penta6":
            alt = ids[3:] + ids[:3]
        elif et == "penta15":
            alt = None      # 15 节点的奇置换会把中点挪到别的棱（几何变了）
        else:
            alt = ids[len(ids) // 2:] + ids[:len(ids) // 2]
        d2 = _ok(XYZ[[i - 1 for i in alt]]) if alt is not None else None
        if d2 is not None and d2 > DET_MIN:
            keep.append(alt); dmins.append(d2); fixed_total += 1
        else:
            nbad += 1
'''
s = s[:i0] + new + s[i1:]
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched，syntax OK")
