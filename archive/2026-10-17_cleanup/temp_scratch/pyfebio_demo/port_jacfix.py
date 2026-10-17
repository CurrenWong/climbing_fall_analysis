"""把 build_feb.py 的 Jacobian 自检/修正/剔除 环节搬进 ankle_lig_feb.py。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

anchor = """lig_blocks = {}
if LIG_MODE == "shell":
    for r in lig_recs:
        lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in q]
                                 for q in r["quads"]]

w('  <Mesh>')"""
assert s.count(anchor) == 1, f"anchor hit {s.count(anchor)}"

FIX = '''lig_blocks = {}
if LIG_MODE == "shell":
    for r in lig_recs:
        lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in q]
                                 for q in r["quads"]]

# ── ★ Jacobian 自检 + 朝向修正 + 剔除 ────────────────────────────────
# 这一步在 build_feb.py 里就有，第一版新脚本**漏搬**了 ⇒ FEBio 直接报
# "Negative jacobian detected during mesh initialization"。
# 教训（LRN-20261003-018）：**复用已验证管线的结构时，不能只搬"看得见"的部分**
# —— 这类"清理/修正"环节最容易被漏掉，而它恰恰是能不能跑起来的关键。
GP2 = 0.577350269189626


def detJ_hex8(P):
    """P:(8,3) 在 8 个 Gauss 点算 det(J)，返回最小值。"""
    xi = np.array([[-GP2, -GP2, -GP2], [GP2, -GP2, -GP2], [GP2, GP2, -GP2],
                   [-GP2, GP2, -GP2], [-GP2, -GP2, GP2], [GP2, -GP2, GP2],
                   [GP2, GP2, GP2], [-GP2, GP2, GP2]])
    ref = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    worst = 1e30
    for sg in xi:
        dN = np.zeros((8, 3))
        for i in range(8):
            a, b, c = ref[i]
            dN[i, 0] = a * (1 + b * sg[1]) * (1 + c * sg[2]) / 8
            dN[i, 1] = b * (1 + a * sg[0]) * (1 + c * sg[2]) / 8
            dN[i, 2] = c * (1 + a * sg[0]) * (1 + b * sg[1]) / 8
        worst = min(worst, float(np.linalg.det(P.T @ dN)))
    return worst


def detJ_tet4(P):
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0],
                                         P[3] - P[0]])))


FN = {"hex8": detJ_hex8, "tet4": detJ_tet4}
DET_MIN = 0.5      # ★ 不能只 >0：近退化薄片会最先翻负，按同类量级设阈
bad_total = fixed_total = 0
print("\\n=== Jacobian 自检 + 朝向修正（det(J) 必须 > 0）===")
for key in list(blocks):
    dom, et = key
    fn = FN[et]
    keep, nbad, dmins = [], 0, []
    for ids in blocks[key]:
        P = XYZ[[i - 1 for i in ids]]
        d = fn(P)
        if d > DET_MIN:
            keep.append(ids); dmins.append(d); continue
        # ★ 修正朝向必须用**奇置换**：
        #   tet4 = 交换两个节点（1 次交换）；hex8 = 交换前后两半。
        #   ✗ 整体倒序是 3 次交换（偶置换），行列式符号**不变**！
        alt = ([ids[1], ids[0]] + ids[2:]) if et == "tet4" else \\
              (ids[len(ids) // 2:] + ids[:len(ids) // 2])
        d2 = fn(XYZ[[i - 1 for i in alt]])
        if d2 > 0:
            keep.append(alt); dmins.append(d2); fixed_total += 1
        else:
            nbad += 1
    blocks[key] = keep
    bad_total += nbad
    mn = min(dmins) if dmins else float("nan")
    print(f"  {'OK ' if not nbad else '!! '} {dom:22s} {et:6s} "
          f"保留 {len(keep):>6,} 剔除 {nbad:>4,}  det(J) min={mn:.4e}")
print(f"  翻转修正 {fixed_total} 个；仍无法修正而剔除 {bad_total} 个")

# 韧带壳：quad4 的退化检查（绕向是任意的，只查退化，不翻转）
if LIG_MODE == "shell":
    nbad_shell = 0
    for nm in list(lig_blocks):
        keep = []
        for ids in lig_blocks[nm]:
            P = XYZ[[i - 1 for i in ids]]
            d = np.cross(P[2] - P[0], P[3] - P[1])
            if np.linalg.norm(d) > 1e-6:
                keep.append(ids)
            else:
                nbad_shell += 1
        lig_blocks[nm] = keep
    print(f"  韧带壳退化剔除 {nbad_shell} 个")

w(\'  <Mesh>\')'''

s = s.replace(anchor, FIX)
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
