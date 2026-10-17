"""加 TRUNC_MM：把胫骨/腓骨截断到踝关节近端 TRUNC_MM 以内。

理由（两层）：
  ① 速度：胫骨一整根带 25,000+ 单元（spon 6651 + 骨髓 14482），
     每步牛顿迭代都在解 32k 方程 —— 一轮要 20 分钟，没法迭代。
  ② 物理：**踝关节模型不需要整根胫骨**。吴恺的模型就是"足踝部"。
     截断后载荷加在截面上，与"胫腓骨上截面 600 N"的表述一致。
截断后新的顶面由「外表面 + 法向对齐」逻辑自动识别，不需额外处理。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = 'FIX_NAVICULAR = os.environ.get("FIX_NAVICULAR", "1") == "1"'
assert s.count(old) == 1
s = s.replace(old, old + '''
# 胫骨/腓骨截断到距踝关节近端 TRUNC_MM 以内（0 = 不截断）
TRUNC_MM = float(os.environ.get("TRUNC_MM", "0"))''')

old2 = 'print(f"[mode] DRIVER={DRIVER} DISP_MM={DISP_MM}")'
assert s.count(old2) == 1
s = s.replace(old2, old2 + '\nprint(f"[mode] TRUNC_MM={TRUNC_MM}")')

# 截断：在构建 blocks 之后、Jacobian 修正之前，按「沿轴投影」过滤胫/腓单元
anchor = "lig_blocks = {}\nLIG_ETYPE ="
assert s.count(anchor) == 1
s = s.replace(anchor, '''# ── 胫/腓骨截断 ──────────────────────────────────────────────
if TRUNC_MM > 0:
    import numpy as _np
    _tal = _np.array([coords[i][:] for i in
                      sorted({c for k, v in blocks.items() for c in v
                              for c in c})[:0]] or [[0, 0, 0]], float)
    # 以距骨质心为参考点（踝关节位置）
    _tal_ids = sorted(dom_nodes["talus_cort"]) if "talus_cort" in dom_nodes \\
        else []
    _ref = _np.mean([[coords[i - 1][0], coords[i - 1][1], coords[i - 1][2]]
                     for i in _tal_ids], axis=0) if _tal_ids else _np.zeros(3)
    _keep = 0
    _drop = 0
    for _k in list(blocks):
        _dom = _k[0]
        if _dom.startswith(("tibia_", "fibula_")):
            _new = []
            for _c in blocks[_k]:
                _pts = _np.array([[coords[i - 1][0], coords[i - 1][1],
                                   coords[i - 1][2]] for i in _c])
                _cen = _pts.mean(0) - _ref
                if float(_cen @ axis) <= TRUNC_MM:
                    _new.append(_c)
                else:
                    _drop += 1
            _keep += len(_new)
            blocks[_k] = _new
    print(f"[trunc] 截断到 {TRUNC_MM} mm 内：保留 {_keep:,} 剔除外侧 {_drop:,}")

lig_blocks = {}
LIG_ETYPE =''')

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
