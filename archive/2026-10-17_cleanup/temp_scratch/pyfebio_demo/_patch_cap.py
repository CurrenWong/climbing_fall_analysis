"""清理 cartilage.py 参数区（有一段被失败 edit 弄坏的重复块），并加对侧面限制。"""
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\cartilage.py")
lines = P.read_text(encoding="utf-8").split("\n")

# 找到两段重复块的边界：第一段从 "# ★ FILL_FRAC = 0.5" 起，第二段从干净的
# "FILL_FRAC = 0.5" 赋值起、到 "drop whole cells with a sliver tet" 行止。
i0 = next(i for i, l in enumerate(lines) if l.startswith("# ★ FILL_FRAC = 0.5"))
i1 = next(i for i, l in enumerate(lines) if l.startswith("FILL_FRAC = 0.5"))
i2 = next(i for i, l in enumerate(lines)
          if l.startswith("MIN_TET_V = 0.05 if _Q else 1e-9")
          and "sliver tet" in l and i > i1)
print(f"删除行 {i0+1}..{i1+1}（坏块），保留 {i1+1}..{i2+1}")

clean = [
    "# ★ FILL_FRAC = 0.5: the two layers exactly fill the gap and are tangent,",
    "#   NOT interpenetrating.",
    "#   (0.6 / 0.52 were tried to make contact active at t=0; initial",
    "#    interpenetration fires a contact force at t=0 and the first Newton step",
    "#    crushes 4-5 elements. The real fix for rigid-body modes is clamping the",
    "#    loaded cross-section -- see the disp driver in build_feb.py.)",
] + lines[i1:i2 + 1] + [
    "",
    "# ---- cap each layer by the ALREADY-BUILT partner offset surface ----------",
    "# Each side used t = 0.5 * (its own local gap) and both offset toward the",
    "# *nearest* point of the opposite bone. Where the two sides are not",
    "# geometrically symmetric, one layer pokes into the other at t=0; penalty",
    "# contact then fires a huge force and crushes elements in step 1",
    "# (reproduced: the 2-joint model dies with '2 negative jacobians' at ~zero",
    "# load). Capping by the partner's offset surface removes it by construction.",
    "# Order matters: the partner must be generated first (joint list is A,B,A,B).",
    'PARTNER = {"tibiotalar_tibia": "tibiotalar_talus",',
    '           "tibiotalar_talus": "tibiotalar_tibia",',
    '           "subtalar_talus": "subtalar_calcaneus",',
    '           "subtalar_calcaneus": "subtalar_talus"}',
    'CAP_SAFETY = float(os.environ.get("CAP_SAFETY", "0.97") or 0.97)',
    "_OFFSET_PTS: dict[str, np.ndarray] = {}",
]
lines = lines[:i0] + clean + lines[i2 + 1:]
t = "\n".join(lines)

# ② thickness 后加 cap
old = """    # 厚度：min(t_target, FILL_FRAC × gap)
    thick = {v: min(T_TARGET, FILL_FRAC * gv) for v, gv in g.items()}"""
assert t.count(old) == 1, f"anchor2 x{t.count(old)}"
t = t.replace(old, """    # 厚度：min(t_target, FILL_FRAC × gap)
    thick = {v: min(T_TARGET, FILL_FRAC * gv) for v, gv in g.items()}
    # ★ 追加「对侧 offset 面」上限，避免初始过盈（见上方 PARTNER 注释）
    _pj = PARTNER.get(j)
    if CAP_SAFETY > 0 and _pj in _OFFSET_PTS:
        from scipy.spatial import cKDTree as _KD
        _pts = np.array([xyz[idx[v]] for v in art])
        _d, _ = _KD(_OFFSET_PTS[_pj]).query(_pts)
        _ncap = 0
        for _i, v in enumerate(art):
            _lim = CAP_SAFETY * float(_d[_i])
            if _lim < thick[v]:
                thick[v] = _lim
                _ncap += 1
        print(f"    对侧({_pj})面限制: {_ncap}/{len(art)} 节点被压薄 (最近 {dmin:.3f} mm)"
              if (dmin := float(_d.min())) is not None else "")""")

# ③ 生成后登记 offset 点
old2 = """    if not base_faces:
        print("    ⚠️ 没有完整落在关节面内的面片，跳过")
        continue"""
assert t.count(old2) == 1, f"anchor3 x{t.count(old2)}"
t = t.replace(old2, old2 + '\n    _OFFSET_PTS[j] = np.asarray(local_pts, float)')

P.write_text(t, encoding="utf-8")
import ast
ast.parse(t)
print("OK：参数区去重 + 对侧面限制已加入")
