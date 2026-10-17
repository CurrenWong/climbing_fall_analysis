"""让 ankle_lig_feb.py 导出「韧性读数所需」的最小数据：
  · fe_id -> thums_id 映射
  · 每条韧带的杆对（FE 编号）+ 名称 + EA
这样读数脚本可以独立于组装脚本工作。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

anchor = 'feb = WORK / f"ankle_lig_{_tag}.feb"'
assert s.count(anchor) == 1
s = s.replace(anchor, '''# ── 导出读数所需数据（供 read_lig.py 用）──
_dump = {}
_dump["n_node"] = np.array(len(coords))
_fe = np.zeros(len(coords) + 1, dtype=np.int64)
for _k, _v in node_id_of.items():
    if isinstance(_k, int):
        _fe[_v] = _k
_dump["fe_to_thums"] = _fe
_dump["axis"] = axis
_names, _eas, _pairs, _npair = [], [], [], []
for r in lig_recs:
    pairs = r.get("trusses") or []
    fp = [[node_id_of[int(a)], node_id_of[int(b)]] for a, b in pairs]
    _names.append(r["name"])
    _eas.append(r["ea"])
    _pairs.append(fp)
    _npair.append(len(fp))
_mx = max(_npair) if _npair else 0
_TR = -np.ones((len(_names), max(_mx, 1), 2), dtype=np.int64)
for i, fp in enumerate(_pairs):
    for j, pr in enumerate(fp):
        _TR[i, j] = pr
_dump["lig_names"] = np.array(_names)
_dump["lig_ea"] = np.array(_eas)
_dump["lig_pairs"] = _TR
_dump["lig_npair"] = np.array(_npair)
np.savez(WORK / f"lig_read_{_tag}.npz", **_dump)
print(f"读数数据 -> {WORK / f'lig_read_{_tag}.npz'}")

feb = WORK / f"ankle_lig_{_tag}.feb"''')

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
