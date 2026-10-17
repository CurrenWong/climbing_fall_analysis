"""让 Control 关键项可由环境变量调节，便于收敛调参。

FEBio 手册 §3.3 的判据：有 dtol(位移)/etol(能量)/rtol(残差) 三个，
满足其一即可。原配置 etol=0.001 偏严 + dtol=0.01 ⇒ 每步上百次迭代。
"""
from pathlib import Path
import ast
import re

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = 'TOPFIX = os.environ.get("TOPFIX", "lateral").strip() or "lateral"'
assert s.count(old) == 1
s = s.replace(old, old + '''
# 收敛容差/迭代上限（env 可调，便于调参）
ETOL = os.environ.get("ETOL", "0.001").strip() or "0.001"
DTOL = os.environ.get("DTOL", "0.01").strip() or "0.01"
MAX_REFS = os.environ.get("MAX_REFS", "50").strip() or "50"
MAX_UPS = os.environ.get("MAX_UPS", "10").strip() or "10"''')

# Control 块里的常量改成变量
for oldv, newv in (("<etol>0.001</etol>", "<etol>{ETOL}</etol>"),
                   ("<dtol>0.01</dtol>", "<dtol>{DTOL}</dtol>"),
                   ("<max_refs>50</max_refs>", "<max_refs>{MAX_REFS}</max_refs>"),
                   ("<max_ups>10</max_ups>", "<max_ups>{MAX_UPS}</max_ups>")):
    if oldv in s:
        s = s.replace(oldv, newv)
        print("  patch", oldv, "->", newv)
    elif newv in s:
        print("  already", newv)
    else:
        print("  !! not found", oldv)

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("syntax OK")

# 确认 Control 块是 f-string（否则 {ETOL} 不会被替换）
i = s.find("<Control>")
seg = s[max(0, i - 400):i + 60]
print("\n--- Control 前的上下文 ---")
print(seg[-500:])
