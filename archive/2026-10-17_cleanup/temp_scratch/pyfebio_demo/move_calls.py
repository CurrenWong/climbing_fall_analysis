"""把 top/bottom 面的**调用**移到 blocks 建好之后。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

CALLS = '''top_tibia = extreme_faces("tibia_cort", "max")
top_fibula = extreme_faces("fibula_cort", "max")
bot_calc = extreme_faces("calcaneus_cort", "min")
print(f"顶面：胫骨 {len(top_tibia)} 面 / 腓骨 {len(top_fibula)} 面；"
      f"跟骨底面 {len(bot_calc)} 面")
'''
assert s.count(CALLS) == 1, f"calls1 {s.count(CALLS)}"
s = s.replace(CALLS, "")

CALLS2 = '''A_tib, A_fib = area_proj(top_tibia), area_proj(top_fibula)
share = A_fib / (A_tib + A_fib) if FIB_SHARE < 0 else FIB_SHARE
print(f"顶面投影面积：胫骨 {A_tib:.1f} mm²，腓骨 {A_fib:.1f} mm² "
      f"=> 腓骨分担 {share*100:.1f}% 载荷")
'''
assert s.count(CALLS2) == 1, f"calls2 {s.count(CALLS2)}"
s = s.replace(CALLS2, "")

LIGCALL = 'tal_faces = extreme_faces("talus_cort", "min")'
assert s.count(LIGCALL) == 1, f"ligcall {s.count(LIGCALL)}"
s = s.replace(LIGCALL, "")

anchor = "w('  <Mesh>')"
assert s.count(anchor) == 1
s = s.replace(anchor, '''# ── 面的选取（必须在 blocks 建好之后）──
top_tibia = extreme_faces("tibia_cort", "max")
top_fibula = extreme_faces("fibula_cort", "max")
bot_calc = extreme_faces("calcaneus_cort", "min")
print(f"外表面顶/底：胫骨顶 {len(top_tibia)} / 腓骨顶 {len(top_fibula)} / "
      f"跟骨底 {len(bot_calc)}")

A_tib, A_fib = area_proj(top_tibia), area_proj(top_fibula)
share = A_fib / (A_tib + A_fib) if FIB_SHARE < 0 else FIB_SHARE
print(f"顶面投影面积：胫骨 {A_tib:.1f} mm²，腓骨 {A_fib:.1f} mm² "
      f"=> 腓骨分担 {share*100:.1f}% 载荷")
tal_faces = extreme_faces("talus_cort", "min")

w('  <Mesh>')''')

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
