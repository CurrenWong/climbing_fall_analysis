"""DET_MIN 改为环境变量（楔形单元对偏移扭曲更敏感，需按需放宽）。"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = "DET_MIN = 0.5      # ★ 不能只 >0：近退化薄片会最先翻负，按同类量级设阈"
assert s.count(old) == 1, s.count(old)
s = s.replace(old, '''DET_MIN = float(os.environ.get("DET_MIN", "0.5"))
# ★ 不能只 >0：近退化薄片会最先翻负。但**闸门太严会把承力单元也挖掉**
#   （实测 DET_MIN=0.5 时软骨楔形被剔 14+29 个，接触面 150/198→138/171，
#     同位移下支反力从 259 N 掉到 117 N）⇒ 可按需放宽。
#   FEBio 自身只要求积分点 det>0，故放宽到 0.02~0.1 通常仍能初始化。''')
p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched，syntax OK")
