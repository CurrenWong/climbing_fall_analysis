"""修正 _fixentry.py 的断言（结尾反引号会误判），然后跑。"""
from pathlib import Path

f = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_fixentry.py")
s = f.read_text(encoding="utf-8")
old = 'assert "`build_feb.py`" not in t and "`cartilage.py`" not in t, "开关表裸名残留"'
new = ('# 注意：不能查 "`build_feb.py`" —— 全路径 `scripts/ankle_fe/build_feb.py` '
       '的结尾\n# 本身就含这个子串，会误判。要查的是「整格」形式。\n'
       'assert ("| `build_feb.py` |" not in t) and '
       '("| `cartilage.py` |" not in t), "开关表仍有裸文件格"')
assert s.count(old) == 1
f.write_text(s.replace(old, new), encoding="utf-8")
print("断言已修正")
