import re
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis")
t = (P / "ANKLE_FE_TODO.md").read_text(encoding="utf-8")
L = t.split("\n")
print("ANKLE_FE_TODO.md 中含「阻塞」的行：")
for i, l in enumerate(L, 1):
    if "阻塞" in l:
        print(f"  {i:3d} (前 12 行内={i <= 12})| {l.strip()[:90]}")
print("\n说明：若不在前 12 行（superseded 提示之外），则是历史正文，不影响接手。")

# 修好的路径正则（允许 /）
f = P / "temp/pyfebio_demo/verify_handover.py"
s = f.read_text(encoding="utf-8")
s = s.replace(r'r"([\w\-.]+\.(?:py|md|json|npz|zip|k|fe?b))"',
              r'r"([\w\-./]+\.(?:py|md|json|npz|zip|k|fe?b))"')
f.write_text(s, encoding="utf-8")
print("已修 verify_handover.py 的路径正则（允许 /）")
