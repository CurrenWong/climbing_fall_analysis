"""读 .feb 的 <Contact> 段与求解器收敛参数，定位 max_aug / rtol。"""
import io
import re
from pathlib import Path

OUT = Path(r"D:\Project\climbing_fall_analysis\temp/pyfebio_demo/a1a2")
f = OUT / "E1_ts40.feb"
t = f.read_text(encoding="utf-8", errors="replace")

print("=== <Contact> 段 ===")
for m in re.finditer(r"<contact[^>]*>.*?</contact>", t, re.S):
    print(m.group(0)[:700])
    print("-" * 60)
print("\n=== 收敛/求解器相关参数（在任意段里搜）===")
for key in ("max_aug", "augtol", "augtoltype", "rtol", "ptol", "etol",
            "dtol", "mactol", "max_refs", "max_ups", "diverge", "line_search",
            "lstol", "ncut"):
    for m in re.finditer(rf"<{key}[^>]*>[^<]*</{key}>", t):
        print("  ", m.group(0)[:120])
print("\n=== <Control> 段 ===")
m = re.search(r"<Control>.*?</Control>", t, re.S)
print(m.group(0)[:600] if m else "无")
