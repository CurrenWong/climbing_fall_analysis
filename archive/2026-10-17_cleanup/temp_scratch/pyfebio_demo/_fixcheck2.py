import io
import re
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis")
f = P / "temp/pyfebio_demo/verify_handover.py"
s = f.read_text(encoding="utf-8")
# 允许反斜杠与正斜杠的完整相对路径
s = re.sub(r"r\"\[\\w[^/]*\]\+", 'r"[\\\\w\\\\-./\\\\\\\\]+"', s)
f.write_text(s, encoding="utf-8")
print("patched:", 'r"[\\\\w' in s or "[\\w\\-./\\\\]+" in s)

# 直接手工确认那 4 个裸名在入口节里的实际上下文
head_t = (P / "FE_PIPELINE_TODO.md").read_text(encoding="utf-8")
head = head_t[head_t.index("# 🧭 接手入口"):head_t.index("## C. 解出 THUMS")]
for name in ("bone_materials.json", "build_feb.py", "cartilage.py",
             "thums_lowerlimb.npz"):
    for line in head.split("\n"):
        if name in line and "`" in line:
            hits = re.findall(r"`([^`]*" + re.escape(name) + r")`", line)
            for h in hits:
                p = P / h.replace("\\", "/")
                print(f"  {'OK ' if p.exists() else 'MISS'} {h}")
            break
