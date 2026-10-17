import io
import re

t = io.open(r"D:\Project\climbing_fall_analysis\temp\thums\febio\ankle_tibiotalar.log",
            encoding="latin-1", errors="replace").read()
i = t.find("CONTACT INTERFACE DATA")
print("=== CONTACT INTERFACE DATA 段 ===")
print(t[i:i + 2200])

print("\n=== 所有被截断的警告（含 'runs better' 之类）===")
for m in re.finditer(r"runs better", t):
    s = max(0, m.start() - 700)
    print(t[s:m.start() + 700])
    print("-" * 70)
