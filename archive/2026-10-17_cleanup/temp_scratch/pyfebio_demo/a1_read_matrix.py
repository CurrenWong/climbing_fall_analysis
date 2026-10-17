import io
import re

for t in ("talus", "calcaneus", "both", "nocontact"):
    s = io.open(rf"temp/pyfebio_demo/a1_matrix/{t}.febio.log",
                encoding="utf-8", errors="replace").read()
    print("=" * 20, t)
    print(" 头部警告:", [x.strip() for x in
                         re.findall(r"\*\s{2,}(.+?)\s{2,}\*", s)][:8])
    i = s.find("beginning time step 1")
    seg = s[i:i + 800].replace("\t", "")
    print(seg)
    j = s.find("beginning time step 2")
    if j < 0:
        print("  >> 无 step2（step1 就失败）")
        k = s.find("Negative jacobian")
        print("  首个 negative jacobian 距 step1 起点:",
              (k - i) if k > 0 else "未出现")
