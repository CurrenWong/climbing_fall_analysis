"""解读 FEBio -output_negative_jacobians 输出：首个翻转发生在哪个时间步/迭代、哪些域。"""
import re
from collections import Counter
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\temp/pyfebio_demo/a1_matrix/calcaneus_negjac.txt")
s = P.read_text(encoding="utf-8", errors="replace")
lines = s.split("\n")

# elem id -> 域
import io
feb = Path(r"D:\Project\climbing_fall_analysis\temp/thums/febio_lig/ankle_lig_WK.feb")
ft = feb.read_text(encoding="latin-1")
rng = []
for m in re.finditer(r'<Elements type="(\w+)" name="([^"]+)">(.*?)</Elements>',
                     ft, re.S):
    et, name, body = m.groups()
    ids = [int(x) for x in re.findall(r'<elem id="(\d+)">', body)]
    rng.append((min(ids), max(ids), name))


def dom(eid):
    for a, b, n in rng:
        if a <= eid <= b:
            return n
    return f"?{eid}"


hits = [(i, l) for i, l in enumerate(lines, 1)
        if "Negative jacobian was detected at element" in l]
print(f"总元件级负 Jacobian 报告 {len(hits)} 条")
first = hits[0][0]
# 向上找最近的 time step / iteration 上下文
ctx_start = max(0, first - 40)
print("\n--- 首次翻转的上下文 ---")
print("\n".join(f"{i}: {lines[i-1]}" for i in range(ctx_start + 1, first + 2)))

ids = [int(re.search(r"element (\d+)", l).group(1)) for _, l in hits]
c = Counter(dom(e) for e in ids)
print("\n翻转元件按域统计（全部报告）:")
for k, v in c.most_common():
    print(f"   {k:34s} {v}")
# 每个域的翻转元件（去重）
for k in c:
    es = sorted({e for e in ids if dom(e) == k})
    print(f"   {k}: 元件号 {es[:20]}{' ...' if len(es) > 20 else ''}（共 {len(es)} 个）")

# 每次 "N negative jacobians detected" 的时间步归属
marks = [(i, l) for i, l in enumerate(lines, 1)
         if re.search(r"\d+ negative jacobians detected|Negative jacobian detected\.", l)]
step = None
print("\n每次批量报告的时间步:")
for i, l in marks:
    st = None
    for j in range(i - 1, 0, -1):
        if "beginning time step" in lines[j - 1]:
            st = lines[j - 1]
            break
    print(f"   L{i}: {l.strip()[:60]}  ← {st}")
