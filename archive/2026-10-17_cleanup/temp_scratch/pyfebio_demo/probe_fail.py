"""看外翻(−5)失败的求解日志：负 Jacobian 在哪、几何 detJ 如何。"""
import io
import re

R = r"D:\Project\climbing_fall_analysis"
for tag in ("IE_500", "IE_-500"):
    p = R + rf"\temp\pyfebio_demo\_solve_{tag}.log"
    try:
        t = io.open(p, encoding="utf-8", errors="replace").read()
    except OSError:
        print(f"✗ 无 {p}")
        continue
    print("=" * 78)
    print(f"=== {tag}   日志 {len(t)} 字符")
    print("=" * 78)
    # 几何诊断行
    for i, l in enumerate(t.splitlines()):
        if re.search(r"det\(?J|detJ|最小.*行列式|负体积|inverted|negative jacobian|退化",
                     l, re.I):
            print(f"  {i+1:>5}: {l.strip()[:132]}")
    # 失败点上下文
    for m in re.finditer(r"[^\n]*negative jacobian[^\n]*", t, re.I):
        s = max(0, m.start() - 250)
        print("\n  --- negative jacobian 上下文 ---")
        print("  " + t[s:m.end() + 150].replace("\n", "\n  ")[:900])
    # 收敛/时间步信息
    for m in re.finditer(r"[^\n]*(Time step|converge|Failed|error)[^\n]*", t[:30000], re.I):
        s = m.group(0).strip()
        if s:
            print(f"  ·  {s[:130]}")
    break_ = None
