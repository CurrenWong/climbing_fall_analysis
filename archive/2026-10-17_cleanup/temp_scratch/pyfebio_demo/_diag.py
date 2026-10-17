"""精确诊断：E1 的 <solver> 段 + 失败消息 + 接触/弹簧的当前状态。"""
import io
import re
from pathlib import Path

OUT = Path(r"D:\Project\climbing_fall_analysis\temp/pyfebio_demo/a1a2")

for tag in ("E1_ts40", "E2_A1A2", "E3_a2_only"):
    f = OUT / f"{tag}.feb"
    if not f.exists():
        print(f"缺 {tag}.feb"); continue
    t = f.read_text(encoding="utf-8", errors="replace")
    print("=" * 74)
    print(f"### {tag}")
    m = re.search(r"<solver>.*?</solver>", t, re.S)
    print("  <solver>:", (m.group(0)[:400].replace("\n", " ") if m else "无"))
    # 失败消息
    lg = OUT / f"{tag}.febio.log"
    if lg.exists():
        lt = lg.read_text(encoding="utf-8", errors="replace")
        tail = lt[-1600:]
        # 找最后的 error / failed 行
        msgs = [l.strip() for l in lt.split("\n")
                if ("failed" in l.lower() or "maximum" in l.lower()
                    or "reached max" in l.lower() or "singular" in l.lower())]
        print("  失败/上限类消息（末尾 6 条）:")
        for s in msgs[-6:]:
            print("   ", s[:140])
        # 最后一个 Nonlinear status 块
        blocks = re.split(r"Nonlinear solution status: time=", lt)
        if len(blocks) > 1:
            last = blocks[-1]
            print("  最后一个求解块（t=" + last.split("\n")[0].strip() + "）:")
            for ln in last.strip().split("\n")[:12]:
                print("   ", ln[:130])
