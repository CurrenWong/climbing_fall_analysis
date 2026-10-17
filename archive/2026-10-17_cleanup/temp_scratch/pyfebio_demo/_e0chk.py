"""对照：E0（成功收敛）的最后一块判据 vs E2p（失败）—— 判定 residual REQ=0 的含义。"""
import io
import re
from pathlib import Path

OUT = Path(r"D:\Project\climbing_fall_analysis\temp/pyfebio_demo")
CANDS = [OUT / "a1a2/E0_base.febio.log", OUT / "E0_base.febio.log"]
for c in CANDS:
    if c.exists():
        f = c
        break
else:
    # 自己找
    hits = list(Path(r"D:\Project\climbing_fall_analysis\temp").rglob("E0_base.febio.log"))
    f = hits[0] if hits else None
print("E0 log =", f)
if f:
    t = f.read_text(encoding="utf-8", errors="replace")
    print("含 NORMAL TERMINATION:", "N O R M A L" in t)
    blocks = re.split(r"Nonlinear solution status: time=", t)
    print(f"总求解块数={len(blocks)-1}")
    for b in blocks[1:]:
        tt = b.split("\n")[0].strip()
        m = re.search(r"residual\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)", b)
        e = re.search(r"energy\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)", b)
        d = re.search(r"displacement\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)", b)
        if m and float(m.group(2)) < 1.0:
            print(f"  t={tt:<12} res CUR={m.group(2)} REQ={m.group(3)}"
                  + (f" | E CUR={e.group(2)} REQ={e.group(3)}" if e else "")
                  + (f" | D CUR={d.group(2)} REQ={d.group(3)}" if d else ""))
    # 收敛汇总
    for l in t.split("\n"):
        if "convergence" in l.lower() or "N O R M A L" in l or "iterations" in l:
            print("   ", l.strip()[:120])
