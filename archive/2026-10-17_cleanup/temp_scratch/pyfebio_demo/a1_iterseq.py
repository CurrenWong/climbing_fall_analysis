"""对比 talus(能解) vs calcaneus(失败) 前几步迭代的收敛序列。"""
import re
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp/pyfebio_demo/a1_matrix")


def seq(tag, n=6):
    s = (W / f"{tag}.febio.log").read_text(encoding="utf-8", errors="replace")
    i = s.find("beginning time step 1")
    seg = s[i:i + 8000]
    lines = [l.strip() for l in seg.split("\n")]
    out, cur = [], None
    for l in lines:
        if re.fullmatch(r"\d+", l):
            cur = {"it": l}
            out.append(cur)
        elif l.startswith("step from line search") and cur is not None:
            cur["ls"] = l.split("=")[1].strip()
        elif l.startswith("residual") and cur is not None:
            cur["r"] = l.split()[1:4]
        elif l.startswith("energy") and cur is not None:
            cur["E"] = l.split()[1:4]
        elif l.startswith("displacement") and cur is not None:
            cur["d"] = l.split()[1:4]
        elif l.startswith("ERROR") or "jacobian" in l.lower():
            out.append({"err": l})
        if len(out) >= n:
            break
    res = []
    for o in out:
        if "err" in o:
            res.append("   >> " + o["err"])
        else:
            res.append(f"   it{o['it']:>2} ls={o.get('ls','?'):<10} "
                       f"r_cur={o.get('r',['?'])[1]:<12} "
                       f"E_cur={o.get('E',['?'])[1]:<12} "
                       f"|du|={o.get('d',['?'])[1]}")
    return res


for tag in ("talus", "calcaneus", "both", "nocontact"):
    print(f"=== {tag}")
    for l in seq(tag):
        print(l)
