"""看 E1p / E2p 的失败签名：是否仍是残差震荡？最后一块的判据是否满足？"""
import io
import re
from pathlib import Path

OUT = Path(r"D:\Project\climbing_fall_analysis\temp/pyfebio_demo/a1a2")
for tag in ("E1p", "E2p", "E1p_bfgs"):
    lg = OUT / f"{tag}.febio.log"
    if not lg.exists():
        print(f"缺 {tag}"); continue
    t = lg.read_text(encoding="utf-8", errors="replace")
    print("=" * 74)
    print(f"### {tag}")
    msgs = re.findall(r"failed to converge at time : ([\d.]+)", t)
    print(f"  失收敛次数={len(msgs)}  最早={msgs[0] if msgs else '—'} 最晚={msgs[-1] if msgs else '—'}")
    negj = re.findall(r"(\d+) negative jacobians?", t)
    print(f"  负J事件={len(negj)}")
    blocks = re.split(r"Nonlinear solution status: time=", t)
    if len(blocks) > 1:
        last = blocks[-1]
        tt = last.split("\n")[0].strip()
        print(f"  最后一块 t={tt}")
        m = re.search(r"residual\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)", last)
        e = re.search(r"energy\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)", last)
        d = re.search(r"displacement\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)", last)
        for nm, mm in (("residual", m), ("energy", e), ("disp", d)):
            if mm:
                ini, cur, req = (float(x) for x in mm.groups())
                ok = "✅满足" if (req == 0 and cur == 0) or (req > 0 and cur < req) else "❌未达"
                print(f"    {nm:<9} INIT={ini:.3e} CUR={cur:.3e} REQ={req:.3e} {ok}")
        # 最后 5 个 CURRENT residual（震荡检查）
        cs = [float(x) for x in re.findall(r"residual\s+[\d.e+-]+\s+([\d.e+-]+)", last)]
        print(f"    末段 residual 序列（{len(cs)} 个）: "
              + " ".join(f"{c:.2e}" for c in cs[-8:]))
    for l in t.split("\n"):
        if "augmentation" in l or "sliding interface" in l:
            pass
    augs = re.findall(r"augmentation # (\d+)", t)
    print(f"  augmentation 最大编号={max(int(a) for a in augs) if augs else '—'}")
