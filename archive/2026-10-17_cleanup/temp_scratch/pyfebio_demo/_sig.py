import glob
import io
import os
import re

for p in sorted(glob.glob(r"temp/pyfebio_demo/a1a2/*.febio.log")):
    t = io.open(p, encoding="utf-8", errors="replace").read()
    tag = os.path.basename(p).split(".febio")[0]
    m = re.findall(r"failed to converge at time : ([\d.]+)", t)
    m2 = re.findall(r"(\d+) negative jacobians?", t)
    r0 = re.search(r"residual\s+([\d.e+-]+)\s+([\d.e+-]+)", t)
    ne = re.search(r"Nr of equations \.+ : (\d+)", t)
    print(f"{tag:<12} 负J事件={len(m2):<3} 首次失收敛 t={m[:1]} "
          f"首残差 INITIAL={r0.group(1) if r0 else '?'} "
          f"当前={r0.group(2) if r0 else '?'} 方程数={ne.group(1) if ne else '?'}")
