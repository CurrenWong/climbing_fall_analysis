import re
import sys
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio\febio.log")
txt = p.read_text(encoding="utf-8", errors="replace")
normal = "N O R M A L" in txt.replace("  ", " ")
steps = re.findall(r"Number of time steps completed \.+: (\d+)", txt)
ft = re.findall(r"failed to converge at time :\s*([0-9.eE+-]+)", txt)
neg = re.findall(r"(\d+) negative jacobians detected", txt)
ts = re.findall(r"^Step = (\d+)", txt, re.M)
dm = re.findall(r"D multiplier :\s*([\d.eE+-]+)\s*([\d.eE+-]+)", txt)
gap = re.findall(r"maximum gap\s*:\s*([\d.eE+-]+)", txt)
print("log size            :", len(txt))
print("NORMAL termination  :", normal)
print("steps completed     :", steps[-1] if steps else None)
print("last Step= line     :", ts[-1] if ts else None)
print("fail_t (all, last4) :", ft[-4:] if ft else None)
print("neg jacobian blocks :", len(neg))
print("last D multiplier   :", dm[-1] if dm else None)
print("last maximum gap    :", gap[-1] if gap else None)
print("diverging count     :", txt.count("Problem is diverging"))