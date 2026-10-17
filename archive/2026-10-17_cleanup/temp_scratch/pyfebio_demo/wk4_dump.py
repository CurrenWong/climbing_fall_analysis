import io
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\paper/_text/吴恺_足踝FE.txt")
t = io.open(p, encoding="utf-8", errors="ignore").read()
compact = [l.strip() for l in t.split("\n") if l.strip()]

out = []
out.append("=== 结果段（韧带位移）===")
out.append(" ".join(compact[158:180]))
out.append("")
out.append("=== 边界条件段 ===")
out.append(" ".join(compact[148:165]))
out.append("")
out.append("=== 材料表 ===")
out.append(" ".join(compact[110:132]))
out.append("")
out.append("=== 验证段 ===")
out.append(" ".join(compact[208:222]))
D = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_wk_out.txt")
D.write_text("\n".join(out), encoding="utf-8")
print("wrote", D.stat().st_size)
