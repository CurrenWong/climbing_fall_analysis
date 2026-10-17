"""打印 FEBio 报错前后的上下文。"""
import io

t = io.open(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_lig_last.log",
            encoding="utf-8", errors="replace").read()
L = t.splitlines()
i = [k for k, l in enumerate(L) if "FEBio rc=" in l]
k = i[-1] if i else 0
print("\n".join(L[max(0, k - 8):k + 44]))
