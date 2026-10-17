"""按行定位，稳健修 §5 限制表 + 下一步；并检查 §4 是否完整。"""
import io

D = r"D:\Project\climbing_fall_analysis\docs\踝FE可用模型使用说明.md"
ls = io.open(D, encoding="utf-8").read().splitlines()
print("--- 含 59,590 的行 ---")
for i, l in enumerate(ls):
    if "59,590" in l:
        print(f"  {i+1}: {repr(l[:130])}")
print("--- 含 '下一步（按优先级）' 的行 ---")
for i, l in enumerate(ls):
    if "下一步" in l:
        print(f"  {i+1}: {repr(l[:130])}")
print("--- 含 '口径' 的行 ---")
for i, l in enumerate(ls):
    if "口径" in l:
        print(f"  {i+1}: {repr(l[:130])}")
