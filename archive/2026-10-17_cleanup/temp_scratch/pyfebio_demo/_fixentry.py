"""入口节：开关表的裸文件名 -> 表头一句 + 去掉文件列。"""
import re
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\FE_PIPELINE_TODO.md")
t = P.read_text(encoding="utf-8")
b = len(t.encode("utf-8"))
done = []


def sub(old, new, tag):
    global t
    assert t.count(old) == 1, f"[{tag}] 命中 {t.count(old)}"
    t = t.replace(old, new)
    done.append(tag)


sub("## 环境开关（全部有默认值，直接跑即可）\n\n"
    "| 开关 | 默认 | 文件 | 作用 |\n|---|---|---|---|\n",
    "## 环境开关（全部有默认值，直接跑即可）\n\n"
    "开关分布在两个文件：`scripts/ankle_fe/build_feb.py`、"
    "`scripts/ankle_fe/cartilage.py`\n\n"
    "| 开关 | 默认 | 作用 |\n|---|---|---|\n", "sw-head")

for cell in ("| `build_feb.py` |", "| `cartilage.py` |"):
    n = t.count(cell)
    t = t.replace(cell, "|")
    done.append(f"去裸名 x{n}")

# 注意：不能查 "`build_feb.py`" —— 全路径 `scripts/ankle_fe/build_feb.py` 的结尾
# 本身就含这个子串，会误判。要查的是「整格」形式。
assert ("| `build_feb.py` |" not in t) and ("| `cartilage.py` |" not in t), "开关表仍有裸文件格"

sub("| `temp/pyfebio_demo/a6_reaction.py` | 读 `.xplt` → 支反力 + 各域 von Mises |",
    "| `temp/thums/thums_lowerlimb.npz` | THUMS 提取结果（骨网格 + 韧带壳，21,608 节点）|\n"
    "| `temp/thums/bone_materials.json` | C 的材料表（皮质/松质/骨髓）|\n"
    "| `temp/pyfebio_demo/a6_reaction.py` | 读 `.xplt` → 支反力 + 各域 von Mises |",
    "artifacts")

P.write_text(t, encoding="utf-8")
print(f"已改: {done}")
print(f"{b:,} -> {len(t.encode('utf-8')):,} bytes")

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
head = t[t.index("# 🧭 接手入口"):t.index("## C. 解出 THUMS")]
paths = sorted(set(re.findall(r"`([\w\-./\\]+\.(?:py|md|json|npz|zip|k|fe?b))`",
                              head)))
miss = [p for p in paths
        if not (P / p.replace("\\", "/")).exists()
        and not (WS / p).exists()]
print(f"\n入口节 {len(paths)} 个路径，{len(miss)} 个解析不了: {miss}")
for p in paths:
    print(f"  {'OK ' if p not in miss else 'MISS'} {p}")
