"""修使用说明最后的致命残留：主命令仍指向拉伸。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
D = os.path.join(R, "docs", "踝FE可用模型使用说明.md")
d = io.open(D, encoding="utf-8").read()
n0 = len(d.encode("utf-8"))
REPS = [
    ("DRIVER=dispx DISP_MM=-1.0 LOAD_SIGN=+1 FIX_TALUS=transverse TRUNC_MM=60 ^",
     "DRIVER=dispx DISP_MM=+1.0 FIX_TALUS=transverse TRUNC_MM=60 ^"),
    ("- **先跑这条** ⇒ 应看到 `NORMAL TERMINATION`、`[load] LOAD_SIGN=+1 ⇒ …压缩…`、",
     "- **先跑这条** ⇒ 应看到 `NORMAL TERMINATION`、顶面**向下朝距骨**移动（压实）、"),
    ("> 🚨 **`DISP_MM` 必须为负、`LOAD_SIGN` 必须为 `+1`** —— 见 §2 与 §5 的符号说明。",
     "> 🚨 **压实 = `DISP_MM` 取正**（`DISP_MM>0`）。`LOAD_SIGN` 对 `dispx` **无效、不要设**。\n"
     "> **依据（不可被约定翻转）**：`.feb` 写入 `z = −DISP_MM`｜顶面 dz 朝距骨｜"
     "距骨固定+无韧带时 `+1.0→149.59 N` vs `−1.0→0.0000 N`。见 §2 与 §5。"),
    ("| `LOAD_SIGN=-1`（默认）/ `DISP_MM>0` | **把关节拉开** ❌ | `+0.25 → +187.77 N`（**拉伸**，韧带在撑）|",
     "| **`DISP_MM < 0`** | **把关节拉开** ❌ | `−0.25 → −132.56 N`；`−1.0 → −608.82 N`（**拉伸**，韧带在撑）|"),
]
n = 0
for old, new in REPS:
    if old in d:
        d = d.replace(old, new)
        n += 1
    else:
        print(f"  ⚠ 未命中：{old[:60]}")
io.open(D, "w", encoding="utf-8").write(d)
print(f"✓ 使用说明改写 {n} 处  {n0} -> {len(d.encode('utf-8'))} B")

# 复查
d2 = io.open(D, encoding="utf-8").read()
bad = [(i + 1, l.strip()) for i, l in enumerate(d2.splitlines())
       if ("必开" in l and "LOAD_SIGN" in l) or "DISP_MM=-1.0" in l
       or ("LOAD_SIGN=-1" in l and "拉开" in l)]
print(f"\n残留：{len(bad)} 处")
for i, l in bad:
    print(f"  {i:>4}: {l[:140]}")
print("\n主命令：")
for i, l in enumerate(d2.splitlines()):
    if "DRIVER=dispx" in l:
        print(f"  {i+1}: {l.strip()[:130]}")
