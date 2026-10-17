"""更新 MEMORY.md 的踝关节段（加入韧带轮状态，等量替换、不超 10KB）。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
s = M.read_text(encoding="utf-8")
lines = s.split("\n")

i0 = next(i for i, l in enumerate(lines) if l.startswith("- **交接文档**"))
i1 = next(i for i in range(i0, len(lines)) if lines[i].startswith("- **详见**"))

NEW = """- **接手入口**：`FE_PIPELINE_TODO.md`「🧭 接手入口」节（自包含）；历史 `docs/踝关节FE建模交接.md`
- **管线状态**：**C ✅ / B ✅ / A ✅**（`rc=0`/16.6 s，支反力 **436.8 N**，软骨 p99 **2.58/4.27 MPa**
  → 外推 600 N 胫侧 ≈3.4 MPa，**与吴恺 3.74 MPa 吻合**）；几何源 THUMS V7.1（mm-s-tonne-N-MPa
  ⇒ E 已是 MPa；`*_center_spon` 是**骨髓腔**非松质骨）
- **韧带轮**（细节→`ANKLE_LIG_TODO.md`）：吴恺 4 条已定位（胫跟=`81101202`、胫舟=`81101203`、
  跟腓=`81100801`、下胫腓=`81101101`；EA 510/510/510/**75.4** MPa、厚 2.0 mm）；
  **模型跑通 `rc=0`/5.4 s** —— 关键 **`TRUNC_MM=60` 截断胫腓骨**（32488→9109 单元，**160 倍**）
- ❌ **排序复现=负结果**：ΔL **0.001~0.08 mm** vs 吴恺 1.92~5.24 mm；**胫跟名次颠倒**
  （吴恺拉伸最大 / 本模型缩短最多）——纯压缩下胫骨朝固定足骨**靠近** ⇒ 胫-足韧带**全缩短**
- **未完成**：subtalar 接触面非合法 2-流形；**最高价值下一步 = 拿到吴恺「图 4」原始数值**
- **详见**：`memory/2026-10-03/ankle-fe-modeling.md` §六–§十、`memory/2026-10-04/ankle-fe-solve-and-handover.md` §四"""

lines[i0:i1 + 1] = NEW.split("\n")
out = "\n".join(lines)
b = len(out.encode("utf-8"))
print(f"新 MEMORY.md = {b:,} / 10240  (余量 {10240 - b:+,})")
if b > 10240:
    print("!! 超限，需要再精简")
else:
    M.write_text(out, encoding="utf-8")
    print("已写入")
