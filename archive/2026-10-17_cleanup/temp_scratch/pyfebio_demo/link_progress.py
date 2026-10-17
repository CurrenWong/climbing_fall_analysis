"""① 确认本页不被 gitignore ② 从 README / 入口文档挂上链接 ③ 校验。"""
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
NEW = ROOT / "docs/总体进度.md"

# ① git 是否忽略
r = subprocess.run(["git", "check-ignore", "-v", str(NEW)], cwd=ROOT,
                   capture_output=True, text=True)
print(f"git check-ignore  rc={r.returncode}  "
      f"{'★ 被忽略（不可追踪）' if r.returncode == 0 else 'OK：未被忽略，可追踪'}")
if r.stdout.strip():
    print("   ", r.stdout.strip())

# ② README 挂链接
P = ROOT / "README.md"
s = P.read_text(encoding="utf-8")
BLOCK = ("\n> 🧭 **总体进度总览** 见 [`docs/总体进度.md`](docs/总体进度.md) ——\n"
         "> 一页看清**两条线 + 一条子线**的进度、阻塞与下一步（**追踪入口**）。\n")
if "总体进度.md" not in s:
    anchor = "> 📄 **阶段性总结与经验记录** 见"
    i = s.find(anchor)
    s = s[:i] + BLOCK.lstrip("\n") + "\n" + s[i:]
    P.write_text(s, encoding="utf-8")
print(f"README.md -> {P.stat().st_size:,} B  "
      f"[{'OK' if '总体进度.md' in P.read_text(encoding='utf-8') else 'XX'}]")

# ③ FE_PIPELINE_TODO 入口挂链接
P = ROOT / "FE_PIPELINE_TODO.md"
s = P.read_text(encoding="utf-8")
if "总体进度.md" not in s:
    s = s.replace("# 🧭 接手入口（先读这一节就够，其余是细节）",
                  "# 🧭 接手入口（先读这一节就够，其余是细节）\n\n"
                  "> 🗺️ **全局定位**：本文件是【子线·踝关节 FE】的入口；"
                  "两条线的总体进度看 `docs/总体进度.md`。", 1)
    P.write_text(s, encoding="utf-8")
print(f"FE_PIPELINE_TODO.md -> {P.stat().st_size:,} B  "
      f"[{'OK' if '总体进度.md' in P.read_text(encoding='utf-8') else 'XX'}]")

# ④ MEMORY.md 指针（用磁盘口径看余量）
M = WS / "MEMORY.md"
s = M.read_text(encoding="utf-8")
before = M.stat().st_size
if "总体进度" not in s:
    OLD = "- **入口**：`FE_PIPELINE_TODO.md`「🧭接手」｜`ANKLE_FIX_TODO.md`｜`docs/Anderson双论文_关键信息.md`"
    NEWL = ("- **🗺️总体进度**：`docs/总体进度.md`（两条线+子线一页总览）\n"
            "- **入口**：`FE_PIPELINE_TODO.md`「🧭接手」｜`ANKLE_FIX_TODO.md`｜"
            "`docs/Anderson双论文_关键信息.md`")
    if OLD in s:
        s = s.replace(OLD, NEWL, 1)
        M.write_text(s, encoding="utf-8")
print(f"MEMORY.md disk = {M.stat().st_size:,} / 10240  "
      f"(before {before:,})  [{'OK' if '总体进度' in M.read_text(encoding='utf-8') else 'XX'}]")
