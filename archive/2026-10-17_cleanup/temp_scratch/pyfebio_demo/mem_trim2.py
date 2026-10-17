"""继续精简主 MEMORY（等长替换 + 长块改指针），细节已在 CONV_TODO.md / references。"""
import io

M = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
s = io.open(M, encoding="utf-8").read()
n0 = len(s.encode("utf-8"))

PAIRS = [
    # 600 N 行收紧
    ("""- ✅ **600 N 已达成（压缩侧，全程 NORMAL）**（`[LRN-20261005-072]`）：`LOAD_SIGN=+1`+`DISP_MM=-1.0`
  ⇒ **−608.82 N**；标定 −0.25→−132.56｜−0.5→−292.22｜−0.8→−479.60｜−1.2→−741.55｜−1.6→−1017.36 N
  （k=530→636 N/mm）。**⚠️ `[LRN-064]`「1.60mm→592.6 N」是拉伸值、已作废**""",
     """- ✅ **600 N 已达（压缩侧，全 NORMAL）**（`[LRN-20261005-072]`）：`LOAD_SIGN=+1`+`DISP_MM=-1.0`
  ⇒ **−608.82 N**；标定 −0.25/−0.5/−0.8/−1.2/−1.6 mm → −132.56/−292.22/−479.60/−741.55/−1017.36 N。
  **⚠️ `[LRN-064]` 的 592.6 N 是拉伸值、已作废**"""),
    # 姿势行收紧
    ("""- 📘 **唯一运行入口**（`[LRN-063]`）= `docs/踝FE可用模型使用说明.md`；**姿势=`POSE_TILT`/`POSE_TILT_ML`**（倾斜式；
  压缩侧 8 组全 NORMAL，剪切 **0→79.7 N**、Fz 570~646 N；几何旋转式 `POSE_DF` 因关节面相切而受阻。
  **⚠️ 该功能曾因 `_outer()` 只认旧软骨键而全崩**，`[ERR-20261005-010]` 已修）""",
     """- 📘 **唯一运行入口**=`docs/踝FE可用模型使用说明.md`（`[LRN-063]`）；**姿势=`POSE_TILT`/`POSE_TILT_ML`**
  （倾斜式；压缩侧 8 组全 NORMAL，剪切 0→79.7 N、Fz 570~646 N）。`[ERR-20261005-010]` 已修 `_outer()` 旧键崩溃"""),
    # 刚度 2 倍：长行 → 指针
    ("""- ⚠️ **刚度 2 倍未归因**：tet4 259.3 vs penta6 116.9 N（恒 2.2~2.35，五假设排除，`CONV_TODO.md`）""",
     """- ⚠️ **刚度 ~2× 未归因**（tet4 259.3 vs penta6 116.9 N，五假设排除）→ `CONV_TODO.md`"""),
    # F 未收敛块：长块 → 指针
    ("""- ⚠️ **F 未收敛但非伪影**（`[LRN-054/056]`）：k=F/δ 趋近 **49/151/451/880** N/mm
  (罚 1e-3→1) vs Anderson 独立估计 **1059**；但 **p≥2 全负 Jacobian ⇒ 无法验收敛**，
  真值可能高 1.5~2× ⇒ `436.8N` **偏""",
     """- ⚠️ **F 未收敛但非伪影**（`[LRN-054/056]`）：k 随罚因子 49/151/451/880 N/mm vs Anderson 估计 1059；
  **p≥2 全负 Jacobian ⇒ 无法验收敛**；`436.8N` **偏低"""),
]
for old, new in PAIRS:
    c = s.count(old)
    if c != 1:
        print(f"⚠️ 未命中({c})：{old.splitlines()[0][:60]}")
        continue
    s = s.replace(old, new)
    print(f"✓ {old.splitlines()[0][:52]}")
io.open(M, "w", encoding="utf-8").write(s)
n1 = len(s.encode("utf-8"))
print(f"\nMEMORY.md {n0} -> {n1} / 10240  余量 {10240 - n1}")
