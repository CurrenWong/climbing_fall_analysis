"""压缩前最后一件事：零回归复验 + 收口日记。"""
import io
import os
import re
import subprocess

R = r"D:\Project\climbing_fall_analysis"
PY = R + r"\.venv\Scripts\python.exe"
ACC = R + r"\scripts\ankle_fe\acceptance.py"
CM = ["DET_MIN=0.0", "TRUNC_MM=60", "AXIS_SRC=full", "ALIGN_AXIS=1",
      "FIX_TALUS=transverse"]

print("=== 零回归复验（默认路径，不得变）===")
for tag, extra in (
        ("默认旧格式 disp 0.1472", ["DRIVER=disp", "DISP_MM=0.1472"]),
        ("TIE+CART_NEW m=1 disp 1.6", ["DRIVER=disp", "DISP_MM=1.6", "TIE=1",
                                       "CART_NEW=1",
                                       "CART_DIR=" + R + r"\temp\thums\cart_m1"]),
        ("压缩 -1.0（新能力）", ["DRIVER=dispx", "DISP_MM=-1.0", "LOAD_SIGN=+1"])):
    r = subprocess.run([PY, "-u", ACC] + CM + extra, capture_output=True,
                       text=True, cwd=R, timeout=1800)
    t = (r.stdout or "") + (r.stderr or "")
    nm = ("NORMAL" if "NORMAL TERMINATION" in t else "ERROR")
    fz = re.search(r"两顶面合力 = \[\s*[-\d.e+]+\s+[-\d.e+]+\s+([-\d.e+]+)", t)
    print(f"  {tag:<26} {nm:>6}  Fz={float(fz.group(1)) if fz else 0:>10.4f} N")

D = r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-05.md"
ADD = """
## 📦 压缩前收口（2026-10-05 末）

**本轮全部落盘并校验 PASS**（AST / 死代码 0 / MEMORY 指针 11 个无缺失 / LRN 真实重号已修 / MEMORY 10212≤10240）。

### 修掉的 4 个真 bug（都会静默出错）
| ID | 内容 |
|---|---|
| `ERR-008` | `<traction lc=…>` **被静默忽略** ⇒ 满载荷瞬加；载荷曲线必须挂 **`<scale>`** |
| `ERR-009` | **加载方向符号反了** ⇒ 一直把关节**拉开（拉伸）**；历史值 187.77/209.55/592.6 全是拉伸 |
| `ERR-010` | `_outer()` 只认旧软骨键 ⇒ **任何 `POSE_TILT` 崩溃**；+ "没删旧 `.feb`"导致 8 姿势同值的假象 |
| `ERR-011` | 几何姿势落位平移**符号反**（越推越开）⇒ 3 法失败后定为**有界一维极小化** |

### 关键能力（压缩后可直接用）
- **压缩侧 600 N**：`LOAD_SIGN=+1` + `DISP_MM=-1.0` ⇒ **−608.82 N**（全程 NORMAL，余量到 −1017 N）
- **姿势**：倾斜式 ±20°（剪切 0→79.7 N）｜**几何式 `POSE_IE` ±5°（真内外翻）**｜`POSE_DF` ±2°
- **⚠️ `POSE_IE` 标签反向**：+IE 实际是**外翻**
- **逐韧带**：内翻拉外侧(CFL 3.513→2.842↓)、外翻拉内侧(MD2 3.708→3.824↑) ⇒ 物理正确
- **判据**：①✅ ②⚠️需重定义 ③✅方向 ④✅ ⑤✅（罚因子 **逐位鲁棒**）

### 压缩后从哪接手
1. `docs/踝FE可用模型使用说明.md`（**唯一运行入口**，§1 命令 / §2 压缩标定 / §3 姿势 / §4 判据 / §5 限制 / §6 坑）
2. `ANKLE_FIX_TODO.md`（顶部**已加更正横幅**，防读到作废数字）
3. `MEMORY.md` 踝段 + `references/fe-modeling-gotchas.md`

### 下一步
① 重定义 C2（有/无姿势两配置对比）② **修韧带插入点几何（ATFL 建模 35~37mm vs 解剖 15~20mm）**
③ 给 C3 定量阈值 ④ 主线 P3（2D 运动学含转动 DOF）仍在等踝线收尾
"""
s = io.open(D, encoding="utf-8").read()
io.open(D, "a", encoding="utf-8").write(ADD)
print(f"\ndaily {len(s)} -> {len(s)+len(ADD)} 字符")
