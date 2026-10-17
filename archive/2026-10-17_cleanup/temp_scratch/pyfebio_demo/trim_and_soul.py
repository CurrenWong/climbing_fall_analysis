"""腾 MEMORY 空间 + 把流程铁律写进 SOUL.md。"""
import io

W = r"C:\Users\curre\.copaw\workspaces\default"
MEM = W + r"\MEMORY.md"
SOUL = W + r"\SOUL.md"

# ── 1) MEMORY 找最长行并压 ──
m = io.open(MEM, encoding="utf-8").read()
n0 = len(m.encode("utf-8"))
print("最长 6 行：")
for nb, i, s in sorted(((len(l.encode("utf-8")), i + 1, l)
                        for i, l in enumerate(m.splitlines())), reverse=True)[:6]:
    print(f"  {nb:>4} B 行{i:>4}: {s[:110]}")

REPS = [
    ("- **venv 无 pip**：`D:\\Program\\uv\\uv.exe pip install --python <venv>\\Scripts\\python.exe <pkg>`",
     "- **venv 无 pip**：`D:\\Program\\uv\\uv.exe pip install --python <venv>\\Scripts\\python.exe <pkg>`"),
    ("- **edit payload 过大/含特殊符号会失败且可能部分应用** ⇒ 改完一律 `ast.parse()`",
     "- **edit payload 过大/含特殊符号会失败且可能**部分应用 ⇒ 改完一律 `ast.parse()`"),
    ("- **cmd 别用 `python -c \"多行\"`**（IndentationError）⇒ 写临时 `.py`；**cmd 无 `tail`** ⇒ `more`",
     "- **cmd 别用 `python -c` 多行**（3 次踩坑）⇒ 写临时 `.py`；**cmd 无 `tail`** ⇒ `more`"),
    ("  压缩侧现全 `NORMAL`。**几何没问题**（互穿仅 6 节点；接触 31→344mm²≈生理 300mm²）",
     "  压缩侧现全 `NORMAL`。**几何没问题**（互穿 6 节点；接触 31→344mm²≈生理 300mm²）"),
    ("- **`uv pip install gmsh` ≠ 官方 SDK**（occ 残缺 ⇒ STL→实体不通）；官方在 "
     "`D:\\Program\\gmsh-sdk\\`（DLL 在 `lib/`）。**同版本号可能是两份二进制，判据=实测 API**",
     "- **`uv pip install gmsh` ≠ 官方 SDK**（occ 残缺 ⇒ STL→实体不通）；官方在 "
     "`D:\\Program\\gmsh-sdk\\`。**同版本号可能是两份二进制，判据=实测 API**"),
]
for old, new in REPS:
    if old in m and old != new:
        m = m.replace(old, new, 1)
io.open(MEM, "w", encoding="utf-8").write(m)
n1 = len(m.encode("utf-8"))
print(f"MEMORY.md {n0} -> {n1} B  {'✅' if n1 <= 10240 else f'❌ 超 {n1-10240}'}")

# ── 2) SOUL.md 加一条方法准则 ──
s = io.open(SOUL, encoding="utf-8").read()
ADD = """
**方向和符号，必须量出来，不能"读"出来。** 涉及方向/正负/标签的结论，先说清"哪条证据是
**不可被约定翻转**的"：产物文件里到底写了什么（结构）、材料点往哪边动（几何）、
有没有一个只有单一解释的算例（物理）、或者一个外部锚定（解剖/独立参照）。
**只凭"反力正负"或"哪个数更像"就下结论，是本项目反复返工的根因**（加载方向连翻两次）。
**另外：检查器本身也要有已知正对照** —— 没跑过对照的判据不算证据；
**更正要替换，不要叠加**（只前置新说明会留下两条相反论断）。
"""
if "方向和符号，必须量出来" in s:
    print("· SOUL.md 已有该条")
else:
    s = s.rstrip() + "\n" + ADD
    io.open(SOUL, "w", encoding="utf-8").write(s)
    print(f"✓ SOUL.md 已加入方法准则 -> {len(s.encode('utf-8'))} B")
