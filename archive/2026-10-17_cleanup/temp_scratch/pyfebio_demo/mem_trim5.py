from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
L = 10240
s = M.read_text(encoding="utf-8")
CUTS = [
    ("**`uv pip install gmsh` ≠ 官方 SDK**：pip wheel 的 occ 绑定残缺 ⇒ STL→实体走不通；\n"
     "  官方 SDK 在 `D:\\Program\\gmsh-sdk\\`，DLL 在 `lib/` 不在 `bin/`。**通用教训：同版本号可能是两份不同二进制 → 判据是实测 API 可调用性**",
     "**`uv pip install gmsh` ≠ 官方 SDK**：pip wheel 的 occ 绑定残缺 ⇒ STL→实体走不通；"
     "官方 SDK 在 `D:\\Program\\gmsh-sdk\\`，DLL 在 `lib/` 不在 `bin/`。"
     "**同版本号可能是两份不同二进制 → 判据是实测 API 可调用性**"),
    ("- **cmd 下别用 `python -c \"多行\"`**（`IndentationError`）⇒ 写临时 `.py` 再跑；\n"
     "  **cmd 没有 `tail`** ⇒ 用 `more` 或重定向后 `read_file`",
     "- **cmd 下别用 `python -c \"多行\"`**（`IndentationError`）⇒ 写临时 `.py`；"
     "**cmd 没有 `tail`** ⇒ 用 `more` 或重定向后 `read_file`"),
]
for a, b in CUTS:
    if a in s:
        s = s.replace(a, b, 1)
    M.write_text(s, encoding="utf-8")
sz = len(s.encode())
print(f"disk = {sz:,} 余量 {L - sz:+,}")
