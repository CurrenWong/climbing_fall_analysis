from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
L = 10240
s = M.read_text(encoding="utf-8")
CUTS = [
    ("\n- **详见**：`memory/2026-10-03/ankle-fe-modeling.md`", ""),
    ("**Python（数据分析）**", "Python（数据分析）"),
    ("- **`uv pip install gmsh` ≠ 官方 SDK**：pip wheel 的 occ 绑定残缺 ⇒ STL→实体走不通；"
     "官方 SDK 在 `D:\\Program\\gmsh-sdk\\`，DLL 在 `lib/` 不在 `bin/`。"
     "**同版本号可能是两份不同二进制 → 判据是实测 API 可调用性**",
     "- **`uv pip install gmsh` ≠ 官方 SDK**：pip wheel 的 occ 绑定残缺 ⇒ STL→实体走不通；"
     "官方 SDK 在 `D:\\Program\\gmsh-sdk\\`，DLL 在 `lib/` 非 `bin/`。"
     "**同版本号可能是两份二进制 → 判据是实测 API 可调用性**"),
    ("- **edit payload 过大/含特殊符号会失败，且可能部分应用留下损坏文件**\n"
     "  ⇒ 改完一律 `ast.parse()` 验语法，别只看工具返回值",
     "- **edit payload 过大/含特殊符号会失败，且可能部分应用留损坏文件**\n"
     "  ⇒ 改完一律 `ast.parse()` 验语法"),
    ("| 4 | **微信 title ≤64 字符**（否则 45003）",
     "| 4 | **微信 title ≤64 字**（否则 45003）"),
]
for a, b in CUTS:
    if a in s and a != b:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
sz = len(s.encode())
print(f"disk = {sz:,} / {L}  余量 {L - sz:+,}")
