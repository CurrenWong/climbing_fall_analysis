"""给 MEMORY.md 留出合理余量（目标 ≥+40）。"""
from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
LIMIT = 10240
s = M.read_text(encoding="utf-8")

CUTS = [
    ("｜详见 `memory/2026-10-03/ankle-fe-modeling.md`", ""),
    ("最尖锐候选=**接触面定向**", "最尖锐=**接触面定向**"),
    ("**载荷-位移关系不可信**", "**载荷-位移不可信**"),
    ("**`taskkill /IM` 会杀 agent 自己** —— 用 PID", "**`taskkill /IM` 会杀自己** —— 用 PID"),
    ("**禁止杀所有 Python**（自杀）", "**禁止杀所有 Python**"),
    ("**删除非临时文件前必问**", "**删非临时文件前必问**"),
    ("**copaw MCP 是死代码** —— 直接用 skill 或 npx server，别等「内置工具」",
     "**copaw MCP 是死代码** —— 直接用 skill / npx，别等内置工具"),
    ("**`--help` 不全信** —— 以 `tools/list` 为准（Flint 多 1 个 tool）",
     "**`--help` 不全信** —— 以 `tools/list` 为准"),
]
for a, b in CUTS:
    if M.stat().st_size < LIMIT - 40:      # ★ 余量已够 ⇒ 停（原来写成 >= 导致一次都不裁）
        break
    if a in s:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
        print(f"  -> {M.stat().st_size:,} (余量 {LIMIT - M.stat().st_size:+,})")
sz = M.stat().st_size
print(f"\n最终 disk = {sz:,} / {LIMIT}   余量 {LIMIT - sz:+,}")

# 复核踝关节段完整
print("\n--- 踝关节段 ---")
for ln in s.split("\n"):
    if any(k in ln for k in ("FE_PIPELINE_TODO", "刚度 2 倍", "根因=接触区",
                             "Anderson 真值", "韧带线已降级")):
        print("  " + ln)
