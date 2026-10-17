"""给 MEMORY.md 腾出 ~45 字节余量。"""
from pathlib import Path

p = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
LIMIT = 10240
s = p.read_text(encoding="utf-8")
FREES = [
    ("**保留原文**（搬运/转载不做「创造性总结」，每日早报除外）",
     "**保留原文**（搬运/转载不「创造性总结」，早报除外）"),
    ("**`taskkill /IM` 会杀自己** —— 用 PID", "**`taskkill /IM` 会杀自己**，用 PID"),
    ("**edit 工具 payload 过大/含特殊符号会失败",
     "**edit payload 过大/含特殊符号会失败"),
    ("**先行动后报告**（能跑就跑）", "**先行动后报告**"),
    ("**emoji 用在合适处**", "**emoji 用在合适处**"),
    ("references/miniprogram-ci.md` | 小程序编译 + 预览二维码（私钥、`?.` 兼容、agent-ui 补丁）|",
     "references/miniprogram-ci.md` | 小程序编译 + 预览二维码（私钥、`?.`、agent-ui 补丁）|"),
]
for a, b in FREES:
    if a != b and a in s:
        s = s.replace(a, b, 1)
p.write_text(s, encoding="utf-8")
sz = len(s.encode())
print(f"disk = {sz:,} / {LIMIT}   余量 {LIMIT - sz:+,}")
