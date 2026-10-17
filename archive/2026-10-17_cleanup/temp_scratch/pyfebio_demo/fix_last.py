"""最终修：腾出空间并补回 CONV 指针。"""
from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
LIMIT = 10240
s = M.read_text(encoding="utf-8")

# ① 腾空间：把两处较长表述压短（原本就在，信息不丢）
FREES = [
    ("**图片引用**：与 md 同目录、**纯文件名**", "**图片引用**：与 md 同目录、纯文件名"),
    ("**微信 title ≤ 64 字符**（否则 45003）", "**微信 title ≤64 字符**（否则 45003）"),
    ("**cmd `set X=Y && Z` 末尾追加 space**（40013）", "**cmd `set X=Y && Z` 尾追 space**（40013）"),
    ("**`--help` 不全信** —— 以 `tools/list` 为准", "**`--help` 不全信** —— 查 `tools/list`"),
    ("**删非临时文件前必问**", "**删非临时文件必问**"),
    ("**FE 接触建模四铁律 + tet 翻转奇偶性 + FEBio `.feb` 语法**",
     "**FE 接触四铁律 + tet 翻转奇偶性 + `.feb` 语法**"),
]
for a, b in FREES:
    if a in s:
        s = s.replace(a, b, 1)

# ② 补回 CONV 指针（紧凑）
if "CONV_TODO" not in s:
    s = s.replace("；材料已排除（比值恒", "；五假设已排除，详见 `CONV_TODO.md`"
                  "（比值恒", 1)

M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"disk = {sz:,} / {LIMIT}   余量 {LIMIT - sz:+,}")
while sz > LIMIT:
    for a, b in (("（比值恒 2.2~2.35）⇒ 剩离散/接触", "⇒ 剩离散/接触"),
                 ("**与实验一致**；吴恺 1.92~5.24mm 隐含 8~21% **不符实验**",
                  "**与实验一致**；吴恺 1.92~5.24mm **不符实验**"),
                 ("**弱弹簧+解除内外翻内旋**", "**弱弹簧+松内外翻内旋**")):
        if a in s:
            s = s.replace(a, b, 1)
            M.write_text(s, encoding="utf-8")
            sz = M.stat().st_size
            print(f"  -> {sz:,}  余量 {LIMIT - sz:+,}")
            break
    else:
        break

print(f"\n最终 disk = {sz:,} / {LIMIT}   余量 {LIMIT - sz:+,}")
print("CONV 指针:", "OK" if "CONV_TODO" in s else "缺")
print("\n--- 踝关节段 ---")
for ln in s.split("\n"):
    if any(k in ln for k in ("FE_PIPELINE_TODO", "刚度 2 倍", "根因=接触区",
                             "Anderson 真值", "韧带线已降级", "详见")):
        print("  " + ln)
