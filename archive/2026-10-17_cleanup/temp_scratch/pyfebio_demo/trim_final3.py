"""最后 47 字节：合并「详见」行到入口行 + 若干缩短。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
LIMIT = 10240
s = M.read_text(encoding="utf-8")

# 合并：把「详见」并入入口行，删掉独立行（省一整行开销）
a = ("- **入口**：接手=`FE_PIPELINE_TODO.md`「🧭 接手入口」｜计划=`ANKLE_FIX_TODO.md`｜"
     "真值=`docs/Anderson双论文_关键信息.md`")
b = ("- **入口**：“`FE_PIPELINE_TODO.md`「🧭 接手入口」｜`ANKLE_FIX_TODO.md`｜"
     "`docs/Anderson双论文_关键信息.md`｜详见 `memory/2026-10-03/ankle-fe-modeling.md`")
if a in s:
    s = s.replace(a, b, 1)
s = s.replace("\n- **详见**：`memory/2026-10-03/ankle-fe-modeling.md`", "", 1)
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"  合并后 {sz:,} (余量 {LIMIT - sz:+,})")

MORE = [
    ("**入口**：“", "**入口**："),          # 去掉误加的中文引号
    ("- **韧带线已降级**：Tochigi 实测站姿下韧带本就松弛（CFL 0.2~2.1%）⇒ 我的 ΔL≈0 "
     "**与实验一致**；吴恺 1.92~5.24mm 隐含 8~21% **不符实验**",
     "- **韧带线已降级**：Tochigi 实测站姿韧带本就松弛（CFL 0.2~2.1%）⇒ 我 ΔL≈0 "
     "**与实验一致**；吴恺 1.92~5.24mm 隐含 8~21% **不符实验**"),
    ("- **Anderson 真值**：all-hex｜距骨**弱弹簧 +解除内外翻内旋**｜骨刚性｜**不建韧带**｜"
     "软骨 E=12 ν=0.42、厚 1.5~1.7mm｜精度 **10~15%**",
     "- **Anderson 真值**：all-hex｜距骨**弱弹簧+解除内外翻内旋**｜骨刚性｜**不建韧带**｜"
     "软骨 E=12 ν=0.42、厚 1.5~1.7mm｜精度 **10~15%**"),
    ("- ⭐ **根因取向=接触区欠分辨**：仅 **73 节点**(~2.5mm)；加密到 692 后软骨体积\n"
     "  涨 **2.3×** ⇒ 粗网格低估面积",
     "- ⭐ **根因取向=接触区欠分辨**：仅 **73 节点**(~2.5mm)；加密到 692 后软骨体积涨\n"
     "  **2.3×** ⇒ 粗网格低估面积"),
]
for x, y in MORE:
    if M.stat().st_size <= LIMIT:
        break
    if x in s:
        s = s.replace(x, y, 1)
        M.write_text(s, encoding="utf-8")
        print(f"  -> {M.stat().st_size:,} (余量 {LIMIT - M.stat().st_size:+,})")
sz = M.stat().st_size
print(f"\n最终 disk = {sz:,} / {LIMIT}  {'OK' if sz <= LIMIT else '超限 ' + str(sz - LIMIT)}")
if sz <= LIMIT:
    print("\n--- 踝关节段最终内容 ---")
    for ln in s.split("\n"):
        if "FE_PIPELINE_TODO" in ln or "刚度 2 倍" in ln or "根因取向" in ln \
           or "Anderson 真值" in ln or "韧带线已降级" in ln or "材料已排除" in ln:
            print("  " + ln)
