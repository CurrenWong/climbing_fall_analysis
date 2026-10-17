"""MEMORY.md：加入踝子线下一步指针 + 等量削减到 ≤10240。"""
from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
s = M.read_text(encoding="utf-8")

# ① 追加指针（替换掉已过时的 "堵点=接触罚因子上限"）
OLD = """- ✅ **几何没问题**（点在-tet 互穿仅 6 节点）；接触 31mm²@t0→**344mm²**@加载 ≈ Anderson
  生理 300mm²。堵点=**接触罚因子上限**。`[LRN-054/056]`"""
NEW = """- ✅ **几何没问题**（点在-tet 互穿仅 6 节点）；接触 31mm²@t0→**344mm²**@加载 ≈ Anderson
  生理 300mm²。`[LRN-054/056]`
- 🎯 **子线下一步（已定）**：目标=**接 OpenSim + 判不同姿势**。卡口=**力控制病态**
  （距骨全固定⇒刚体模态；`FIX_TALUS=all`⇒双端夹持⇒踝不能内外翻）。解=**Anderson §2.2
  弱弹簧 500/50 N/mm + 解除内外翻内旋**；验收=**姿势排序对罚因子鲁棒**。
  ⚠️ 韧带线结案是**循环论证**（ΔL≈0 是双端夹持强制的）⇒ 待重开。`[LRN-057]`"""

# ② 等量削减
TRIMS = [
    ("""- **venv 里没有 pip**：`D:\\Program\\uv\\uv.exe pip install --python <venv>\\Scripts\\python.exe <pkg>`
- **`uv pip install gmsh` ≠ 官方 SDK**：pip wheel 的 occ 绑定残缺 ⇒ STL→实体走不通；官方 SDK 在 `D:\\Program\\gmsh-sdk\\`，DLL 在 `lib/` 非 `bin/`。**同版本号可能是两份二进制；判据=实测 API 可调用性**
- **cmd 下别用 `python -c "多行"`**（`IndentationError`）⇒ 写临时 `.py`；**cmd 没有 `tail`** ⇒ 用 `more` 或重定向后 `read_file`
- **edit payload 过大/含特殊符号会失败，且可能部分应用留损坏文件**
  ⇒ 改完一律 `ast.parse()` 验语法""",
     """- **venv 无 pip**：`D:\\Program\\uv\\uv.exe pip install --python <venv>\\Scripts\\python.exe <pkg>`
- **`uv pip install gmsh` ≠ 官方 SDK**（occ 绑定残缺 ⇒ STL→实体不通）；官方在 `D:\\Program\\gmsh-sdk\\`，DLL 在 `lib/`。**同版本号可能是两份二进制；判据=实测 API 可调用性**
- **cmd 别用 `python -c "多行"`**（IndentationError）⇒ 写临时 `.py`；**cmd 无 `tail`** ⇒ `more`
- **edit payload 过大/含特殊符号会失败且可能部分应用** ⇒ 改完一律 `ast.parse()`"""),
    ("""- **中文优先**（专有名词保留英文）；**先行动后报告**；
  **宁缺毋滥**（搜索硬 3-5 条）；**保留原文**（搬运/转载不「创造性总结」，早报除外）
- **emoji 用在合适处**：Discord/Slack 表情回应自然；正文中不用
- **重要数字全列**：人名/航班号/ID 这种 → 显式标注，不要「之前说过」
- **负面结果也要如实报**：试过但假设错了的，标明"假设错误 + 实测数据\"""",
     """- **中文优先**（专有名词留英文）；**先行动后报告**；**宁缺毋滥**（搜索硬 3-5 条）；
  **保留原文**（搬运不「创造性总结」，早报除外）
- **emoji 用在合适处**（Discord/Slack 表情回应）；正文中不用
- **重要数字全列**（人名/航班号/ID）；**负面结果也如实报**（标"假设错误 + 实测数据"）"""),
    ("""1. **MEMORY.md ≤ 10 KB 硬上限**
2. 新内容默认进 `references/`，主 MEMORY 只留指针；项目专属细节归对应 references
3. **搬，不是删** —— 迁走前先确认新文件里有全部信息
4. 整理完必验三项：**尺寸达标 + 无死链（每个路径都存在）+ 信息未丢**
5. 更新时同步改 `last_updated`；重大重构前先备份到 `memory_archive/`""",
     """1. **MEMORY.md ≤ 10 KB 硬上限**
2. 新内容默认进 `references/`，主 MEMORY 只留指针；项目细节归对应 references
3. **搬，不是删**（迁走前确认新文件有全部信息）；整理完必验：**尺寸 + 无死链 + 信息未丢**
4. 更新同步改 `last_updated`；重大重构前备份到 `memory_archive/`"""),
]
if OLD in s:
    s = s.replace(OLD, NEW, 1)
else:
    print("⚠ 未匹配指针锚点")
for a, b in TRIMS:
    if a in s:
        s = s.replace(a, b, 1)
    else:
        print(f"⚠ 未匹配削减: {a[:40]}...")
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
t = M.read_text(encoding="utf-8")
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'✓ OK' if sz<=10240 else '✗ 超'}")
print(f"  LRN-057 在档: {'Y' if 'LRN-057' in t else 'n'}")
