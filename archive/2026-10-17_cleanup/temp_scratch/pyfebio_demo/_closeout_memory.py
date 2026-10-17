"""当日记忆追加 §十一：会话收口 + 交接指针（压缩上下文后靠它恢复）。"""
from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default"
         r"\memory\2026-10-03\ankle-fe-modeling.md")
t = M.read_text(encoding="utf-8")
b = len(t.encode("utf-8"))

add = """

---

## 十一、会话收口（用户："更新文档，我准备压缩上下文"）

### 🎯 压缩后从哪接手 —— 一句话

**读 `D:\\Project\\climbing_fall_analysis\\FE_PIPELINE_TODO.md` 的「🧭 接手入口」节**
（自包含：状态 / 一键复现 / 下一步排序 / 负面结果清单 / 环境开关 / 别再踩的坑 /
关键文件）。其余文档都是它的展开。

### 本次收口做了什么

| 文件 | 动作 |
|---|---|
| `FE_PIPELINE_TODO.md` | **新增「🧭 接手入口」节**（自包含，压缩后只读这一节就够）；修正过时状态行；入口标注"命令已实跑验证"；开关表去掉裸文件列改表头一句 + 全路径 |
| `ANKLE_FE_TODO.md` | 顶部加 **⚠️ 已取代** 提示 —— IITD 路线已停止，接手请去 THUMS 线（否则新会话可能又去走那条死路）|
| `docs/踝关节FE建模交接.md` | §13 顶部互指入口节 |
| `.learnings/LEARNINGS.md` | LRN-019（接触面必须合法 2-流形）/ LRN-020（优化前先确认是不是瓶颈）|
| `MEMORY.md` | 19.5 KB → **10.2 KB**（SecLog 12 条 gotcha 详情搬进新 `references/seclog-miniprogram-gotchas.md`；修掉 1 处原有死链；维护规则加"整理完必验三项"）|

### 交接前的实跑验证（不是"应该能跑"，是跑过了）

从零重跑整条链，4 步全 rc=0：
```
thums_materials.py -> cartilage.py -> build_feb.py(FEBio rc=0 / 17.3 s) -> a6_reaction.py
读回：支反力 436.817 N、软骨 p99 2.5799 MPa
```
与文档各处记录**逐字一致** ✅。入口节里 13 个路径引用**全部可解析** ✅

### 验证脚本（留在 `temp/pyfebio_demo/`，压缩后可直接复用）

- `verify_memory.py` —— MEMORY 三验：尺寸 / **每个本地路径都存在** / 信息未丢
- `verify_handover2.py` —— 交接九项：入口存在、命令已验证、状态一致、路径全解析…
- `final_handover.py` —— ASCII 输出的最终版（cmd 下中文/emoji 会吞行，验证脚本要避免）
- `diag_facets.py` —— **判定接触面是否合法 2-流形**（数每条边被几个面用）

### 当前状态一句话

**踝关节（胫距）FE 已跑通并通过吴恺量级验证**；剩余全是**加料**（韧带排序、
subtalar 接触面修复、force 驱动 600 N），不是**修地基**。下一步建议**先加 8 组踝韧带**
（不经过接触路径，能拿到吴恺的定性排序）。

### 本次会话我犯的错（都已记入 LRN）

1. `edit_file` **部分应用**留下损坏/重复代码块 ⇒ 改完必须 `ast.parse()` + **通读**
2. 诊断脚本把骨节点**原始 node ID** 当 coords 下标用 ⇒ 面识别全失败
3. 补丁脚本锚点**命中 2 次/0 次**（锚点文本在注释里也出现）⇒ 锚点要选只出现一次的
4. `assert "`build_feb.py`" not in t` **误判** —— 全路径结尾就含这个子串
5. `verify_handover.py` 里 `P` 误设成**文件**路径 ⇒ 13 个路径全报 MISS（假警报）
6. cmd 下 `more` + 中文/emoji 会**吞输出行** ⇒ 验证脚本输出改纯 ASCII"""

M.write_text(t + add, encoding="utf-8")
print(f"memory {b:,} -> {len(M.read_text(encoding='utf-8').encode('utf-8')):,} bytes")
print("章节:", [l.strip()[:30] for l in
                M.read_text(encoding="utf-8").split("\n")
                if l.startswith("## ")])
