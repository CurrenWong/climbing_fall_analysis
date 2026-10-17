"""落盘：子线协调机制 + A1 任务 ID。"""
import io
import os

W = r"C:\Users\curre\.copaw\workspaces\default"
DAY = os.path.join(W, "memory", "2026-10-06.md")

ADD = """
## 🤝 子线协调机制（用户 2026-10-06 指定）

**用户原话**：「保存这个计划，然后用 subagent 逐步推进，你作为主 agent 协调计划实现和目标」

**分工**：
- **我（主 agent）**：规划 / 写自包含任务书 / **验收**（跑校验、核数字、查回归）/ 落盘 / 决定继续或退回
- **subagent**：执行单个步骤，串行不并行（后步依赖前步数字）

**计划文件**：`ANKLE_ROADMAP_TODO.md`（项目根，3643B，**已被 `.gitignore:13: *_TODO.md` 覆盖**，
经 `git check-ignore -v` 确认）
内含：接距下关节的真正意义（一条链串三件事）／A→D 四阶段+验收判据／
A1 后立刻检验的假设／5 项并行修复债／终点交付物／诚实不确定性／协调协议／8 步状态表

**四阶段一句话**：
A1 接距下关节 → A2 解除 `FIX_TALUS=transverse`（**唯一不确定步骤**）→
A3 摩擦+tie 罚因子（验收=绝对刚度进 Anderson `1059 N/mm` 的 ±20%，现 49/151/451/880）→
B1 二分找根 + B2 双步序（接 OpenSim 的闸门）→ C1/C2 (载荷,角度) 失效边界 → D runner+风险表

**📌 当前任务 ID（压缩后据此续）**：
| 步骤 | task_id | session | 状态 |
|---|---|---|---|
| **A1** 接距下关节 | **`task-16deba15fcd7`** | `sub-76983873` | running（3600s） |

**A1 任务书要点**（已在 spawn 参数里，简报回来须核对）：
几何预检**只测不改**（间隙 >1mm 或穿透 >0.1mm ⇒ 先报告不硬接）→
接接触对进 `.feb` → 四项验收带数字（含接触对 / rc=0+NORMAL /
载荷路径从纯韧带剪切变有接触 / 与基线 592.6020 N 的 Δ% 及解释）→
**明确 PASS 判断**；已声明 592.6020 会变是**预期非回归**；禁碰 A2/A3。

**硬规则**（每个任务书都写）：方向结论须有不可被约定翻转的证据｜检测与处置分离｜
失败如实报不外推｜删旧 `.feb` 查 rc｜不用 `python -c` 多行｜细网格须 `LAUGON=AUGLAG`

---

## 🔧 git（另一 agent 负责，我只管自己的）

- ✅ **我已提交 `c58e7a4`**：`scripts/ankle_fe/` 16 文件 / 6,385 行（踝子线代码入库）
  - cmd **不支持 heredoc** ⇒ 写临时文件再 `git commit -F <file>`
  - `.gitignore:13` 的内容是 `*_TODO.md`（grep 输出的 `> ` 是显示前缀，不是文件内容）
- ⚠️ **发现 82 个 .py / 41,069 行未进 git**（`scripts/opensim_fe/` 49 个、
  `src/climbing/coupling/` 18 个、15 个 tests、30 个 docs 含 `总体进度.md` 与
  `踝FE可用模型使用说明.md`；另有 `悬崖攀登坠落伤害仿真真方案.md` 被删非我所为）
  —— **用户 2026-10-06：「其他 git 有别的 agent 跟进，你关注你负责的文件夹即可」⇒ 不再跟进**
"""
d = io.open(DAY, encoding="utf-8").read()
if "子线协调机制（用户 2026-10-06 指定）" in d:
    print("· 已有")
else:
    io.open(DAY, "a", encoding="utf-8").write(ADD)
    print(f"✓ 日记 -> {len(io.open(DAY, encoding='utf-8').read().encode('utf-8')):,} B")
