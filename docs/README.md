# docs/

抱石落地骨与韧带损伤建模项目的设计文档、阶段报告、论文稿与文献笔记。

> **本页自 2026-10-10 起是「权威索引」** —— 解决"文档太多、不知道信哪份"。
> **规则：一个主题只有一份权威。** 其余是过程件或归档件，**不要从里面取结论**。

---

## 🎯 权威索引（一个主题 → 只看这一份）

| 我想知道… | 只看这一份 |
|---|---|
| **现在整体到哪了 / 下一步做什么** | [`reports/总体进度.md`](reports/总体进度.md) —— 一页总览，**§四 = 下一步**，§五 = V/A 靶子总表 |
| **论文写什么 / 数字从哪来** | [`paper/论文v2_数学建模版.md`](paper/论文v2_数学建模版.md) ＋ [`paper/论文v2_出处索引.md`](paper/论文v2_出处索引.md)（逐条溯源） |
| **当前结论与可信边界** | [`reports/最终结论_2026-10-07.md`](reports/最终结论_2026-10-07.md) |
| **线 1 路线 P0–P5** | [`methods/方案更新_v2.md`](methods/方案更新_v2.md) |
| **主线 P2 踝（解析模型）结论** | [`reports/P2_踝结论_2026-10-07.md`](reports/P2_踝结论_2026-10-07.md) |
| **工况错配到底占多少（P5 前置侦察）** | [`notes/P5_CONDITION_ALIGN.md`](notes/P5_CONDITION_ALIGN.md) |
| **踝 FE 子线怎么用 / 到哪了** | [`notes/ANKLE_FIX_TODO.md`](notes/ANKLE_FIX_TODO.md) **顶部「🧭 接手入口」** |
| **踝 FE 接触根因 + 罚标定** | [`notes/CASEB_CONTACT_ROOTCAUSE.md`](notes/CASEB_CONTACT_ROOTCAUSE.md) §1–§14 |
| **踝子线的定位（供不供论文数）** | [`reports/踝子线_scope决策_2026-10-07.md`](reports/踝子线_scope决策_2026-10-07.md)（§八 = 在论文里的角色） |
| **线 2 OpenSim → FE 怎么做** | [`reports/OpenSim_FE交接.md`](reports/OpenSim_FE交接.md) |
| **S5 唯一可引用数字** | [`reports/S5_QUOTABLE.md`](reports/S5_QUOTABLE.md)（正文结论必须引它） |
| **S5 脚-垫接触契约 / 门禁** | [`s5-contact/S5_contact_contract.md`](s5-contact/S5_contact_contract.md) |
| **踝 FE 运行入口（命令）** | [`methods/踝FE可用模型使用说明.md`](methods/踝FE可用模型使用说明.md) |

---

## 🗂 目录

| 目录 | 作用 | 权威件 |
| --- | --- | --- |
| `paper/` | **论文成品**（v2 稿 ＋ 出处索引 ＋ 评审 ＋ 汇总报告） | `论文v2_*.md` |
| `reports/` | **阶段结论与进度**（这是取结论的地方） | `总体进度.md` |
| `methods/` | 仿真方案、方法学说明、建模交接 | `方案更新_v2.md` |
| `notes/` | 阶段性技术记录（根因、研究、结论） | `CASEB_CONTACT_ROOTCAUSE.md`、`P5_CONDITION_ALIGN.md` |
| `s5-contact/` | S5 脚-垫接触收敛系列 17 篇，**成套阅读** | `S5_contact_contract.md` |
| `archive/` | **已归档（只读历史，不再是权威）** | 见 [`archive/README.md`](archive/README.md) |
| （根目录） | `README.md`、`rigid_vtarget_2026-10-08.md` | — |

---

## ⚠️ 三层文档，别混（这是"搞不清楚"的主因）

| 层 | 特征 | 数量 | 能不能取结论 |
| --- | --- | --- | --- |
| **① 持久权威** | `docs/**` 且**不以 `_TODO.md` 结尾**，进了 git | **120** | ✅ **唯一可以取结论的层** |
| **② 会话草稿** | **`*_TODO.md`**（根目录 **13** ＋ `docs/notes/` **23**） | **36** | ❌ **被 `.gitignore` 排除，只在本机**，换机器/清库即消失 |
| **③ 归档** | `docs/archive/`（3 份文档 ＋ 1 份 README） | 3 | ❌ 只读历史 |

> 口径：项目自有 markdown **156** 份（不含 `.venv/`、根 `paper/` 第三方文献、`model/`）。
> 归档/待判/绝不动的处置判断 → [`reports/文档处置说明.md`](reports/文档处置说明.md)。

⇒ **看到 `*_TODO.md` 就知道那是"当时的待办/接手笔记"，不是结论。**
⇒ 真正的结论只在 ① 里；② 过期后应把结论**搬进 ①** 再归档或删除。

> 另：项目文档普遍用**短路径**引用（写 `docs/X.md`、实际在 `docs/reports/X.md`），
> 这是全库惯例，**不是死链**。

---

## 📖 推荐阅读顺序

1. [`reports/总体进度.md`](reports/总体进度.md) —— 先看一页总览与 §四 下一步
2. [`reports/最终结论_2026-10-07.md`](reports/最终结论_2026-10-07.md) —— 结论与可信边界
3. [`paper/论文v2_数学建模版.md`](paper/论文v2_数学建模版.md) —— 论文本体（数字溯源走 `论文v2_出处索引.md`）
4. [`methods/方案更新_v2.md`](methods/方案更新_v2.md) —— 路线与 P0–P5 定义
5. 需要细节再进对应子目录

## 📦 已归档

[`archive/`](archive/) 存放已被取代的稿件（当前：**论文初稿** → 由 v2 取代）。
**归档 ≠ 删除**，内容完整可查；⚠️ 但**章节号不可与 v2 互换**，详见 [`archive/README.md`](archive/README.md)。

## 第三方模型与文献

`model/`（AM50、THUMS、足踝 FE 模型等）与 `paper/` 下的 `MinerU_markdown_*.md`
（吴恺论文、Anderson 双论文、FEBio 手册等第三方资料）**不随本仓库分发**，因各自许可约束。
复现所需的获取方式、版本与放置约定见仓库根目录 `README.md`。
