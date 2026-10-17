"""补当日记忆：① 三处文档收口 ② 1→3 bug 修复 ③ 新建总览。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
DAILY = WS / "memory/2026-10-04.md"

SEC = """

---

## 十二、三处收口 + 自查 bug + 新建总览（2026-10-04 晚）

### 12.1 tied-elastic 结果同步进 6 处文档
`FE_PIPELINE_TODO.md`（入口「三件必须知道」表新增 ✅ 新通路行）·
`ANKLE_FIX_TODO.md`（P1d2a ✅ / P1d2b 🔴）· `CONV_TODO.md`（交叉引用节）·
`.learnings/LEARNINGS.md`（+LRN-046~050）· 本文件 §十一 · `MEMORY.md`。

### 12.2 ★ 自查出并修复一个「静默丢 25% 体积」的 bug（LRN-051）
用户追问「**加密什么意思**」时回看代码，发现 `cart_refine.py` 的**三角形红加密写成了 1→3**：

```
正确 1→4:  (v0,m01,m20) (m01,v1,m12) (m20,m12,v2) + ★(m01,m12,m20)
我写的 1→3: 只有 3 个角的 —— 中央那个丢了
```

后果：每个三角底面少 1/4 面积 ⇒ 对应楔形少 1/4 体积。
**且 3 个角楔形各自都是合法正体积 ⇒ 负体积检查完全抓不到、FEBio 也不报错，完全静默。**

| 关节 | 修前偏差 | **修后** |
|---|---|---|
| tibiotalar_tibia | −14.1% | **−0.77%** |
| tibiotalar_talus | −17.6% | **−0.77%** |
| subtalar_talus | −15.9% | **−1.50%** |
| subtalar_calcaneus | −17.0% | **−2.05%** |

单元数 93 → **372**（= 61×4 + 32×4 ✓；修前 340，**正好少 32 = 三角底面个数**）。

**更值得记的是我的错误归因**：当时把 −14% 判断成「非凸单元的 6-四面体 fan 分解
近似偏差」—— 一个听起来很专业的数值/物理解释，**于是我没去查代码**。

> **教训（LRN-051）**
> 1. 拓扑重划分（加密/粗化）必须把**体积守恒当硬闸门**，偏差 >1% **先查代码**，
>    不要先找物理解释。「非凸」「近似误差」这类解释最容易骗过自己。
> 2. **子单元数必须手算核对**（母单元数 × 每类子单元数）。
>    我看到 93→340 没去算 61×4+32×4=372，**漏检一次**；
>    事后 372−340=32 直接指认了 bug。
>
> **不受影响**：hex8 侧 **28%** 子单元翻负仍然成立（与三角漏面独立），P1d2b 阻塞依旧。

### 12.3 新建可追踪的总览文档 `docs/总体进度.md`
计划层面盘点后固化为**一页总览**（两条线 + 一条子线），8 节：
〇目标与两线关系 / 一线1 P0–P5 / 二线2 S0–S5 / 三子线踝FE / 四下一步 /
五验证靶子 V1–V9 + A1–A6 / 六待决策 3 件 / 七专项文档索引 / 八总结。
带状态标记约定（✅⏳🚧⛔⚪）、「最后更新」字段、「更新本页的时机」。

**★ 命名决策（重要）**：`.gitignore` 含 `*_TODO.md` ⇒
`FE_PIPELINE_TODO.md` / `ANKLE_FIX_TODO.md` / `CONV_TODO.md` **都不进 git**（会话本地笔记）。
所以新总览**刻意不叫 `_TODO`**，实测 `git check-ignore` 确认**可追踪**。
已挂链到 `README.md`（顶部第 3 条，标「追踪入口」）、
`FE_PIPELINE_TODO.md`（🧭 接手入口加「全局定位」）、`MEMORY.md`（踝段首行 🗺️ 总览指针）。

> **推论（值得记住）**：要长期留存的结论必须落到
> `阶段总结.md` / `docs/*.md` / `.learnings/` —— `*_TODO.md` 系列不进版本库。

### 12.4 计划层面的结论（用户要求「从总体层面汇报」）
- **线 1 主线**：v1 ✅ / P0 ✅ / P1 ✅，但 P1 硬结果 = **真实抱石下一例骨折都报不出来**
  （最高利用率 0.725，可信度 13.8%）⇒ **根因结构性：无转动自由度，算不出占 26–28% 的踝伤**
- **线 2 耦合链**：S0–S3.0/S3.C/S3.D 全 ✅，但**绝对骨折高度不可移植**
  （粗体素"7 m"=伪影 / 平滑 gmsh gauge 26.5 m / 方向 C ≈2–3 m / 最佳 BC 5.43 m vs [Y25] 7–9 m）
- **下一步首选 = 主线 P3 二维运动学**（唯一能同时动摇 V5+V8 的动作）；线 2 建议转入收尾
- **待用户决策 3 件**：P3 vs P5 先做哪个 / 线 2 绝对高度是否封档 / 踝子线投入程度
"""

if SESS.exists():
    t = SESS.read_text(encoding="utf-8")
    if "十二、三处收口" not in t:
        SESS.write_text(t.rstrip() + SEC, encoding="utf-8")
    print(f"session note -> {SESS.stat().st_size:,} B")
else:
    print("MISS session note")

# 索引
if DAILY.exists():
    t = DAILY.read_text(encoding="utf-8")
    add = ("\n- `2026-10-04/ankle-fe-solve-and-handover.md` → §十一 tied-elastic 通路已验证；"
           "§十二 三处收口 + LRN-051（三角 1→3 静默丢 25% 体积）+ 新建 `docs/总体进度.md`\n")
    if "十二 三处收口" not in t:
        DAILY.write_text(t.rstrip() + "\n" + add, encoding="utf-8")
    print(f"daily note -> {DAILY.stat().st_size:,} B")
else:
    print("MISS daily note")
