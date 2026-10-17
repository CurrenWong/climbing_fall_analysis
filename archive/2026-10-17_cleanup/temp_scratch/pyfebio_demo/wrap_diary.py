"""收口④：日记 append（压缩前收口记录）。"""
import io
import os

W = r"C:\Users\curre\.copaw\workspaces\default"
D = os.path.join(W, "memory", "2026-10-06.md")
ADD = """

## 📦 压缩前收口（2026-10-06 末）

**改了三处文档**（净字节≈0，MEMORY 从 10,222 → **10,175 B**，余量 65）：
1. **MEMORY.md 踝段重写**：换掉三条过时条目 ——
   ① 「`[LRN-080]` 韧带 77% / 接触 28%」→ **改为已撤回**（实体应变能 147.9 vs 外力功 474.1 N·mm
   ⇒ ~70% 功在罚弹簧）
   ② 「刚度 ~2× 未归因」→ **改为三条建模伪影已定量**（缺距下关节 44.8% / 无摩擦切向滑移 39.5% /
   tie 罚滑移 13.9%）
   ③ 新增最重要的一条：**相对结论可用、绝对数值不可用**（姿势排序跨两套软骨网格一致）
   ④ 另加：骨折裕度首版、`area_proj` 修复、11 个翻转单元不可剔除、P1 应力独立验证、
   h5 nodesets ≠ .feb NodeSet、patch test 须用畸变单元
2. **`ANKLE_FIX_TODO.md` 顶部新增「🧭 接手入口」**（12 行：一句话状态 / 三件必知 /
   下一步按序 / 一键复现 / 章节索引 / 关键文件）→ 62,109 B
3. **`FE_PIPELINE_TODO.md` 迁移指针更新**（原指向 §第六轮 → 改指新入口）→ 31,409 B

**校验**：MEMORY 10,175/10,240 ✅｜LEARNINGS 102 / ERRORS 23 **标题级重号 0**｜
过时句 2 条已清 ✓｜关键内容 0 缺失｜**ALL PASS ✅**

**压缩后接手的一句话**：读 `ANKLE_FIX_TODO.md` 顶部「🧭 接手入口」即可；
**下一步 = 接入距下关节**（最大单一伪影 44.8%，软骨已生成 ⇒ 纯接线）。
"""
t = io.open(D, encoding="utf-8").read()
if "压缩前收口（2026-10-06 末）" in t:
    print("· 已有")
else:
    io.open(D, "a", encoding="utf-8").write(ADD)
    print(f"✓ 日记 -> {len(io.open(D, encoding='utf-8').read().encode('utf-8')):,} B")
