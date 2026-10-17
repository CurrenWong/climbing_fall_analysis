"""收尾：修横幅残留 + MEMORY 压到 ≤10240（删已陈旧内容）。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
TOD = os.path.join(R, "ANKLE_FIX_TODO.md")
MEM = os.path.join(W, "MEMORY.md")
USE = os.path.join(R, "docs", "踝FE可用模型使用说明.md")

# ── 1) 横幅里的错说法（我自己加的，最该先改）──
t = io.open(TOD, encoding="utf-8").read()
OLD_B = ("> **⚠️ `POSE_IE` 打印的\"内翻+\"是反的：+IE 实际是外翻。**\n")
NEW_B = ("> **⚠️ `POSE_IE` 符号（2026-10-06 已再更正）**：`check_c2.py` 的两个**解剖锚定**判据\n"
         "> （足底法向 / 足缘高度，\"外侧\"由**腓骨**定义）证明 **+IE = 内翻**，**代码注释正确**。\n"
         "> 先前写的\"+IE 实际是外翻\"**已作废**（`[ERR-20261006-002]`，根因=用数据反推标签的循环论证）。\n")
if OLD_B in t:
    t = t.replace(OLD_B, NEW_B)
    io.open(TOD, "w", encoding="utf-8").write(t)
    print(f"✓ 横幅已更正 -> {len(t.encode('utf-8'))} B")
else:
    print("⚠ 横幅目标串未命中")

# ── 2) MEMORY 压缩：删已被取代的旧结论 ──
m = io.open(MEM, encoding="utf-8").read()
OLD1 = """- ⚠️ **F 未收敛但非伪影**（`[LRN-054/056]`）：k 随罚因子 49/151/451/880 N/mm vs Anderson 估计 1059；
  **p≥2 全负 Jacobian ⇒ 无法验收敛**；`436.8N` **偏低低**。`|u|max` 非判据（悬臂弯 0.17°）；
  **几何没问题**（互穿仅 6 节点；接触 31→**344mm²**≈生理 300mm²）
"""
NEW1 = """- ⚠️ 旧「F 未收敛/罚因子伪影/`|u|max` 判据/几何质疑」**均已被 2026-10-05 的符号修正取代**：
  压缩侧现全 `NORMAL`。**几何没问题**（互穿仅 6 节点；接触 31→344mm²≈生理 300mm²）。详见 `references/`
"""
OLD2 = """**最后更新**: 2026-10-05　🚨 **踝线重大更正（`[ERR-20261005-009]`）**：
**加载方向符号一直反了** —— `-LAX` 把关节**拉开**（拉伸）而非压实，历次「力控制崩在 ~12 N」的真因；
历史"成功值" 187.77/**592.6 N 全是拉伸**。
**正解 `LOAD_SIGN=+1` + `DISP_MM<0`（压缩）⇒ −1.0 mm = −608.82 N，全程 NORMAL**。
另修：`lc` 挂 `<scale>`、`sym_stiffness` 在 `<solver>`、`_outer()` 新格式崩溃。
**整理者**: Spark ⚡　- **详见**：`docs/踝FE可用模型使用说明.md`"""
NEW2 = """**最后更新**: 2026-10-06　✅ **C2 重定义完成并 PASS**（`+IE=内翻`，代码注释本来就对；
我先前"标签反了"**已作废**）；🎯 **韧带应变阈值只能用 (c) 总应变**（(b) 增量与风险反相关）；
⚠️ **外翻不收敛 ⇒ 目前只有内翻侧可用**。详见 `references/ankle-fe-criteria.md`。
（上一轮：加载方向符号修正 `[ERR-20261005-009]` ⇒ 压缩侧 `DISP_MM=-1.0` = **−608.82 N** 全 NORMAL。）
**整理者**: Spark ⚡"""
for old, new, tag in ((OLD1, NEW1, "陈旧F结论"), (OLD2, NEW2, "最后更新块")):
    if old in m:
        m = m.replace(old, new)
        print(f"· MEMORY 替换【{tag}】")
    else:
        print(f"⚠ MEMORY 未命中【{tag}】")
io.open(MEM, "w", encoding="utf-8").write(m)
n = len(m.encode("utf-8"))
print(f"MEMORY.md -> {n} B  {'✅ ≤10240' if n <= 10240 else f'❌ 超 {n-10240}'}")

# ── 3) 全程残留复查 ──
print("\n--- 残留复查 ---")
for f, nm in ((MEM, "MEMORY"), (TOD, "ANKLE_FIX_TODO"), (USE, "使用说明")):
    x = io.open(f, encoding="utf-8").read()
    bad = [p for p in ("实际是外翻", "−IE 才是内翻", "-IE 才是内翻", "是反的")
           if p in x]
    print(f"  {nm:<16} {'✓ 干净' if not bad else '⚠ ' + str(bad)}")
