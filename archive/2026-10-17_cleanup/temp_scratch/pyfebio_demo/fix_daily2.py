"""修 daily 索引（守卫改用 LRN ID，不用节号）+ 报告。"""
from pathlib import Path

D = Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-04.md")
t = D.read_text(encoding="utf-8")
MARK = "LRN-20261005-059"
LINE = ("\n- `2026-10-04/ankle-fe-solve-and-handover.md` → **§十七(10-05)** 执行 ANKLE_GOAL "
        "Step1~2（`[LRN-20261005-059]` / `[ERR-20261005-007]`）：`<Loads>` **恒为空是真 bug**"
        "（`continue` 后死代码 ⇒ `DRIVER=force` 载荷被静默丢弃）；力控制失败真因=**关节初始张开**"
        "（**0.6 N 也不收敛 463 s** ⇒ 与载荷大小无关，非"刚体模态"；`FIX_TALUS=all`/`FIB_SHARE=0`/"
        "`NO_CONTACT=1` 全无效）；找到**稳定的姿势可用配置 `DRIVER=dispx` + `FIX_TALUS=transverse`**"
        "（NORMAL/5 s/0 负 Jacobian，与 `all` 响应显著不同 90 vs 162 N）；补**文本读数出口 `LOG_DATA`**"
        "（语法出 FEBio 手册 §3.19.1 + 附录 E.2.1）+ 验收脚本 `scripts/ankle_fe/acceptance.py`；"
        "⚠️ 未解释：transverse 反力在 x、all 反力在 z（顶面 BC 同为 prescribe x）\n")
if MARK not in t:
    D.write_text(t.rstrip() + LINE, encoding="utf-8")
    print(f"✓ daily 已追加  -> {D.stat().st_size:,} B")
else:
    print("daily 已含该索引（跳过）")
tt = D.read_text(encoding="utf-8")
print(f"  [{'OK' if MARK in tt else 'XX'}] daily 含 {MARK}")
print(f"  daily 复检节号出现次数: §十七 = {tt.count('§十七')}（含自动笔记的 1 处，属正常）")
