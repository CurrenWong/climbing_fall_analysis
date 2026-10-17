"""门（gate）的**上下文管理器** —— 把「改全局 dict → 用完恢复」变成**机制**，而不是靠人写 `try/finally`。

背景（`results/STRUCT_AUDIT_2026-10-17.md` §③）：`src/climbing/pad.py` 有 **4 个模块级可变 dict**
（`RIGID_STANDING_CHAIN` / `RIGID_CONTACT` / `STROKE_ANCHOR` / `TWO_SEG_LEG`），**8 个脚本在改它们**。
只要有一处忘恢复，**后续所有运行被静默污染** —— 本会话已经踩过同类坑
（`TWO_SEG_LEG` 未清零、`m_theta_table` 传错格式等）。

用法（**推荐**）::

    from climbing.gates import gate
    import climbing.pad as P

    with gate(P.TWO_SEG_LEG, k_theta_nm_per_rad=350.0, segment_l_m=0.425,
              m_theta_table=TABLE):
        r = P.simulate_boulder_fall(...)     # 退出时**自动**按原键值恢复（含"原本不存在"⇒ 删除）

多门同时::

    with gate(P.TWO_SEG_LEG, ...), gate(P.RIGID_CONTACT, k_n_per_m=2.1e6):
        ...

实现要点：
  * 只恢复**本上下文改过的键**（不整表回滚）⇒ 可与其它作用域嵌套；
  * `None` 值 ⇒ **删除该键**（把"门关闭"表达为键不存在），比 `=None` 更干净；
  * 恢复在 `finally` 里 ⇒ 异常路径同样恢复。
"""
from __future__ import annotations

import contextlib

__all__ = ["gate"]

_SENTINEL = object()


@contextlib.contextmanager
def gate(d: dict, **overrides):
    """临时覆盖模块级门 dict `d` 的若干键；退出时按原状恢复。"""
    saved = {k: d.get(k, _SENTINEL) for k in overrides}
    try:
        for k, v in overrides.items():
            if v is None:
                d.pop(k, None)
            else:
                d[k] = v
        yield d
    finally:
        for k, v in saved.items():
            if v is _SENTINEL:
                d.pop(k, None)
            else:
                d[k] = v
