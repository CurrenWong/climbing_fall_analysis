# -*- coding: utf-8 -*-
"""载荷链**指纹** —— 让「拿旧基线比新结果」这件事变得可判定（闸门 10 的工具）。

**背景**（`results/opensim_fe/STALE_LOADCHAIN_AUDIT.md`）：归档的
``s2_thums_height_sweep{,_anat}.json`` 是 **2026-10-08 刚性载荷链标定之前**的产物；
`gauge_max` 与 `subtalar_kN` **同比例**变化（差 <0.3%）⇒ σ 的差异 **100% 由载荷解释**。
更麻烦的是**变化不是常数缩放**：``current/archived`` 随高度从 **1.1202（h=1）**
走到 **0.9400（h=10）**、再到 1.0117（h=50），**两次穿过 1.0**。

**⇒ 结论时的规矩（闸门 10）**：比任何两个数之前，先证它们**出自同一条链**。
本模块提供该判定的**最小工具**：把「载荷-高度表」压成一个短指纹。

设计要点
--------
* **纯函数、无 I/O**（除 ``table_from_json``）⇒ 可单测、无副作用。
* **只对载荷表取值**（不掺时间戳/文件名）⇒ 同链必得同指纹，跨链必不同。
* 指纹取 **sha256 前 16 位 hex**，对 ``(h, F_N)`` 排序后按固定格式序列化；
  数值按 **6 位有效数字** 归一，避免浮点末位噪声造成假不匹配。

用法::

    from climbing.coupling.loadchain import fingerprint, table_from_json, compare

    fp_now = fingerprint({1: 18029.4, 2: 22723.6, ...})
    tab = table_from_json("results/opensim_fe/s2_thums_height_sweep.json")
    print(compare(tab, current_table))
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path

__all__ = ["fingerprint", "table_from_json", "compare", "FP_LEN",
           "DEFAULT_LOAD_KEYS"]

#: 指纹长度（hex 字符数）
FP_LEN = 16

#: 从 sweep JSON 的每行里按顺序尝试的「载荷」字段名（单位 N 或 kN 由 :key 后缀决定）
DEFAULT_LOAD_KEYS = ("subtalar_kN", "subtalar_N", "sub_kN")


def _norm(x: float) -> str:
    """把载荷值归一成 6 位有效数字的字符串（吸收浮点末位噪声）。"""
    return "%.6g" % float(x)


def fingerprint(loads: dict[float, float], *, scale: float = 1.0) -> str:
    """载荷-高度表的短指纹。

    Parameters
    ----------
    loads : dict[float, float]
        ``{高度(m): 载荷}``。单位无所谓 —— 只要**同一链的所有调用方用同一单位**
        （``scale`` 可把 kN 统一成 N）。
    scale : float
        乘在载荷上的常数（如 kN→N 用 1000）。

    Returns
    -------
    str
        ``sha256`` 前 :data:`FP_LEN` 位 hex。
    """
    items = sorted((float(h), float(v) * scale) for h, v in loads.items())
    payload = "|".join("%s=%s" % (("%g" % h), _norm(v)) for h, v in items)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()[:FP_LEN]


def table_from_json(path: str | Path, *, act: float = 0.0,
                    keys: tuple[str, ...] = DEFAULT_LOAD_KEYS) -> dict[float, float]:
    """从 sweep JSON 里抽 ``{h: 载荷}``（只取 ``act`` 匹配的行）。

    兼容两种 JSON 形态：**行列表**（本仓库现行）与 ``{"rows": [...]}``。
    载荷字段按 ``keys`` 顺序找第一个存在的；字段名以 ``_kN`` 结尾者自动 ×1000 归一为 N。
    """
    obj = json.load(io.open(path, encoding="utf-8"))
    rows = obj.get("rows", []) if isinstance(obj, dict) else obj
    out: dict[float, float] = {}
    for r in rows:
        if float(r.get("act", 0.0)) != float(act):
            continue
        for k in keys:
            if k in r:
                v = float(r[k]) * (1000.0 if k.endswith("_kN") else 1.0)
                out[float(r["h"])] = v
                break
    if not out:
        raise ValueError("在 %s 里找不到载荷列（试过 %s）" % (path, list(keys)))
    return out


def compare(a: dict[float, float], b: dict[float, float]) -> dict:
    """比较两张表：**先在共有高度上比指纹**，再给出逐高度比值。

    ⚠️ 指纹必须建在 **``set(a) & set(b)``** 上：两张表的高度集合往往不同
    （一份只跑了 h=1..6、另一份 h=1..50），直接整表取指纹会把「同链」误判成「异链」。

    ``ratio`` 在共有高度上应为常数（±1e-3）才算**同链**；否则**形状变了**。
    """
    common = sorted(set(a) & set(b))
    aa = {h: a[h] for h in common}
    bb = {h: b[h] for h in common}
    fa, fb = fingerprint(aa), fingerprint(bb)
    ratios = {h: (bb[h] / aa[h]) for h in common if aa[h]}
    lo = min(ratios.values()) if ratios else None
    hi = max(ratios.values()) if ratios else None
    spread = (hi - lo) if (lo is not None) else None
    return {
        "fp_a": fa, "fp_b": fb, "same_fingerprint": fa == fb,
        "n_a": len(a), "n_b": len(b), "n_common": len(common),
        "common_heights": common,
        "ratio_min": lo, "ratio_max": hi, "ratio_spread": spread,
        "ratio_span_pct": (100.0 * spread / ((lo + hi) / 2.0)
                           if spread is not None and (lo + hi) else None),
        "is_constant_scale": (spread is not None and spread < 1e-3),
        "ratios": ratios,
    }
