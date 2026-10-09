"""跨进程并行 map —— 批量实验（probe / 参数扫描 / 蒙特卡洛）专用。

约定
----
* ``fn`` 必须是**模块级**函数，``items`` 元素须可 pickle；
  Windows 默认 spawn：调用方脚本必须有 ``if __name__ == "__main__"`` guard。
* 返回值严格按 ``items`` 顺序 —— 调用方先在主进程里按 RNG 流顺序
  生成好全部任务，再并行执行，即可保持输出与串行**逐位一致**。
* ``workers<=1`` 或单任务时退化为普通串行 for，方便调试。

用法::

    from climbing.batch import parallel_map
    rows = parallel_map(run_cell, tasks, progress=lambda i, n: print(f"{i}/{n}"))
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from typing import Any, Callable, Iterable


def default_workers(n_items: int) -> int:
    """默认并行度：留 2 核给系统，上限 14（16 逻辑核的常见超配点）。"""
    cpu = os.cpu_count() or 2
    return max(1, min(n_items, cpu - 2, 14))


def parallel_map(
    fn: Callable[[Any], Any],
    items: Iterable[Any],
    *,
    workers: int | None = None,
    chunksize: int = 1,
    progress: Callable[[int, int], None] | None = None,
) -> list:
    """进程池 map，保序返回。串行退化条件：workers<=1 或 items 长度 <2。"""
    items = list(items)
    if not items:
        return []
    w = workers or default_workers(len(items))
    if w <= 1 or len(items) < 2:
        out = []
        for i, it in enumerate(items):
            out.append(fn(it))
            if progress:
                progress(i + 1, len(items))
        return out

    ctx = mp.get_context("spawn")
    out: list = []
    with ProcessPoolExecutor(max_workers=w, mp_context=ctx) as ex:
        for i, res in enumerate(ex.map(fn, items, chunksize=chunksize)):
            out.append(res)
            if progress:
                progress(i + 1, len(items))
    return out
