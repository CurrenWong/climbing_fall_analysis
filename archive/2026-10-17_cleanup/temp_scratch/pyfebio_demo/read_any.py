"""读 xplt 的末态支反力（沿加载轴求和）+ 软骨 p99 应力。可用于任意文件。"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from a6_reaction import *      # noqa: F403  复用已验证的读取逻辑
import a6_reaction as A6


def report(fp):
    fp = Path(fp)
    print(f"\n=== {fp.name} ===")
    A6.main(fp) if hasattr(A6, "main") else None


if __name__ == "__main__":
    for a in sys.argv[1:]:
        print(a)
