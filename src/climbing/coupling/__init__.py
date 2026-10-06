"""OpenSim ↔ FEBio 耦合子包（方案见 ``docs/OpenSim_FE复现方案.md``）。

模块
----
febio_run   FEBio 可执行文件定位与安全调用（S0.2）。

设计原则沿用仓库约定：**让"没找到 / 没跑成"不可能静默通过**。
"""

from __future__ import annotations
