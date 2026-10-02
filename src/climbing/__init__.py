"""攀岩坠落仿真包。

子模块
------
rope      Phase A-1  绳索坠落 1D 模型（非线性绳 + 保护器两态）
belay     Phase A-1  保护器（ATC / Grigri / 8 字结 / Reverso）摩擦与保持力
pad       Phase B-1  抱石软垫非线性 + 渐进接触面积模型
metrics   Phase A-4/B-4 伤害指标（HIC、峰值 g、阈值判定）
"""

__version__ = "0.1.0"

G = 9.80665  # 标准重力加速度 m/s^2
