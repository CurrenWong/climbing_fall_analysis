"""重写 Control 段为手册 §3.3 的推荐配置（含 BFGS 拟牛顿、显式容差）。

原来缺的东西：
  · <qn_method type="BFGS"><max_ups>10</max_ups>  —— 拟牛顿加速，FEBio 的标准做法
  · 显式 <dtol>/<etol>/<rtol>/<lstol>/<min_residual>/<rhoi>
实测：缺这些时每步要 ~90 次刚度重生成、3800 次 RHS 求值仍在磨
（残差判据 REQUIRED=0 永不满足，只能靠位移判据，而它又过严）。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
L = p.read_text(encoding="utf-8").split("\n")

i0 = next(i for i, l in enumerate(L) if l.strip() == "w('  <Control>')")
i1 = next(i for i, l in enumerate(L) if l.strip() == "w('  </Control>')")

NEW = [
    "w('  <Control>')",
    "w(f'    <analysis>STATIC</analysis><time_steps>{TIME_STEPS}</time_steps>'",
    "  f'<step_size>{1.0/TIME_STEPS:.6f}</step_size>')",
    "w('    <solver type=\"solid\">')",
    "w('      <linear_solver type=\"pardiso\"/>')",
    "w('      <max_refs>50</max_refs>')",
    "w('      <diverge_reform>1</diverge_reform>')",
    "w('      <reform_each_time_step>1</reform_each_time_step>')",
    "w('      <dtol>0.01</dtol>')",
    "w('      <etol>0.001</etol>')",
    "w('      <rtol>0</rtol>')",
    "w('      <lstol>0.9</lstol>')",
    "w('      <min_residual>1e-20</min_residual>')",
    "w('      <rhoi>-2</rhoi>')",
    "w('      <qn_method type=\"BFGS\"><max_ups>10</max_ups></qn_method>')",
    "w('    </solver>')",
    "w('    <time_stepper>')",
    "w('      <dtmin>0.0001</dtmin><dtmax>0.05</dtmax>')",
    "w('      <max_retries>10</max_retries><opt_iter>10</opt_iter>')",
    "w('    </time_stepper>')",
    "w('  </Control>')",
]
L[i0:i1 + 1] = NEW
p.write_text("\n".join(L), encoding="utf-8")
s = p.read_text(encoding="utf-8")
ast.parse(s)
print("Control 段已重写 + syntax OK")
