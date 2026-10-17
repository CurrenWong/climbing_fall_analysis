"""放宽求解控制：max_refs 25->100、opt_iter 15->25、加自动步长控制。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = """w('    <time_stepper type="default"><max_retries>10</max_retries>'
  '<opt_iter>15</opt_iter><cutback>0.25</cutback></time_stepper>')
w('    <solver type="solid"><linear_solver type="pardiso"/>'
  '<max_refs>25</max_refs><reform_each_time_step>1</reform_each_time_step></solver>')"""
assert s.count(old) == 1, f"hit {s.count(old)}"
s = s.replace(old, """w('    <time_stepper type="default"><max_retries>20</max_retries>'
  '<opt_iter>25</opt_iter><cutback>0.25</cutback>'
  '<dtmin>1e-6</dtmin><dtmax>0.05</dtmax></time_stepper>')
# ★ max_refs 必须够大：实测每步要 ~97 次牛顿迭代，25 次会在第 26 次
#   "Max nr of reformations reached" 直接终止（这曾被误当成奇异/机械问题）。
w('    <solver type="solid"><linear_solver type="pardiso"/>'
  '<max_refs>100</max_refs><reform_each_time_step>1</reform_each_time_step>'
  '<diverge_reform>1</diverge_reform></solver>')""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
