"""加 TOPFIX 开关：force 驱动下顶面是否加横向约束。

吴恺的机理（论文 §3 讨论原文）：
  「外踝较内踝较长并靠后 ⇒ 踝轴向后下方成角 ⇒ 负重时小腿外翻 + 轻微内旋
    ⇒ **外侧跟腓位移 < 内侧，胫跟 > 胫舟**」
⇒ 排序来自**耦合的外翻+内旋**，不是纯轴向压缩。
而 DRIVER=disp 把顶面三自由度全指定（夹持）⇒ **禁掉了这个机理本身**。
⇒ force 驱动 + 顶面**不加横向约束**（只靠接触 + 韧带约束）才是吴恺的配置。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = 'LOADDIR = os.environ.get("LOADDIR", "axis").strip() or "axis"'
assert s.count(old) == 1
s = s.replace(old, old + '''
# force 模式下顶面横向是否约束：none=不加（吴恺配置，允许外翻/内旋）
TOPFIX = os.environ.get("TOPFIX", "lateral").strip() or "lateral"''')

old2 = """    if not _faces:      # 空 NodeSet 不能写 BC（FEBio: invalid value for node_set）
        continue"""
assert s.count(old2) == 1
s = s.replace(old2, """    if not _faces:      # 空 NodeSet 不能写 BC（FEBio: invalid value for node_set）
        continue
    if DRIVER == "force" and TOPFIX == "none":
        continue        # ★ 顶面完全不约束：靠接触 + 韧带定运动学（允许外翻/内旋）""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
