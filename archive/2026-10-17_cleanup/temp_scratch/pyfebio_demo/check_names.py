"""从 parts_index.json 查 8110xxxx 的权威名称 —— 验证 ATFL/CFL 映射是否错。"""
import io
import json
import os

R = r"D:\Project\climbing_fall_analysis"
p = os.path.join(R, "temp", "opensim_fe", "parts", "parts_index.json")
d = json.load(io.open(p, encoding="utf-8"))
print(f"parts_index.json 记录 {len(d)} 条")
by_pid = {}
for r in d:
    by_pid[int(r["pid"])] = r

LIG_DEF = [
    (81100701, "PTFL_talofibular_posterior"),
    (81100801, "CFL_calcaneofibular"),
    (81100901, "ATFL_talofibular_anterior"),
    (81101001, "PTCL_talocalcaneal_posterior"),
    (81101101, "SYND_tibiofibular"),
    (81101201, "MD1_tibiotalar_post"),
    (81101202, "MD2_tibiocalcaneal"),
    (81101203, "MD3_tibionavicular"),
    (81101204, "MD4_tibiotalar_ant"),
    (81101409, "SPRING_calcaneonavicular"),
]
print(f"\n{'我们的 pid':>10} {'我们用的名':<30} {'THUMS 权威 title':<40} 一致?")
print("─" * 104)
for pid, our in LIG_DEF:
    r = by_pid.get(pid)
    t = r["title"] if r else "**（该 pid 不在表中）**"
    ok = "?"
    if r:
        tl = t.lower()
        key = our.split("_")[0].lower()
        ok = "✓" if key in tl else "⚠ 看内容"
    print(f"{pid:>10} {our:<30} {t:<40} {ok}")

# 额外：列出所有 810 与 811 的 pid+title（脚踝区）
print("\n=== 81100xxx–81101xxx 全部（韧带区）===")
for pid in sorted(by_pid):
    if 81100000 <= pid < 81102000:
        print(f"  {pid:>10}  {by_pid[pid]['title']}")
print("\n=== 81000xxx–81002xxx 全部（骨区）===")
for pid in sorted(by_pid):
    if 81000000 <= pid < 81003000:
        print(f"  {pid:>10}  {by_pid[pid]['title']}")
