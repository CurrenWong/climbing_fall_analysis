import io

p = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
t = io.open(p, encoding="utf-8").read()
b = len(t.encode("utf-8"))
for old, new in (
    ("**② 踝关节 FE 第二轮**：", "**② 踝关节 FE 二轮**："),
    ("12 条 gotcha 详情搬进新 ", "12 条 gotcha 搬进 "),
    ("（主表改「症状→一句话根因」速查）；", "（主表改「症状→根因」速查）；"),
    ("；**先行动后报告**（能跑就跑，跑出问题再回看）；", "；**先行动后报告**（能跑就跑）；"),
    ("**宁缺毋滥**（搜索硬性 3-5 条，凑数扣分）", "**宁缺毋滥**（搜索硬 3-5 条）"),
    ("- **负面结果也要如实报**：试过但假设错了的，标明\"假设错误 + 实测数据\"，\n  别包装成进展",
     "- **负面结果也要如实报**：试过但假设错了的，标明\"假设错误 + 实测数据\""),
):
    t = t.replace(old, new)
io.open(p, "w", encoding="utf-8").write(t)
a = len(t.encode("utf-8"))
ok = "OK, 余量 " + str(10240 - a) if a <= 10240 else "FAIL 差 " + str(a - 10240)
print(f"{b} -> {a} bytes  {ok}")
