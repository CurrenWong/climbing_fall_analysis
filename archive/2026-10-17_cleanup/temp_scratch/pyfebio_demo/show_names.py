"""看已有的 THUMS 部件名元数据。"""
import io
import json
import os

R = r"D:\Project\climbing_fall_analysis"
P = os.path.join(R, "temp", "opensim_fe")

for f in ("parts\\_all_names.txt", "parts\\parts_index.json",
          "shells\\shell_pid_stats.csv", "shells\\shell_sections.json"):
    p = os.path.join(P, f)
    print("=" * 90)
    print(f"=== {f}   {os.path.getsize(p) if os.path.exists(p) else '不存在'}")
    print("=" * 90)
    if not os.path.exists(p):
        continue
    if f.endswith(".json"):
        d = json.load(io.open(p, encoding="utf-8"))
        if isinstance(d, dict):
            items = list(d.items())[:60]
            for k, v in items:
                s = json.dumps(v, ensure_ascii=False)[:110]
                print(f"  {k}: {s}")
        else:
            for v in d[:60]:
                print(f"  {json.dumps(v, ensure_ascii=False)[:130]}")
    else:
        t = io.open(p, encoding="utf-8", errors="ignore").read()
        ls = t.splitlines()
        print(f"  （{len(ls)} 行）")
        for l in ls[:70]:
            print(f"  {l[:126]}")
