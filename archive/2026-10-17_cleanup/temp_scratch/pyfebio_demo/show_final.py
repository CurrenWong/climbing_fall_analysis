from pathlib import Path

t = Path(r"temp/pyfebio_demo/_v.txt").read_text(encoding="utf-8")
print(t[-1300:])
print("=== MEMORY ankle section ===")
s = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md").read_text(encoding="utf-8")
KEYS = ("FE_PIPELINE_TODO", "tied-elastic", "刚度 2 倍", "根因=接触区",
        "Anderson 真值", "韧带线", "LRN-046", "cart_refine")
for ln in s.split("\n"):
    if any(k in ln for k in KEYS):
        print("  " + ln)
