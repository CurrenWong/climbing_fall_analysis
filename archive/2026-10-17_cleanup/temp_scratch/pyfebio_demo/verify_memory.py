"""MEMORY.md 整理后的验证：尺寸 / 死链 / 信息是否丢失。"""
import re
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
t = M.read_text(encoding="utf-8")
size = len(t.encode("utf-8"))

print(f"=== 尺寸 ===\nMEMORY.md {size:,} bytes / 上限 10240 "
      f"→ {'✅ 合规' if size <= 10240 else '❌ 超标'}"
      f"（整理前 {Path(WS/'memory_archive'/'MEMORY.md.bak-20261004').stat().st_size:,}）")

print("\n=== 死链检查（MEMORY.md 里提到的每个本地路径）===")
paths = sorted(set(re.findall(r"`([A-Za-z0-9_\-./\\:]+\.(?:md|py|json|npz|k|pdf))`", t)))
bad = 0
for p in paths:
    cands = [WS / p, Path(r"D:\Project\climbing_fall_analysis") / p,
             Path(r"D:\Project\SecLog") / p, Path(p)]
    hit = next((c for c in cands if c.exists()), None)
    if hit is None:
        print(f"  ❌ 找不到: {p}")
        bad += 1
print(f"  {len(paths)} 个路径，{bad} 个失效" + ("  ✅" if bad == 0 else ""))

print("\n=== references/ 目录实况 vs 索引 ===")
refs = {p.name for p in (WS / "references").glob("*.md")}
listed = set(re.findall(r"`references/([\w\-.]+\.md)`", t))
for r in sorted(refs):
    print(f"  {'索引中' if r in listed else '⚠️ 未列入索引'}  {r}")
extra = listed - refs
if extra:
    print(f"  ❌ 索引指向不存在的文件: {sorted(extra)}")

print("\n=== 信息是否丢失：12 条 SecLog gotcha 是否都在新 reference 里 ===")
seclog = (WS / "references/seclog-miniprogram-gotchas.md").read_text(encoding="utf-8")
KEYS = {
    "冷启动缓存/二次invoke": "二次 invoke",
    "compile.js upload/假发布": "compile.js upload",
    "tcb CLI fn deploy 弹窗": "fn deploy",
    "hy3 纯文本/幻觉": "纯文本",
    "DeepSeek 迁移": "deepseek-flash",
    "TokenHub 额度包 403": "pkg_tcb_tokencredits_",
    "glm-5v-turbo 未启用": "UpdateAIModel",
    "fileID 真 bucket": "真 bucket",
    "fileID 不抛异常": "不抛异常",
    "UTC 月边界": "TZ_OFFSET",
    "list days 跨月": "month:'this'",
    "round2/Number.EPSILON": "Number.EPSILON",
    "WXML 不能调方法": "xxxText",
    "miniprogram-ci ?. ": "Unexpected token",
    "input min-height": "min-height: 88rpx",
    "pyfebio/xplt": "",
}
for k, probe in KEYS.items():
    if not probe:
        continue
    ok = probe in seclog
    print(f"  {'✅' if ok else '❌'} {k}")

print("\n=== gotcha 20/21 细节仍在 references ===")
fe = (WS / "references/fe-modeling-gotchas.md").read_text(encoding="utf-8")
for probe in ("奇置换", "2-流形" if "2-流形" in fe else "刚体模态", "NORMAL TERMINATION"):
    print(f"  {'✅' if probe in fe else '❌'} fe: {probe}")
digest = WS / "digest/procedure/git-history-secret-purge.md"
print(f"  {'✅' if digest.exists() else '❌'} git-history-secret-purge.md "
      f"({digest.stat().st_size if digest.exists() else 0:,} bytes)")
