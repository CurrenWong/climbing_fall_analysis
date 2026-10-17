"""tied-elastic 参数快扫：**不重建网格**，只在已写好的 .feb 上就地改参数再跑 FEBio。

为什么：网格（含 Jacobian 自检 + 写 1.9 MB 文件）约 15-20 s，而参数扫描
根本不需要动网格。分离后每轮只剩 FEBio 求解，可把 max_retries 压到 0 让它
一失败就退出（快速判生死），不必等满 20 步。

用法:
    .venv\\Scripts\\python.exe temp\\pyfebio_demo\\fasttie.py sweeps.json
或
    .venv\\Scripts\\python.exe temp\\pyfebio_demo\\fasttie.py   # 内置默认扫描
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

FEB = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
SRC = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio\ankle_tibiotalar.feb")
WORK = SRC.parent / "_fast"
WORK.mkdir(parents=True, exist_ok=True)


def patch(text: str, ov: dict) -> str:
    """把 ov 里的键写进 .feb；tied-elastic 段的 <penalty> 等需定位到 contact 块。"""
    for k, v in ov.items():
        if k in ("time_steps", "step_size", "max_retries", "opt_iter", "max_refs"):
            text = re.sub(rf"<{k}>[^<]*</{k}>", f"<{k}>{v}</{k}>", text)
        elif k in ("penalty", "tolerance", "search_tol", "search_radius",
                   "auto_penalty", "two_pass", "laugon", "symmetric_stiffness"):
            # 只改 tied-elastic 接触块里的参数，别碰 sliding-elastic
            def rep(m, k=k, v=v):
                blk = m.group(0)
                if re.search(rf"<{k}>[^<]*</{k}>", blk):
                    return re.sub(rf"<{k}>[^<]*</{k}>", f"<{k}>{v}</{k}>", blk)
                # 缺这个标签就插在块尾
                return blk.replace("</contact>", f"<{k}>{v}</{k}>\n      </contact>")
            text = re.sub(r"<contact name=\"tie_.*?</contact>", rep, text,
                          flags=re.S)
    return text


def run(tag: str, ov: dict) -> dict:
    dst = WORK / f"t_{tag}.feb"
    dst.write_text(patch(SRC.read_text(encoding="utf-8"), ov), encoding="utf-8")
    log = WORK / f"t_{tag}.log"
    t0 = time.time()
    with log.open("w", encoding="utf-8") as fh:
        p = subprocess.run([str(FEB), str(dst)], stdout=fh,
                           stderr=subprocess.STDOUT, timeout=1800)
    txt = log.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"(\d+)\s+negative jacobians", txt)
    f = re.search(r"failed to converge at time\s*:\s*([\d.]+)", txt)
    inv = re.findall(r'surface "([^"]+)" has (\d+) invalid facets', txt)
    inv = {a: b for a, b in inv if not a.startswith("surf_")}
    print(f"  rc={p.returncode}  {time.time()-t0:6.1f}s  "
          f"negJac={m.group(1) if m else '-':>4}  "
          f"failT={f.group(1) if f else '-':>8}"
          + (f"  invalid={inv}" if inv else ""), flush=True)
    return {"tag": tag, "rc": p.returncode, "neg": m.group(1) if m else None,
            "fail_t": f.group(1) if f else None, "invalid": inv,
            "log": str(log), "feb": str(dst)}


BASE = {"time_steps": 20, "max_retries": 0, "opt_iter": 8,
        "penalty": 10, "auto_penalty": 1, "tolerance": 1.0,
        "search_tol": 0.01, "search_radius": 1.0}

SWEEPS = [
    ("base",        {}),
    ("pen1",        {"penalty": 1}),
    ("pen10",       {"penalty": 10}),
    ("pen100",      {"penalty": 100}),
    ("tol0.02",     {"tolerance": 0.02}),
    ("tol0.1",      {"tolerance": 0.1}),
    ("srad0.3",     {"search_radius": 0.3}),
    ("srad3",       {"search_radius": 3.0}),
    ("stol0.001",   {"search_tol": 0.001}),
    ("laugon1",     {"laugon": 1}),
    ("pen1000",     {"penalty": 1000}),
    ("pen10_tol0.02_srad0.3", {"penalty": 10, "tolerance": 0.02,
                               "search_radius": 0.3}),
]

if __name__ == "__main__":
    if len(sys.argv) > 1:
        cases = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    else:
        cases = SWEEPS
    print(f"源 .feb: {SRC.stat().st_size/1e6:.2f} MB   "
          f"max_retries={BASE['max_retries']}（一失败就退）")
    out = []
    for tag, ov in cases:
        o = dict(BASE) | ov
        print(f"[{tag}] {' '.join(f'{k}={v}' for k, v in ov.items()) or '(基准)'}")
        try:
            out.append(run(tag, o))
        except subprocess.TimeoutExpired:
            print("  TIMEOUT")
            out.append({"tag": tag, "rc": -9})
    (WORK / "results.json").write_text(json.dumps(out, ensure_ascii=False, indent=1),
                                       encoding="utf-8")
    ok = [r for r in out if r.get("rc") == 0]
    print(f"\n=== {len(ok)}/{len(out)} 通过 ===")
    for r in out:
        print(f"  {'OK ' if r.get('rc') == 0 else 'xx '} {r['tag']}")
