"""罚参数收敛扫描：跑 FEBio 后直接读支反力，看是否收敛到共形模型的 436.8 N。

这是 tied-elastic 的**唯一判据**：L=0（底面复制、几何与共形模型完全相同）
时，penalty→∞ 必须收敛到共形解 436.817 N；收敛不到就说明 mortar 有系统偏差。
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

FEB = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
SRC = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio\ankle_tibiotalar.feb")
WORK = SRC.parent / "_fast"
AXIS = None


def patch(text: str, ov: dict) -> str:
    for k, v in ov.items():
        if k in ("time_steps", "step_size", "max_retries", "opt_iter", "max_refs"):
            text = re.sub(rf"<{k}>[^<]*</{k}>", f"<{k}>{v}</{k}>", text)
        else:
            def rep(m, k=k, v=v):
                blk = m.group(0)
                if re.search(rf"<{k}>[^<]*</{k}>", blk):
                    return re.sub(rf"<{k}>[^<]*</{k}>", f"<{k}>{v}</{k}>", blk)
                return blk.replace("</contact>", f"<{k}>{v}</{k}>\n      </contact>")
            text = re.sub(r"<contact name=\"tie_.*?</contact>", rep, text, flags=re.S)
    return text


def reaction(xplt: Path):
    import h5py
    import numpy as np
    import pyfebio.xplt as X
    import a6_reaction as A6          # 复用已验证的 tolerant parser
    X.parse_state = A6.tolerant_parse_state
    h5 = WORK / "_r.h5"
    X.to_hdf5(str(xplt), str(h5))
    ax = np.array([0.6766, 0.1363, 0.7237])
    ax = ax / np.linalg.norm(ax)
    with h5py.File(h5, "r") as f:
        sts = sorted(f["states"].keys(), key=int)
        st = sts[-1]
        R = np.asarray(f[f"states/{st}/node_data/reaction forces/all"])
        d = np.asarray(f[f"states/{st}/node_data/displacement/all"])
        grp = f[f"states/{st}/element_data/stress"]
        cart = None
        for k in grp:
            a = np.asarray(grp[k])
            if a.shape[0] == 458:
                s = a
                v = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2
                                   + (s[:, 1] - s[:, 2]) ** 2
                                   + (s[:, 2] - s[:, 0]) ** 2)
                            + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))
                cart = (float(np.median(v)), float(np.percentile(v, 99)),
                        float(v.max()))
        return float((R @ ax).sum()), len(sts), float((d @ ax).min()), cart


CASES = [50, 100, 300, 1000, 3000, 10000]

if __name__ == "__main__":
    base = SRC.read_text(encoding="utf-8")
    print(f"基准（共形模型，无 tie）: 436.817 N\n")
    allr = []
    for p in CASES:
        dst = WORK / f"p{p}.feb"
        dst.write_text(patch(base, {"penalty": p, "auto_penalty": 1,
                                    "max_retries": 0, "opt_iter": 8,
                                    "time_steps": 20}), encoding="utf-8")
        log = WORK / f"p{p}.log"
        t0 = time.time()
        with log.open("w", encoding="utf-8") as fh:
            rc = subprocess.run([str(FEB), str(dst)], stdout=fh,
                                stderr=subprocess.STDOUT,
                                timeout=1800).returncode
        dt = time.time() - t0
        rec = None
        if rc == 0:
            try:
                rec = reaction(WORK / f"p{p}.xplt")
            except Exception as e:      # noqa: BLE001
                rec = f"read-fail {e}"
        print(f"  penalty={p:<6} rc={rc}  {dt:5.1f}s  " +
              (f"R={rec[0]:8.3f} N  状态{rec[1]}  uz_min={rec[2]:+.4f}  "
               f"软骨 p99={rec[3][1]:.3f}" if isinstance(rec, tuple) else str(rec)),
              flush=True)
        allr.append({"penalty": p, "rc": rc, "R": rec[0] if isinstance(rec, tuple) else None,
                     "cart": rec[3] if isinstance(rec, tuple) else None})
    (WORK / "penalty_sweep.json").write_text(
        json.dumps(allr, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n目标 436.817 N；已扫 {[a['penalty'] for a in allr]}")
