# -*- coding: utf-8 -*-
"""pyfebio ``xplt`` 解析器补丁：容忍**无名域**。

**问题（本仓库实测，2026-10-16）**

``pyfebio/xplt.py:675`` 的 ``parse_mesh`` 里::

    mesh_dict["domains"][domain["id"][0]] = domain["name"]

**没有空名兜底**，而同一函数里的 **nodeset 分支有**::

    if nodeset["name"] == "":
        nodeset["name"] = nodeset["id"][0]

**触发条件**：FEBio 在**含接触的模型**里会多写域（实测 ``plantar_bc="contact"`` 的
跟骨模型：``calcaneus``/``trabecular``/**``pad_floor``**/一个 **``id=4, etype=5,
nelems=573`` 的无名域**）⇒ ``parse_mesh`` 抛 ``KeyError: 'name'`` ⇒
``to_hdf5`` 只写出 ``meshes``、**没有 ``states``** ⇒
下游 ``fe_post`` / ``s2_thums_sweep._read_domain_vm`` **全部失效**。

**影响面**：**所有** ``plantar_bc="contact"``（以及任何带接触的 FEBio 模型）的结果
后处理。默认 ``plantar_bc="fixed"`` 的 deck 无无名域 ⇒ 逐位不受影响。

**为什么用 shim 而不改 site-packages**：改 ``.venv`` 里的第三方文件会在重装/升级时
静默丢失，且遮蔽真实上游缺陷。本 shim 显式、可测、可撤。

用法::

    from climbing.coupling.xplt_shim import patch_unnamed_domains
    patch_unnamed_domains()          # 幂等
    from pyfebio import xplt as xplt_mod
    xplt_mod.to_hdf5("m.xplt", "m.hdf5")
"""
from __future__ import annotations

__all__ = ["patch_unnamed_domains", "is_patched"]

_PATCHED = False


def patch_unnamed_domains() -> bool:
    """给 ``pyfebio.xplt._parse_domain_section`` 加空名兜底。幂等。

    Returns
    -------
    bool
        本次调用是否**新做了**补丁（``False`` 表示之前已打过）。
    """
    global _PATCHED
    if _PATCHED:
        return False

    from pyfebio import xplt as X

    orig = X._parse_domain_section

    def _patched(buffer):
        domains = orig(buffer)
        for k, d in enumerate(domains):
            if d.get("name") in (None, "", b""):
                # 与 nodeset 分支同款兜底：用域 id；id 也没有时用序号。
                did = d.get("id")
                try:
                    did = int(did[0]) if hasattr(did, "__len__") else int(did)
                except Exception:                      # noqa: BLE001
                    did = k + 1
                d["name"] = f"__unnamed_{did}"
        return domains

    X._parse_domain_section = _patched
    _PATCHED = True
    return True


def is_patched() -> bool:
    return _PATCHED
