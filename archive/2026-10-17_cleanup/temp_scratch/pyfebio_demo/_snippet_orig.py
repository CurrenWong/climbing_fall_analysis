for surf, F, _faces in (("tibia_top", F_t, top_tibia),):
    if not _use_force or not _faces or F <= 0:
        continue
        # traction = 力 / 面积
        a = area_proj(top_tibia)
        tr = F / a if a > 0 else 0.0
        w("    <surface_load ...>")
