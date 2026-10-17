import h5py, numpy as np

H5 = r"temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = r"temp/thums/febio_lig/ankle_lig_WK.feb"

# --- feb head
with open(FEB, "r", errors="replace") as f:
    lines = f.readlines()
print("feb lines:", len(lines))
print("---- head 25 ----")
for l in lines[:25]:
    print(repr(l[:160]))
print("---- Elements headers ----")
for i, l in enumerate(lines):
    s = l.strip()
    if s.startswith("<Elements") or s.startswith("</Elements") or s.startswith("<NodeSet") or s.startswith("<Solid") or s.startswith("<Surface") or s.startswith("<Material") or s.startswith("<contact") or s.startswith("<rigid") or s.startswith("<bc") or s.startswith("<DiscreteElement"):
        print(i, s[:170])
