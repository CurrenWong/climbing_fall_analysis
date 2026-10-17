@echo off
setlocal
set "PATH=D:\Program\FEBioStudio\bin;%PATH%"
set "DRIVER=disp"
set "DISP_MM=0.1472"
set "TIME_STEPS=20"
set "TIE=1"
set "CART_MESH=cartilage_mesh_t0d.npz"
cd /d D:\Project\climbing_fall_analysis

for %%T in (0.05 0.2 0.5 1.5) do (
  set "TIE_SEC_TOL=%%T"
  .venv\Scripts\python.exe -u scripts\ankle_fe\build_feb.py > temp\pyfebio_demo\_sec_%%T.txt 2>&1
  echo === TOL=%%T ===
  findstr /c:"tie" /c:"FEBio rc=" temp\pyfebio_demo\_sec_%%T.txt
)
