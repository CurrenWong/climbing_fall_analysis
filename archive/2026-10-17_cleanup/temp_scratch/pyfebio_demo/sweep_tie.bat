@echo off
setlocal
set "PATH=D:\Program\FEBioStudio\bin;%PATH%"
set "DRIVER=disp"
set "DISP_MM=0.1472"
set "TIME_STEPS=20"
set "TIE=1"
set "CART_MESH=cartilage_mesh_t0d.npz"
cd /d D:\Project\climbing_fall_analysis

for %%P in (1 10 100) do (
  for %%A in (1 0) do (
    set "TIE_PENALTY=%%P"
    set "TIE_AUTO_PENALTY=%%A"
    .venv\Scripts\python.exe -u scripts\ankle_fe\build_feb.py > temp\pyfebio_demo\_ts_%%P_%%A.txt 2>&1
    echo PEN=%%P AUTO=%%A
    findstr /c:"FEBio rc=" temp\pyfebio_demo\_ts_%%P_%%A.txt
  )
)
