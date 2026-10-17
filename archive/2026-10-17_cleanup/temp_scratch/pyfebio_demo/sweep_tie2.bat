@echo off
setlocal
set "PATH=D:\Program\FEBioStudio\bin;%PATH%"
set "DRIVER=disp"
set "DISP_MM=0.1472"
set "TIME_STEPS=20"
set "TIE=1"
set "TIE_AUTO_PENALTY=1"
set "CART_MESH=cartilage_mesh_t0d.npz"
cd /d D:\Project\climbing_fall_analysis

for %%P in (10 50 100 300) do (
  set "TIE_PENALTY=%%P"
  .venv\Scripts\python.exe -u scripts\ankle_fe\build_feb.py > temp\pyfebio_demo\_tp_%%P.txt 2>&1
  echo === PEN=%%P ===
  findstr /c:"FEBio rc=" temp\pyfebio_demo\_tp_%%P.txt
  findstr /c:"negative jacobians" /c="failed to converge" temp\pyfebio_demo\_tp_%%P.txt
)
