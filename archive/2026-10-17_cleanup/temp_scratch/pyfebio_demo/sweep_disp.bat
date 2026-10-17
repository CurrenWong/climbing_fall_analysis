@echo off
setlocal
set "PATH=D:\Program\FEBioStudio\bin;%PATH%"
set "DRIVER=disp"
set "TIME_STEPS=20"
set "TIE=1"
set "TIE_AUTO_PENALTY=1"
set "TIE_PENALTY=10"
set "CART_MESH=cartilage_mesh_t0d.npz"
cd /d D:\Project\climbing_fall_analysis

for %%D in (0.02 0.05 0.10) do (
  set "DISP_MM=%%D"
  .venv\Scripts\python.exe -u scripts\ankle_fe\build_feb.py > temp\pyfebio_demo\_td_%%D.txt 2>&1
  echo === DISP=%%D ===
  findstr /c:"FEBio rc=" temp\pyfebio_demo\_td_%%D.txt
  findstr /c:"negative jacobians" temp\pyfebio_demo\_td_%%D.txt
  findstr /c:"failed to converge" temp\pyfebio_demo\_td_%%D.txt
)
