@echo off
setlocal
set "PATH=D:\Program\FEBioStudio\bin;%PATH%"
set "DRIVER=disp"
set "DISP_MM=0.1472"
set "TIME_STEPS=20"
set "XPLT=D:\Project\climbing_fall_analysis\temp\thums\febio\ankle_tibiotalar.xplt"
cd /d D:\Project\climbing_fall_analysis

echo ===== 1) BASELINE regression (TIE=0, original npz) =====
set "TIE=0"
set "CART_MESH=cartilage_mesh.npz"
.venv\Scripts\python.exe -u scripts\ankle_fe\build_feb.py > temp\pyfebio_demo\_fin_base.txt 2>&1
findstr /c:"外层接触面" /c:"FEBio rc=" temp\pyfebio_demo\_fin_base.txt
.venv\Scripts\python.exe temp\pyfebio_demo\steps_time.py %XPLT% 2>&1 | findstr /c:"1.00000"

echo.
echo ===== 2) TIED-ELASTIC best point (penalty=4) =====
set "TIE=1"
set "CART_MESH=cartilage_mesh_t0d.npz"
set "TIE_PENALTY=4"
.venv\Scripts\python.exe -u scripts\ankle_fe\build_feb.py > temp\pyfebio_demo\_fin_tie.txt 2>&1
findstr /c:"tie" /c:"FEBio rc=" temp\pyfebio_demo\_fin_tie.txt
.venv\Scripts\python.exe temp\pyfebio_demo\steps_time.py %XPLT% 2>&1 | findstr /c:"1.00000"
