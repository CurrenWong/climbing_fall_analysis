@echo off
setlocal
set "PATH=D:\Program\FEBioStudio\bin;%PATH%"
set "DRIVER=disp"
set "DISP_MM=0.1472"
set "TIME_STEPS=20"
set "XPLT=D:\Project\climbing_fall_analysis\temp\thums\febio\ankle_tibiotalar.xplt"
cd /d D:\Project\climbing_fall_analysis

echo ============ (a) original npz, TIE=0  (expect 436.8) ============
set "TIE=0"
set "CART_MESH=cartilage_mesh.npz"
.venv\Scripts\python.exe -u scripts\ankle_fe\build_feb.py > temp\pyfebio_demo\_iso_a.txt 2>&1
findstr /c:"FEBio rc=" temp\pyfebio_demo\_iso_a.txt
.venv\Scripts\python.exe temp\pyfebio_demo\steps_time.py %XPLT% 2>&1 | findstr /c:"1.00000"

echo ============ (b) t0d npz, TIE=0  (base free) ============
set "CART_MESH=cartilage_mesh_t0d.npz"
.venv\Scripts\python.exe -u scripts\ankle_fe\build_feb.py > temp\pyfebio_demo\_iso_b.txt 2>&1
findstr /c:"FEBio rc=" temp\pyfebio_demo\_iso_b.txt
.venv\Scripts\python.exe temp\pyfebio_demo\steps_time.py %XPLT% 2>&1 | findstr /c:"1.00000"

echo ============ (c) t0d npz, TIE=1 penalty=1e6 ============
set "TIE=1"
set "TIE_PENALTY=1000000"
.venv\Scripts\python.exe -u scripts\ankle_fe\build_feb.py > temp\pyfebio_demo\_iso_c.txt 2>&1
findstr /c:"tie" /c:"FEBio rc=" temp\pyfebio_demo\_iso_c.txt
.venv\Scripts\python.exe temp\pyfebio_demo\steps_time.py %XPLT% 2>&1 | findstr /c:"1.00000"
