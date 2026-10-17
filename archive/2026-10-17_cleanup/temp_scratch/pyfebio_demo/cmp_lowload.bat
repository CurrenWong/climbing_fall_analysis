@echo off
cd /d D:\Project\climbing_fall_analysis
echo ==== A) prism @ 0.1472mm (compare with tet4 baseline 140N) ====
set LIG_MODE=shell
set DRIVER=disp
set TOPFIX=lateral
set TRUNC_MM=60
set FIX_TALUS=all
set CART_ELEM=prism
set TIME_STEPS=20
set DISP_MM=0.1472
.venv\Scripts\python.exe -u scripts\ankle_fe\ankle_lig_feb.py > temp\pyfebio_demo\_a1.txt 2>&1
findstr /C:"FEBio rc" /C:"tibiotalar" /C:"接触面" temp\pyfebio_demo\_a1.txt
.venv\Scripts\python.exe temp\pyfebio_demo\read_last.py
