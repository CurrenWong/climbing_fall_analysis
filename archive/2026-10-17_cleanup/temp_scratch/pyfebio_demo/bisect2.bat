@echo off
cd /d D:\Project\climbing_fall_analysis
echo ==== 1) prism + NO ligaments ====
set LIG_MODE=off
set DRIVER=disp
set TRUNC_MM=60
set FIX_TALUS=all
set CART_ELEM=prism
set TIME_STEPS=5
set DISP_MM=0.10
.venv\Scripts\python.exe -u scripts\ankle_fe\ankle_lig_feb.py > temp\pyfebio_demo\_b1.txt 2>&1
findstr /C:"FEBio rc" /C:"Negative" /C:"initialization" temp\pyfebio_demo\_b1.txt

echo ==== 2) tet4 baseline (expect OK) ====
set LIG_MODE=shell
set CART_ELEM=tet4
.venv\Scripts\python.exe -u scripts\ankle_fe\ankle_lig_feb.py > temp\pyfebio_demo\_b2.txt 2>&1
findstr /C:"FEBio rc" /C:"Negative" /C:"initialization" temp\pyfebio_demo\_b2.txt

echo ==== 3) prism + ligaments ====
set CART_ELEM=prism
.venv\Scripts\python.exe -u scripts\ankle_fe\ankle_lig_feb.py > temp\pyfebio_demo\_b3.txt 2>&1
findstr /C:"FEBio rc" /C:"Negative" /C:"initialization" temp\pyfebio_demo\_b3.txt
