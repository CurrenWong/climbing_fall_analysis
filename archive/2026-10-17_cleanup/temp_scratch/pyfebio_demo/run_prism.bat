@echo off
cd /d D:\Project\climbing_fall_analysis
set LIG_MODE=shell
set DRIVER=disp
set TOPFIX=lateral
set TRUNC_MM=60
set FIX_TALUS=all
set CART_ELEM=prism
set TIME_STEPS=25
set DISP_MM=0.65
.venv\Scripts\python.exe -u scripts\ankle_fe\ankle_lig_feb.py
echo FEBIO_EXIT=%ERRORLEVEL%
