@echo off
REM disp 控制加到 ~600N（TRUNC_MM=60 后只要几秒），用于核对刚度 + 接触面积
cd /d D:\Project\climbing_fall_analysis
set LIG_MODE=shell
set DRIVER=disp
set TOPFIX=lateral
set TRUNC_MM=60
set FIX_TALUS=all
set TIME_STEPS=20
set DISP_MM=0.60
.venv\Scripts\python.exe -u scripts\ankle_fe\ankle_lig_feb.py
echo FEBIO_EXIT=%ERRORLEVEL%
