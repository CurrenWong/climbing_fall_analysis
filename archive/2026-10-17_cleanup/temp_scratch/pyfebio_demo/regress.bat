@echo off
cd /d D:\Project\climbing_fall_analysis
set DRIVER=disp
set DISP_MM=0.1472
set TIME_STEPS=20
.venv\Scripts\python.exe -u scripts\ankle_fe\build_feb.py
echo EXIT=%ERRORLEVEL%
