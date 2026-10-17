@echo off
rem Ankle ligament runner. Edit the CONFIG block below, not the command line.
cd /d D:\Project\climbing_fall_analysis

rem --- CONFIG ---
set "LIG_MODE=shell"
set "DRIVER=force"
set "TRUNC_MM=60"
set "FIX_TALUS=transverse"
set "TOPFIX=lateral"
set "TIME_STEPS=20"
set "MATSET=WK"
rem --------------

.venv\Scripts\python.exe -u scripts\ankle_fe\ankle_lig_feb.py
