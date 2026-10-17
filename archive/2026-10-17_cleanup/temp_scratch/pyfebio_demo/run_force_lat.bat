@echo off
REM 力控制 + 顶面横向约束（允许绕轴扭转）+ 细步长 + 放宽 etol
cd /d D:\Project\climbing_fall_analysis
set LIG_MODE=shell
set DRIVER=force
set TOPFIX=lateral
set TRUNC_MM=60
set FIX_TALUS=all
set TIME_STEPS=40
set FORCE_N=600
set ETOL=0.01
set MAX_REFS=30
.venv\Scripts\python.exe -u scripts\ankle_fe\ankle_lig_feb.py
echo FEBIO_EXIT=%ERRORLEVEL%
