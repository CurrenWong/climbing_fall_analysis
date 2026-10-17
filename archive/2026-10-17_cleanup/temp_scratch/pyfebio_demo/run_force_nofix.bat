@echo off
REM force 驱动 + 顶面不约束（吴恺配置）——让外翻/内旋自由发生
cd /d D:\Project\climbing_fall_analysis
set LIG_MODE=shell
set DRIVER=force
set TOPFIX=none
set TRUNC_MM=60
set FIX_TALUS=all
set TIME_STEPS=20
set FORCE_N=600
.venv\Scripts\python.exe -u scripts\ankle_fe\ankle_lig_feb.py
echo FEBIO_EXIT=%ERRORLEVEL%
