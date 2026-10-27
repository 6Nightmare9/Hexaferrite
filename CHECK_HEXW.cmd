@echo off
setlocal
cd /d "%~dp0"
if not exist "tools\venv\Scripts\python.exe" goto setup
"%~dp0tools\venv\Scripts\python.exe" -B tools\check_installation.py
set "HEXW_CHECK_EXIT=%errorlevel%"
pause
exit /b %HEXW_CHECK_EXIT%
:setup
echo First run SETUP_HEXW.cmd.
pause
exit /b 1
