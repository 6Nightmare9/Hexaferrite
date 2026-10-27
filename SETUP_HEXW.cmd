@echo off
setlocal
cd /d "%~dp0"
where git >nul 2>nul
if errorlevel 1 goto no_git
git lfs install --local
if errorlevel 1 goto failed
git lfs pull
if errorlevel 1 goto failed
if exist "tools\venv\Scripts\python.exe" goto install
where py >nul 2>nul
if errorlevel 1 goto no_python
py -3.14 -m venv tools\venv
if errorlevel 1 goto no_python
:install
cd /d "%~dp0project"
"%~dp0tools\venv\Scripts\python.exe" -m pip install -r requirements-test.txt
if errorlevel 1 goto failed
"%~dp0tools\venv\Scripts\python.exe" -m pip check
if errorlevel 1 goto failed
cd /d "%~dp0"
if not exist "%~dp0registry\FROZEN_MODEL.json" copy /y "%~dp0audit\initial_FROZEN_MODEL.json" "%~dp0registry\FROZEN_MODEL.json" >nul
"%~dp0tools\venv\Scripts\python.exe" -B tools\check_installation.py
if errorlevel 1 goto failed
echo.
echo Setup complete. Double-click START_HEXW.cmd, then open http://127.0.0.1:8765/
pause
exit /b 0
:no_git
echo Install Git for Windows with Git LFS, then clone this repository and retry.
goto failed
:no_python
echo Install Python 3.14 with the Python launcher, then retry.
:failed
echo Setup did not finish. Read the message above. Saved calculations were not deleted.
pause
exit /b 1
