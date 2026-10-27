@echo off
set "HEXW_ROOT=%~dp0"
set OMP_NUM_THREADS=2
set OPENBLAS_NUM_THREADS=2
set MKL_NUM_THREADS=2
set "MPLCONFIGDIR=%~dp0tools\venv\mpl_config"
cd /d "%~dp0project"
echo Open http://127.0.0.1:8765 in your browser. Keep this window open.
"%~dp0tools\venv\Scripts\python.exe" -B -m hexw serve --port 8765
pause
