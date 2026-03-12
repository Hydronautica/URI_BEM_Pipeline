@echo off
echo === Generating mesh and input files ===
python generate_case.py
if %errorlevel% neq 0 goto :error

echo.
echo === Running preProcessor ===
..\bin_windows\preProc.exe .
if %errorlevel% neq 0 goto :error

echo.
echo === Running solver ===
..\bin_windows\solver.exe .
if %errorlevel% neq 0 goto :error

echo.
echo === Running postProcessor ===
..\bin_windows\postProc.exe .
if %errorlevel% neq 0 goto :error

echo.
echo === ALL DONE ===
pause
goto :eof

:error
echo.
echo !!! ERROR at previous step !!!
pause
