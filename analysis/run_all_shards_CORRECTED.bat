@echo off
setlocal
REM ====================================================================
REM  Launch the 4 shards of the CORRECTED run (fit_raen_mm now passes
REM  lambda_fixed * n to cd_sweep).
REM
REM  Differences from run_all_shards.bat:
REM    1. Sets GSE14520_PATH. The original relied on the series matrix
REM       sitting beside the script; it is in Downloads, so every shard
REM       would have stopped at startup with "GEO series matrix not
REM       found".
REM    2. Sets RAENR_OUT explicitly and cd's to the script directory, so
REM       the shards and combine_shards.R agree on where the .rds files
REM       go regardless of where this file is launched from.
REM    3. Moves any previous shard output aside rather than mixing old
REM       and new results in one directory.
REM ====================================================================

set RSCRIPT=C:\Program Files\R\R-4.3.1\bin\Rscript.exe
set SCRIPT=%~dp0run_shard.R
set GSE14520_PATH=C:\Users\alabi\Downloads\GSE14520-GPL3921_series_matrix.txt.gz
set RAENR_OUT=%~dp0raenr_out

cd /d "%~dp0"

if not exist "%RSCRIPT%" (
  echo Rscript.exe not found at:
  echo   %RSCRIPT%
  echo Installed R versions:
  dir /b "C:\Program Files\R" 2^>nul
  pause
  exit /b 1
)

if not exist "%SCRIPT%" (
  echo run_shard.R not found next to this .bat  ^(%SCRIPT%^)
  pause
  exit /b 1
)

if not exist "%GSE14520_PATH%" (
  echo.
  echo GEO series matrix not found at:
  echo   %GSE14520_PATH%
  echo Edit the 'set GSE14520_PATH=' line in this file to match.
  pause
  exit /b 1
)

REM --- confirm the patch is in the file being run ---------------------
findstr /C:"lambda_fixed * n" "%SCRIPT%" >nul
if errorlevel 1 (
  echo.
  echo *** run_shard.R does not contain the corrected penalty scaling. ***
  echo Expected: cd_sweep^(Xw, yw, beta, wp, lambda_fixed * n, al^)
  echo Refusing to launch - this would reproduce the defective results.
  pause
  exit /b 1
)

REM --- move any earlier shard output aside ----------------------------
if exist "%RAENR_OUT%\shard_01.rds" (
  set STAMP=%DATE:/=-%_%TIME::=-%
  echo Existing shard output found. Moving it to raenr_out_previous.
  if not exist "%~dp0raenr_out_previous" mkdir "%~dp0raenr_out_previous"
  move /Y "%RAENR_OUT%\shard_*.rds" "%~dp0raenr_out_previous\" >nul
)

echo Rscript      : %RSCRIPT%
echo Script       : %SCRIPT%
echo GEO matrix   : %GSE14520_PATH%
echo Output       : %RAENR_OUT%
echo Patch check  : PASSED
echo.
echo Launching 4 shards. Expect roughly 3 to 4 hours.
echo Leave all four windows open until each says SHARD COMPLETE.
echo.

start "RAENR shard 1 of 4" cmd /k ""%RSCRIPT%" "%SCRIPT%" 1 4"
start "RAENR shard 2 of 4" cmd /k ""%RSCRIPT%" "%SCRIPT%" 2 4"
start "RAENR shard 3 of 4" cmd /k ""%RSCRIPT%" "%SCRIPT%" 3 4"
start "RAENR shard 4 of 4" cmd /k ""%RSCRIPT%" "%SCRIPT%" 4 4"

echo All four launched.
echo.
echo When every window says SHARD COMPLETE, run in RStudio:
echo   setwd("%~dp0")
echo   Sys.setenv(RAENR_OUT = "%RAENR_OUT%")
echo   source("combine_shards.R")
pause
endlocal
