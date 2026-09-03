@echo off
setlocal
REM ====================================================================
REM  Rebuild ONLY p = 20,000 for shards 1, 2 and 3, one at a time.
REM
REM  Shard 4 already has all five dimensionalities and is not touched.
REM
REM  SEQUENTIAL BY DESIGN. Three of four shards died at p = 20,000 while
REM  running concurrently, and the survivor ran faster once they had
REM  gone. Running these one at a time avoids repeating that.
REM
REM  Expect roughly 2 hours per shard on the evidence of shard 4
REM  (130 minutes for 250 replications at p = 20,000), so about 6 hours
REM  in total. Leave this window open.
REM ====================================================================

set RSCRIPT=C:\Program Files\R\R-4.3.1\bin\Rscript.exe
set SCRIPT=%~dp0topup_p20000.R
set GSE14520_PATH=C:\Users\alabi\Downloads\GSE14520-GPL3921_series_matrix.txt.gz
set RAENR_OUT=%~dp0raenr_out

cd /d "%~dp0"

if not exist "%RSCRIPT%" (
  echo Rscript.exe not found at %RSCRIPT%
  dir /b "C:\Program Files\R" 2^>nul
  pause
  exit /b 1
)
if not exist "%SCRIPT%"  ( echo topup_p20000.R not found & pause & exit /b 1 )
if not exist "%GSE14520_PATH%" (
  echo GEO series matrix not found at %GSE14520_PATH%
  pause
  exit /b 1
)

findstr /C:"lambda_fixed * n" "%~dp0run_shard.R" >nul
if errorlevel 1 (
  echo *** run_shard.R is not patched. Refusing to launch. ***
  pause
  exit /b 1
)

echo Patch check : PASSED
echo Output      : %RAENR_OUT%
echo.

for %%S in (1 2 3) do (
  echo ============================================================
  echo  Shard %%S  -  p = 20000 only
  echo ============================================================
  "%RSCRIPT%" "%SCRIPT%" %%S
  if errorlevel 1 (
    echo.
    echo Shard %%S FAILED. Stopping so the error is visible.
    pause
    exit /b 1
  )
  echo.
)

echo All three shards topped up.
echo.
echo Now run in RStudio:
echo   setwd("%~dp0")
echo   Sys.setenv(RAENR_OUT = "%RAENR_OUT%")
echo   source("combine_shards.R")
pause
endlocal
