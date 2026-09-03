@echo off
setlocal
REM ====================================================================
REM  Launch 4 independent R processes, one per shard.
REM
REM  NOTE ON QUOTING: `cmd /k "C:\Program Files\..." args` fails, because
REM  cmd consumes the first quoted string as its own delimiter. The whole
REM  command must be wrapped in an EXTRA pair of quotes:
REM      cmd /k ""prog" "script" args"
REM  That is why the previous version opened four idle prompts.
REM ====================================================================

REM paths WITHOUT surrounding quotes - they are added at the call site
set RSCRIPT=C:\Program Files\R\R-4.3.1\bin\Rscript.exe
set SCRIPT=%~dp0run_shard.R

if not exist "%RSCRIPT%" (
  echo.
  echo Rscript.exe not found at:
  echo   %RSCRIPT%
  echo.
  echo Installed R versions:
  dir /b "C:\Program Files\R" 2^>nul
  echo.
  echo Edit the 'set RSCRIPT=' line in this file to match, then rerun.
  pause
  exit /b 1
)

if not exist "%SCRIPT%" (
  echo run_shard.R not found next to this .bat  ^(%SCRIPT%^)
  pause
  exit /b 1
)

echo Rscript : %RSCRIPT%
echo Script  : %SCRIPT%
echo.
echo Launching 4 shards. Each writes C:\raenr_out\shard_0N.rds
echo Leave all four windows open until they say SHARD COMPLETE.
echo.

start "RAENR shard 1 of 4" cmd /k ""%RSCRIPT%" "%SCRIPT%" 1 4"
start "RAENR shard 2 of 4" cmd /k ""%RSCRIPT%" "%SCRIPT%" 2 4"
start "RAENR shard 3 of 4" cmd /k ""%RSCRIPT%" "%SCRIPT%" 3 4"
start "RAENR shard 4 of 4" cmd /k ""%RSCRIPT%" "%SCRIPT%" 4 4"

echo All four launched. You can close this window.
pause
endlocal
