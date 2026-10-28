@echo off
REM Sets RSCRIPT to Rscript.exe: from the PATH, else the newest R in Program Files
REM or in the per-user install folder. Used by FHInumerics.bat (one folder up).
set "RSCRIPT="
for /f "delims=" %%i in ('where Rscript.exe 2^>nul') do if not defined RSCRIPT set "RSCRIPT=%%i"
if not defined RSCRIPT for /d %%d in ("%ProgramFiles%\R\R-*" "%LOCALAPPDATA%\Programs\R\R-*") do if exist "%%d\bin\Rscript.exe" set "RSCRIPT=%%d\bin\Rscript.exe"
if not defined RSCRIPT (
  echo.
  echo ========================================================================
  echo FHInumerics could not start: R is not installed on this computer.
  echo.
  echo FHInumerics needs R, version 4.4 or newer. Ask IT to install it, or
  echo install it from https://cran.r-project.org - then start FHInumerics again.
  echo ========================================================================
  echo.
  pause
  exit /b 1
)
exit /b 0
