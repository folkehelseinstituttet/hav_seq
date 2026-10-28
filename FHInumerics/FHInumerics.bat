@echo off
REM Double-click to start FHInumerics. The first start installs what it needs
REM (5-15 minutes, needs the internet), then the app opens in the browser.
REM To put it on the desktop, right-click this file > Send to > Desktop (create shortcut).
REM If you copy the file somewhere else instead, set APP_DIR below to the FHInumerics folder.
setlocal
set "APP_DIR=%~dp0."
title FHInumerics - close this window to stop
if not exist "%APP_DIR%\launcher\launch.R" (
  echo.
  echo FHInumerics cannot find its files.
  echo If you opened it from a .zip file: right-click the .zip file, choose
  echo "Extract All...", and start FHInumerics from the new folder.
  echo.
  pause
  exit /b 1
)
call "%APP_DIR%\launcher\find-r.bat" || exit /b 1
echo Starting FHInumerics... a browser tab will open.
"%RSCRIPT%" "%APP_DIR%\launcher\launch.R"
if errorlevel 1 pause
