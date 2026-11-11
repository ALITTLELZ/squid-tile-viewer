@echo off
REM Double-click to start Tile Viewer (Windows)
cd /d "%~dp0"
set "UV="
where uv >nul 2>nul && set "UV=uv"
if not defined UV if exist "%USERPROFILE%\.local\bin\uv.exe" set "UV=%USERPROFILE%\.local\bin\uv.exe"
if defined UV goto run
echo Tile Viewer uses "uv" (a small Python installer from astral.sh) to set up its dependencies.
choice /M "Install uv now"
if errorlevel 2 exit /b 0
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
set "UV=%USERPROFILE%\.local\bin\uv.exe"
:run
echo Starting Tile Viewer ... (first launch downloads ~300 MB of packages and may take a few minutes)
"%UV%" run --python 3.12 tile_viewer.py %*
if errorlevel 1 pause
