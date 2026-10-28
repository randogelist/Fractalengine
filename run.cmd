@echo off
setlocal EnableExtensions DisableDelayedExpansion
set "fraktal_run_dir=%~dp0"
pushd "%fraktal_run_dir%"
if errorlevel 1 goto folder_error

echo LAUNCHER=RH_LAB_GPU_V1_8_6
echo DEFAULT_SCREEN=CANDIDATE_CATALOGUE
if not exist "fraktalfaktur.exe" goto build
if not exist "fraktalfaktur-lab-gpu-v1.8.6.ok" goto build
if not exist "fraktal-lab-gpu-sources.sig" goto build
rem CMD-only dependency signature. build.cmd can always force a full rebuild.
(
for /r "src" %%F in (*.cpp *.hpp *.h *.inl) do @echo "%%~fF" %%~zF "%%~tF"
for %%F in ("build.cmd") do @echo "%%~fF" %%~zF "%%~tF"
) >"fraktal-lab-gpu-sources.check.tmp"
if errorlevel 1 goto cache_error
fc /b "fraktal-lab-gpu-sources.sig" "fraktal-lab-gpu-sources.check.tmp" >nul 2>nul
set "fraktal_cache_exit=%errorlevel%"
del /q "fraktal-lab-gpu-sources.check.tmp" >nul 2>nul
if not "%fraktal_cache_exit%"=="0" goto build
goto launch

:build
if not exist ".\build.cmd" goto missing_build
call ".\build.cmd"
set "fraktal_exit=%errorlevel%"
if not "%fraktal_exit%"=="0" goto build_error
if not exist "fraktalfaktur.exe" goto missing_exe
if not exist "fraktalfaktur-lab-gpu-v1.8.6.ok" goto missing_exe

:launch
rem Explicit arguments come last and can override launcher defaults.
if defined FRAKTAL_THREADS goto launch_threads
".\fraktalfaktur.exe" --refine fast %*
set "fraktal_exit=%errorlevel%"
goto launch_done

:launch_threads
".\fraktalfaktur.exe" --refine fast --threads "%FRAKTAL_THREADS%" %*
set "fraktal_exit=%errorlevel%"

:launch_done
popd
exit /b %fraktal_exit%

:build_error
echo ERROR: Build or mathematical self-test failed. The runner was not launched.
echo Read build.log for the compiler or test result.
if exist "build.log" echo BUILD_LOG=build.log
popd
exit /b %fraktal_exit%

:missing_build
echo ERROR: build.cmd is missing. Extract the complete ZIP and retry.
popd
exit /b 1

:missing_exe
echo ERROR: Build success marker or executable is missing.
popd
exit /b 1

:cache_error
echo ERROR: Cannot check source changes. Check this folder is writable.
if exist "fraktal-lab-gpu-sources.check.tmp" del /q "fraktal-lab-gpu-sources.check.tmp" >nul 2>nul
popd
exit /b 1

:folder_error
echo ERROR: Cannot open the extracted runner folder.
echo Folder: "%fraktal_run_dir%"
exit /b 1
