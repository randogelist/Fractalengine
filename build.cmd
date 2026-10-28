@echo off
setlocal EnableExtensions DisableDelayedExpansion
set "fraktal_build_dir=%~dp0"
pushd "%fraktal_build_dir%"
if errorlevel 1 goto folder_error

echo BUILD=ENTER VERSION=1.8.6
set "fraktal_exit=1"
set "fraktal_cxx="
rem Invalidate the receipt before any action that can fail.
if exist "fraktalfaktur-lab-gpu-v1.8.6.ok" del /q "fraktalfaktur-lab-gpu-v1.8.6.ok"
if exist "fraktalfaktur-lab-gpu-v1.8.6.ok" goto marker_error
if exist "fraktalfaktur-build.exe" del /q "fraktalfaktur-build.exe"
if exist "fraktalfaktur-build.exe" goto temp_error
if not exist "src\amp_native.cpp" goto source_error

rem Use only a compiler already installed on this Windows machine.
rem FRAKTAL_CXX accepts one executable path, without flags.
if defined FRAKTAL_CXX goto override
if exist "C:\msys64\ucrt64\bin\g++.exe" set "fraktal_cxx=C:\msys64\ucrt64\bin\g++.exe"
if defined fraktal_cxx goto compiler_ready
if exist "D:\msys64\ucrt64\bin\g++.exe" set "fraktal_cxx=D:\msys64\ucrt64\bin\g++.exe"
if defined fraktal_cxx goto compiler_ready
for /f "delims=" %%G in ('where g++.exe 2^>nul') do if not defined fraktal_cxx set "fraktal_cxx=%%G"
if defined fraktal_cxx goto compiler_ready
if exist "C:\msys64\mingw64\bin\g++.exe" set "fraktal_cxx=C:\msys64\mingw64\bin\g++.exe"
if defined fraktal_cxx goto compiler_ready
if exist "D:\msys64\mingw64\bin\g++.exe" set "fraktal_cxx=D:\msys64\mingw64\bin\g++.exe"
if defined fraktal_cxx goto compiler_ready
goto compiler_missing

:override
set "fraktal_cxx=%FRAKTAL_CXX:"=%"
if not exist "%fraktal_cxx%" goto override_missing

:compiler_ready
for %%G in ("%fraktal_cxx%") do set "fraktal_compiler_dir=%%~dpG"
set "PATH=%fraktal_compiler_dir%;%PATH%"
>"build.log" echo PACKAGE=FRAKTALFAKTUR_RH_LAB_GPU_v1_8_6
if errorlevel 1 goto log_error
>>"build.log" echo COMPILER="%fraktal_cxx%"
"%fraktal_cxx%" -dumpmachine >"fraktal-lab-gpu-target.tmp" 2>>"build.log"
set "fraktal_exit=%errorlevel%"
if not "%fraktal_exit%"=="0" goto compiler_failed
findstr /i /c:"x86_64-w64-mingw32" "fraktal-lab-gpu-target.tmp" >nul
if errorlevel 1 goto compiler_wrong_target
type "fraktal-lab-gpu-target.tmp" >>"build.log"
del /q "fraktal-lab-gpu-target.tmp" >nul 2>nul

echo BUILD=START VERSION=1.8.6
echo COMPILER="%fraktal_cxx%"
echo PROFILE=NATIVE_CPU_STATIC_RUNTIME
echo Build output is being written to build.log. Please wait.
rem No fast-math: preserve finite-value and pole handling.
"%fraktal_cxx%" -std=c++17 -O3 -march=native -mtune=native -static -static-libgcc -static-libstdc++ -municode -o "fraktalfaktur-build.exe" "src\amp_native.cpp" -ld3d11 -ld3dcompiler -ldxgi -lgdi32 -luser32 -lcomdlg32 -lshell32 -lole32 -luuid >>"build.log" 2>&1
set "fraktal_exit=%errorlevel%"
if not "%fraktal_exit%"=="0" goto compile_error
if not exist "fraktalfaktur-build.exe" goto output_error
echo SELF_TEST=START
".\fraktalfaktur-build.exe" --self-test >>"build.log" 2>&1
set "fraktal_exit=%errorlevel%"
if not "%fraktal_exit%"=="0" goto test_error
echo SELF_TEST=PASS
echo GPU_SELF_TEST=START
".\fraktalfaktur-build.exe" --gpu-self-test >>"build.log" 2>&1
set "fraktal_exit=%errorlevel%"
if "%fraktal_exit%"=="2" goto gpu_unavailable
if not "%fraktal_exit%"=="0" goto gpu_test_error
echo GPU_SELF_TEST=PASS
goto tests_done
:gpu_unavailable
echo GPU_SELF_TEST=UNAVAILABLE - CPU fallback remains available.
set "fraktal_exit=0"
:tests_done
move /y "fraktalfaktur-build.exe" "fraktalfaktur.exe" >nul
if errorlevel 1 goto replace_error
(
for /r "src" %%F in (*.cpp *.hpp *.h *.inl) do @echo "%%~fF" %%~zF "%%~tF"
for %%F in ("build.cmd") do @echo "%%~fF" %%~zF "%%~tF"
) >"fraktal-lab-gpu-sources.sig"
if errorlevel 1 goto marker_error
>"fraktalfaktur-lab-gpu-v1.8.6.ok" echo FRAKTALFAKTUR_RH_LAB_GPU_v1_8_6
if errorlevel 1 goto marker_error
if not exist "fraktalfaktur-lab-gpu-v1.8.6.ok" goto marker_error
>>"build.log" echo BUILD=PASS VERSION=1.8.6
echo BUILD=PASS VERSION=1.8.6
popd
exit /b 0

:compile_error
echo BUILD=FAIL
echo ERROR: Your compiler returned exit code %fraktal_exit%.
type "build.log"
goto fail

:test_error
echo SELF_TEST=FAIL
echo ERROR: Bundled mathematical checks returned exit code %fraktal_exit%.
type "build.log"
goto fail

:gpu_test_error
echo GPU_SELF_TEST=FAIL
echo ERROR: GPU shader or short-orbit comparison failed. Read the details below.
type "build.log"
goto fail

:compiler_failed
echo ERROR: Your compiler could not report its Windows target.
type "build.log"
echo Check the existing compiler installation and FRAKTAL_CXX path.
goto fail

:compiler_missing
echo ERROR: Your native Windows x64 MinGW-w64 g++ was not found.
echo Set the path to your existing compiler, then run again:
echo   set "FRAKTAL_CXX=C:\msys64\ucrt64\bin\g++.exe"
echo   call run.cmd
goto fail_one

:override_missing
echo ERROR: FRAKTAL_CXX does not point to an existing compiler file.
echo Path: "%fraktal_cxx%"
echo Set FRAKTAL_CXX to the full path of your native Windows g++.exe.
goto fail_one

:compiler_wrong_target
echo ERROR: The selected compiler does not target Windows x64 MinGW-w64.
echo Path: "%fraktal_cxx%"
type "fraktal-lab-gpu-target.tmp"
echo Set FRAKTAL_CXX to your MSYS2 ucrt64 or mingw64 g++.exe and retry.
goto fail_one

:source_error
echo ERROR: src\amp_native.cpp is missing. Extract the complete ZIP and retry.
goto fail_one

:temp_error
echo ERROR: Cannot remove a previous temporary executable.
echo Close the runner and check this folder is writable, then retry.
goto fail_one

:output_error
echo ERROR: The compiler returned success but produced no executable.
goto fail_one

:replace_error
echo ERROR: Cannot replace fraktalfaktur.exe.
echo Close every runner window and check this folder is writable, then retry.
goto fail_one

:marker_error
echo ERROR: Cannot update the build cache or success marker.
echo Check this folder is writable, then retry.
goto fail_one

:log_error
echo ERROR: Cannot write build.log. Check this folder is writable, then retry.
goto fail_one

:fail_one
set "fraktal_exit=1"
:fail
if exist "fraktalfaktur-build.exe" del /q "fraktalfaktur-build.exe" >nul 2>nul
if exist "fraktal-lab-gpu-target.tmp" del /q "fraktal-lab-gpu-target.tmp" >nul 2>nul
if exist "fraktalfaktur-lab-gpu-v1.8.6.ok" del /q "fraktalfaktur-lab-gpu-v1.8.6.ok" >nul 2>nul
popd
exit /b %fraktal_exit%

:folder_error
echo ERROR: Cannot open the extracted runner folder.
echo Folder: "%fraktal_build_dir%"
exit /b 1
