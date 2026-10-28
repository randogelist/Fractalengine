@echo off
setlocal EnableExtensions DisableDelayedExpansion
echo ENGINE_RUN=CPU_GPU_PARALLELISM
echo This runs native measurements and writes engine-report.txt and engine-samples.csv.
call "%~dp0run.cmd" --engine-probe --probe-full --probe-seconds 120 %*
set "fraktal_probe_exit=%errorlevel%"
if not "%fraktal_probe_exit%"=="0" goto failed
echo.
echo Paste the console output or send this report:
echo "%~dp0engine-report.txt"
exit /b 0
:failed
echo ENGINE_RUN=FAILED
echo Read build.log or engine-report.txt in the extracted folder.
exit /b %fraktal_probe_exit%
