FRAKTALFAKTUR RH LAB GPU — MERGED CATALOGUE v1.8.4
===================================================

This is an additive SOURCE MERGE of the user's original v1.8.0 D3D11
C++17 renderer and the verified v1.8.3 eight-candidate extension.
One launcher: run.cmd; build.cmd is its internal native compiler script.
No Python, web access, Lake, WSL, or Linux needed to RUN this package.

WINDOWS CMD (from Downloads)
  cd /d "%USERPROFILE%\Downloads"
  "%SystemRoot%\System32\tar.exe" -xf "FRAKTALFAKTUR_RH_LAB_GPU_v1_8_4.zip"
  call "FRAKTALFAKTUR_RH_LAB_GPU_v1_8_4\run.cmd"

WHAT WAS PRESERVED
  All 18 original RH / angular / spectral / Cayley / Tricorn /
  finite-Taylor candidate drivers (original registry IDs 0..17)
  Original native Direct3D11 compute rendering and CPU worker scheduler
  Original responsive zoom/pan, live reprojection, progressive refinement,
  cache, GPU fallbacks, screen/layout, presets, save/load, BMP and engine probe.
  Original mathematical symbols/models are unchanged except registration,
  UI routing and the original 18-count self-test (now validates all 26).

EIGHT NEW SELECTABLE DRIVERS (catalogue positions 19..26)
  19  extra_18   Cayley commutator      Mobius-group commutator -> square map
  20  extra_19   Tricorn Cayley        Antiholomorphic Cayley-square map
  21  extra_20   GP echo memory        Two-step complex orbit, four real coords
  22  extra_21   Spectral Newton web   Newton map of five numerical zero heights
  23  extra_22   Halley quintic web    Quintic Halley iteration
  24  extra_23   Entire Collatz parity Entire complex interpolation (lambda=0)
  25  extra_24   Dirichlet eta Newton  Newton for finite eta_9 Dirichlet series
  26  extra_25   Schwarzian pole web   Complex rational Schwarzian correction

Use the catalogue (Tab, arrows, Enter) or launch a candidate directly:
  call run.cmd --renderer extra_19
  call run.cmd --renderer extra_25

CONTROLS
  Wheel zoom; drag pan; right-click select a parameter or show an orbit.
  C/B selects the parameter plane when supported; J selects z-start.
  1-4 presets; P cycles three new color palettes; R resets view;
  S exports a BMP; Ctrl+S / Ctrl+O save/load .view files (new format 17);
  Tab returns to the original catalogue; Esc closes the application.
  Root-start-only maps 22–25 have no variable c parameter plane.

BACKENDS
  GPU: Original D3D11 / compute-shader-5 pipeline. For family 25 the
  additional HLSL compute shader is runtime-compiled (no GLSL or new window).
  CPU: Original parallel C++17 renderer, floating-point double evaluation
  for the added candidates. At deep zoom, the original engine chooses CPU.
  Compiler: existing MSYS2 UCRT64 x86_64-w64-mingw32 g++.

SELF-TEST & BUILD STATUS
  Local Linux C++17 mathematical and integration tests PASS: 189 checks,
  including original 18 entries, all 8 additions, hash parity at 1 and
  4 threads, camera defaults and eight view-file roundtrips.
  Prior Windows v1.8.3 standalone addon GPU HLSL test PASS (ps_5_0).
  The combined v1.8.4 Windows binary and its new compute shader cs_5_0
  have NOT been compiled or run here. run.cmd compiles/tests on Windows.
  If BUILD=FAIL, build.log prints the native diagnostics.

Original mathematical caution: All maps are exploratory finite-iteration
complex dynamics, not proofs of the Riemann hypothesis or Collatz conjecture;
finite-height spectral roots are parameters of an artificial polynomial.

Original documentation: README_BASE_v1_8_0.txt, RH_BRANCH_MAPS.txt,
NEW_STRUCTURES.txt, PERFORMANCE.txt, FINDINGS.txt, VALIDATION.txt.
Added documentation: MERGE_VALIDATION.txt, NEW_GENERATION_MATH.txt.
