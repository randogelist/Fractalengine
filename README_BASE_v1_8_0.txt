FRAKTALFAKTUR RH LAB GPU 1.8.0 - SPECTRAL ECHO GENERATION

NATIVE ENGINE
Windows x64, native C++17, Direct3D 11 / HLSL Compute Shader 5,
CPU parallel scheduler, guarded CPU precision repair and progressive
reprojection/refinement. Based directly on FRAKTALFAKTUR_RH_LAB_GPU_v1_7_2.
One build/launch interface: run.cmd. Existing compiler is used (the current
MSYS2 UCRT64 GCC is detected). No Python files, Python runtime, installation
or network access required. The archive contains source, not binaries.

START (Windows CMD, from the parent folder containing this ZIP)
  tar -xf "FRAKTALFAKTUR_RH_LAB_GPU_v1_8_0.zip"
  call "FRAKTALFAKTUR_RH_LAB_GPU_v1_8_0\run.cmd"

18 REGISTERED RH DRIVERS; FIRST 13 RETAINED WITHOUT FORMULA CHANGES.
Catalogue opens with driver 14 highlighted (Spectral echo determinant).
Select by Enter/double-click. Use Tab for the catalogue; F11 fullscreen;
mouse wheel zooms, drag pans, S exports image, Ctrl+S/O saves/loads views.

ADDED MATHEMATICAL DRIVERS (each has an independent exact Newton target)
14 spectral_echo:            finite-shift correlation determinant
15 laguerre_gram:            infinitesimal Laguerre/log-curvature determinant
16 dyadic_holonomy:         nonlocal multiplicative dilation defect
17 antipodal_wronskian:    odd reflected spectral current
18 projective_schwarzian:  numerator of projective Schwarzian curvature

Curvature/Angular/Eq branches 1-13 and the three focus banks are retained.
Focus bank T switches between Original Lab / Currents / Weaves. 1-4 applies
one of four presets. Native Cayley currents, Tricorn weave and all twelve
feedback carriers apply to the five new analytic targets as before.

MODEL
  T_N(y) = sum(a_k*y^k,k=0..N) = Theta_N(200*y)/Theta(0)
  y = b*z+h,  b>0
  feedback = -lambda * F(y)/(b*F'(y))

The default host is z_next=z^2+c+feedback. Other carriers include
conjugation/Tricorn, Cayley multiplication, memory, twist and blends.
Degree N = 8/16/32. 14's real-shift R uses the existing radius slider.
For 15-18 the R slider is irrelevant to the exact target; b/h, c and
complex lambda still control the images. The Level slider is irrelevant
to all five new drivers. There are no hidden degree-dependent levels.

Use --renderer spectral_echo, laguerre_gram, dyadic_holonomy,
antipodal_wronskian or projective_schwarzian to launch directly.

ENGINE / PRECISION
Same native CPU/GPU scheduling, progressive zoom, live responsive
reprojection, shader precompile and cached orbit metrics as v1.7.2.
The C++ evaluator uses long-double, shared-normalization Horner jets.
The GPU uses float32 guarded Horner and requests CPU repair for poorly
conditioned values; it does not clamp or insert an epsilon into poles.
New radius upload occupies an unused coefficient-atlas metadata texel;
constant-buffer ABI remains unchanged. Renderer cache is radius-aware.

STARTUP VALIDATION
build.cmd performs the full mathematical --self-test and --gpu-self-test
before installing the compiled executable. With no D3D11 device the CPU
fallback is permitted. The new suite checks exact finite-polynomial target
formulas, analytic derivatives, symmetry, echo limit, Möbius-invariant
Schwarzian, GPU/CPU parity, Tricorn/Cayley, charts, recolour cache and
unrepaired GPU pixel coverage. Read build.log on failures.

This is exploratory finite-Taylor complex dynamics, NOT an RH/GRH proof,
a root certificate, an assertion of positive Gram matrices, or a claim
of new mathematical theorems. Radius-16 termination is a visualization
cutoff, not a certified escape radius. Hardware compilation and native
GPU tests require your Windows machine and have not been run here.

OTHER FILES
NEW_STRUCTURES.txt: derivations of all five targets
RH_BRANCH_MAPS.txt: previous mathematical drivers
PATCH_NOTES.txt: changed implementation surfaces
VALIDATION.txt: separate local vs on-device verification receipts
PREVIEW_ATLAS.png: portable C++ Julia-boundary preview of five new modes
run_engine.cmd: retained engine CPU/GPU performance audit command
