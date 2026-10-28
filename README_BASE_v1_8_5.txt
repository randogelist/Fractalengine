FRAKTALFAKTUR RH LAB GPU 1.8.5 - QQRH HILBERT GEOMETRY

Built on the actual integrated v1.8.4 source. All 26 existing entries remain.
Three additions in the same native D3D11/C++ renderer:
  27  extra_26  QQRH Hilbert overlap       Direct finite Hardy-space geometry
  28  extra_27  QQRH Hardy shift fractal  Correlation-driven quadratic dynamics
  29  extra_28  QQRH angular kernel Newton  Finite angular-derivative basins

RUN
Extract the full ZIP and call run.cmd. It uses your existing Windows C++
compiler, runs CPU tests and native GPU checks, and opens the catalogue.
No Python, PowerShell, WSL, download, or separate renderer installation.

  call run.cmd
  call run.cmd --renderer extra_26
  call run.cmd --renderer extra_27
  call run.cmd --renderer extra_28

CONTROL
Wheel zoom; drag pan; Tab catalogue; H formulas/help; P palette; R reset.
For new entries, the radius slider R ranges from .1 to 1.25; r=200R.
The overlap field has no iteration budget. Right-click reads coherence.
Shift fractal: C/B parameter portrait (seed0), J Julia at selected c.
Angular Newton: each pixel is an initial state, with no parameter-c plane.
1-4 presets; S BMP; Ctrl+S / Ctrl+O save/load views.
N=12 means 13 components. Existing controls and saved views are retained.

See QQRH_HILBERT_GEOMETRY.txt for exact formulas, scaling, and provenance.
See QQRH_HILBERT_PREVIEW.png for actual CPU-rendered default views.
See VALIDATION_1_8_5.txt and CPU_VERIFICATION_1_8_5.log for current checks.
Reports whose titles/filenames mention earlier versions are inherited history.

Portable C++ checks passed. Full native Windows compilation and live D3D11
execution are pending your run; no local Lean compilation was performed.
The displayed finite model has no full-Theta tail certification or RH claim.
