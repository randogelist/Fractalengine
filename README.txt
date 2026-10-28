FRAKTALFAKTUR RH LAB GPU 1.8.6 - QQRH HILBERT DEVELOPMENT

Extends the actual v1.8.5 source. All 29 previous entries remain available.
New entries in the original native C++/D3D11 renderer:
  30  extra_29  QQRH projective metric
  31  extra_30  QQRH overlap kernel Newton
  32  extra_31  QQRH kernel feedback fractal

Entry 27 (Hilbert overlap) now has a point inspector showing squared overlap,
projective ray angle (radians), local metric g, and 13 component probabilities.
Entry 30 has the same inspector. Right-click selects a point; I hides/shows it.

RUN IN WINDOWS CMD
Extract the full ZIP, then call run.cmd. It builds with your existing native
Windows C++ compiler, runs CPU and native GPU checks, and opens the catalogue.
No Python, PowerShell, WSL, Lean, download, or separate renderer installation.

  call run.cmd
  call run.cmd --renderer extra_26
  call run.cmd --renderer extra_29
  call run.cmd --renderer extra_30
  call run.cmd --renderer extra_31

CONTROLS
Wheel: zoom; drag: pan; Tab: catalogue; H: formulas/help; P: palette; R: reset.
R slider: scaled radius in [0.1,1.25]; physical Theta radius r=200R.
I: inspector in the two direct geometry fields; right-click: select point.
Kernel feedback and Hardy shift: C/B for parameter c portrait, J for Julia.
Right-click a c portrait to select c and open its Julia view.
Newton views: each pixel is the initial z. Fields have no orbit iteration.
1-4: named presets. Field presets: default, R=.25, R=1.10, z=1 close-up.
Newton presets: default, R=.35, R=1.10, boundary close-up.
Feedback/shift presets: default Julia, parameter c, selected Julia, complex Julia.
S: BMP image; Ctrl+S / Ctrl+O: save/load views.
The inspector selection is a UI readout; it is not stored in saved views/BMPs.

FINITE MODEL
N=12 means 13 coefficients, unchanged from v1.8.5. Both new fractal maps use
H(z)=sum v_n^2 z^n, the coherent-state overlap polynomial. The metric-only
field is radial geometry; the new Newton and feedback maps produce fractals.
See QQRH_HILBERT_DEVELOPMENT.txt for exact formulas and interpretation.
See QQRH_HILBERT_DEVELOPMENT_PREVIEW.png for actual CPU-rendered examples.

VALIDATION
See VALIDATION_1_8_6.txt and CPU_VERIFICATION_1_8_6.log.
Portable C++ checks passed. Native Windows compilation, GUI interaction,
HLSL compilation and D3D11 execution of this exact package remain unverified.
The finite numerical model does not certify full QQRH positivity or RH.
Older-version reports are inherited history, not current build receipts.
